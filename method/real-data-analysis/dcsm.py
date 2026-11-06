"""DCSM comparison with three Weibull experts, fitted on the primary cohort."""
import os, sys, copy
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.special import gamma
from sklearn.model_selection import GroupShuffleSplit
upstream = Path(os.environ.get('DCSM_DIR', ''))
if not (upstream / 'models/dcsm_torch.py').is_file():
    raise RuntimeError('Set DCSM_DIR to a checkout of https://github.com/BojianHou/DCSM')
sys.path.insert(0, str(upstream))
from models.dcsm_torch import DeepClusteringSurvivalMachinesTorch
from utils.losses import conditional_loss
torch.set_num_threads(1)
torch.use_deterministic_algorithms(True)
SEEDS = [42, 73, 666, 777, 1009]
LABELS = ['Early-onset', 'Mid-onset', 'Late-onset']
class CPUDCSM(DeepClusteringSurvivalMachinesTorch):
    def forward(self, x, risk='1'):
        # Same fixed/shared experts as upstream fix=True, without forced .cuda().
        return (self.shape[risk].expand(len(x),-1), self.scale[risk].expand(len(x),-1),
                self.gate[risk](self.embedding(x)) / self.temp)

def make_model(p, seed, times):
    torch.manual_seed(seed)
    model = CPUDCSM(p, 3, layers=[16], temp=1, discount=1, fix=True).double()
    # Distinct training-time quantile initialization prevents identical experts.
    with torch.no_grad():
        model.shape['1'].fill_(np.log(2.0))
        model.scale['1'].copy_(torch.tensor(-np.log(np.quantile(times,[.2,.5,.8]))))
    return model

def fit(x, t, e, seed, valid=None, epochs=2500):
    model = make_model(x.shape[1], seed, t.numpy())
    opt = torch.optim.Adam(model.parameters(),lr=.005,weight_decay=.001)
    best, best_epoch, state, stale = np.inf, 0, None, 0
    for epoch in range(1,epochs+1):
        opt.zero_grad()
        loss = conditional_loss(model,x,t,e,elbo=False)
        if not torch.isfinite(loss): raise RuntimeError('Nonfinite training loss')
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),10)
        opt.step()
        with torch.no_grad():
            model.shape['1'].clamp_(-2,3)
            model.scale['1'].clamp_(-5,5)
        if valid is not None:
            with torch.no_grad(): value = float(conditional_loss(model,*valid,elbo=False))
            if value < best - 1e-6:
                best, best_epoch, state, stale = value, epoch, copy.deepcopy(model.state_dict()), 0
            else: stale += 1
            if epoch >= 200 and stale >= 200: break
    if valid is not None: model.load_state_dict(state)
    else: best_epoch = epochs
    return model, best_epoch, best

def predict(model,x,time_scale):
    with torch.no_grad():
        shape, scale, logits = model(x)
        a, b = np.exp(shape[0].numpy()), np.exp(scale[0].numpy())
        means = gamma(1+1/a)/b*time_scale
        order = np.argsort(means)
        prob = torch.softmax(logits,dim=1).numpy()[:,order]
    return prob, means[order], a[order], b[order]


df = pd.read_csv(sys.argv[1])
mutations=sorted(df.ProtDescription.unique());X=np.column_stack([(df.Sex=='M').astype(float),np.array([[float(m==v) for v in mutations] for m in df.ProtDescription])])
groups=df.Family_ID.to_numpy();dev,test=next(GroupShuffleSplit(n_splits=1,test_size=.25,random_state=20260909).split(X,groups=groups));itr,ival=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=20260910).split(X[dev],groups=groups[dev]));train,val=dev[itr],dev[ival]
assert not set(groups[train])&set(groups[val]) and not set(groups[dev])&set(groups[test])
def tensors(idx,scale):return tuple(torch.tensor(a,dtype=torch.float64) for a in (X[idx],df.time.to_numpy()[idx]/scale,df.status.to_numpy()[idx]))
scale=float(df.time.iloc[train].median());runs=[]
for seed in SEEDS:
 model,epoch,nll=fit(*tensors(train,scale),seed,valid=tensors(val,scale));runs.append({'seed':seed,'epochs':epoch,'validation_nll_scaled_time':nll});print(runs[-1],flush=True)
chosen=min(runs,key=lambda a:a['validation_nll_scaled_time'])
idx=np.arange(len(df));fs=float(df.time.median());full,_,_=fit(*tensors(idx,fs),chosen['seed'],epochs=chosen['epochs']);prob,means,shape,rate=predict(full,tensors(idx,fs)[0],fs)
df['DCSM']=[LABELS[i] for i in prob.argmax(1)]
for i,label in enumerate(LABELS):df['prob_'+label]=prob[:,i]
df['assignment_probability']=prob.max(1);df['expected_age']=prob@means
ann=df.drop_duplicates('ProtDescription').set_index('ProtDescription')
counts=pd.crosstab(df.ProtDescription,df.DCSM).reindex(columns=LABELS,fill_value=0);s=counts.copy();s['patients']=counts.sum(axis=1);s['majority_DCSM']=counts.idxmax(axis=1);tie=counts.eq(counts.max(axis=1),axis=0).sum(axis=1)>1;s.loc[tie,'majority_DCSM']='Tied';s['n_patient_clusters']=(counts>0).sum(axis=1);s['YTAClass']=ann.YTAClass;s.to_csv(sys.argv[2])
