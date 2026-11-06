# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/figures.R'))
results <- readRDS(file.path(folder,'Supplement-Figure-A7.rds'))
keys<-c('n_S','D_target','rep');metrics<-results$metrics;diag<-results$diagnostics
stopifnot(nrow(diag)==600,!anyDuplicated(diag[keys]),nrow(metrics)==1800,!anyDuplicated(metrics[c(keys,'method')]))
d<-merge(metrics,diag,by=keys,all=TRUE);stopifnot(nrow(d)==1800,all(complete.cases(d)),all(table(d$n_S,d$D_target,d$method)==50))
d$method<-factor(d$method,levels=c('Cox-K-means','SCP','TL-SCP'));colors<-c('Cox-K-means'='black',SCP='red2','TL-SCP'='purple2')
write.csv(d,file.path(folder,'Supplement-Figure-A7-performance.csv'),row.names=FALSE);write.csv(diag,file.path(folder,'Supplement-Figure-A7-ranking-diagnostics.csv'),row.names=FALSE)
# Use the same fitted curve to draw the ranking margin and locate its zero.
fit_loess<-function(a,y,span=.8)loess(as.formula(paste(y,'~D_hat')),data=a,span=span,degree=2,family='gaussian',control=loess.control(surface='interpolate',statistics='approximate'))
find_roots<-function(fit,xmin,xmax,n=10001){
 grid<-seq(xmin,xmax,length.out=n);val<-as.numeric(predict(fit,newdata=data.frame(D_hat=grid)))
 stopifnot(all(is.finite(val)))
 ix<-which(val[-length(val)]*val[-1]<0)
 roots<-vapply(ix,function(i)uniroot(function(x)as.numeric(predict(fit,newdata=data.frame(D_hat=x))),c(grid[i],grid[i+1]),tol=1e-12)$root,numeric(1))
 if(any(val==0))roots<-c(roots,grid[val==0]);sort(unique(roots))
}
curve_rows<-list();crossings<-list();sensitivity<-list()
for(ns in c(500,1000)){
 a<-diag[diag$n_S==ns,];fit<-fit_loess(a,'m_TL')
 roots<-find_roots(fit,min(a$D_hat),max(a$D_hat));r<-if(length(roots)==1)roots[1] else NA_real_
 below<-if(is.finite(r))sum(a$D_hat>=r-.025 & a$D_hat<r) else 0L
 above<-if(is.finite(r))sum(a$D_hat>r & a$D_hat<=r+.025) else 0L
 reliable<-length(roots)==1 && below>=10 && above>=10 && r>min(a$D_hat) && r<max(a$D_hat)
 if(!reliable)r<-NA_real_
 root_residual<-if(is.finite(r))as.numeric(predict(fit,newdata=data.frame(D_hat=r))) else NA_real_
 coarse<-find_roots(fit,min(a$D_hat),max(a$D_hat),n=501)
 crossings[[length(crossings)+1]]<-data.frame(n_S=ns,d_0=r,observed_Dhat_min=min(a$D_hat),observed_Dhat_max=max(a$D_hat),n_replicates=nrow(a),span=.8,degree=2,n_crossings=length(roots),interior_supported_crossing=reliable,n_below_within_025=below,n_above_within_025=above,fitted_margin_at_d0=root_residual,coarse_grid_root=if(length(coarse)==1)coarse else NA_real_)
 for(span in c(.7,.8,.9)){
  alt<-fit_loess(a,'m_TL',span);rr<-find_roots(alt,min(a$D_hat),max(a$D_hat))
  sensitivity[[length(sensitivity)+1]]<-data.frame(n_S=ns,span=span,n_crossings=length(rr),d_0=if(length(rr)==1)rr else NA_real_)
 }
 grid<-sort(unique(c(seq(min(a$D_hat),max(a$D_hat),length.out=801),r[is.finite(r)])))
 curve_rows[[length(curve_rows)+1]]<-data.frame(n_S=ns,metric='m_TL',method='TL-SCP',D_hat=grid,y=as.numeric(predict(fit,newdata=data.frame(D_hat=grid))))
 for(metric in c('relErr','NMI'))for(method in levels(d$method)){
  z<-d[d$n_S==ns & d$method==method,];f<-fit_loess(z,metric);xx<-seq(min(z$D_hat),max(z$D_hat),length.out=801)
  curve_rows[[length(curve_rows)+1]]<-data.frame(n_S=ns,metric=metric,method=method,D_hat=xx,y=as.numeric(predict(f,newdata=data.frame(D_hat=xx))))
 }
}
curves<-do.call(rbind,curve_rows);crossings<-do.call(rbind,crossings)
write.csv(crossings,file.path(folder,'Supplement-Figure-A7-crossovers.csv'),row.names=FALSE);write.csv(do.call(rbind,sensitivity),file.path(folder,'Supplement-Figure-A7-span-sensitivity.csv'),row.names=FALSE);write.csv(curves,file.path(folder,'Supplement-Figure-A7-loess-curves.csv'),row.names=FALSE)
print(crossings);print(do.call(rbind,sensitivity))
xlimits<-range(diag$D_hat);error_top<-max(d$relErr)*1.025;margin_limits<-range(c(diag$m_TL,curves$y[curves$metric=='m_TL']))+c(-3,3)
panel<-function(ns,metric){
 is_margin<-metric=='m_TL';a<-if(is_margin)transform(diag[diag$n_S==ns,],method='TL-SCP') else d[d$n_S==ns,]
 z<-curves[curves$n_S==ns & curves$metric==metric,];root<-crossings$d_0[crossings$n_S==ns]
 title<-if(ns==500)switch(metric,relErr='(A) relErr',NMI='(B) NMI',m_TL='(C) Ranking margin') else NULL
 p<-ggplot(a,aes(D_hat,.data[[metric]],color=method))
 if(is_margin)p<-p+geom_hline(yintercept=0,linetype='dashed',color='gray40',linewidth=.4)
 if(is.finite(root))p<-p+geom_vline(xintercept=root,linetype='dotted',color='gray40',linewidth=.4)
 p<-p+geom_point(size=.75,alpha=.2,stroke=0,show.legend=!is_margin)+geom_line(data=z,aes(D_hat,y,color=method),inherit.aes=FALSE,linewidth=.85,show.legend=!is_margin)+
  scale_color_manual(values=colors,breaks=names(colors),labels=c('Cox $K$-means','SCP','TL-SCP'),name=NULL,guide=if(is_margin)'none' else 'legend')+
  scale_x_continuous(breaks=c(.2,.3,.4,.5))+
  labs(x=if(ns==1000)'Empirical discordance $\\widehat D$' else NULL,y=if(is_margin)'$m_{\\rm TL}$' else metric,title=title,subtitle=paste0('$n_S=',ns,'$'))+
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(hjust=.5,size=13),plot.subtitle=element_text(hjust=.5,size=12),axis.title=element_text(size=12),axis.text=element_text(size=12),legend.text=element_text(size=12),legend.position='top',legend.key.width=grid::unit(1.4,'lines'),plot.margin=margin(5,5,5,5))
 if(is_margin){
  p<-p+coord_cartesian(xlim=xlimits,ylim=margin_limits)+scale_y_continuous(breaks=c(-40,-20,0,20,40))
  if(is.finite(root))p<-p+annotate('text',x=root+.008,y=27,label=sprintf('$d_0=%.3f$',root),hjust=0,size=4.2,color='gray20')
 }else if(metric=='NMI')p<-p+scale_y_continuous(breaks=seq(0,1,.25))+coord_cartesian(xlim=xlimits,ylim=c(0,1)) else p<-p+coord_cartesian(xlim=xlimits,ylim=c(0,error_top))
 p
}
combined<-(panel(500,'relErr')+panel(500,'NMI')+panel(500,'m_TL'))/(panel(1000,'relErr')+panel(1000,'NMI')+panel(1000,'m_TL'))+plot_layout(guides='collect')&theme(legend.position='top')
save_figure(combined,file.path(folder,'Supplement-Figure-A7.pdf'),11.8,6.6)
