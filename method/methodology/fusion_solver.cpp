#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// Algebraically identical distance-penalized iterations to the frozen simulation.
// [[Rcpp::export]]
Rcpp::List pcox_distance_fast(const Eigen::Map<Eigen::VectorXd> y,
 const Eigen::Map<Eigen::VectorXd> status, const Eigen::Map<Eigen::MatrixXd> X,
 const Eigen::Map<Eigen::MatrixXd> D, int k, Eigen::VectorXd beta0,
 int maxiters=1000, double tol=1e-4) {
 const int n=X.rows(),p=X.cols(),q=D.rows();
 Eigen::VectorXd eta=X*beta0, beta=beta0;
 int outer=0;
 for(outer=1;outer<=maxiters;outer++) {
  Eigen::VectorXd a=eta.array().exp(), w=Eigen::VectorXd::Zero(n), score=status;
  Eigen::VectorXi has=Eigen::VectorXi::Zero(n);
  for(int i=0;i<n;i++) if(status[i]==1) {
   double den=0; for(int l=0;l<n;l++)if(y[l]>=y[i])den+=a[l];
   for(int l=0;l<n;l++)if(y[l]>=y[i]) {
    has[l]=1; double vv=(a[l]*den-a[l]*a[l])/(den*den);
    if(std::isfinite(den)&&std::isfinite(vv))w[l]+=vv;
    score[l]-=a[l]/den;
   }
  }
  Eigen::VectorXd z=eta;
  for(int l=0;l<n;l++)if(has[l])z[l]+=score[l]/w[l];
  Eigen::VectorXd wn=w/n;
  Eigen::MatrixXd XtW=X.transpose()*wn.asDiagonal(), XtWX=XtW*X;
  Eigen::VectorXd XtWy=XtW*z;
  beta=beta0;double rho=1e-3;
  for(int inner=1;inner<=maxiters;inner++) {
   Eigen::VectorXd db=D*beta,proj=Eigen::VectorXd::Zero(q);
   std::vector<int> ix(q);std::iota(ix.begin(),ix.end(),0);
   std::stable_sort(ix.begin(),ix.end(),[&](int i,int j){return std::abs(db[i])<std::abs(db[j]);});
   for(int j=q-k;j<q;j++)if(j>=0)proj[ix[j]]=db[ix[j]];
   Eigen::VectorXd grad=XtWX*beta-XtWy+rho*D.transpose()*(db-proj);
   Eigen::VectorXd dg=D*grad;
   double den=grad.dot(XtWX*grad)+rho*dg.squaredNorm();
   double step=grad.squaredNorm()/den;
   Eigen::VectorXd next=beta-step*grad;
   bool done=(next-beta).norm()<tol*(1+beta.norm());beta=next;
   if(done)break;rho*=2;
  }
  if(!beta.allFinite())Rcpp::stop("Nonfinite coefficient in distance fit");
  double change=(beta-beta0).norm()/(1+beta0.norm());
  if(change<tol)break;beta0=beta;eta=X*beta0;
 }
 return Rcpp::List::create(Rcpp::Named("beta")=beta,Rcpp::Named("iters")=std::min(outer,maxiters));
}
