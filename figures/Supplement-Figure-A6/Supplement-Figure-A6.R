# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/simulation/helpers.R'))
source(file.path(root,'method/methodology/scp_tlscp.R'))
source(file.path(root,'method/simulation/rank_discordance.R'))
use_fast_distance(root)
repetitions <- as.integer(Sys.getenv('TLSCP_REPS','50'))
controls <- readRDS(file.path(folder,'Supplement-Figure-A6-target-only.rds'))
stopifnot(repetitions<=max(controls$metrics$rep))
source_sizes <- c(500,1000)
# TLSCP_LIMIT is useful for checking the installation before the full experiment.
limit <- as.integer(Sys.getenv('TLSCP_LIMIT','0'))
rows <- list(); diagnostics <- list(); iterations <- list(); counter <- 0L

shifts <- c(.05,.10,.25,.35,.50)
for(ns in source_sizes) for(shift in shifts) {
  set.seed(2024)
  beta <- sample(c(rep(2,30),rep(.5,30),rep(-1,30)))
  X <- matrix(runif(200*90,-.5,.5),200,90)
  beta_source <- log(rank(beta,ties.method='first'))
  XS <- matrix(runif(ns*90,-.5,.5),ns,90)+matrix(rnorm(ns*90,0,shift),ns,90)
  for(j in seq_len(repetitions)) {
    if(limit>0 && counter>=limit)next
    set.seed(2024+j);Y <- simulate_outcomes(X,beta);YS <- simulate_outcomes(XS,beta_source)
    f <- trans_SCP(Y,X,YS,XS,penalty='distance',nclusters=1:5,maxiters=1000,tol=1e-4)
    k <- which.min(f$IC_path);b <- f$Beta_path[,k]
    result <- data.frame(rep=j,method='TL-SCP',K=k,relErr=sum((b-beta)^2)/sum(beta^2),
                         NMI=nmi(cluster_label(beta),cluster_label_dist(b,k)))
    comparison <- rbind(controls$metrics[controls$metrics$rep==j,],result)
    comparison$n_S <- ns;comparison$setting <- shift
    rows[[length(rows)+1]] <- comparison
    counter <- counter+1L
    saveRDS(do.call(rbind,rows),file.path(folder,'Supplement-Figure-A6.rds'))
    message('Source n=',ns,', covariate shift=',shift,', replicate ',j)
  }
}
if(limit==0 && repetitions==50)source(file.path(folder,'Supplement-Figure-A6-from-rds.R'))
