# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/simulation/helpers.R'))
source(file.path(root,'method/methodology/scp_tlscp.R'))
use_fast_distance(root)
save_results <- function() {
  methods <- c('SCP'='SCP','TL-SCP'='tSCP','Cox-K-means'='cox')
  settings <- nS_Seq
  rows <- lapply(names(methods),function(method) {
    code <- methods[[method]]
    error <- get(paste0('Err_',code)); accuracy <- get(paste0('NMI_M_',code))
    runtime <- get(paste0('time_',code))
    do.call(rbind,lapply(seq_along(settings),function(s) {
      data.frame(setting=settings[s],rep=seq_len(nReps),method=method,
                 relErr=error[s,],NMI=accuracy[s,],Runtime=runtime[s,])
    }))
  })
  d <- do.call(rbind,rows);d <- d[is.finite(d$relErr)&is.finite(d$NMI),]
  saveRDS(d,file.path(folder,'Supplement-Figure-A5.rds'))
}



n_T <- 200
p <- 90
nS_Seq <- c(500, 1000, 1500, 2000)
nnS <- length(nS_Seq)

sig <- 0.5
censor_rate <- 0.3

nReps <- as.integer(Sys.getenv('TLSCP_REPS','50'))

kSeq <- 1:5
maxiters <- 1e3
tol <- 1e-4

Err_SCP <- matrix(NA, nrow = nnS, ncol = nReps)
Err_tSCP <- matrix(NA, nrow = nnS, ncol = nReps)
Err_cox <- matrix(NA, nrow = nnS, ncol = nReps)


NMI_M_SCP <- matrix(NA, nrow = nnS, ncol = nReps)
NMI_M_tSCP <- matrix(NA, nrow = nnS, ncol = nReps)
NMI_M_cox <- matrix(NA, nrow = nnS, ncol = nReps)



Beta_SCP <- list()
Beta_tSCP <- list()
Beta_cox <- list()

k_SCP <- matrix(NA, nrow = nnS, ncol = nReps)
k_tSCP <- matrix(NA, nrow = nnS, ncol = nReps)
k_cox <- matrix(NA, nrow = nnS, ncol = nReps)

time_SCP <- matrix(NA, nrow = nnS, ncol = nReps)
time_tSCP <- matrix(NA, nrow = nnS, ncol = nReps)
time_cox <- matrix(NA, nrow = nnS, ncol = nReps)

for (i in 1:nnS) {
  
  set.seed(2024)
  
  n_S <- nS_Seq[i]
  n <- n_T+n_S
  beta_T <- c(rep(2, p/3), rep(0.5, p/3), rep(-1, p/3))
  beta_T <- permute(beta_T)
  cluster.true <- cluster_label(beta_T)
  X_T <- matrix(runif(n_T*p, -1/2, 1/2), nrow = n_T, ncol = p)

  beta_S <- log(rank(beta_T, ties.method = "first"))
  X_S <- matrix(runif(n_S*p, -1/2, 1/2), nrow = n_S, ncol = p)
  
  beta_SCP <- matrix(NA, nrow = p, ncol = nReps)  
  beta_tSCP <- matrix(NA, nrow = p, ncol = nReps)
  beta_cox <- matrix(NA, nrow = p, ncol = nReps) 
  
  for (j in 1:nReps) {
    
    cat(paste('Current stage:', i, 'case', j, 'reps '))
    
    set.seed(2024+j)
    eta_T <- as.vector(X_T%*%beta_T+sig*rnorm(n_T)) 
    y_T <- rep(NA, n_T)
    for (k in 1:n_T) {
      y_T[k] <- rexp(1, rate=exp(eta_T[k]))
    }
    status_T <- rep(1, n_T)
    status_T[sample(1:n_T, size=censor_rate*n_T)] <- 0
    Y_T <- cbind(time=y_T, status=status_T)
 
    eta_S <- as.vector(X_S%*%beta_S+sig*rnorm(n_S)) 
    y_S <- rep(NA, n_S)
    for (k in 1:n_S) {
      y_S[k] <- rexp(1, rate=exp(eta_S[k]))
    }
    status_S <- rep(1, n_S)
    status_S[sample(1:n_S, size=censor_rate*n_S)] <- 0
    Y_S <- cbind(time=y_S, status=status_S)

    t0 <- proc.time()
    fit_SCP <- SCP(y=Y_T[, 1], status=Y_T[, 2], X_T, penalty="distance", nclusters=kSeq,
                          ebic.gm=0, maxiters=maxiters, tol=tol)
    
    t1 <-proc.time()-t0
    time_SCP[i, j] <- t1[[3]]
    idx1 <- which.min(fit_SCP$IC_path)
    beta_SCP[, j] <- fit_SCP$Beta_path[, idx1]
    Err_SCP[i, j] <- norm(beta_SCP[, j] - beta_T, "2")^2/norm(beta_T, "2")^2 
    k_SCP[i, j] <- kSeq[idx1]
    cluster.SCP <- cluster_label_dist(beta_SCP[, j], k_SCP[i, j])
    NMI_M_SCP[i, j] <- nmi(cluster.true, cluster.SCP)
    cat('SCP Done \n')

    t0 <- proc.time()
    fit_tSCP <- trans_SCP(Y_T, X_T, Y_S, X_S, penalty="distance", nclusters=kSeq,
                          ebic.gm=0,  maxiters=maxiters, tol=tol)
    
    t2 <-proc.time()-t0
    time_tSCP[i, j] <- t2[[3]]
    idx2 <- which.min(fit_tSCP$IC_path)
    beta_tSCP[, j] <- fit_tSCP$Beta_path[, idx2]
    Err_tSCP[i, j] <- norm(beta_tSCP[, j] - beta_T, "2")^2/norm(beta_T, "2")^2 
    k_tSCP[i, j] <- kSeq[idx2]
    cluster.tSCP <- cluster_label_dist(beta_tSCP[, j], k_tSCP[i, j])
    NMI_M_tSCP[i, j] <- nmi(cluster.true, cluster.tSCP)
    cat('SCP Done \n')

    t0 <- proc.time()
    fit_cox <- cox_initial(X_T, Y_T, family = "cox", lambda = 0)
    beta0 <- as.vector(fit_cox$beta) 
    if(norm(beta0, "2")==0){
      fit_cox <- cox_initial(X_T, Y_T, family = "cox", lambda = 1e-2)
      beta0 <- as.vector(fit_cox$beta)
    }
    t3 <-proc.time()-t0
    
    time_cox[i, j] <- t3[[3]]
    beta_cox[, j] <- beta0
    Err_cox[i, j] <- norm(beta_cox[, j] - beta_T, "2")^2/norm(beta_T, "2")^2
    k_select <- silhouette_path(as.matrix(beta_cox[, j]), kmeans, method = "silhouette")
    k_cox[i, j] <- which.max(k_select$data$y)
    cluster.cox <- cluster_label_kmeans(beta0, which.max(k_select$data$y))
    NMI_M_cox[i, j] <- nmi(cluster.true, cluster.cox)
    
    cat('cox-kmeans Done \n')
    save_results()
  }
  
  cat('case', i, 'finish \n')
  
  Beta_SCP[[i]] <- beta_SCP
  Beta_tSCP[[i]] <- beta_tSCP
  Beta_cox[[i]] <- beta_cox
    save_results()
}
save_results()
source(file.path(folder,'Supplement-Figure-A5-from-rds.R'))
