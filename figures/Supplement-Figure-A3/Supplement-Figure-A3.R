# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/simulation/helpers.R'))
source(file.path(root,'method/methodology/scp_penalties.R'))
if(!requireNamespace('gurobi',quietly=TRUE)) stop('The penalty comparison requires the Gurobi R interface and a license.')
library(gurobi)
save_results <- function() {
  methods <- c('SCP-Lasso'='Lasso','SCP-SCAD'='SCAD','SCP-MCP'='MCP','SCP-Distance'='Dist')
  settings <- cs_rate
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
  saveRDS(d,file.path(folder,'Supplement-Figure-A3.rds'))
}



n <- 300
p <- 30
sig <- 0.5

cs_rate <- c(0.1, 0.3, 0.5, 0.7)
nRate <- length(cs_rate)

nReps <- as.integer(Sys.getenv('TLSCP_REPS','50'))

lambdaSeq <- 10^seq(0, -4, length.out=30)
lambdaSeq2 <- 10^seq(0, -4, length.out=30)
kSeq <- 1:5
rho <- 1
maxiters <- 1e3
tol <- 1e-5

Err_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
Err_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
Err_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
Err_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
Err_cox <- matrix(NA, nrow = nRate, ncol = nReps)


NMI_M_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
NMI_M_cox <- matrix(NA, nrow = nRate, ncol = nReps)

Beta_Lasso <- list()
Beta_SCAD <- list()
Beta_MCP <- list()
Beta_cox <- list()
Beta_Dist <- list()

k_Dist <- matrix(NA, nrow = nRate, ncol = nReps)

time_Lasso <- matrix(NA, nrow = nRate, ncol = nReps)
time_SCAD <- matrix(NA, nrow = nRate, ncol = nReps)
time_MCP <- matrix(NA, nrow = nRate, ncol = nReps)
time_Dist <- matrix(NA, nrow = nRate, ncol = nReps)
time_cox <- matrix(NA, nrow = nRate, ncol = nReps)

for (i in 1:nRate) {
  
  set.seed(2025)
  beta <- c(rep(2, p/3), rep(0, p/3), rep(-1, p/3))
  beta <- permute(beta)
  cluster.true <- cluster_label(beta)
  X <- matrix(runif(n*p, -1/2, 1/2), nrow = n, ncol = p)
  
  censor_rate <- cs_rate[i]
  
  beta_Lasso <- matrix(NA, nrow = p, ncol = nReps)
  beta_SCAD <- matrix(NA, nrow = p, ncol = nReps)
  beta_MCP <- matrix(NA, nrow = p, ncol = nReps)
  beta_cox <- matrix(NA, nrow = p, ncol = nReps)  
  beta_Dist <- matrix(NA, nrow = p, ncol = nReps)

  for (j in 1:nReps) {
    
    cat(paste('Current stage:', i, 'case', j, 'reps '))
    
    set.seed(2025+j)
    eta <- as.vector(X%*%beta+sig*rnorm(n)) 
    y <- rep(NA, n)
    for (k in 1:n) {
      y[k] <- rexp(1, rate=exp(eta[k]))
    }
    status <- rep(1, n)
    status[sample(1:n, size=censor_rate*n)] <- 0
    
    Y <- cbind(time=y, status=status)

    t0 <- proc.time()
    fit_bcards_lasso <- SCP(y, status, X, penalty="Lasso", lambdaSeq = lambdaSeq,
                                   ebic.gm=0, maxiters=maxiters, tol=tol)
    t1 <-proc.time()-t0
    time_Lasso[i, j] <- t1[[3]]
    idx1 <- which.min(fit_bcards_lasso$IC_path)
    beta_Lasso[, j] <- fit_bcards_lasso$Beta_path[, idx1]
    Err_Lasso[i, j] <- norm(beta_Lasso[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- silhouette_path(as.matrix(beta_Lasso[, j]), kmeans, method = "silhouette")
    cluster.Lasso <- cluster_label_kmeans(beta_Lasso[, j], which.max(k_select$data$y))
    NMI_M_Lasso[i, j] <- nmi(cluster.true, cluster.Lasso)
    cat('Lasso Done ')

    t0 <- proc.time()
    fit_bcards_scad <- SCP(y, status, X, penalty="SCAD", lambdaSeq = lambdaSeq,
                                  ebic.gm=0, maxiters=maxiters, tol=tol)
    t3 <-proc.time()-t0
    time_SCAD[i, j] <- t3[[3]]
    idx3 <- which.min(fit_bcards_scad$IC_path)
    beta_SCAD[, j] <- fit_bcards_scad$Beta_path[, idx3]
    Err_SCAD[i, j] <- norm(beta_SCAD[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- silhouette_path(as.matrix(beta_SCAD[, j]), kmeans, method = "silhouette")
    cluster.SCAD <- cluster_label_kmeans(beta_SCAD[, j], which.max(k_select$data$y))
    NMI_M_SCAD[i, j] <- nmi(cluster.true, cluster.SCAD)
    
    cat('SCAD Done ')

    t0 <- proc.time()
    fit_bcards_mcp <- SCP(y, status, X, penalty="MCP", lambdaSeq = lambdaSeq,
                                 ebic.gm=0, maxiters=maxiters, tol=tol)
    t4 <-proc.time()-t0
    time_MCP[i, j] <- t4[[3]]
    idx4 <- which.min(fit_bcards_mcp$IC_path)
    beta_MCP[, j] <- fit_bcards_mcp$Beta_path[, idx4]
    Err_MCP[i, j] <- norm(beta_MCP[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- silhouette_path(as.matrix(beta_MCP[, j]), kmeans, method = "silhouette")
    cluster.MCP <- cluster_label_kmeans(beta_MCP[, j], which.max(k_select$data$y))
    NMI_M_MCP[i, j] <- nmi(cluster.true, cluster.MCP)
    
    cat('MCP Done ')

    t0 <- proc.time()
    fit_cox <- cox_initial(X, Y, family = "cox", lambda = 0)
    beta0 <- as.vector(fit_cox$beta) 
    
    if(norm(beta0, "2")==0){
      fit_cox <- cox_initial(X, Y, family = "cox", lambda = 1e-2)
      beta0 <- as.vector(fit_cox$beta)
    }
    t5 <-proc.time()-t0
    time_cox[i, j] <- t5[[3]]
    beta_cox[, j] <- beta0
    Err_cox[i, j] <- norm(beta_cox[, j] - beta, "2")^2/norm(beta, "2")^2
    k_select <- silhouette_path(as.matrix(beta_cox[, j]), kmeans, method = "silhouette")
    cluster.cox <- cluster_label_kmeans(beta0, which.max(k_select$data$y))
    NMI_M_cox[i, j] <- nmi(cluster.true, cluster.cox)
    
    cat('Cox Done ')

    t0 <- proc.time()
    fit_bcards_dist <- SCP(y, status, X, penalty="distance", nclusters=kSeq,
                                  ebic.gm=0, maxiters=maxiters, tol=tol)
    
    t6 <-proc.time()-t0
    time_Dist[i, j] <- t6[[3]]
    idx6 <- which.min(fit_bcards_dist$IC_path)
    beta_Dist[, j] <- fit_bcards_dist$Beta_path[, idx6]
    Err_Dist[i, j] <- norm(beta_Dist[, j] - beta, "2")^2/norm(beta, "2")^2 
    k_Dist[i, j] <- kSeq[idx6]
    cluster.Dist <- cluster_label_dist(beta_Dist[, j], k_Dist[i, j])
    NMI_M_Dist[i, j] <- nmi(cluster.true, cluster.Dist)
    
    cat('Dist Done \n')
    save_results()
  }
  
  cat('case', i, 'finish \n')
  
  Beta_Lasso[[i]] <- beta_Lasso
  Beta_SCAD[[i]] <- beta_SCAD
  Beta_MCP[[i]] <- beta_MCP
  Beta_cox[[i]] <- beta_cox
  Beta_Dist[[i]] <- beta_Dist
    save_results()
}
 
save_results()
source(file.path(folder,'Supplement-Figure-A3-from-rds.R'))
