# SCP penalty comparisons: implementation saved with Supplementary Figures
# A.1--A.4.

cluster_label <- function(beta) {
  threshold_value <- max(abs(beta)) / 100
  ss <- sort(beta, index.return = TRUE)
  sorted_beta <- ss$x
  ind_set <- which(diff(sorted_beta) > threshold_value)
  label_vec <- rep(NA, length(beta))
  ind_set <- c(0, ind_set, length(beta))
  for (i in 1:(length(ind_set) - 1)) {
    sorted_idx_set <- (ind_set[i] + 1):(ind_set[i + 1])
    original_idx_set <- ss$ix[sorted_idx_set]
    label_vec[original_idx_set] <- i
  }
  return(label_vec)
}

cluster_label_dist <- function(beta, k) {
  ss <- sort(beta, index.return = TRUE)
  sorted_beta <- ss$x
  top_diff <- head(sort(diff(sorted_beta), decreasing = TRUE), k - 1)
  top_loc <- sort(match(top_diff, diff(sorted_beta)))
  label_vec <- rep(NA, length(beta))
  ind_set <- c(0, top_loc, length(beta))
  for (i in 1:(length(ind_set) - 1)) {
    sorted_idx_set <- (ind_set[i] + 1):(ind_set[i + 1])
    original_idx_set <- ss$ix[sorted_idx_set]
    label_vec[original_idx_set] <- i
  }
  return(label_vec)
}

cluster_label_kmeans <- function(beta, k) {
  if (length(k) == 1) {
    kmeans_results <- kmeans(beta, k)
    label_vec <- kmeans_results$cluster
  } else {
    numK <- length(k)
    loss <- Inf
    for (j in 1:numK) {
      kmeans_results <- kmeans(beta, k[j])
      if (kmeans_results$totss < loss) {
        label_vec <- kmeans_results$cluster
      }
    }
  }
  return(label_vec)
}

fused_Dist <- function(y, X, w = rep(1, length(y)), D, k, rho = 0.001, beta0 = rep(
                         0,
                         ncol(X)
                       ), maxiters = 100, tol = 1e-05) {
  n <- length(w)
  W <- Diagonal(n = n, x = w)
  XtW <- t(X) %*% W
  XtWy <- XtW %*% y
  XtWX <- XtW %*% X
  beta <- beta0
  for (i in 1:maxiters) {
    beta_last <- beta
    a <- XtWX %*% beta_last
    Dbeta <- as.vector(D %*% beta_last)
    Dbeta_proj <- Proj_sparse(Dbeta, k)
    gradient <- as.vector(a - XtWy + rho * t(D) %*% (Dbeta - Dbeta_proj))
    Av <- XtWX %*% gradient
    vAv <- t(gradient) %*% Av
    Dv <- as.vector(D %*% gradient)
    stepsize <- as.numeric(norm(gradient, "2")^2 / (vAv + rho * norm(Dv, "2")^2))
    beta <- beta_last - stepsize * gradient
    if (norm(as.matrix(beta_last - beta), "f") < tol * (1 + norm(
      as.matrix(beta_last),
      "f"
    ))) {
      break
    }
    rho <- 2 * rho
  }
  return(list(beta = beta))
}

Proj_sparse <- function(x, k) {
  res <- sort(abs(x), method = "quick", index.return = TRUE)
  ind <- tail(res$ix, k)
  xnew <- rep(0, length(x))
  xnew[ind] <- x[ind]
  return(xnew)
}

fused_Lasso <- function(y, X, w = rep(1, length(y)), D, lambda) {
  n <- length(y)
  p <- ncol(X)
  r <- nrow(D)
  W <- Diagonal(n = n, x = w)
  XtW <- t(X) %*% W
  XtWy <- XtW %*% y
  XtWX <- XtW %*% X
  P <- Matrix(0, nrow = p + 2 * r, ncol = p + 2 * r)
  P[1:p, 1:p] <- 0.5 * XtWX
  q <- c(as.vector(-XtWy), rep(lambda, 2 * r))
  A <- cbind(D, Diagonal(n = r, x = -1), Diagonal(n = r, x = 1))
  b <- rep(0, r)
  model <- list()
  model$A <- A
  model$Q <- P
  model$obj <- as.vector(q)
  model$modelsense <- "min"
  model$rhs <- b
  model$sense <- rep("=", r)
  model$lb <- c(rep(-Inf, p), rep(0, 2 * r))
  params <- list(OutputFlag = 0)
  result <- gurobi(model, params)
  beta <- result$x[1:p]
  v <- result$x[(p + 1):(p + r)] - result$x[(p + r + 1):(p + 2 * r)]
  return(list(beta = beta, v = v))
}

fused_SCAD <- function(y, X, w = rep(1, length(y)), D, lambda, a = 3.7, beta0 = rep(
                         0,
                         ncol(X)
                       ), maxiters = 100, tol = 1e-05) {
  for (k in 1:maxiters) {
    Dbeta <- D %*% beta0
    L <- Diagonal(n = nrow(D), x = scad_prime(Dbeta, lambda = lambda, a = a))
    D1 <- L %*% D
    beta <- fused_Lasso(y, X, w, D = D1, lambda = 1)$beta
    if (norm(beta - beta0, "2") < tol * (1 + norm(beta0, "2"))) {
      break
    }
    beta0 <- beta
  }
  return(list(beta = beta, iters = k))
}

fused_MCP <- function(
  y, X, w, D, lambda, a = 3, beta0 = rep(0, ncol(X)), maxiters = 1000,
  tol = 1e-05
) {
  cur_beta <- beta0
  for (i in 1:maxiters) {
    Dbeta <- D %*% cur_beta
    prime_mat <- diag(as.vector(MCP_prime(Dbeta, lambda, a)))
    new_D <- prime_mat %*% D
    next_beta <- fused_Lasso(y, X, w, D = new_D, lambda = 1)$beta
    A <- norm(next_beta - cur_beta, "2") / (1 + norm(cur_beta, "2"))
    if (A < tol) {
      break
    }
    cur_beta <- next_beta
  }
  return(list(beta = cur_beta, iters = i))
}

scad_prime <- function(beta, lambda, a = 3.7) {
  prime <- rep(0, length(beta))
  for (j in 1:length(beta)) {
    prime[j] <- ifelse(abs(beta[j]) <= lambda, lambda, max(
      a * lambda - abs(beta[j]),
      0
    ) / (a - 1))
  }
  return(prime)
}

MCP_prime <- function(beta, lambda, a) {
  size <- length(beta)
  prime_result <- matrix(0, nrow = 1, ncol = size)
  thresh <- lambda * a
  for (i in 1:size) {
    prime_result[1, i] <- ifelse(abs(beta[i, 1]) <= thresh, (lambda - abs(beta[
      i,
      1
    ]) / a) * sign(beta[i, 1]), 0)
  }
  return(prime_result)
}

myGetDkn <- function(k, n) {
  D1n <- function(n) {
    D1n <- matrix(0, nrow = (n - 1), ncol = n)
    for (i in 1:(n - 1)) {
      D1n[i, i] <- -1
      D1n[i, i + 1] <- 1
    }
    return(D1n)
  }
  if (k == 1) {
    return(D1n(n))
  } else {
    Dkn <- D1n(n)
    for (i in 1:(k - 1)) {
      m <- n - i
      Dkn <- D1n(m) %*% Dkn
    }
    A <- Matrix(Dkn, sparse = TRUE)
    return(A)
  }
}

pcox <- function(y, status, X, D, penalty = "Lasso", lambda = 0, k = 2, beta0 = rep(
                   0,
                   ncol(X)
                 ), maxiters = 1000, tol = 1e-04) {
  n <- length(y)
  eta0 <- X %*% beta0
  for (i in 1:maxiters) {
    wz <- update_wz(eta0, y, status)
    w <- wz$w
    z <- wz$z
    if (penalty == "Lasso") {
      beta <- fused_Lasso(y = z, X = X, w = w / n, D = D, lambda = lambda)$beta
    }
    if (penalty == "SCAD") {
      beta <- fused_SCAD(
        y = z, X = X, w = w / n, D = D, lambda = lambda, beta0 = beta0,
        maxiters = maxiters, tol = tol
      )$beta
    }
    if (penalty == "MCP") {
      beta <- fused_MCP(
        y = z, X = X, w = w / n, D = D, lambda = lambda, beta0 = beta0,
        maxiters = maxiters, tol = tol
      )$beta
    }
    if (penalty == "distance") {
      beta <- fused_Dist(
        y = z, X = X, w = w / n, D = D, k = k, beta0 = beta0,
        maxiters = maxiters, tol = tol
      )$beta
    }
    A <- norm(beta - beta0, "2") / (1 + norm(beta0, "2"))
    if (A < tol) {
      break
    } else {
      beta0 <- beta
      eta0 <- X %*% beta0
    }
  }
  return(list(beta = beta, iters = i))
}

pcox_path <- function(y, status, X, D, penalty = "Lasso", lambdaSeq = seq(1, 0.001,
                        length.out = 10
                      ), kSeq = 1:3, ebic.gm = 0, maxiters = 1000, tol = 1e-04) {
  n <- nrow(X)
  p <- ncol(X)
  if (penalty == "Lasso" || penalty == "SCAD" || penalty == "MCP") {
    nlam <- length(lambdaSeq)
    Beta <- matrix(NA, nrow = p, ncol = nlam)
    Iters <- rep(NA, nlam)
    IC <- rep(NA, nlam)
    Y <- cbind(time = y, status = status)
    fit_cox <- cox_initial(X, Y, family = "cox", alpha = 0, lambda = 1e-04)
    beta0 <- runmed(fit_cox$beta, k = 5)
    for (i in 1:nlam) {
      lambda <- lambdaSeq[i]
      fit <- pcox(y, status, X, D,
        penalty = penalty, lambda = lambda, beta0 = beta0,
        maxiters = maxiters, tol = tol
      )
      beta0 <- Beta[, i] <- fit$beta
      Iters[i] <- fit$iters
      nu <- length(which(abs(D %*% beta0) > 0.001))
      tau <- choose(nrow(D), nu)
      IC[i] <- -2 * pll_cox(y, status, X, beta0) + log(n) * nu + 2 * ebic.gm *
        log(tau)
    }
  }
  if (penalty == "distance") {
    nk <- length(kSeq)
    Beta <- matrix(NA, nrow = p, ncol = nk)
    Iters <- rep(NA, nk)
    IC <- rep(NA, nk)
    Y <- cbind(time = y, status = status)
    fit_cox <- cox_initial(X, Y, family = "cox", alpha = 0, lambda = 1e-04)
    beta0 <- runmed(fit_cox$beta, k = 5)
    for (i in 1:nk) {
      fit <- pcox(y, status, X, D,
        penalty = "distance", k = kSeq[i], beta0 = beta0,
        maxiters = maxiters, tol = tol
      )
      Beta[, i] <- as.vector(fit$beta)
      Iters[i] <- fit$iters
      nu <- length(which(abs(D %*% fit$beta) > 0.001))
      tau <- choose(nrow(D), nu)
      IC[i] <- -2 * pll_cox(y, status, X, fit$beta) + log(n) * nu + 2 * ebic.gm *
        log(tau)
    }
  }
  return(list(Beta = Beta, Iters = Iters, IC = IC))
}

pll_cox <- function(y, status, X, beta, weight = rep(1, length(y))) {
  fail_times <- unique(y[status == 1])
  m <- length(fail_times)
  likelihood_vec <- rep(NA, m)
  for (j in 1:m) {
    D <- which(y == fail_times[j] & status == 1)
    R <- which(y >= fail_times[j])
    Xb1 <- X[D, ] %*% beta
    part1 <- sum(diag(weight[D], nrow = length(D)) %*% Xb1)
    Xb2 <- X[R, ] %*% beta
    part2 <- sum(weight[D]) * log_sum_exp_weighted(x = Xb2, w = weight[R])
    likelihood_vec[j] <- part1 - part2
  }
  return(sum(likelihood_vec))
}

log_sum_exp_weighted <- function(x, w) {
  a <- max(x)
  sum_exp <- sum(w * exp(x - a))
  result <- a + log(sum_exp)
  return(result)
}

update_wz <- function(eta, y, status) {
  n <- length(y)
  w <- rep(0, n)
  z <- rep(NA, n)
  for (k in 1:n) {
    indset <- which(y <= y[k] & status == 1)
    if (length(indset) == 0) {
      w[k] <- 0
      z[k] <- eta[k]
    } else {
      s <- 0
      for (i in indset) {
        R <- which(y >= y[i])
        c <- sum(exp(eta[R]))
        a <- exp(eta[k])
        if (is.infinite(c) || is.na((a * c - a^2) / c^2)) {
          w[k] <- w[k] + 0
        } else {
          w[k] <- w[k] + (a * c - a^2) / c^2
        }
        s <- s + a / c
      }
      z[k] <- eta[k] + (status[k] - s) / w[k]
    }
  }
  return(list(w = w, z = z))
}

SCP <- function(
  y, status, X, penalty = "Lasso", lambdaSeq = seq(1, 0.001, length.out = 10),
  MM = FALSE, nclusters = 2:4, ebic.gm = 0, maxiters = 1000, tol = 1e-04
) {
  n <- nrow(X)
  p <- ncol(X)
  Y <- cbind(time = y, status = status)
  fit0 <- cox_initial(X, Y, family = "cox", alpha = 0, lambda = 1e-04)
  beta0 <- as.vector(fit0$beta)
  ss <- sort(beta0, index.return = TRUE)
  tau <- ss$ix
  X_order <- X[, tau]
  D <- myGetDkn(k = 1, n = ncol(X_order))
  fit <- pcox_path(y, status, X_order, D,
    penalty = penalty, lambdaSeq = lambdaSeq,
    kSeq = nclusters - 1, ebic.gm = ebic.gm, maxiters = maxiters, tol = tol
  )
  beta_path <- fit$Beta
  for (i in 1:ncol(beta_path)) {
    beta_order <- fit$Beta[, i]
    for (j in 1:length(beta_order)) {
      beta_path[tau[j], i] <- beta_order[j]
    }
  }
  return(list(Beta_path = beta_path, Iter_path = fit$Iters, IC_path = fit$IC))
}

SCP_allpairs <- function(
  y, status, X, penalty = "Lasso", lambdaSeq = seq(1, 0.001,
    length.out = 10
  ), nclusters = 2:4, ebic.gm = 0, screening = FALSE, maxiters = 1000,
  tol = 1e-04
) {
  n <- nrow(X)
  p <- ncol(X)
  if (screening) {
    Y <- survival::Surv(y, status)
    sis_fit <- SIS(X, Y, family = "cox", penalty = "lasso", tune = "bic", varISIS = "vanilla")
    nnz_id <- sis_fit$ix
    if (length(nnz_id) > n) {
      nnz_num <- floor(p / (4 * log(p)))
      sis_fit <- SIS(X, Y,
        family = "cox", penalty = "lasso", tune = "bic", varISIS = "vanilla",
        nsis = nnz_num
      )
      nnz_id <- sis_fit$ix
    }
    X0 <- X[, nnz_id]
  } else {
    X0 <- X
  }
  p0 <- ncol(X0)
  D <- allpair_D(p0)
  fit <- pcox_path(y, status, X0, D, penalty = penalty, lambdaSeq = lambdaSeq, kSeq = nclusters -
    1, ebic.gm = ebic.gm, maxiters = maxiters, tol = tol)
  if (screening) {
    Beta_path <- matrix(0, nrow = p, ncol = ncol(fit$Beta))
    for (i in 1:ncol(fit$Beta)) {
      Beta_path[nnz_id, i] <- fit$Beta[, i]
    }
  } else {
    Beta_path <- fit$Beta
  }
  return(list(Beta_path = Beta_path, Iter_path = fit$Iters, IC_path = fit$IC))
}

allpair_D <- function(n) {
  A <- combn(n, 2)
  m <- ncol(A)
  D <- Matrix(0, nrow = m, ncol = n, sparse = TRUE)
  for (i in 1:m) {
    D[i, A[, i]] <- c(-1, 1)
  }
  return(D)
}
