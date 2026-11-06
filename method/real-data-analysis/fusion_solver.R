#' Create cluster labels for a given vector of coefficients in a data-driven manner.
#'
#' @param beta The vector of coefficients.
#' @export
cluster_label <- function(beta) {
  threshold_value <- max(abs(beta)) / 100
  ss <- sort(beta, index.return = TRUE)
  sorted_beta <- ss$x
  # threshold_value<- quantile(diff(sorted_beta), 0.9)
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


#' Create cluster labels for a given vector of coefficients from the distance penalization
#'
#' @param beta The vector of coefficients.
#' @param k The number of clusters
#' @export
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


#' Create cluster labels for a given vector of coefficients in a data-driven manner.
#'
#' @param beta The vector of coefficients.
#' @param k The number of clusters
#' @export
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


#' Compute k-th order fusion matrix.
#'
#' @param k The order of the fusion matrix.
#' @param n The number of time points.
#' @export
myGetDkn <- function(k, n) {
  D1n <- function(n) {
    # a function generates the first order difference matrix with n columns
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


#' Fit the penalized Cox regression model at a fixed value of the penalization parameter.
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param D The fusion matrix in the penalty term.
#' @param penalty The penalty function, options include 'Lasso'(default), 'SCAD', 'MCP', and 'distance'.
#' @param lambda The regularization parameter lambda.
#' @param k The regularization parameter k for the distance penalty (equal to the number of nonzeros in Dbeta).
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm.
#' @export
pcox <- function(
  y, status, X, D, penalty = "Lasso", lambda = 0, k = 2, rho = 0.001,
  beta0 = rep(0, ncol(X)), maxiters = 1000, tol = 1e-04
) {
  n <- length(y)
  # get the initial value of eta
  eta0 <- X %*% beta0

  for (i in 1:maxiters) {
    # compute the Hessian matrix W and z through a function
    wz <- update_wz(eta0, y, status)
    w <- wz$w
    z <- wz$z


    # update beta through a fused penalized LS problem
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
        y = z, X = X, w = w / n, D = D, k = k, rho = rho, beta0 = beta0,
        maxiters = maxiters, tol = tol
      )$beta
    }


    # check convergence
    A <- norm(beta - beta0, "2") / (1 + norm(beta0, "2"))
    if (A < tol) {
      # terminate the loop
      break
    } else {
      # update for the next iteration
      beta0 <- beta
      eta0 <- X %*% beta0
    }
  }

  # return the results
  return(list(beta = beta, iters = i))
}


#' Compute the solution of weighted LS + fusion-distance penalty.
#'
#' @param y The vector of response.
#' @param X The design matrix.
#' @param w The vector of weights in the least squares term.
#' @param D The fusion matrix.
#' @param k The regularization parameter k for the distance penalty (equal to the number of nonzeros in Dbeta).
#' @param rho The initial value of rho, increases at each iteration to reach to a large value.
#' @param beta0 The initial vector of coefficients.
#' @param maxiters The maximum number of iterations.
#' @param tol The tolerance value for early stopping the algorithm.
#' @export
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
#' Compute the partial log-likelihood for Cox regression.
#'
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @param X The design matrix.
#' @param beta The vector of coefficients.
#' @param weight The vector of weights associate with each observation, default is 1.
#' @export
pll_cox <- function(y, status, X, beta, weight = rep(1, length(y))) {
  fail_times <- unique(y[status == 1])
  m <- length(fail_times)
  # Define a vector to record likelihood at each specific failure time
  likelihood_vec <- rep(NA, m)

  for (j in 1:m) {
    D <- which(y == fail_times[j] & status == 1) # consider ties
    R <- which(y >= fail_times[j])
    # Compute the current stage likelihood with the given weights
    Xb1 <- X[D, ] %*% beta
    part1 <- sum(diag(weight[D], nrow = length(D)) %*% Xb1)
    Xb2 <- X[R, ] %*% beta
    part2 <- sum(weight[D]) * log_sum_exp_weighted(x = Xb2, w = weight[R])

    likelihood_vec[j] <- part1 - part2
  }

  # Return the partial log-likelihood
  return(sum(likelihood_vec))
}


#' A function to apply the log-sum-exp trick
#'
#' @param X the parameter we wants to compute the log exp sum
#' @param w The weights of the exp
#' @export
log_sum_exp_weighted <- function(x, w) {
  a <- max(x) # Find the maximum value to scale the computation
  sum_exp <- sum(w * exp(x - a)) # Compute the weighted sum of exponentials in a numerically stable way
  result <- a + log(sum_exp) # Adjust back by adding 'a' and applying log
  return(result)
}


#' Compute w and z in the transformed weighted LS problem for Cox regression
#'
#' @param eta The eta vector in the definition of w and z, eta=Xbeta.
#' @param y The vector of survival time.
#' @param status The vector of status of each object, 1 for failure and 0 for right-censored.
#' @export
update_wz <- function(eta, y, status) {
  # get the number of sample size, which is the dimension of w and z
  n <- length(y)
  w <- rep(0, n)
  z <- rep(NA, n)

  for (k in 1:n) {
    # compute C_k, the set of i with t_i < y_k
    indset <- which(y <= y[k] & status == 1)

    if (length(indset) == 0) {
      w[k] <- 0 # before is 1.
      z[k] <- eta[k] #+ (status[k] - 1)
    } else {
      # compute w and a part in z
      s <- 0
      for (i in indset) {
        # compute R_i, the set of indices j with y_j >= t_i
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

      # compute z
      z[k] <- eta[k] + (status[k] - s) / w[k]
    }
  }

  return(list(w = w, z = z))
}
#' Weighted rank averaging
#'
#' @param v1 The first ranking
#' @param v2 The second ranking
#' @param w1 The weight vector associated with v1
#' @param w2 The weight vector associated with v2
#' @export
wRankAveg <- function(v1, v2, w1, w2) {
  r1 <- rank(v1, ties.method = "first")
  r2 <- rank(v2, ties.method = "first")

  r <- w1 * r1 + w2 * r2
  tau <- sort(r, index.return = TRUE)$ix
}
