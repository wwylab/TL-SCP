# Retain the Breslow tie convention across glmnet releases.
cox_initial <- function(...) {
  args <- list(...)
  if ("cox.ties" %in% names(formals(glmnet::glmnet))) {
    args$cox.ties <- "breslow"
  }
  do.call(glmnet::glmnet, args)
}

library(Matrix)
library(glmnet)

# Match the maximum-entropy normalization used by aricode::NMI.
nmi <- function(truth, labels) {
  p <- table(truth, labels) / length(truth)
  a <- rowSums(p)
  b <- colSums(p)
  expected <- outer(a, b)
  positive <- p > 0
  information <- sum(p[positive] * log(p[positive] / expected[positive]))
  entropy <- function(q) -sum(q[q > 0] * log(q[q > 0]))
  denominator <- max(entropy(a), entropy(b))
  if (denominator == 0) {
    1
  } else {
    information / denominator
  }
}

# Evaluate the same one-dimensional K-means silhouette path as fviz_nbclust.
silhouette_path <- function(x, FUNcluster = kmeans, method = "silhouette", k.max = 10) {
  # Match the random starts in factoextra 1.0.7 used for the reported
  # simulations.
  set.seed(123)
  d <- as.matrix(dist(x))
  score <- numeric(k.max)
  for (k in 2:k.max) {
    fit <- FUNcluster(x, k)
    groups <- split(seq_len(nrow(d)), fit$cluster)
    widths <- vapply(seq_len(nrow(d)), function(i) {
      own <- groups[[as.character(fit$cluster[i])]]
      if (length(own) == 1L) {
        return(0)
      }
      a <- mean(d[i, setdiff(own, i)])
      others <- groups[names(groups) != as.character(fit$cluster[i])]
      b <- min(vapply(others, function(j) mean(d[i, j]), numeric(1)))
      if (max(a, b) == 0) {
        0
      } else {
        (b - a) / max(a, b)
      }
    }, numeric(1))
    score[k] <- mean(widths)
  }
  list(data = data.frame(clusters = seq_len(k.max), y = score))
}
permute <- function(x) sample(x, length(x), replace = FALSE)
runmed <- function(x, ...) stats::runmed(as.numeric(x), ...)

use_fast_distance <- function(root, r_sort = FALSE) {
  Rcpp::sourceCpp(file.path(root, "method/methodology/fusion_solver.cpp"), cacheDir = file.path(
    tempdir(),
    "tlscp-cpp"
  ), showOutput = FALSE)
  original <- pcox
  assign("pcox", function(
    y, status, X, D, penalty = "Lasso", lambda = 0, k = 2,
    beta0 = rep(0, ncol(X)), maxiters = 1000, tol = 1e-04
  ) {
    if (penalty != "distance") {
      return(original(y, status, X, D, penalty, lambda, k, beta0, maxiters, tol))
    }
    pcox_distance_fast(
      as.numeric(y), as.numeric(status), as.matrix(X), as.matrix(D),
      k, as.numeric(beta0), maxiters, tol, r_sort
    )
  }, envir = .GlobalEnv)
}

simulate_outcomes <- function(X, beta) {
  n <- nrow(X)
  eta <- as.vector(X %*% beta + 0.5 * rnorm(n))
  time <- vapply(seq_len(n), function(i) rexp(1, rate = exp(eta[i])), numeric(1))
  status <- rep(1, n)
  status[sample(seq_len(n), 0.3 * n)] <- 0
  cbind(time = time, status = status)
}

generate_target_controls <- function(nsource, repetitions, permuted_source = FALSE) {
  set.seed(2024)
  beta <- sample(c(rep(2, 30), rep(0.5, 30), rep(-1, 30)))
  X <- matrix(runif(200 * 90, -0.5, 0.5), 200, 90)
  rank <- rank(beta, ties.method = "first")
  if (permuted_source) {
    rank[1:5] <- sample(rank[1:5])
  }
  XS <- matrix(runif(nsource * 90, -0.5, 0.5), nsource, 90) + matrix(rnorm(nsource *
    90, 0, 0.25), nsource, 90)
  truth <- cluster_label(beta)
  rows <- list()
  coefficients <- list()
  for (j in seq_len(repetitions)) {
    set.seed(2024 + j)
    Y <- simulate_outcomes(X, beta)
    # Preserve the original random stream used before K-means initialization.
    unused_source <- simulate_outcomes(XS, log(rank))
    f <- SCP(Y[, 1], Y[, 2], X,
      penalty = "distance", nclusters = 1:5, maxiters = 1000,
      tol = 1e-04
    )
    k <- which.min(f$IC_path)
    b <- f$Beta_path[, k]
    rows[[length(rows) + 1]] <- data.frame(rep = j, method = "SCP", K = k, relErr = sum((b -
      beta)^2) / sum(beta^2), NMI = nmi(truth, cluster_label_dist(b, k)))
    coefficients[[paste0("SCP_", j)]] <- b
    b <- as.numeric(cox_initial(X, Y, family = "cox", lambda = 0)$beta)
    if (sum(b^2) == 0) {
      b <- as.numeric(cox_initial(X, Y, family = "cox", lambda = 0.01)$beta)
    }
    k <- which.max(silhouette_path(as.matrix(b))$data$y)
    cl <- cluster_label_kmeans(b, k)
    rows[[length(rows) + 1]] <- data.frame(
      rep = j, method = "Cox-K-means", K = k,
      relErr = sum((b - beta)^2) / sum(beta^2), NMI = nmi(truth, cl)
    )
    coefficients[[paste0("Cox_", j)]] <- b
    message("Target controls: replicate ", j)
  }
  list(
    metrics = do.call(rbind, rows), coefficients = coefficients, beta_T = beta,
    X_T = X
  )
}
