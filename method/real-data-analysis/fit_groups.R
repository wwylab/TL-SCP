OPT_MAXITERS <- 150L
OPT_TOL <- 1e-04

build_design_matrix <- function(df, mut_levels) {
  sex <- ifelse(df$Sex == "M", 1, 0)
  mut <- factor(df$ProtDescription, levels = mut_levels)
  X_mut <- Matrix::sparse.model.matrix(~ 0 + mut)
  colnames(X_mut) <- sub("^mut", "", colnames(X_mut))
  cbind(gender = sex, X_mut)
}

fit_grouped_method <- function(
  X_T, Y_T, X_order, D, tau, nc_seq, covar_count = 1,
  method_name = "SCP"
) {
  nC <- length(nc_seq)
  llk <- bic <- aic <- nu <- rep(NA_real_, nC)
  beta_mat <- matrix(NA_real_, nrow = ncol(X_T), ncol = nC)
  rownames(beta_mat) <- colnames(X_T)

  init_fit <- cox_initial(X_order, Y_T, family = "cox", alpha = 0, lambda = 1e-04)
  beta_init <- as.numeric(init_fit$beta)
  previous_beta_order <- NULL
  for (j in seq_along(nc_seq)) {
    nclusters <- nc_seq[j]
    cat(sprintf(
      "[%s] fitting candidate with %d clusters (%d/%d)\n", method_name,
      nclusters, j, length(nc_seq)
    ))

    if (nclusters == 1) {
      if (covar_count > 0) {
        covar_df <- as.data.frame(as.matrix(X_order[, seq_len(covar_count),
          drop = FALSE
        ]))
        names(covar_df) <- paste0("x", seq_len(covar_count))
        covar_df$time <- Y_T[, 1]
        covar_df$status <- Y_T[, 2]
        cox_fit <- survival::coxph(survival::Surv(time, status) ~ .,
          data = covar_df,
          ties = "breslow"
        )
        beta_cov <- as.numeric(stats::coef(cox_fit))
        beta_cov[!is.finite(beta_cov)] <- 0
      } else {
        beta_cov <- numeric(0)
      }
      beta_order_hat <- c(beta_cov, rep(0, ncol(X_order) - covar_count))
      beta_hat <- c(beta_cov, rep(0, ncol(X_T) - covar_count))
      beta_mat[, j] <- beta_hat
      previous_beta_order <- beta_order_hat
      llk[j] <- pll_cox(Y_T[, 1], Y_T[, 2], X_T, beta_hat)
      nu[j] <- 0
      bic[j] <- -2 * llk[j]
      aic[j] <- -2 * llk[j]
      next
    }

    starts <- list(beta_init)
    if (!is.null(previous_beta_order)) {
      starts[[length(starts) + 1]] <- previous_beta_order
    }
    best_out <- NULL
    best_llk <- -Inf
    best_beta_hat <- NULL
    for (s in seq_along(starts)) {
      out_s <- pcox(Y_T[, 1], Y_T[, 2], X_order, D,
        penalty = "distance", lambda = 0,
        k = nclusters - 1, rho = 1, beta0 = starts[[s]], maxiters = OPT_MAXITERS,
        tol = OPT_TOL
      )

      beta_order <- out_s$beta[-seq_len(covar_count)]
      beta_mut <- beta_order
      beta_mut[tau] <- beta_order
      beta_hat_s <- c(out_s$beta[seq_len(covar_count)], beta_mut)
      llk_s <- pll_cox(Y_T[, 1], Y_T[, 2], X_T, beta_hat_s)
      if (is.finite(llk_s) && llk_s > best_llk) {
        best_out <- out_s
        best_llk <- llk_s
        best_beta_hat <- beta_hat_s
      }
    }

    if (is.null(best_out)) {
      stop("No finite solution found for ", method_name, " with ", nclusters,
        " clusters.",
        call. = FALSE
      )
    }

    beta_hat <- best_beta_hat
    beta_mat[, j] <- beta_hat
    previous_beta_order <- best_out$beta

    llk[j] <- best_llk
    nu[j] <- length(which(abs(D %*% best_out$beta) > 0.01))
    bic[j] <- -2 * llk[j] + log(nrow(X_T)) * nu[j]
    aic[j] <- -2 * llk[j] + 2 * nu[j]
  }

  idx <- which.min(bic)
  list(method = method_name, idx = idx, nclusters = nc_seq[idx], beta = beta_mat[
    ,
    idx
  ], beta_matrix = beta_mat, bic = bic, aic = aic, llk = llk, nu = nu)
}
