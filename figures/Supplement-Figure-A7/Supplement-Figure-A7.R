# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/simulation/helpers.R"))
source(file.path(root, "method/methodology/scp_tlscp.R"))
source(file.path(root, "method/simulation/rank_discordance.R"))
use_fast_distance(root)
repetitions <- as.integer(Sys.getenv("TLSCP_REPS", "50"))
controls <- readRDS(file.path(folder, "Supplement-Figure-A7-target-only.rds"))
stopifnot(repetitions <= max(controls$metrics$rep))
source_sizes <- c(500, 1000)
# TLSCP_LIMIT is useful for checking the installation before the full experiment.
limit <- as.integer(Sys.getenv("TLSCP_LIMIT", "0"))
rows <- list()
diagnostics <- list()
iterations <- list()
counter <- 0L

discordance <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4)
original_rank <- wRankAveg
preliminary <- NULL
wRankAveg <- function(v1, v2, w1, w2) {
  preliminary <<- list(source = as.numeric(v1), target = as.numeric(v2))
  original_rank(v1, v2, w1, w2)
}
for (ns in source_sizes) {
  set.seed(2024)
  beta <- sample(c(rep(2, 30), rep(0.5, 30), rep(-1, 30)))
  X <- matrix(runif(200 * 90, -0.5, 0.5), 200, 90)
  XS <- matrix(runif(ns * 90, -0.5, 0.5), ns, 90) + matrix(
    rnorm(ns * 90, 0, 0.25),
    ns, 90
  )
  for (D in discordance) {
    for (j in seq_len(repetitions)) {
      if (limit > 0 && counter >= limit) {
        next
      }
      seed <- 104729L + ((match(ns, source_sizes) - 1L) * length(discordance) + match(
        D,
        discordance
      )) * 1000L + j
      design <- calibrate_independent_ranking(beta, D, seed)
      beta_source <- log(design$ranks)
      margin <- population_margin(beta, beta_source, 200, ns)$m_TL
      set.seed(2024 + j)
      Y <- simulate_outcomes(X, beta)
      YS <- simulate_outcomes(XS, beta_source)
      fit <- trans_SCP(Y, X, YS, XS,
        penalty = "distance", nclusters = 1:5, maxiters = 1000,
        tol = 1e-04
      )
      k <- which.min(fit$IC_path)
      b <- fit$Beta_path[, k]
      result <- data.frame(
        rep = j, method = "TL-SCP", K = k, relErr = sum((b - beta)^2) / sum(beta^2),
        NMI = nmi(cluster_label(beta), cluster_label_dist(b, k))
      )
      comparison <- rbind(controls$metrics[controls$metrics$rep == j, ], result)
      comparison$n_S <- ns
      comparison$D_target <- D
      rows[[length(rows) + 1]] <- comparison
      diagnostics[[length(diagnostics) + 1]] <- data.frame(
        n_S = ns, D_target = D,
        rep = j, D_realized = design$discordance$rate, D_hat = empirical_discordance(
          preliminary$target,
          preliminary$source
        ), m_TL = margin
      )
      iterations[[length(iterations) + 1]] <- data.frame(
        n_S = ns, D_target = D,
        rep = j, K = 1:5, iterations = fit$Iter_path, selected = 1:5 == k
      )
      counter <- counter + 1L
      saveRDS(
        list(metrics = do.call(rbind, rows), diagnostics = do.call(rbind, diagnostics)),
        file.path(folder, "Supplement-Figure-A7.rds")
      )
      write.csv(do.call(rbind, iterations), file.path(folder, "Supplement-Figure-A7-iterations.csv"),
        row.names = FALSE
      )
      message("Source n=", ns, ", discordance=", D, ", replicate ", j)
    }
  }
}
if (limit == 0 && repetitions == 50) source(file.path(folder, "Supplement-Figure-A7-from-rds.R"))
