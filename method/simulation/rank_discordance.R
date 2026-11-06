# Population discordance excludes ties in both numerator and denominator.
population_discordance <- function(bt, bs) {
  pairs <- which(upper.tri(outer(bt, bt, "-")) & outer(bt, bt, "-") != 0, arr.ind = TRUE)
  reversed <- (bt[pairs[, 1]] - bt[pairs[, 2]]) * (bs[pairs[, 1]] - bs[pairs[, 2]]) <
    0
  list(rate = mean(reversed), count = sum(reversed), denominator = nrow(pairs))
}
# Theory Section B: counts of lower-ranked predictors OUTSIDE one's target group.
population_margin <- function(bt, bs, nt, ns) {
  outside <- outer(bt, bt, "!=")
  BT <- rowSums(outer(bt, bt, ">") & outside)
  BS <- rowSums(outer(bs, bs, ">") & outside)
  target_pairs <- which(outer(bt, bt, "<"), arr.ind = TRUE)
  wt <- nt / (nt + ns)
  ws <- ns / (nt + ns)
  margins <- wt * (BT[target_pairs[, 2]] - BT[target_pairs[, 1]]) + ws * (BS[target_pairs[
    ,
    2
  ]] - BS[target_pairs[, 1]])
  worst <- which.min(margins)
  list(m_TL = min(margins), B_T = BT, B_S = BS, w_T = wt, w_S = ws, worst_j = target_pairs[
    worst,
    1
  ], worst_k = target_pairs[worst, 2])
}
empirical_discordance <- function(bt, bs) {
  delta <- outer(bt, bt, "-") * outer(bs, bs, "-")
  mean(delta[upper.tri(delta)] < 0)
}
# Independent Gaussian displacement of an aligned ranking. The exact cross-group
# crossing times give the amplitude needed for the requested inversion count.
# Within-target-group tie breaking is randomized independently in every draw.
calibrate_independent_ranking <- function(bt, D, seed) {
  set.seed(seed)
  p <- length(bt)
  pairs <- which(outer(bt, bt, "<"), arr.ind = TRUE)
  count <- round(D * nrow(pairs))
  stopifnot(abs(count / nrow(pairs) - D) < 1e-12)
  ordered <- order(bt, runif(p))
  base <- integer(p)
  base[ordered] <- seq_len(p)
  attempts <- 0L
  repeat {
    attempts <- attempts + 1L
    offset <- rnorm(p)
    den <- offset[pairs[, 1]] - offset[pairs[, 2]]
    eligible <- den > 0
    times <- sort((base[pairs[eligible, 2]] - base[pairs[eligible, 1]]) / den[eligible])
    if (length(times) > count) {
      break
    }
  }
  amplitude <- if (count == 0) {
    0
  } else {
    mean(times[c(count, count + 1)])
  }
  score <- base + amplitude * offset
  ranks <- rank(score, ties.method = "first")
  disc <- population_discordance(bt, log(ranks))
  stopifnot(disc$count == count)
  list(
    ranks = as.integer(ranks), base = base, offset = offset, amplitude = amplitude,
    attempts = attempts, seed = seed, discordance = disc
  )
}
