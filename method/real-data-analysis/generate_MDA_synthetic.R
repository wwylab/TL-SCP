# Generate a synthetic MD Anderson target cohort for testing the analysis scripts.
# Run with Rscript from any directory. No clinical data are read.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))

set.seed(2024)
# Mutation names only, from saved Table F.2; counts and effects are synthetic.
mutations <- c(
  "p.D281N", "p.G245S", "p.M237I", "p.M246V", "p.P177R", "p.R110L", "p.R158H",
  "p.R175H", "p.R181C", "p.R181H", "p.R196*", "p.R248Q", "p.R248W", "p.R267Q",
  "p.R267W", "p.R273H", "p.R280T", "p.R282W", "p.R283C", "p.R290H", "p.R306*",
  "p.R337C", "p.R337H", "p.R342*", "p.S241fs", "p.T125T", "p.W146*", "p.Y220C"
)
# Figure A.5 effects, with nearly equal groups for 28 mutation categories.
beta_target <- sample(rep(c(2, 0.5, -1), times = c(10, 9, 9)))

n <- 202L
mutation <- sample(rep(mutations, length.out = n))
male <- rbinom(n, 1, 0.5)
X <- cbind(gender = male, outer(mutation, mutations, "==") * 1)
colnames(X) <- c("gender", mutations)
stopifnot(all(X %in% c(0, 1)), all(rowSums(X[, -1, drop = FALSE]) == 1))

# Cox event times in years, with a sex coefficient of 0.5 and no frailty term.
baseline_rate <- 0.05
event_rate <- baseline_rate * exp(as.vector(X %*% c(0.5, beta_target)))

# Choose a common censoring rate for 30% expected censoring, before drawing times.
censor_rate <- uniroot(function(rate) {
  mean(rate / (rate + event_rate)) - 0.3
}, interval = c(0, max(event_rate)), tol = 1e-12)$root
event_time <- rexp(n, rate = event_rate)
censor_time <- rexp(n, rate = censor_rate)

# Pair carriers within each mutation; an odd last carrier is alone.
within_mutation <- ave(seq_len(n), mutation, FUN = seq_along)
family <- (within_mutation + 1L) %/% 2L
result <- data.frame(
  Individual_code = sprintf("SYN-T-%04d", seq_len(n)),
  Family_ID = sprintf("SYN-T-M%03d-F%02d", match(mutation, mutations), family),
  Sex = ifelse(male == 1, "M", "F"),
  ProtDescription = mutation,
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time)
)
stopifnot(
  nrow(result) == n, !anyNA(result), !anyDuplicated(result$Individual_code),
  all(is.finite(result$time) & result$time > 0),
  all(result$status %in% c(0L, 1L)),
  length(unique(result$status)) == 2L,
  length(unique(result$ProtDescription)) == 28L,
  min(table(result$ProtDescription)) >= 2L
)
write.csv(result, file.path(folder, "MDA-synthetic.csv"), row.names = FALSE)
message("Wrote MDA-synthetic.csv: ", n, " synthetic carriers and 28 mutation categories.")
