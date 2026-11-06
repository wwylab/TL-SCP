# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))

source(file.path(root, "method/real-data-analysis/analysis.R"))
a <- fit_clinical(root)
b <- fit_clinical(root, mda = TRUE)
shared <- merge(a$annotations, b$annotations, by = "ProtDescription", suffixes = c(
  "_Database",
  "_MDA"
))
changed <- shared[shared[["TL-SCP_Database"]] != shared[["TL-SCP_MDA"]], ]
# Refit the family bootstrap from authorized records; retain only standard
# errors.
bootstrap <- function() {
  tmp <- tempfile("family-bootstrap-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  for (tag in c("NCI", "MDA")) {
    d <- if (tag == "NCI") {
      a$target
    } else {
      b$target
    }
    if (!"Family_ID" %in% names(d)) {
      stop("The family bootstrap requires Family_ID in the authorized MDA input.")
    }
    write.csv(d, file.path(tmp, paste0(tag, ".csv")), row.names = FALSE)
  }
  code <- system2(Sys.getenv("PYTHON", "python3"), c(
    shQuote(file.path(root, "method/real-data-analysis/family_bootstrap.py")),
    shQuote(file.path(tmp, "NCI.csv")), shQuote(file.path(tmp, "MDA.csv")), shQuote(file.path(
      tmp,
      "summary.csv"
    ))
  ))
  if (code != 0) {
    stop("Family bootstrap did not complete.")
  }
  read.csv(file.path(tmp, "summary.csv"))
}
se <- bootstrap()
ranks <- do.call(rbind, lapply(c("NCI", "MDA"), function(tag) {
  fit <- if (tag == "NCI") {
    a
  } else {
    b
  }
  d <- fit$ranks
  d$group <- sub("-onset", "", fit$annotations[["TL-SCP"]][match(d$Mutation, fit$annotations$ProtDescription)])
  d$cohort <- tag
  d$mutation <- sub("^p\\.", "", d$Mutation)
  d$highlight <- ifelse(d$Mutation %in% changed$ProtDescription, "Disconcordant",
    "Other"
  )
  d$y <- c(Early = 3, Mid = 2, Late = 1)[d$group]
  for (g in unique(d$group)) {
    ix <- which(d$group == g)
    ix <- ix[order(d$weighted_rank[ix], d$Mutation[ix])]
    runs <- cumsum(c(TRUE, diff(d$weighted_rank[ix]) > 0.45))
    for (run in unique(runs)) {
      ii <- ix[runs == run]
      n <- length(ii)
      if (n > 1) {
        d$y[ii] <- d$y[ii] + seq(if (n == 2) {
          -0.12
        } else {
          -0.22
        }, if (n == 2) {
          0.12
        } else {
          0.22
        }, length.out = n)
      }
    }
  }
  d
}))
forest <- data.frame(Mutation = changed$ProtDescription)
for (tag in c("NCI", "MDA")) {
  z <- ranks[ranks$cohort == tag, ]
  s <- se[se$cohort == tag, ]
  centered <- z$target_beta_raw - z$target_beta_raw[match("p.R175H", z$Mutation)]
  forest[[paste0("beta_", tag)]] <- centered[match(forest$Mutation, z$Mutation)]
  forest[[paste0("SE_", tag)]] <- s$bootstrap_SE[match(forest$Mutation, s$Mutation)]
}
forest$delta <- forest$beta_MDA - forest$beta_NCI
forest$SE_delta <- sqrt(forest$SE_NCI^2 + forest$SE_MDA^2)
forest$lower <- forest$delta - 1.96 * forest$SE_delta
forest$upper <- forest$delta + 1.96 * forest$SE_delta
forest$shift <- paste0("$", substr(changed[["TL-SCP_Database"]], 1, 1), "\\to ", substr(
  changed[["TL-SCP_MDA"]],
  1, 1
), "$")
forest <- forest[order(forest$delta), ]
forest$row <- nrow(forest):1
saveRDS(list(ranks = ranks, forest = forest), file.path(folder, "Main-Figure-7.rds"))
source(file.path(folder, "Main-Figure-7-from-rds.R"))
