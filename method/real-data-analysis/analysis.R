# Retain the Breslow tie convention across glmnet releases.
cox_initial <- function(...) {
  args <- list(...)
  if ("cox.ties" %in% names(formals(glmnet::glmnet))) {
    args$cox.ties <- "breslow"
  }
  do.call(glmnet::glmnet, args)
}

# Prepare the two clinical analyses used in the manuscript.
library(survival)
library(Matrix)
library(glmnet)
library(dplyr)

prepare_database <- function(file, classes) {
  read.csv(file, stringsAsFactors = FALSE) |>
    filter(
      Class %in% classes, Germline_carrier %in% c("Confirmed", "Obligatory"),
      !is.na(Age) | !is.na(Age_at_diagnosis), !is.na(Sex), !is.na(ProtDescription),
      ProtDescription != "p.?"
    ) |>
    mutate(
      time = ifelse(is.na(Age_at_diagnosis), Age, Age_at_diagnosis),
      status = as.integer(!is.na(Age_at_diagnosis))
    ) |>
    filter(!is.na(time), time > 0) |>
    group_by(Individual_code) |>
    arrange(is.na(Age_at_diagnosis), Age_at_diagnosis, .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    as.data.frame()
}

prepare_mda <- function(file) {
  # Authorized users may supply the reviewed analysis CSV directly.
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    d <- read.csv(file, check.names = FALSE)
    required <- c("Individual_code", "Sex", "ProtDescription", "time", "status")
    stopifnot(all(required %in% names(d)), !anyDuplicated(d$Individual_code))
    return(d)
  }
  e <- new.env(parent = emptyenv())
  load(file, e)
  frames <- Filter(function(x) {
    is.data.frame(x) && all(c(
      "pid", "sex", "aminoacidchange",
      "diagAge"
    ) %in% names(x))
  }, as.list(e))
  if (length(frames) != 1) {
    stop("Supply the reviewed MDA data frame or the analysis CSV described in the README.")
  }
  x <- frames[[1]]
  age <- if ("age" %in% names(x)) {
    as.numeric(x$age)
  } else {
    rep(NA_real_, nrow(x))
  }
  death <- if ("ageDeath" %in% names(x)) {
    as.numeric(x$ageDeath)
  } else {
    rep(NA_real_, nrow(x))
  }
  d <- data.frame(
    Individual_code = x$pid, Sex = x$sex, ProtDescription = x$aminoacidchange,
    time = ifelse(!is.na(x$diagAge), x$diagAge, ifelse(!is.na(age), age, death)),
    status = as.integer(!is.na(x$diagAge))
  )
  if ("Family_ID" %in% names(x)) {
    d$Family_ID <- x$Family_ID
  }
  d |>
    filter(
      !is.na(Sex), !is.na(ProtDescription), ProtDescription != "p.?", !is.na(time),
      time > 0
    ) |>
    group_by(Individual_code) |>
    arrange(desc(status), time, .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    as.data.frame()
}

retain_mutations <- function(target, source) {
  nt <- table(target$ProtDescription)
  ns <- table(source$ProtDescription)
  mutations <- sort(intersect(names(nt)[nt >= 2], names(ns)))
  mutations <- mutations[!(nt[mutations] %in% c(2, 3) & ns[mutations] == 1)]
  list(target = target[target$ProtDescription %in% mutations, ], source = source[source$ProtDescription %in%
    mutations, ], mutations = mutations)
}

clinical_inputs <- function(root, mda = FALSE) {
  raw <- file.path(root, "tables/Main-Table-1/raw-data")
  source <- prepare_database(file.path(raw, "GermlineDownload_r21.csv"), c(
    "LFL",
    "FH", "noFH", "Other"
  ))
  if (mda) {
    file <- Sys.getenv("MDA_DATA")
    if (!nzchar(file) || !file.exists(file)) {
      stop("Set MDA_DATA to your authorized MDA analysis file. See the README for its required columns.")
    }
    target <- prepare_mda(file)
  } else {
    target <- prepare_database(file.path(raw, "GermlineDownload_r21.csv"), c(
      "LFS",
      "TP53_Chompret"
    ))
  }
  retain_mutations(target, source)
}

yta_classes <- function(mutations, root) {
  catalog <- read.csv(file.path(root, "tables/Main-Table-1/raw-data/yta-catalog.csv"))
  y <- as.character(catalog$cluster_rec[match(mutations, catalog$ProtDescription)])
  truncation <- grepl("*", mutations, fixed = TRUE) | grepl("fs", mutations, fixed = TRUE) |
    grepl("del", mutations, fixed = TRUE) | grepl("insIQ", mutations, fixed = TRUE)
  y[truncation] <- "0"
  y[is.na(y)] <- "Unclassified"
  y
}

fit_clinical <- function(root, mda = FALSE) {
  source(file.path(root, "method/real-data-analysis/fusion_solver.R"), local = TRUE)
  source(file.path(root, "method/real-data-analysis/fit_groups.R"), local = TRUE)
  # Vectorized working weights are algebraically identical to update_wz.
  scalar_update <- update_wz
  update_wz <- function(eta, y, status) {
    eta <- as.numeric(eta)
    a <- exp(eta)
    event <- which(status == 1)
    risk <- outer(y, y[event], ">=")
    den <- colSums(risk * a)
    pr <- outer(a, den, "/")
    if (any(!is.finite(pr))) {
      return(scalar_update(eta, y, status))
    }
    w <- rowSums(risk * pr * (1 - pr))
    score <- status - rowSums(risk * pr)
    z <- eta + score / w
    z[rowSums(risk) == 0] <- eta[rowSums(risk) == 0]
    list(w = w, z = z)
  }
  q <- clinical_inputs(root, mda)
  T <- q$target
  S <- q$source
  m <- q$mutations
  X <- build_design_matrix(T, m)
  XS <- build_design_matrix(S, m)
  Y <- cbind(time = T$time, status = T$status)
  YS <- cbind(time = S$time, status = S$status)
  bt <- as.numeric(cox_initial(X, Y, family = "cox", alpha = 0, lambda = 1e-04)$beta[-1])
  bs <- as.numeric(cox_initial(XS, YS, family = "cox", alpha = 0, lambda = 1e-04)$beta[-1])
  nt <- as.numeric(table(factor(T$ProtDescription, levels = m)))
  ns <- as.numeric(table(factor(S$ProtDescription, levels = m)))
  ranks <- data.frame(
    Mutation = m, target_n = nt, source_n = ns, target_beta_raw = bt,
    source_beta_raw = bs, target_rank = rank(bt, ties.method = "first"), source_rank = rank(bs,
      ties.method = "first"
    ), target_weight = nt / (nt + ns), weighted_rank = (nt *
      rank(bt, ties.method = "first") + ns * rank(bs, ties.method = "first")) / (nt +
      ns)
  )
  ann <- data.frame(
    ProtDescription = m, YTAClass = yta_classes(m, root), target_n = nt,
    source_n = ns
  )
  fits <- list()
  for (method in c("SCP", "TL-SCP")) {
    tau <- if (method == "SCP") {
      order(bt)
    } else {
      wRankAveg(bs, bt, ns / (nt + ns), nt / (nt + ns))
    }
    XO <- cbind(X[, 1, drop = FALSE], X[, -1, drop = FALSE][, tau, drop = FALSE])
    k <- if (mda) {
      4L
    } else {
      1:5
    }
    f <- fit_grouped_method(X, Y, XO, cbind(0, myGetDkn(1, length(m))), tau, k,
      covar_count = 1, method_name = method
    )
    b <- f$beta[-1]
    cl <- cluster_label_dist(b, f$nclusters)
    means <- tapply(b, cl, mean)
    order <- rank(-means, ties.method = "first")[as.character(cl)]
    labels <- c("Early-onset", "Mid-onset", "Late-onset", "Group 4")[order]
    if (mda && method == "TL-SCP") {
      # Reported grouping: the fourth fitted group is Late; R337H joins Mid by
      # rank.
      labels[labels == "Group 4"] <- "Late-onset"
      ann$raw_TLSCP_group <- labels
      labels[m == "p.R337H"] <- "Mid-onset"
    }
    ann[[method]] <- labels
    fits[[method]] <- f
  }
  for (col in c("YTAClass", "SCP", "TL-SCP")) {
    T[[col]] <- ann[[col]][match(
      T$ProtDescription,
      m
    )]
  }
  list(target = T, source = S, annotations = ann, ranks = ranks, fits = fits)
}

km_data <- function(d, column, levels, confidence = 0.9, title = column) {
  d <- d[d[[column]] %in% levels, ]
  d$Group <- factor(d[[column]], levels = levels)
  fit <- survfit(Surv(time, status) ~ Group, d, conf.int = confidence)
  s <- summary(fit)
  curve <- data.frame(
    time = s$time, surv = s$surv, lower = s$lower, upper = s$upper,
    Group = sub("^Group=", "", as.character(s$strata))
  )
  curve <- rbind(
    data.frame(time = 0, surv = 1, lower = 1, upper = 1, Group = levels),
    curve
  )
  curve <- curve[order(curve$Group, curve$time), ]
  curve$lower[is.na(curve$lower)] <- curve$surv[is.na(curve$lower)]
  curve$upper[is.na(curve$upper)] <- curve$surv[is.na(curve$upper)]
  counts <- do.call(rbind, lapply(levels, function(g) {
    z <- d[d$Group == g, ]
    data.frame(Group = g, MutationN = length(unique(z$ProtDescription)), PatientN = nrow(z))
  }))
  lr <- survdiff(Surv(time, status) ~ Group, d)
  list(
    curve = curve, counts = counts, p = pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE),
    confidence = confidence, title = title
  )
}
