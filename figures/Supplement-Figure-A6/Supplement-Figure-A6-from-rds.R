# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/visualization/figures.R"))
d <- readRDS(file.path(folder, "Supplement-Figure-A6.rds"))
settings <- sort(unique(d$setting))
label <- "Covariate shift level $\\nu$"
methods <- c("Cox-K-means", "SCP", "TL-SCP")
reference <- "TL-SCP"
controls <- c("SCP", "Cox-K-means")
metrics <- c("relErr", "NMI")
rows <- list()
tests <- list()
for (ns in sort(unique(d$n_S))) {
  p <- simulation_boxplots(d[d$n_S == ns, ], settings, label, reference, controls,
    methods, metrics,
    scale_data = d
  )
  tests[[as.character(ns)]] <- transform(attr(p, "tests"), n_S = ns)
  p <- lapply(p, function(panel) panel + labs(title = paste0("$n_S=", ns, "$")))
  rows[[as.character(ns)]] <- wrap_plots(p, ncol = 2)
}
p <- (rows[[1]] / rows[[2]]) + plot_layout(guides = "collect") & theme(legend.position = "top")
write.csv(do.call(rbind, tests), file.path(folder, "Supplement-Figure-A6-paired-comparisons.csv"),
  row.names = FALSE
)
save_figure(p, file.path(folder, "Supplement-Figure-A6.pdf"), 9, 6.6)
