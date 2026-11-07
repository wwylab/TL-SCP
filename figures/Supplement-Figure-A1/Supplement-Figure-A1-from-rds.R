# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/visualization/figures.R"))
d <- readRDS(file.path(folder, "Supplement-Figure-A1.rds"))
settings <- c(0.1, 0.3, 0.7, 1)
label <- "Noise level $\\sigma$"
methods <- c("SCP-Lasso", "SCP-Allpair")
reference <- "SCP-Lasso"
controls <- "SCP-Allpair"
metrics <- c("relErr", "NMI", "Runtime")
panels <- simulation_boxplots(d, settings, label, reference, controls, methods, metrics)
write.csv(attr(panels, "tests"), file.path(folder, "Supplement-Figure-A1-paired-comparisons.csv"),
  row.names = FALSE
)
p <- wrap_plots(panels, ncol = length(metrics), guides = "collect") & theme(legend.position = "top")
save_figure(p, file.path(folder, "Supplement-Figure-A1.pdf"), 10.2, 3.6)
