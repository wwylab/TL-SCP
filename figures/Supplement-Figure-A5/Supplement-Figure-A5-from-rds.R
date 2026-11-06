# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/visualization/figures.R"))
d <- readRDS(file.path(folder, "Supplement-Figure-A5.rds"))
settings <- sort(unique(d$setting))
label <- "Source sample size $n_S$"
methods <- c("Cox-K-means", "SCP", "TL-SCP")
reference <- "TL-SCP"
controls <- c("SCP", "Cox-K-means")
metrics <- c("relErr", "NMI")
panels <- simulation_boxplots(d, settings, label, reference, controls, methods, metrics)
write.csv(attr(panels, "tests"), file.path(folder, "Supplement-Figure-A5-paired-comparisons.csv"),
  row.names = FALSE
)
p <- wrap_plots(panels, ncol = length(metrics), guides = "collect") & theme(legend.position = "top")
save_figure(p, file.path(folder, "Supplement-Figure-A5.pdf"), 8.3, 3.6)
