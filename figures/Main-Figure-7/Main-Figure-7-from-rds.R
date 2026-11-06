# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/figures.R'))
d <- readRDS(file.path(folder,'Main-Figure-7.rds'))
source(file.path(root,'method/visualization/ranks.R'))
p <- rank_comparison_plot(d)
save_figure(p,file.path(folder,'Main-Figure-7.pdf'),12,5.6)
