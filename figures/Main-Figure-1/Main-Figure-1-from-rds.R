# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/figures.R'))
d <- readRDS(file.path(folder,'Main-Figure-1.rds'))
p <- km_panel(d)
save_figure(p, file.path(folder,'Main-Figure-1.pdf'), 5.2, 3.6)
