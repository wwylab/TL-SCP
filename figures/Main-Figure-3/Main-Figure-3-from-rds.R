# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/figures.R'))
d <- readRDS(file.path(folder,'Main-Figure-3.rds'))
p <- km_panel(d$yta) + km_panel(d$tlscp, y_label=FALSE)
save_figure(p,file.path(folder,'Main-Figure-3.pdf'),8.2,2.65)
