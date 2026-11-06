# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/visualization/figures.R"))
d <- readRDS(file.path(folder, "Main-Figure-6.rds"))
panels <- list(alluvial_panel(d, "YTAClass", "TL-SCP", "YTA", "TL-SCP"))
save_alluvial(panels, file.path(folder, "Main-Figure-6.pdf"), root)
