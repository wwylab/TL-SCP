# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/visualization/figures.R"))
d <- readRDS(file.path(folder, "Main-Figure-2.rds"))
d$cox$ProtDescription <- d$cox$Mutation
panels <- list(
  alluvial_panel(d$dcsm, "YTAClass", "majority_DCSM", "YTA", "DCSM", "(a) DCSM"),
  alluvial_panel(
    d$cox, "group_source", "group_target", "Non-LFS/\\\\non-Chompret",
    "LFS/Chompret", "(b) Cox-K-means"
  )
)
save_alluvial(panels, file.path(folder, "Main-Figure-2.pdf"), root)
