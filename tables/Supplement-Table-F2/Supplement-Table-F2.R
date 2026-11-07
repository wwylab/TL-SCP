# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/real-data-analysis/analysis.R"))
result <- estimate_groups(root, mda = TRUE)$annotations
saveRDS(result, file.path(folder, "Supplement-Table-F2.rds"))
source(file.path(folder, "Supplement-Table-F2-from-rds.R"))
