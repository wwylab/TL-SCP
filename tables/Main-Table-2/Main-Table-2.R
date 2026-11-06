# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root, "method/real-data-analysis/analysis.R"))
a <- fit_clinical(root)$annotations
b <- fit_clinical(root, mda = TRUE)$annotations
result <- merge(a, b, by = "ProtDescription", suffixes = c("_Database", "_MDA"))
saveRDS(result, file.path(folder, "Main-Table-2.rds"))
source(file.path(folder, "Main-Table-2-from-rds.R"))
