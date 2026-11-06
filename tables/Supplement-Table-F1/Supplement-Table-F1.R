# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/real-data-analysis/analysis.R'))
result <- fit_clinical(root,mda=FALSE)$annotations
saveRDS(result,file.path(folder,'Supplement-Table-F1.rds'))
source(file.path(folder,'Supplement-Table-F1-from-rds.R'))
