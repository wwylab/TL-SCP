# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/real-data-analysis/analysis.R'))
fit <- fit_clinical(root,mda=FALSE)
result <- fit$annotations
saveRDS(result,file.path(folder,'Main-Figure-4.rds'))
source(file.path(folder,'Main-Figure-4-from-rds.R'))
