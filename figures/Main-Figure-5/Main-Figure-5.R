# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/real-data-analysis/analysis.R'))
fit <- fit_clinical(root,mda=TRUE)
result <- list(yta=km_data(fit$target,'YTAClass',c('0','A','B','C','D'),.90,'YTA classes'),
               tlscp=km_data(fit$target,'TL-SCP',c('Early-onset','Mid-onset','Late-onset'),.90,'TL-SCP groups'))
saveRDS(result,file.path(folder,'Main-Figure-5.rds'))
source(file.path(folder,'Main-Figure-5-from-rds.R'))
