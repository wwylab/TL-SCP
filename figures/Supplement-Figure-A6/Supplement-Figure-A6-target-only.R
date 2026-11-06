# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/simulation/helpers.R'))
source(file.path(root,'method/methodology/scp_tlscp.R'))
use_fast_distance(root)
controls <- generate_target_controls(1000,as.integer(Sys.getenv('TLSCP_REPS','50')),FALSE)
saveRDS(controls,file.path(folder,'Supplement-Figure-A6-target-only.rds'))
