# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/tables.R'))
d <- readRDS(file.path(folder,'Main-Table-2.rds'))
write_concordance_table(d,file.path(folder,'Main-Table-2.tex'))
