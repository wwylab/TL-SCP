# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))
source(file.path(root,'method/visualization/figures.R'))
d <- readRDS(file.path(folder,'Supplement-Figure-A4.rds'))
settings <- c(30,60,90,120)
label <- 'Number of predictors $p$'
methods <- c('SCP-Lasso','SCP-SCAD','SCP-MCP','SCP-Distance')
reference <- 'SCP-Distance'
controls <- c('SCP-MCP','SCP-SCAD','SCP-Lasso')
metrics <- c('relErr','NMI')
panels <- simulation_boxplots(d,settings,label,reference,controls,methods,metrics)
write.csv(attr(panels,'tests'),file.path(folder,'Supplement-Figure-A4-paired-comparisons.csv'),row.names=FALSE)
p <- wrap_plots(panels,ncol=length(metrics),guides='collect') & theme(legend.position='top')
save_figure(p,file.path(folder,'Supplement-Figure-A4.pdf'),8.3,3.9)
