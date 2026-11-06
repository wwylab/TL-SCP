# Run with Rscript from any directory.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
folder <- dirname(normalizePath(sub("^--file=", "", arg[1])))
root <- normalizePath(file.path(folder, "../.."))

source(file.path(root,'method/real-data-analysis/analysis.R'))
inputs <- clinical_inputs(root)
annotations <- yta_classes(inputs$mutations,root)
target <- inputs$target
target$YTAClass <- annotations[match(target$ProtDescription,inputs$mutations)]
# DCSM is a patient-clustering comparison. Only mutation-level counts are exported.
tmp <- tempfile('dcsm-');dir.create(tmp)
write.csv(target,file.path(tmp,'target.csv'),row.names=FALSE)
python <- Sys.getenv('PYTHON','python3')
status <- system2(python,c(shQuote(file.path(root,'method/real-data-analysis/dcsm.py')),
                 shQuote(file.path(tmp,'target.csv')),shQuote(file.path(tmp,'mutation-counts.csv'))))
if(status!=0) { unlink(tmp,recursive=TRUE);stop('DCSM fitting failed.') }
dcsm <- read.csv(file.path(tmp,'mutation-counts.csv'),check.names=FALSE)
unlink(tmp,recursive=TRUE)
assignments <- list()
for(role in c('source','target')) {
  q <- inputs[[role]]; m <- sort(unique(q$ProtDescription))
  Z <- outer(q$ProtDescription,m,'==')*1;male <- as.integer(q$Sex=='M')
  Y <- Surv(q$time,q$status)
  model <- coxph(Y~male+ridge(Z,theta=nrow(q)*1e-4,scale=FALSE),ties='breslow',
                control=coxph.control(iter.max=200,eps=1e-10))
  beta <- as.numeric(coef(model))[-1]
  set.seed(20260928);km <- kmeans(matrix(beta,ncol=1),3,nstart=100,iter.max=100)
  group <- c('Early','Mid','Late')[match(km$cluster,order(as.numeric(km$centers),decreasing=TRUE))]
  assignments[[role]] <- data.frame(Mutation=m,beta=beta,group=group)
}
cox <- merge(assignments$source,assignments$target,by='Mutation',suffixes=c('_source','_target'))
saveRDS(list(dcsm=dcsm,cox=cox),file.path(folder,'Main-Figure-2.rds'))
source(file.path(folder,'Main-Figure-2-from-rds.R'))
