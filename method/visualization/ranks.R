rank_comparison_plot <- function(d) {
  ranks <- d$ranks; forest <- d$forest
  colors <- c(Early='#E76F65',Mid='#258A47',Late='#5289D4')
  panel <- function(tag,title) {
    z <- ranks[ranks$cohort==tag,]
    ggplot(z,aes(weighted_rank,y)) +
      geom_hline(yintercept=1:3,color='grey92') +
      geom_point(aes(fill=group,size=highlight),shape=21,color='grey35') +
      ggrepel::geom_text_repel(data=z[z$mutation %in% sub('^p\\.','',forest$Mutation),],
             aes(label=mutation),size=3.5,seed=20261019,max.overlaps=Inf,max.iter=20000,max.time=10) +
      scale_fill_manual(values=colors,breaks=c('Early','Mid','Late')) + scale_size_manual(values=c(Other=2.2,Disconcordant=4.2),guide='none') +
      scale_x_reverse(limits=c(nrow(z)+.5,.5)) +
      scale_y_continuous(breaks=1:3,labels=c('Late','Mid','Early'),limits=c(.55,3.8)) +
      labs(title=title,x='Weighted average rank',y=NULL,fill=NULL) +
      theme_classic(base_size=11,base_family='serif') + theme(legend.position='bottom')
  }
  a <- panel('NCI','(A) TP53 Database LFS/Chompret')
  b <- panel('MDA','(B) MD Anderson prospective cohort')
  width <- diff(range(c(forest$lower,forest$upper,0)))
  shiftx <- max(forest$upper)+.24*width
  c <- ggplot(forest,aes(delta,row)) +
    geom_vline(xintercept=0,linetype='dashed',color='grey40') +
    geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.18,linewidth=.5) +
    geom_point(size=2) + geom_text(aes(x=shiftx,label=shift),size=3.1) +
    scale_y_continuous(breaks=forest$row,labels=sub('^p\\.','',forest$Mutation)) +
    labs(title='(C) Between-cohort Cox-effect differences',x='Change in initial log hazard ratio',y=NULL) +
    theme_classic(base_size=11,base_family='serif')
  ((a/b)+plot_layout(guides='collect') | c) + plot_layout(widths=c(1.1,1))
}
