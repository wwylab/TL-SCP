# LaTeX tables are built from the same mutation summaries used in the figures.
escape_tex <- function(x) {
  x <- gsub('_','\\\\_',as.character(x),fixed=TRUE)
  x <- gsub('%','\\\\%',x,fixed=TRUE)
  x
}
write_annotation_table <- function(d,file) {
  order <- match(d[['TL-SCP']],c('Early-onset','Mid-onset','Late-onset'))
  d <- d[order(order,d$ProtDescription),]
  rows <- vapply(seq_len(nrow(d)),function(i) {
    x <- d[i,]; yta <- if(x$YTAClass=='Unclassified')'Unknown' else x$YTAClass
    paste0('\\texttt{',sub('^p\\.','',x$ProtDescription),'} & ',x$target_n,' & ',x$source_n,
           ' & ',x[['TL-SCP']],' & ',x$SCP,' & ',yta,' \\\\')
  },character(1))
  writeLines(c('\\begin{tabular}{lrrlll}','\\toprule',
               'Mutation & Target carriers & Source carriers & TL-SCP & SCP & YTA \\\\',
               '\\midrule',rows,'\\midrule',
               paste0('Total & ',sum(d$target_n),' & ',sum(d$source_n),' & & & \\\\'),
               '\\bottomrule','\\end{tabular}'),file)
}
write_composition_table <- function(d,file) {
  lev <- c('0','A','B','C','D','Unclassified')
  counts <- function(a,col,g)c(sum(a$YTAClass==g),sum(a[[col]][a$YTAClass==g]))
  rows <- vapply(lev,function(g)paste(c(g,counts(d$database,'target_n',g),
                    counts(d$database,'source_n',g),counts(d$mda,'target_n',g)),collapse=' & '),character(1))
  totals <- c('Total',nrow(d$database),sum(d$database$target_n),nrow(d$database),
              sum(d$database$source_n),nrow(d$mda),sum(d$mda$target_n))
  writeLines(c('\\begin{tabular}{lrrrrrr}','\\toprule',
    '& \\multicolumn{2}{c}{Database LFS/Chompret} & \\multicolumn{2}{c}{Database non-LFS/non-Chompret} & \\multicolumn{2}{c}{MD Anderson} \\\\',
    'YTA class & Mutations & Carriers & Mutations & Carriers & Mutations & Carriers \\\\',
    '\\midrule',paste0(rows,' \\\\'),'\\midrule',paste0(paste(totals,collapse=' & '),' \\\\'),
    '\\bottomrule','\\end{tabular}'),file)
}
write_concordance_table <- function(d,file) {
  lev <- c('Early-onset','Mid-onset','Late-onset')
  cell <- function(a,b) {
    z <- d[d[['TL-SCP_Database']]==a & d[['TL-SCP_MDA']]==b,]
    if(!nrow(z))return('--')
    z <- z[order(z$ProtDescription),]
    text <- paste0('\\texttt{',sub('^p\\.','',z$ProtDescription),'}')
    if(a!=b)text <- paste0(text,' (',z$target_n_Database,', ',z$target_n_MDA,')')
    paste0('\\begin{tabular}{@{}l@{}}',paste(text,collapse='\\\\'),'\\end{tabular}')
  }
  rows <- vapply(lev,function(a)paste(c(a,vapply(lev,function(b)cell(a,b),character(1))),collapse=' & '),character(1))
  writeLines(c('\\begin{tabular}{lccc}','\\toprule',
    '& \\multicolumn{3}{c}{MD Anderson TL-SCP label} \\\\',
    'Database TL-SCP label & Early-onset & Mid-onset & Late-onset \\\\',
    '\\midrule',paste0(rows,' \\\\'),'\\bottomrule','\\end{tabular}'),file)
}
