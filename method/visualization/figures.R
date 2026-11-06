# Plotting functions shared by the manuscript figures.
library(ggplot2)
library(patchwork)

onset <- c('Early-onset', 'Mid-onset', 'Late-onset')
onset_colors <- setNames(c('#F8766D', '#00BA38', '#619CFF'), onset)
yta_colors <- c('0' = '#999999', A = '#F8766D', B = '#00BA38',
                C = '#619CFF', D = '#B279A2', Unclassified = '#aaaaaa')

save_figure <- function(plot, file, width, height) {
  # TikZ retains the Palatino text used in the manuscript.
  compiler <- Sys.which('pdflatex')
  if (!nzchar(compiler)) stop('Install a TeX distribution and put pdflatex on PATH.')
  tex <- sub('\\.pdf$', '.tex', file)
  options(tikzLatex = compiler)
  packages <- c('\\usepackage{mathpazo}', '\\usepackage[T1]{fontenc}',
                '\\usepackage{amsmath}', '\\usepackage{tikz}',
                '\\usepackage[active,tightpage,psfixbb]{preview}',
                '\\PreviewEnvironment{pgfpicture}', '\\setlength\\PreviewBorder{0pt}')
  tikzDevice::tikz(tex, width = width, height = height, standAlone = TRUE,
                  packages = packages, sanitize = FALSE)
  print(plot)
  dev.off()
  code <- system2(compiler, c('-interaction=nonstopmode', '-halt-on-error',
                    paste0('-output-directory=', shQuote(dirname(tex))), shQuote(tex)),
                  stdout = sub('\\.pdf$', '.log.txt', file), stderr = sub('\\.pdf$', '.log.txt', file))
  if (code != 0) stop('Figure compilation failed: ', tex)
  for (ext in c('.aux', '.log', '.log.txt')) unlink(sub('\\.pdf$', ext, file))
}

km_panel <- function(x, title = NULL, y_label = TRUE) {
  lev <- as.character(x$counts$Group)
  colors <- if (all(lev %in% names(yta_colors))) yta_colors else onset_colors
  labels <- setNames(paste0(sub('-onset', '', lev), ': ', x$counts$MutationN,
                           ' (', x$counts$PatientN, ')'), lev)
  ptext <- if (x$p < .001) 'Log-rank $p < 0.001$' else sprintf('Log-rank $p = %.3f$', x$p)
  ggplot(x$curve, aes(time, surv, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .18, color = NA) +
    geom_step(linewidth = .8) +
    scale_color_manual(values = colors, breaks = lev, labels = labels) +
    scale_fill_manual(values = colors, breaks = lev, labels = labels) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(xlim = c(0, 105), ylim = c(0, 1)) +
    labs(title = title, subtitle = ptext, x = 'Age',
         y = if (y_label) 'Cancer-free probability' else NULL,
         color = x$title, fill = x$title) +
    theme_bw(base_size = 13.5, base_family = 'serif') +
    theme(legend.position = 'inside', legend.position.inside = c(.98, .98),
          legend.justification = c(1, 1), legend.background = element_rect(fill = 'white', color = 'grey75'),
          legend.title = element_text(size = 9), legend.text = element_text(size = 8.5),
          legend.key.height = grid::unit(.45, 'lines'), panel.grid.minor = element_blank())
}

alluvial_panel <- function(d, left, right, left_title, right_title, title = NULL) {
  # Each ribbon represents one mutation, with its name on both sides.
  clean <- function(z) sub('-onset', '', as.character(z))
  d$left <- clean(d[[left]]); d$right <- clean(d[[right]])
  d$mutation <- sub('^p\\.', '', d$ProtDescription)
  order_groups <- function(x) {
    ref <- c('0', 'A', 'B', 'C', 'D', 'Unclassified', 'Early', 'Mid', 'Late', 'Group 4', 'Tied')
    c(intersect(ref, unique(x)), setdiff(unique(x), ref))
  }
  pos <- function(group, other) {
    lev <- order_groups(group)
    idx <- order(match(group, lev), match(other, order_groups(other)), d$mutation)
    y <- numeric(length(idx)); gaps <- cumsum(c(0, diff(match(group[idx], lev)) != 0))
    y[idx] <- length(idx) - seq_along(idx) + gaps[length(gaps)] - gaps
    y
  }
  d$y1 <- pos(d$left, d$right); d$y2 <- pos(d$right, d$left)
  colors <- c(Early = '#F8766D', Mid = '#00BA38', Late = '#619CFF',
              'Group 4' = '#888888', Tied = '#888888', yta_colors)
  ribbons <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    t <- seq(0, 1, length.out = 45); ease <- 3*t^2 - 2*t^3
    yy <- d$y1[i] + (d$y2[i] - d$y1[i])*ease
    data.frame(x = c(t, rev(t)), y = c(yy+.38, rev(yy-.38)),
               mutation = d$mutation[i], group = d$right[i])
  }))
  ymax <- max(d$y1, d$y2) + 2.4
  ggplot() + geom_polygon(data = ribbons, aes(x, y, group = mutation, fill = group), alpha = .4) +
    geom_segment(data = d, aes(x = -.025, xend = .025, y = y1, yend = y1, color = right), linewidth = 2) +
    geom_segment(data = d, aes(x = .975, xend = 1.025, y = y2, yend = y2, color = right), linewidth = 2) +
    geom_text(data = d, aes(x = -.04, y = y1, label = paste0(mutation, ' [', left, ']'), color = right),
              hjust = 1, size = 2.7, family = 'serif') +
    geom_text(data = d, aes(x = 1.04, y = y2, label = paste0(mutation, ' [', right, ']'), color = right),
              hjust = 0, size = 2.7, family = 'serif') +
    annotate('text', x = c(0, 1), y = ymax, label = c(left_title, right_title), size = 3.8, family = 'serif') +
    scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
    coord_cartesian(xlim = c(-.62, 1.65), ylim = c(-1, ymax+1), clip = 'off') +
    labs(title = title) + theme_void(base_family = 'serif') +
    theme(legend.position = 'none', plot.title = element_text(size = 13, hjust = .5),
          plot.margin = margin(5, 2, 5, 2))
}

simulation_boxplots <- function(d, settings, setting_label, reference, controls,
                               method_levels, metrics = c('relErr', 'NMI')) {
  d$method <- factor(d$method, levels = method_levels)
  d$setting <- factor(d$setting, levels = settings)
  colors <- c('Cox-K-means' = 'black', SCP = 'red2', 'TL-SCP' = 'purple2',
              'SCP-Lasso' = 'darkcyan', 'SCP-SCAD' = 'blue', 'SCP-MCP' = 'purple',
              'SCP-Distance' = 'red', 'SCP-Allpair' = 'black')
  panels <- list(); tests <- list()
  for (metric in metrics) {
    dd <- d[is.finite(d[[metric]]), ]
    span <- diff(range(dd[[metric]])); if (span == 0) span <- 1
    brackets <- list()
    for (setting in settings) {
      pvals <- numeric(length(controls)); ns <- integer(length(controls))
      for (j in seq_along(controls)) {
        a <- dd[dd$setting == setting & dd$method == reference, c('rep', metric)]
        b <- dd[dd$setting == setting & dd$method == controls[j], c('rep', metric)]
        pairs <- merge(a, b, by = 'rep'); delta <- pairs[[2]] - pairs[[3]]
        ns[j] <- length(delta)
        pvals[j] <- if (all(delta == 0)) 1 else wilcox.test(delta, exact = FALSE, correct = TRUE)$p.value
      }
      adjusted <- p.adjust(pvals, 'holm')
      for (j in seq_along(controls)) {
        i <- match(setting, settings); offset <- function(m) .75*(match(m, method_levels)-(length(method_levels)+1)/2)/length(method_levels)
        y <- max(dd[dd$setting == setting, metric]) + span*(.08 + .13*(j-1))
        brackets[[length(brackets)+1]] <- data.frame(x=i+offset(controls[j]), xend=i+offset(reference), y=y,
                                  label=if(adjusted[j]<.001)'$p<0.001$' else sprintf('$p=%.3f$',adjusted[j]))
        tests[[length(tests)+1]] <- data.frame(metric, setting, reference, control=controls[j],
                                              n=ns[j], p=pvals[j], p_Holm=adjusted[j])
      }
    }
    b <- do.call(rbind, brackets)
    b$label_y <- b$y + span*.035
    panels[[metric]] <- ggplot(dd, aes(setting, .data[[metric]], fill=method)) +
      geom_boxplot(width=.7, outlier.shape=NA, position=position_dodge(.75), linewidth=.3) +
      geom_segment(data=b,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,linewidth=.3) +
      geom_text(data=b,aes(x=(x+xend)/2,y=label_y,label=label),inherit.aes=FALSE,size=3) +
      scale_fill_manual(values=colors,drop=FALSE) + labs(x=setting_label,y=metric,fill=NULL) +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),legend.position='top')
  }
  attr(panels,'tests') <- do.call(rbind,tests)
  panels
}
