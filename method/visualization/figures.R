# Plotting functions shared by the manuscript figures.
library(ggplot2)
library(patchwork)

onset <- c("Early-onset", "Mid-onset", "Late-onset")
onset_colors <- setNames(c("#F8766D", "#00BA38", "#619CFF"), onset)
yta_colors <- c(
  `0` = "#999999", A = "#F8766D", B = "#00BA38", C = "#619CFF", D = "#B279A2",
  Unclassified = "#aaaaaa"
)

save_figure <- function(plot, file, width, height) {
  # TikZ retains the Palatino text used in the manuscript.
  compiler <- Sys.which("pdflatex")
  if (!nzchar(compiler)) {
    stop("Install a TeX distribution and put pdflatex on PATH.")
  }
  tex <- sub("\\.pdf$", ".tex", file)
  options(tikzLatex = compiler)
  packages <- c(
    "\\usepackage{mathpazo}", "\\usepackage[T1]{fontenc}", "\\usepackage{amsmath}",
    "\\usepackage{tikz}", "\\usepackage[active,tightpage,psfixbb]{preview}", "\\PreviewEnvironment{pgfpicture}",
    "\\setlength\\PreviewBorder{0pt}"
  )
  tikzDevice::tikz(tex,
    width = width, height = height, standAlone = TRUE, packages = packages,
    sanitize = FALSE
  )
  print(plot)
  dev.off()
  compile_figure(tex)
}

km_panel <- function(x, title = NULL, y_label = TRUE) {
  lev <- as.character(x$counts$Group)
  x$curve$Group <- factor(x$curve$Group, levels = lev)
  colors <- if (all(lev %in% names(yta_colors))) {
    yta_colors
  } else {
    onset_colors
  }
  labels <- setNames(paste0(
    sub("-onset", "", lev), ": ", x$counts$MutationN, " (",
    x$counts$PatientN, ")"
  ), lev)
  ptext <- if (x$p < 0.001) {
    "Log-rank $p < 0.001$"
  } else {
    sprintf("Log-rank $p = %.3f$", x$p)
  }
  ggplot(x$curve, aes(time, surv, color = Group, fill = Group)) +
    geom_ribbon(aes(
      ymin = lower,
      ymax = upper
    ), alpha = 0.18, color = NA) +
    geom_step(linewidth = 0.8) +
    scale_color_manual(
      values = colors,
      breaks = lev, labels = labels
    ) +
    scale_fill_manual(
      values = colors, breaks = lev,
      labels = labels
    ) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(xlim = c(
      0,
      105
    ), ylim = c(0, 1)) +
    labs(title = title, subtitle = ptext, x = "Age", y = if (y_label) {
      "Cancer-free probability"
    } else {
      NULL
    }, color = x$title, fill = x$title) +
    theme_bw(base_size = 13.5, base_family = "serif") +
    guides(
      color = guide_legend(ncol = 1, byrow = TRUE, title.position = "top"),
      fill = guide_legend(ncol = 1, byrow = TRUE, title.position = "top")
    ) +
    theme(
      legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(
        1,
        1
      ), legend.background = element_rect(
        fill = "white", color = "grey75",
        linewidth = 0.2
      ), legend.key.width = grid::unit(0.8, "lines"), legend.key.spacing.y = grid::unit(
        0,
        "pt"
      ), legend.margin = margin(2, 3, 2, 3), axis.title = element_text(size = 13.5),
      axis.text = element_text(size = 11.5), legend.title = element_text(size = 9),
      legend.text = element_text(size = 8.5), legend.key.height = grid::unit(
        0.45,
        "lines"
      ), panel.grid.minor = element_blank()
    )
}

# Assemble aggregate mutation annotations for the TikZ alluvial renderer.
alluvial_panel <- function(d, left, right, left_title, right_title, title = "") {
  clean <- function(x) sub("-onset", "", as.character(x))
  left_group <- clean(d[[left]])
  right_group <- clean(d[[right]])
  ordering <- c(
    "0", "A", "B", "C", "D", "Unclassified", "Early", "Mid", "Late",
    "Group 4", "Tied"
  )
  data.frame(
    mutation = sub("^p\\.", "", d$ProtDescription), left = left_group, right = right_group,
    left_order = paste(intersect(ordering, unique(left_group)), collapse = "|"),
    right_order = paste(intersect(ordering, unique(right_group)), collapse = "|"),
    left_name = left_title, right_name = right_title, title = title
  )
}

save_alluvial <- function(panels, file, root) {
  rows <- do.call(rbind, lapply(seq_along(panels), function(i) {
    transform(panels[[i]], panel = i)
  }))
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv))
  write.csv(rows, csv, row.names = FALSE)
  tex <- sub("\\.pdf$", ".tex", file)
  renderer <- file.path(root, "method/visualization/alluvials.py")
  status <- system2(Sys.getenv("PYTHON", "python3"), shQuote(c(renderer, csv, tex)))
  if (status != 0) {
    stop("Alluvial rendering failed: ", file)
  }
  compile_figure(tex)
}

compile_figure <- function(tex) {
  compiler <- Sys.which("pdflatex")
  if (!nzchar(compiler)) {
    stop("Install a TeX distribution and put pdflatex on PATH.")
  }
  log <- sub("\\.tex$", ".log.txt", tex)
  status <- system2(compiler, c("-interaction=nonstopmode", "-halt-on-error", paste0(
    "-output-directory=",
    shQuote(dirname(tex))
  ), shQuote(tex)), stdout = log, stderr = log)
  if (status != 0) {
    stop("Figure compilation failed: ", tex)
  }
  for (ext in c(".aux", ".log", ".log.txt")) unlink(sub("\\.tex$", ext, tex))
}

simulation_boxplots <- function(
  d, settings, setting_label, reference, controls, method_levels,
  metrics = c("relErr", "NMI"), scale_data = d
) {
  d$method <- factor(d$method, levels = method_levels)
  d$setting <- factor(d$setting, levels = settings)
  colors <- c(
    `Cox-K-means` = "black", SCP = "red2", `TL-SCP` = "purple2", `SCP-Lasso` = "darkcyan",
    `SCP-SCAD` = "blue", `SCP-MCP` = "purple", `SCP-Distance` = "red", `SCP-Allpair` = "black"
  )
  panels <- list()
  tests <- list()
  source_size <- if ("n_S" %in% names(scale_data)) {
    scale_data$n_S
  } else {
    rep(1, nrow(scale_data))
  }
  for (metric in metrics) {
    dd <- d[is.finite(d[[metric]]), ]
    whiskers <- unlist(lapply(split(scale_data[[metric]], interaction(source_size,
      scale_data$setting, scale_data$method,
      drop = TRUE
    )), function(x) boxplot.stats(x)$stats))
    top <- max(whiskers)
    transfer <- reference == "TL-SCP"
    if (metric == "NMI") {
      base <- if (transfer) {
        max(0.85, top + 0.06)
      } else {
        max(1, top) + 0.035
      }
      step <- if (transfer) {
        0.11
      } else {
        0.105
      }
      tip <- 0.015
      ceiling <- if (transfer) {
        base + 0.21
      } else {
        base + (length(controls) - 1) * step + 4 * tip
      }
    } else if (transfer) {
      top <- top * 1.04
      base <- top * 1.04
      step <- top * 0.13
      tip <- top * 0.015
      ceiling <- top * 1.3
    } else {
      base <- max(top * 1.06, 1e-05)
      step <- base * 0.13
      tip <- base * 0.018
      ceiling <- base + (length(controls) - 1) * step + 4 * tip
    }
    brackets <- list()
    for (setting in settings) {
      pvals <- numeric(length(controls))
      ns <- integer(length(controls))
      for (j in seq_along(controls)) {
        a <- dd[dd$setting == setting & dd$method == reference, c("rep", metric)]
        b <- dd[dd$setting == setting & dd$method == controls[j], c(
          "rep",
          metric
        )]
        pairs <- merge(a, b, by = "rep")
        delta <- pairs[[2]] - pairs[[3]]
        ns[j] <- length(delta)
        pvals[j] <- if (all(delta == 0)) {
          1
        } else {
          wilcox.test(delta, exact = FALSE, correct = TRUE)$p.value
        }
      }
      adjusted <- p.adjust(pvals, "holm")
      for (j in seq_along(controls)) {
        i <- match(setting, settings)
        offset <- function(m) {
          0.72 * (match(m, method_levels) - (length(method_levels) +
            1) / 2) / length(method_levels)
        }
        y <- base + step * (j - 1)
        brackets[[length(brackets) + 1]] <- data.frame(
          x = i + offset(controls[j]),
          xend = i + offset(reference), y = y, label = if (adjusted[j] < 0.001) {
            "$p<0.001$"
          } else {
            sprintf("$p=%.3f$", adjusted[j])
          }
        )
        tests[[length(tests) + 1]] <- data.frame(metric, setting, reference,
          control = controls[j], n = ns[j], p = pvals[j], p_Holm = adjusted[j]
        )
      }
    }
    b <- do.call(rbind, brackets)
    b$label_y <- b$y + tip
    b$bracket_bottom <- b$y - tip
    p <- ggplot(dd, aes(setting, .data[[metric]], color = method)) +
      geom_boxplot(
        width = 0.72,
        outlier.shape = NA, position = if (transfer) {
          position_dodge2()
        } else {
          position_dodge(0.72)
        }, linewidth = 0.45, fill = "white"
      ) +
      geom_segment(
        data = b,
        aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3,
        color = "black"
      ) +
      geom_segment(data = b, aes(
        x = x, xend = x, y = bracket_bottom,
        yend = y
      ), inherit.aes = FALSE, linewidth = 0.3, color = "black") +
      geom_segment(
        data = b,
        aes(x = xend, xend = xend, y = bracket_bottom, yend = y), inherit.aes = FALSE,
        linewidth = 0.3, color = "black"
      ) +
      geom_text(
        data = b, aes(
          x = (x + xend) / 2,
          y = label_y, label = label
        ), inherit.aes = FALSE, size = 3.5, vjust = 0,
        color = "black"
      ) +
      scale_color_manual(
        values = colors[method_levels], drop = FALSE,
        name = NULL
      ) +
      labs(x = setting_label, y = metric) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(), axis.title = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 12), legend.position = "top",
        legend.key.width = grid::unit(1.2, "lines")
      ) +
      coord_cartesian(ylim = c(
        0,
        ceiling
      ))
    if (metric == "NMI") {
      p <- p + scale_y_continuous(breaks = seq(0, 1, 0.25))
    }
    panels[[metric]] <- p
  }
  attr(panels, "tests") <- do.call(rbind, tests)
  panels
}

yta_zero_a_plot <- function(results) {
  d <- results$curve
  counts <- results$counts
  p <- results$p
  labels <- setNames(paste0(
    counts$Group, " (", counts$PatientN, " patients)\n",
    counts$MutationN, " mutations"
  ), counts$Group)
  plot <- ggplot(d, aes(time, surv, color = Group, fill = Group)) +
    geom_ribbon(aes(
      ymin = lower,
      ymax = upper
    ), alpha = 0.13, colour = NA) +
    geom_step(linewidth = 0.7) +
    scale_color_manual(values = c(
      `0` = "black",
      A = "#FF00FF"
    ), labels = labels) +
    scale_fill_manual(values = c(
      `0` = "black",
      A = "#FF00FF"
    ), labels = labels) +
    annotate("text", x = 2, y = 0.08, label = sprintf(
      "Log-rank $p=%.3f$",
      p
    ), hjust = 0, size = 5.3) +
    coord_cartesian(xlim = c(0, 70), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 70)) +
    scale_y_continuous(breaks = seq(
      0,
      1, 0.25
    )) +
    labs(
      x = "Age (years)", y = "Cancer-free survival\nprobability",
      color = "YTA group", fill = "YTA group"
    ) +
    theme_classic(base_size = 17, base_family = "serif") +
    theme(
      legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(
        1,
        1
      ), axis.text = element_text(size = 15), legend.title = element_text(size = 14),
      legend.text = element_text(size = 13), legend.key.height = grid::unit(
        22,
        "pt"
      ), legend.background = element_rect(fill = "white", color = NA)
    )

  plot
}
