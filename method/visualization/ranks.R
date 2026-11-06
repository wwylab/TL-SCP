rank_comparison_plot <- function(data) {
  d <- data$ranks
  forest <- data$forest
  pal <- c(Early = "#E76F65", Mid = "#258A47", Late = "#5289D4")
  discordant <- forest$Mutation
  d$group <- factor(d$group, levels = c("Early", "Mid", "Late"))
  d$highlight <- factor(d$highlight, levels = c("Other", "Disconcordant"))
  panel <- function(tag, title) {
    z <- d[d$cohort == tag, ]
    z <- z[order(z$highlight), ]
    ggplot(z, aes(weighted_rank, y)) +
      geom_hline(
        yintercept = 1:3, color = "grey92",
        linewidth = 0.3
      ) +
      geom_point(aes(
        fill = group, size = highlight, alpha = highlight,
        colour = highlight, stroke = highlight
      ), shape = 21) +
      ggrepel::geom_text_repel(
        data = z[z$Mutation %in%
          discordant, ], aes(label = mutation), family = "serif", size = 3.5, nudge_y = 0.32,
        box.padding = 0.32, point.padding = 0.28, min.segment.length = 0, segment.color = "grey45",
        segment.size = 0.3, seed = 20261019, max.overlaps = Inf, max.iter = 20000,
        max.time = 10, ylim = c(0.78, 3.7), force = 3
      ) +
      scale_fill_manual(
        values = pal,
        drop = FALSE
      ) +
      scale_colour_manual(
        values = c(Other = "transparent", Disconcordant = "black"),
        breaks = "Disconcordant"
      ) +
      scale_size_manual(
        values = c(Other = 2.2, Disconcordant = 4.2),
        guide = "none"
      ) +
      scale_alpha_manual(
        values = c(Other = 0.3, Disconcordant = 1),
        guide = "none"
      ) +
      scale_discrete_manual(aesthetics = "stroke", values = c(
        Other = 0,
        Disconcordant = 1.5
      ), guide = "none") +
      scale_x_reverse(limits = c(nrow(z) +
        0.5, 0.5), breaks = seq(5, nrow(z), 5), expand = expansion(mult = c(
        0.01,
        0.015
      ))) +
      scale_y_continuous(limits = c(0.55, 3.8), breaks = 1:3, labels = c(
        "Late",
        "Mid", "Early"
      ), expand = c(0, 0)) +
      labs(title = title, x = expression("Weighted average rank" ~
        italic(r) * "(" * italic(j) * ")"), y = NULL) +
      theme_classic(
        base_size = 11,
        base_family = "serif"
      ) +
      theme(
        plot.title = element_text(
          size = 11, hjust = 0,
          margin = margin(b = 6)
        ), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9),
        plot.margin = margin(8, 10, 7, 8)
      ) +
      guides(fill = guide_legend(
        order = 1,
        override.aes = list(
          shape = 21, size = 3, alpha = 1, colour = "transparent",
          stroke = 0
        )
      ), colour = guide_legend(order = 2, override.aes = list(
        shape = 21,
        fill = "grey75", size = 4.2, alpha = 1, stroke = 1.5
      )))
  }
  A <- panel("NCI", expression("(A)" ~ italic("TP53") ~ "Database LFS/Chompret"))
  B <- panel("MDA", "(B) MD Anderson prospective cohort")
  left <- (A / B) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  if (nrow(forest)) {
    span <- range(c(forest$lower, forest$upper, 0))
    width <- diff(span)
    span <- span + c(-0.08, 0.45) * width
    shiftx <- max(forest$upper) + 0.24 * width
    C <- ggplot(forest, aes(delta, row)) +
      geom_vline(
        xintercept = 0, linetype = "dashed",
        color = "grey40", linewidth = 0.4
      ) +
      geom_errorbar(aes(xmin = lower, xmax = upper),
        orientation = "y", width = 0.18, linewidth = 0.5
      ) +
      geom_point(size = 2) +
      geom_text(aes(x = shiftx, label = shift), size = 3.1) +
      annotate("text",
        x = shiftx, y = nrow(forest) + 0.7, label = "Shift", size = 3.1
      ) +
      scale_y_continuous(
        breaks = forest$row,
        labels = sub("^p\\.", "", forest$Mutation), limits = c(0.5, nrow(forest) +
          1), expand = c(0, 0)
      ) +
      scale_x_continuous(limits = span, expand = c(
        0,
        0
      )) +
      labs(
        title = "(C) Between-cohort Cox-effect differences", x = "Change in initial log hazard ratio",
        y = NULL
      ) +
      theme_classic(base_size = 11, base_family = "serif") +
      theme(
        plot.title = element_text(size = 11),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.grid.major.x = element_line(
          color = "grey92",
          linewidth = 0.3
        ), plot.margin = margin(8, 10, 7, 8)
      )
  } else {
    C <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No discordant mutation labels") +
      theme_void()
  }

  (left | C) + plot_layout(widths = c(1.1, 1))
}
