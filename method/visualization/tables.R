# LaTeX tables are built from the same mutation summaries used in the figures.
escape_tex <- function(x) {
  x <- gsub("_", "\\\\_", as.character(x), fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}
write_annotation_table <- function(d, file) {
  order <- match(d[["TL-SCP"]], c("Early-onset", "Mid-onset", "Late-onset"))
  d <- d[order(order, d$ProtDescription), ]
  rows <- vapply(seq_len(nrow(d)), function(i) {
    x <- d[i, ]
    yta <- if (x$YTAClass == "Unclassified") {
      "Unknown"
    } else {
      x$YTAClass
    }
    paste0(
      "\\texttt{", sub("^p\\.", "", x$ProtDescription), "} & ", x$target_n,
      " & ", x$source_n, " & ", x[["TL-SCP"]], " & ", x$SCP, " & ", yta, " \\\\"
    )
  }, character(1))
  # Separate onset groups with the same small gap used in the supplement.
  rows <- unlist(lapply(seq_along(rows), function(i) {
    c(rows[i], if (i %in% cumsum(table(factor(d[["TL-SCP"]], levels = c(
      "Early-onset",
      "Mid-onset", "Late-onset"
    ))))[1:2]) {
      "\\addlinespace[2pt]"
    })
  }))
  writeLines(
    c(
      "\\begingroup", "\\fontsize{9.5}{11.4}\\selectfont", "\\renewcommand{\\arraystretch}{1}",
      "\\setlength{\\tabcolsep}{6pt}", "\\begin{tabular*}{0.98\\textwidth}{@{\\extracolsep{\\fill}}lrrlll@{}}",
      "\\toprule", "& \\multicolumn{2}{c}{Carriers} & \\multicolumn{3}{c}{Annotation} \\\\",
      "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}", "Mutation & Target & Source & TL-SCP & SCP & YTA \\\\",
      "\\midrule", rows, "\\midrule", paste0(
        "Total & ", sum(d$target_n), " & ",
        sum(d$source_n), " & & & \\\\"
      ), "\\bottomrule", "\\end{tabular*}", "\\endgroup"
    ),
    file
  )
}
write_composition_table <- function(d, file) {
  lev <- c("0", "A", "B", "C", "D", "Unclassified")
  counts <- function(a, col, g) {
    c(sum(a$YTAClass == g), sum(a[[col]][a$YTAClass ==
      g]))
  }
  rows <- vapply(lev, function(g) {
    paste(c(g, counts(d$database, "target_n", g), counts(
      d$database,
      "source_n", g
    ), counts(d$mda, "target_n", g)), collapse = " & ")
  }, character(1))
  totals <- c(
    "Total", nrow(d$database), sum(d$database$target_n), nrow(d$database),
    sum(d$database$source_n), nrow(d$mda), sum(d$mda$target_n)
  )
  writeLines(c(
    "\\begingroup", "\\small", "\\setlength{\\tabcolsep}{4pt}", "\\begin{tabular}{@{}lrrrrrr@{}}",
    "\\toprule", "& \\multicolumn{2}{c}{\\shortstack{\\emph{TP53} Database\\\\LFS/Chompret}} & \\multicolumn{2}{c}{\\shortstack{\\emph{TP53} Database\\\\non-LFS/non-Chompret}} & \\multicolumn{2}{c}{\\shortstack{MD Anderson\\\\prospective cohort}} \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}", "YTA class & Mutations & Patients & Mutations & Patients & Mutations & Patients \\\\",
    "\\midrule", paste0(rows, " \\\\"), "\\midrule", paste0(
      paste(totals, collapse = " & "),
      " \\\\"
    ), "\\bottomrule", "\\end{tabular}", "\\endgroup"
  ), file)
}
write_concordance_table <- function(d, file) {
  lev <- c("Early-onset", "Mid-onset", "Late-onset")
  cell <- function(a, b) {
    z <- d[d[["TL-SCP_Database"]] == a & d[["TL-SCP_MDA"]] == b, ]
    if (!nrow(z)) {
      return("--")
    }
    z <- z[order(z$ProtDescription), ]
    text <- paste0("\\texttt{", sub("^p\\.", "", z$ProtDescription), "}")
    if (a != b) {
      text <- paste0(text, " (", z$target_n_Database, ", ", z$target_n_MDA, ")")
    }
    if (a == b) {
      text <- vapply(split(text, ceiling(seq_along(text) / 2)), paste,
        collapse = ", ",
        FUN.VALUE = character(1)
      )
    }
    if (length(text) == 1) {
      return(text)
    }
    paste0("\\begin{tabular}[c]{@{}c@{}}", paste(text, collapse = "\\\\"), "\\end{tabular}")
  }
  rows <- vapply(lev, function(a) {
    paste(c(a, vapply(
      lev, function(b) cell(a, b),
      character(1)
    )), collapse = " & ")
  }, character(1))
  rows <- c(
    paste0(rows[1], " \\\\"), "\\addlinespace[0.4em]", paste0(rows[2], " \\\\"),
    "\\addlinespace[0.4em]", paste0(rows[3], " \\\\")
  )
  writeLines(c(
    "\\begingroup", "\\small", "\\renewcommand{\\arraystretch}{1.18}",
    "\\begin{tabular*}{0.95\\textwidth}{@{\\extracolsep{\\fill}}lccc@{}}", "\\toprule",
    "& \\multicolumn{3}{c}{MD Anderson TL-SCP label} \\\\", "\\cmidrule(lr){2-4}",
    "Primary TL-SCP label & Early-onset & Mid-onset & Late-onset \\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular*}", "\\endgroup"
  ), file)
}
