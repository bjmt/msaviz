#' Per-position percent identity track
#'
#' Builds a ggplot showing the fraction of non-gap cells at each alignment
#' position that match the reference sequence (`Aln == "Ref"`). The x-axis
#' is configured to line up exactly with [msaHeatmap()] so the two plots
#' can be stacked via [composeMSA()] (or any patchwork composition) with
#' shared coordinates.
#'
#' `y.axis.pos` defaults to `"right"` to match [msaHeatmap()]'s default
#' `names.pos = "right"`. This keeps the panel's left edge flush with the
#' heatmap's panel's left edge, which is what `composeMSA()` needs for
#' clean panel alignment when a `left` companion (e.g. [msaDendro()]) is
#' also present.
#'
#' @param alnDF A tibble produced by [msa2DF()].
#' @param aln.size Total alignment width in positions. If `NULL` (default)
#'   it is taken from `attr(alnDF, "aln.size")` and falls back to
#'   `max(alnDF$Position)`.
#' @param style One of `"bar"` (default, `geom_col`) or `"line"`
#'   (`geom_step`).
#' @param fill,line.colour Aesthetics for the bar fill or step line.
#' @param text.size Size in pt of axis text.
#' @param hide.x.labels Hide the x-axis text + ticks. Default `TRUE` so
#'   the track plays nicely as a top strip over a heatmap.
#' @param y.axis.pos Side of the panel to draw the y-axis on. Defaults to
#'   `"right"`, matching `msaHeatmap()`'s default `names.pos`. Set to
#'   `"left"` if you called `msaHeatmap(names.pos = "left")` or are
#'   showing the PID track standalone.
#'
#' @return A ggplot object.
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACGAACGT",
#'   seq4 = "ACGTACGT"
#' )
#' msaPID(msa2DF(aln, reference = "seq1"))
#'
#' @export
msaPID <- function(alnDF, aln.size = NULL,
  style = c("bar", "line"),
  fill = "grey40", line.colour = "black",
  text.size = 6, hide.x.labels = TRUE,
  y.axis.pos = c("right", "left")) {

  style <- match.arg(style)
  y.axis.pos <- match.arg(y.axis.pos)

  if (is.null(aln.size)) {
    aln.size <- attributes(alnDF)$aln.size
    if (is.null(aln.size)) aln.size <- max(alnDF$Position)
  }

  stats <- posStats(alnDF)

  p <- ggplot(stats, aes(x = .data[["Position"]], y = .data[["pctIdentity"]]))
  if (style == "bar") {
    p <- p + geom_col(fill = fill, width = 1, na.rm = TRUE)
  } else {
    p <- p + geom_step(colour = line.colour, direction = "mid", na.rm = TRUE)
  }

  p <- p +
    scale_x_continuous(limits = c(0.5, aln.size + 0.5), expand = c(0, 0),
      breaks = integer_breaks(c(1, aln.size)), name = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
      breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
      name = "% identity",
      position = y.axis.pos) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(),
      axis.text.x = element_text(size = text.size, colour = "black"),
      axis.text.y = element_text(size = text.size, colour = "black"),
      axis.title.y = element_text(size = text.size, colour = "black"),
      axis.ticks.y = element_line()
    )

  if (hide.x.labels) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt")
    )
  }

  p
}
