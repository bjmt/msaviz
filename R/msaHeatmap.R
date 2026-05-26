#' Save a multiple sequence alignment heatmap
#'
#' Thin wrapper around [ggplot2::ggsave()] that auto-derives a sensible height
#' (scaling with the number of sequences) and width (at least twice the height
#' or 100 mm, whichever is larger) when none are supplied. Pass `height` and
#' `width` explicitly to override.
#'
#' @param msaHeatmap A ggplot object as produced by [msaHeatmap()].
#' @param file Output file path. Format is inferred from the extension by
#'   [ggplot2::ggsave()].
#' @param height,width Numeric image dimensions. When `NULL`, computed from
#'   the number of sequences and the axis text size.
#' @param unit Length unit used for `height` and `width`. Defaults to `"mm"`.
#' @param ... Additional arguments forwarded to [ggplot2::ggsave()].
#'
#' @return Invisibly returns the path of the saved file.
#'
#' @examples
#' \dontrun{
#' aln <- c(seq1 = "ACGTACGT", seq2 = "ACGTTCGT", seq3 = "ACGAACGT")
#' p <- msaHeatmap(msa2DF(aln, reference = "seq1"))
#' saveHeatmap(p, tempfile(fileext = ".pdf"))
#' }
#'
#' @export
saveHeatmap <- function(msaHeatmap, file, height = NULL, width = NULL,
  unit = "mm", ...) {

  if (is.null(height)) {
    nseqs <- levels(msaHeatmap$data$Sequence)
    if (is.null(nseqs)) {
      nseqs <- length(unique(as.character(msaHeatmap$data$Sequence)))
    } else {
      nseqs <- length(nseqs)
    }
    tsize <- msaHeatmap$theme$axis.text.y$size
    if (is.null(tsize)) {
      tsize <- 8.8
    }
    extra <- 6
    if (!is.null(msaHeatmap$theme$legend.position) &&
      msaHeatmap$theme$legend.position == "none") {
      extra <- extra - 3
    } else if (!is.null(msaHeatmap$guides$nrow) && msaHeatmap$guides$nrow > 1) {
      extra <- extra + 2 * (msaHeatmap$guides$nrow - 1)
    }
    height <- ((nseqs + extra) * tsize) / ggplot2::.pt
    height <- convertUnit(unit(height, "mm"), unit, valueOnly = TRUE)
  }

  if (is.null(width)) {
    height2 <- convertUnit(unit(height, unit), "mm", valueOnly = TRUE)
    width <- max(c(2 * height2, 100))
    width <- convertUnit(unit(width, "mm"), unit, valueOnly = TRUE)
  }

  message("Saving ", round(width, 3), " x ", round(height, 3), " ", unit, " image")

  ggsave(file, msaHeatmap, height = height, width = width, units = unit, ...)

  invisible(file)
}

#' Draw a multiple sequence alignment as a heatmap
#'
#' Builds a ggplot heatmap from a long-format alignment data frame (as
#' produced by [msa2DF()]). Each tile is one (sequence, position) cell. Colour
#' can encode the per-letter identity of each cell (`column = "Letter"`) or
#' the relationship to a reference sequence (`column = "Aln"`, default).
#'
#' For large alignments, `raster = TRUE` (the default) rasterizes the tile
#' layer at the supplied DPI, producing a much smaller and faster output than
#' a vector layer of millions of rectangles.
#'
#' @param alnDF A tibble produced by [msa2DF()]. Must contain at minimum the
#'   `Sequence` and `Position` columns plus whichever of `Aln` / `Letter` is
#'   chosen by `column`.
#' @param column Which column of `alnDF` to map to the fill aesthetic. One of
#'   `"Aln"` (default) or `"Letter"`.
#' @param gap.colour Fill colour for `NA` (gap) tiles. Defaults to `NA`,
#'   leaving them transparent.
#' @param aln.colours Named character vector of fill colours, used when
#'   `column = "Aln"`. Defaults to `c(Alt = "red", Ref = "black")`.
#' @param letter.colours Named character vector of fill colours, used when
#'   `column = "Letter"`. When `NULL` ggplot's default discrete palette is
#'   used.
#' @param row.order Character vector giving the sequence display order, top to
#'   bottom. Defaults to the factor levels in `alnDF$Sequence`.
#' @param names.pos Side of the plot where sequence names are drawn.
#' @param legend.pos Where to place the fill legend, or `"none"` to suppress
#'   it.
#' @param legend.nrow Number of rows in the legend.
#' @param groups,group.box.width,group.colours,group.labels,group.pos
#'   Reserved for upcoming row-group annotation support. Currently accepted
#'   for forward compatibility but have no effect on the returned plot.
#' @param raster Logical. Rasterize the tile layer via [rasterGeom()] to
#'   shrink vector output and speed up rendering. Defaults to `TRUE`.
#' @param raster.dpi DPI for the rasterized layer.
#' @param raster.type Graphics device type passed through to [grDevices::png()].
#' @param text.size,text.size.x,text.size.y,text.size.legend Sizes (in pt) for
#'   the various text elements. The three specific ones default to the value
#'   of `text.size`.
#' @param trim.names Logical. Trim long sequence names with an ellipsis.
#' @param trim.names.nchar Max number of characters in a trimmed name.
#' @param x.breaks,x.labels,y.breaks,y.labels Override the default x/y scale
#'   breaks and labels.
#' @param hide.x.labels,hide.y.labels Hide axis text and ticks on the
#'   respective axis.
#'
#' @return A ggplot object.
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACGAACGT"
#' )
#' alnDF <- msa2DF(aln, reference = "seq1")
#' msaHeatmap(alnDF)
#'
#' @export
msaHeatmap <- function(alnDF, column = c("Aln", "Letter"), gap.colour = NA,
  aln.colours = c(Alt = "red", Ref = "black"), letter.colours = NULL,
  row.order = levels(alnDF$Sequence), names.pos = c("right", "left"),
  legend.pos = c("top", "bottom", "left", "right", "none"),
  legend.nrow = 1,
  groups = NULL, group.box.width = grid::unit(2, "mm"),
  group.colours = colorRampPalette(palette.colors())(length(groups)),
  group.labels = if (is.factor(groups)) levels(groups) else unique(groups),
  group.pos = c("left", "right"),
  raster = TRUE, raster.dpi = 150, raster.type = getOption("bitmapType"),
  text.size = 6, text.size.x = text.size, text.size.y = text.size,
  text.size.legend = text.size,
  trim.names = TRUE, trim.names.nchar = 25,
  x.breaks = NULL, x.labels = NULL, y.breaks = NULL, y.labels = NULL,
  hide.x.labels = FALSE, hide.y.labels = FALSE
) {

  names.pos <- match.arg(names.pos)
  legend.pos <- match.arg(legend.pos)
  column <- match.arg(column)
  group.pos <- match.arg(group.pos)

  if (any(!c("Sequence", "Position") %in% colnames(alnDF))) {
    stop("alnDF must contain columns Sequence and Position")
  }

  if (column == "Letter" && !"Letter" %in% colnames(alnDF)) {
    stop("alnDF does not contain Letter column")
  }

  if (column == "Aln" && !"Aln" %in% colnames(alnDF)) {
    stop("alnDF does not contain Aln column")
  }

  aln.size <- attributes(alnDF)$aln.size
  if (is.null(aln.size)) aln.size <- max(alnDF$Position)

  if (is.null(x.labels)) x.labels <- waiver()
  if (is.null(y.breaks)) y.breaks <- waiver()
  if (is.null(y.labels)) y.labels <- waiver()
  if (is.null(x.breaks)) {
    x.breaks <- scale_x_continuous(limits = c(0.5, aln.size + 0.5))$get_breaks()
    if (is.na(x.breaks[1])) {
      x.breaks <- c(1, x.breaks[-1])
    }
  }

  if (is.null(row.order)) {
    row.order <- levels(alnDF$Sequence)
  }
  if (is.null(row.order)) {
    row.order <- unique(as.character(alnDF$Sequence))
  }

  tileGeom <- geom_tile(colour = NA)
  if (raster) {
    tileGeom <- rasterGeom(tileGeom, png.dpi = raster.dpi, png.bg = gap.colour,
      png.type = raster.type)
  }

  if (trim.names) {
    if (is(y.labels, "waiver") || is.null(y.labels)) {
      y.labels <- rev(row.order)
    }
    tofix <- nchar(y.labels) > trim.names.nchar
    y.labels[tofix] <- paste0(strtrim(y.labels[tofix], max(c(1, trim.names.nchar - 3))), "...")
  }

  p <- ggplot(alnDF, aes(.data[["Position"]],
      factor(as.character(.data$Sequence), levels = rev(row.order)),
      fill = .data[[column]])) +
    tileGeom +
    scale_y_discrete(position = names.pos, breaks = y.breaks, labels = y.labels,
      name = NULL) +
    scale_x_continuous(limits = c(0.5, aln.size + 0.5), expand = c(0, 0),
      breaks = x.breaks, labels = x.labels, name = NULL) +
    xlab(NULL) +
    ylab(NULL) +
    theme(
      legend.text = element_text(size = text.size.legend, colour = "black"),
      legend.position = legend.pos,
      legend.margin = margin(),
      plot.margin = margin(),
      legend.key.size = unit(2, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.x = element_text(size = text.size.x, colour = "black"),
      axis.text.y = element_text(size = text.size.y, colour = "black")
    )

  if (column == "Aln" && !is.null(aln.colours)) {
    p <- p +
      scale_fill_manual(values = aln.colours, na.value = gap.colour, name = NULL)
  } else if (column == "Letter" && !is.null(letter.colours)) {
    p <- p +
      scale_fill_manual(values = letter.colours, na.value = gap.colour, name = NULL)
  } else {
    p <- p +
      scale_fill_discrete(na.value = gap.colour, name = NULL)
  }

  p <- p + guides(fill = guide_legend(nrow = legend.nrow, byrow = TRUE))

  if (hide.x.labels) {
    p <- p +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.length.x = unit(0, "pt"))
  }
  if (hide.y.labels) {
    p <- p +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"))
  }

  p

}
