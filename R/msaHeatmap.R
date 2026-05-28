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
    # Read the legend hints stashed by msaHeatmap() rather than poking into
    # ggplot2's internals, whose shape varies between versions.
    legend_pos <- attr(msaHeatmap, "legend.pos")
    legend_nrow <- attr(msaHeatmap, "legend.nrow")
    extra <- 6
    if (identical(legend_pos, "none")) {
      extra <- extra - 3
    } else if (is.numeric(legend_nrow) && length(legend_nrow) == 1L &&
      legend_nrow > 1) {
      extra <- extra + 2 * (legend_nrow - 1)
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
#' @param x.labels.rotate,y.labels.rotate Rotation angle in degrees applied
#'   to the axis text on each axis. Defaults to `0` (no rotation). Useful
#'   for fitting every position label on long alignments
#'   (`x.labels.rotate = 90` plus a denser `x.breaks =`).
#' @param emphasize Character vector of values to emphasize by drawing the
#'   matching cells as enlarged tiles overlaid on the base layer. Defaults
#'   to `"Alt"` (the natural "make SNPs pop" case). Set to `NULL` or
#'   `character(0)` to disable. The overlay only renders when at least one
#'   of `emphasize.size.x` or `emphasize.size.y` is greater than `1`.
#' @param emphasize.by Name of the column in `alnDF` to match `emphasize`
#'   against. Defaults to `"Aln"` so the SNP case works regardless of
#'   `column`. Set to `"Letter"` to enlarge specific residues.
#' @param emphasize.size Master expansion factor for emphasized cells.
#'   Defaults to `1` (no overlay).
#' @param emphasize.size.x,emphasize.size.y Per-axis expansion factors.
#'   Default to `emphasize.size`. Values greater than `1` make emphasized
#'   tiles overlap their neighbours on the corresponding axis. Cells at
#'   the edge of the alignment (position `1`, position `aln.size`, first
#'   or last sequence row) are clipped to the panel so they never extend
#'   past the plot limits.
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
  raster = TRUE, raster.dpi = 150,
  raster.type = getOption("bitmapType"),
  text.size = 6, text.size.x = text.size, text.size.y = text.size,
  text.size.legend = text.size,
  trim.names = TRUE, trim.names.nchar = 25,
  x.breaks = NULL, x.labels = NULL, y.breaks = NULL, y.labels = NULL,
  hide.x.labels = FALSE, hide.y.labels = FALSE,
  x.labels.rotate = 0, y.labels.rotate = 0,
  emphasize = "Alt", emphasize.by = "Aln",
  emphasize.size = 1,
  emphasize.size.x = emphasize.size,
  emphasize.size.y = emphasize.size
) {

  names.pos <- match.arg(names.pos)
  legend.pos <- match.arg(legend.pos)
  column <- match.arg(column)

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
    x.breaks <- integer_breaks(c(1, aln.size))
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

  alnDF$Sequence <- factor(as.character(alnDF$Sequence), levels = rev(row.order))

  emphGeom <- NULL
  if ((emphasize.size.x > 1 || emphasize.size.y > 1) &&
      length(emphasize) > 0) {
    if (!emphasize.by %in% colnames(alnDF)) {
      stop("emphasize.by = '", emphasize.by, "' is not a column in alnDF",
        call. = FALSE)
    }
    emph_rows <- as.character(alnDF[[emphasize.by]]) %in% emphasize
    emph_data <- alnDF[emph_rows, , drop = FALSE]
    if (nrow(emph_data) > 0) {
      # Cap the enlarged tile at the panel limits so it does not bleed past
      # the first/last position or first/last sequence row.
      nrows <- nlevels(alnDF$Sequence)
      y_idx <- as.integer(emph_data$Sequence)
      emph_data$emph_xmin <- pmax(emph_data$Position - emphasize.size.x / 2, 0.5)
      emph_data$emph_xmax <- pmin(emph_data$Position + emphasize.size.x / 2,
        aln.size + 0.5)
      emph_data$emph_ymin <- pmax(y_idx - emphasize.size.y / 2, 0.5)
      emph_data$emph_ymax <- pmin(y_idx + emphasize.size.y / 2, nrows + 0.5)
      emphGeom <- geom_rect(data = emph_data, colour = NA,
        aes(xmin = .data[["emph_xmin"]], xmax = .data[["emph_xmax"]],
            ymin = .data[["emph_ymin"]], ymax = .data[["emph_ymax"]]))
      if (raster) {
        emphGeom <- rasterGeom(emphGeom, png.dpi = raster.dpi,
          png.bg = gap.colour, png.type = raster.type)
      }
    }
  }

  p <- ggplot(alnDF, aes(.data[["Position"]], .data[["Sequence"]],
      fill = .data[[column]])) +
    tileGeom
  if (!is.null(emphGeom)) p <- p + emphGeom
  p <- p +
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

  if (x.labels.rotate != 0) {
    x_hjust <- if (x.labels.rotate %% 180 == 0) 0.5 else 1
    x_vjust <- if (abs(x.labels.rotate) == 90) 0.5 else 1
    p <- p + theme(
      axis.text.x = element_text(angle = x.labels.rotate,
        hjust = x_hjust, vjust = x_vjust,
        size = text.size.x, colour = "black")
    )
  }
  if (y.labels.rotate != 0) {
    y_hjust <- if (y.labels.rotate %% 180 == 0) 1 else 0.5
    y_vjust <- 0.5
    p <- p + theme(
      axis.text.y = element_text(angle = y.labels.rotate,
        hjust = y_hjust, vjust = y_vjust,
        size = text.size.y, colour = "black")
    )
  }

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

  attr(p, "row.order") <- row.order
  attr(p, "legend.pos") <- legend.pos
  attr(p, "legend.nrow") <- legend.nrow
  p

}
