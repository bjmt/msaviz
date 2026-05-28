#' @noRd
.lollipop_ref_letters <- function(alnDF, aln_size) {
  ref_name <- attr(alnDF, "reference")
  if (is.null(ref_name)) return(NULL)
  ref_rows <- alnDF[as.character(alnDF$Sequence) == ref_name, , drop = FALSE]
  if (nrow(ref_rows) == 0L) return(NULL)
  out <- rep(NA_character_, aln_size)
  out[ref_rows$Position] <- as.character(ref_rows$Letter)
  out
}

#' Per-sequence lollipop chart of alignment differences
#'
#' Builds a ggplot lollipop chart from a long-format alignment data frame
#' (as produced by [msa2DF()]). One row of the discrete y-axis is allocated
#' per sequence, and each highlighted cell within a row gets a vertical
#' stem topped with a circular head and an optional `RefLetter Position
#' AltLetter` label (e.g. `A123T`). A per-position horizontal baseline
#' runs along each row at `y = 0` so that gaps in the alignment show up
#' as visible breaks in the baseline.
#'
#' This is the function that drew Supplementary Figure 20 of Hodgins et
#' al. (2023), *Nature Communications* 14:5475, generalised so it can be
#' used as a drop-in substitute for [msaHeatmap()] inside [composeMSA()]
#' (dendrogram, row + column group strips, and a PID track all line up
#' with lollipop rows).
#'
#' @param alnDF A tibble produced by [msa2DF()]. Must contain `Sequence`,
#'   `Position`, `Letter`, and `Aln`.
#' @param highlight Character vector of values in `highlight.by` to draw
#'   lollipops for. Defaults to `"Alt"` (every non-reference cell).
#' @param highlight.by Column of `alnDF` to match `highlight` against.
#'   Defaults to `"Aln"`, so the natural "lollipop every SNP" case works
#'   regardless of `column =` choices upstream. Set to `"Letter"` to
#'   lollipop specific residues.
#' @param y Optional column name in `alnDF` giving a numeric stem height
#'   per cell. When `NULL` (default), every highlighted cell gets a stem
#'   of height 1. When supplied, heights are normalised against
#'   `max(alnDF[[y]], na.rm = TRUE)` so the tallest stem still fits
#'   within `stem.height`.
#' @param labels Logical. Draw `RefPosAlt` labels above each head. When
#'   the reference sequence is not present in `alnDF` (e.g. the
#'   `msa2DF(..., keep.consensus = FALSE)` default consensus path),
#'   labels fall back to position-only (`"123T"`) and a one-time message
#'   points at `keep.consensus = TRUE`.
#' @param baseline Logical. Draw per-position horizontal baseline
#'   segments. Gaps in the alignment show as visible breaks.
#' @param drop.empty Logical. Drop sequences with no highlighted cells
#'   before plotting.
#' @param row.order Character vector giving the sequence display order,
#'   top to bottom. Names not present in `alnDF$Sequence` are kept as
#'   empty rows, which reproduces the "spacer" trick from the paper.
#' @param names.pos Side of the plot where sequence names are drawn.
#' @param trim.names Logical. Trim long sequence names with an ellipsis.
#' @param trim.names.nchar Max number of characters in a trimmed name.
#' @param stem.height Maximum stem height within a row, in row-height
#'   units. Default `0.7` leaves 30% of each row free for the label.
#' @param point.size,point.shape,point.fill,point.colour Aesthetics for
#'   the lollipop heads. `point.fill` is the constant fill colour used
#'   when `head.fill.by` is `NULL` (the default).
#' @param head.fill.by Optional column name in `alnDF`. When set, each
#'   head is filled according to that column's value, useful for
#'   colouring SNPs by residue identity (`head.fill.by = "Letter"`) or
#'   by SNP category (`head.fill.by = "Aln"`). Defaults to `NULL`
#'   (constant `point.fill`).
#' @param head.fill.colours Named character vector of fill colours,
#'   used with `head.fill.by`. When `NULL`, ggplot's default discrete
#'   palette is used. The bundled palettes [msa_palette_DNA],
#'   [msa_palette_RNA], and [msa_palette_AA] are well-suited here.
#' @param segment.lwd,segment.colour Aesthetics for the stems.
#' @param baseline.lwd,baseline.colour Aesthetics for the baseline.
#' @param label.size Size of the `geom_text` labels above the heads.
#' @param label.nudge.y Vertical offset of labels above each head, in
#'   units of `stem.height` (so labels track the configured stem
#'   height).
#' @param text.size Size of the sequence-name axis text.
#' @param hide.x.labels Hide the x-axis text and ticks. Default `TRUE` so
#'   the plot reads clean as in the paper.
#' @param x.breaks,x.labels Override the default x-axis breaks and
#'   labels.
#'
#' @return A ggplot object. Carries `row.order`, `legend.pos`,
#'   `legend.nrow`, and `aln.size` attributes so that [composeMSA()] can
#'   stack it with companion plots the same way it stacks [msaHeatmap()]
#'   output.
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACGAACGT"
#' )
#' alnDF <- msa2DF(aln, reference = "seq1", drop.gaps = FALSE,
#'   keep.consensus = TRUE)
#' msaLollipop(alnDF)
#'
#' @export
msaLollipop <- function(alnDF,
  highlight = "Alt", highlight.by = "Aln",
  y = NULL,
  labels = TRUE, baseline = TRUE, drop.empty = TRUE,
  row.order = levels(alnDF$Sequence),
  names.pos = c("right", "left"),
  trim.names = TRUE, trim.names.nchar = 25,
  stem.height = 0.7,
  point.size = 2, point.shape = 21,
  point.fill = "white", point.colour = "black",
  head.fill.by = NULL, head.fill.colours = NULL,
  segment.lwd = 0.5, segment.colour = "black",
  baseline.lwd = 0.5, baseline.colour = "black",
  label.size = 2, label.nudge.y = 0.25,
  text.size = 6,
  hide.x.labels = TRUE,
  x.breaks = NULL, x.labels = NULL
) {

  names.pos <- match.arg(names.pos)

  if (any(!c("Sequence", "Position", "Letter", "Aln") %in% colnames(alnDF))) {
    stop("alnDF must contain columns Sequence, Position, Letter, Aln",
      call. = FALSE)
  }
  if (!highlight.by %in% colnames(alnDF)) {
    stop("highlight.by = '", highlight.by, "' is not a column in alnDF",
      call. = FALSE)
  }
  if (!is.null(y) && !y %in% colnames(alnDF)) {
    stop("y = '", y, "' is not a column in alnDF", call. = FALSE)
  }
  if (!is.null(head.fill.by) && !head.fill.by %in% colnames(alnDF)) {
    stop("head.fill.by = '", head.fill.by, "' is not a column in alnDF",
      call. = FALSE)
  }
  if (!is.numeric(stem.height) || length(stem.height) != 1L ||
      stem.height <= 0 || stem.height >= 1) {
    stop("stem.height must be a single number in (0, 1)", call. = FALSE)
  }

  aln_size <- attr(alnDF, "aln.size")
  if (is.null(aln_size)) aln_size <- max(alnDF$Position)

  if (is.null(row.order)) {
    row.order <- levels(alnDF$Sequence)
  }
  if (is.null(row.order)) {
    row.order <- unique(as.character(alnDF$Sequence))
  }

  hit_rows <- as.character(alnDF[[highlight.by]]) %in% highlight
  hits <- alnDF[hit_rows, , drop = FALSE]

  if (drop.empty) {
    non_empty <- unique(as.character(hits$Sequence))
    row.order <- row.order[row.order %in% non_empty]
    if (length(row.order) == 0L) {
      row.order <- non_empty
    }
  }

  ref_letters <- if (labels) .lollipop_ref_letters(alnDF, aln_size) else NULL
  fallback_msg <- labels && is.null(ref_letters)

  if (nrow(hits) > 0L) {
    hits$y_norm <- if (is.null(y)) {
      rep(1, nrow(hits))
    } else {
      vals <- as.numeric(hits[[y]])
      if (any(!is.na(vals) & vals < 0)) {
        stop("y column '", y, "' must be non-negative", call. = FALSE)
      }
      m <- suppressWarnings(max(vals, na.rm = TRUE))
      if (!is.finite(m) || m <= 0) m <- 1
      vals / m
    }
    if (!is.null(y)) hits <- hits[!is.na(hits$y_norm), , drop = FALSE]
  } else {
    hits$y_norm <- numeric(0)
  }

  if (labels && nrow(hits) > 0L) {
    alt <- as.character(hits$Letter)
    if (is.null(ref_letters)) {
      hits$lab <- paste0(hits$Position, alt)
    } else {
      ref <- ref_letters[hits$Position]
      ref[is.na(ref)] <- "?"
      hits$lab <- paste0(ref, hits$Position, alt)
    }
  }

  base_rows <- alnDF[as.character(alnDF$Aln) != "Gap" &
    !is.na(alnDF$Letter) &
    as.character(alnDF$Sequence) %in% row.order, , drop = FALSE]

  hits <- hits[as.character(hits$Sequence) %in% row.order, , drop = FALSE]

  seq_levels <- rev(row.order)

  hits$Sequence <- factor(as.character(hits$Sequence), levels = seq_levels)
  base_rows$Sequence <- factor(as.character(base_rows$Sequence),
    levels = seq_levels)

  hits$y_int <- as.integer(hits$Sequence)
  base_rows$y_int <- as.integer(base_rows$Sequence)
  hits$y_top <- hits$y_int + stem.height * hits$y_norm
  hits$y_lab <- hits$y_top + label.nudge.y * stem.height

  if (is.null(x.breaks)) x.breaks <- integer_breaks(c(1, aln_size))
  if (is.null(x.labels)) x.labels <- waiver()

  p <- ggplot(hits, aes(x = .data[["Position"]], y = .data[["y_int"]]))

  if (baseline && nrow(base_rows) > 0L) {
    p <- p + geom_segment(
      data = base_rows,
      aes(x = .data[["Position"]] - 0.5, xend = .data[["Position"]] + 0.5,
        y = .data[["y_int"]], yend = .data[["y_int"]]),
      colour = baseline.colour, linewidth = baseline.lwd,
      inherit.aes = FALSE
    )
  }

  if (nrow(hits) > 0L) {
    p <- p +
      geom_segment(
        aes(x = .data[["Position"]], xend = .data[["Position"]],
          y = .data[["y_int"]], yend = .data[["y_top"]]),
        colour = segment.colour, linewidth = segment.lwd,
        inherit.aes = FALSE
      ) +
      (if (is.null(head.fill.by)) {
        geom_point(
          aes(x = .data[["Position"]], y = .data[["y_top"]]),
          shape = point.shape, size = point.size,
          fill = point.fill, colour = point.colour,
          inherit.aes = FALSE
        )
      } else {
        geom_point(
          aes(x = .data[["Position"]], y = .data[["y_top"]],
              fill = .data[[head.fill.by]]),
          shape = point.shape, size = point.size,
          colour = point.colour,
          inherit.aes = FALSE
        )
      })
    if (labels) {
      p <- p + geom_text(
        aes(x = .data[["Position"]], y = .data[["y_lab"]],
          label = .data[["lab"]]),
        size = label.size, vjust = 0, inherit.aes = FALSE
      )
    }
  }

  y_axis_labels <- seq_levels
  if (trim.names) {
    tofix <- nchar(y_axis_labels) > trim.names.nchar
    y_axis_labels[tofix] <- paste0(
      strtrim(y_axis_labels[tofix], max(c(1, trim.names.nchar - 3))), "...")
  }

  p <- p +
    scale_x_continuous(limits = c(0.5, aln_size + 0.5), expand = c(0, 0),
      breaks = x.breaks, labels = x.labels, name = NULL) +
    scale_y_continuous(
      breaks = seq_along(seq_levels),
      labels = y_axis_labels,
      limits = c(0.5, length(seq_levels) + stem.height +
        label.nudge.y * stem.height + 0.5),
      expand = c(0, 0),
      position = names.pos,
      name = NULL
    ) +
    theme(
      legend.position = "none",
      plot.margin = margin(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.ticks.y.right = element_blank(),
      axis.text.x = element_text(size = text.size, colour = "black"),
      axis.text.y = element_text(size = text.size, colour = "black")
    )

  if (!is.null(head.fill.by) && !is.null(head.fill.colours)) {
    p <- p + scale_fill_manual(values = head.fill.colours, name = NULL)
  }

  if (hide.x.labels) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.length.x = unit(0, "pt")
    )
  }

  if (fallback_msg) {
    message("msaLollipop(): reference sequence not present in alnDF; ",
      "labels fall back to position-only. Pass keep.consensus = TRUE ",
      "to msa2DF() (or pick an explicit reference =) for full ",
      "RefPosAlt labels.")
  }

  attr(p, "row.order") <- row.order
  attr(p, "legend.pos") <- "none"
  attr(p, "legend.nrow") <- 1L
  attr(p$data, "aln.size") <- aln_size
  p
}
