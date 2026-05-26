#' @noRd
make_group_strip <- function(heatmap, groups, group.colours = NULL) {
  row_order <- attr(heatmap, "row.order")
  if (is.null(row_order)) {
    hm_data <- heatmap$data
    if (is.null(hm_data) || !"Sequence" %in% colnames(hm_data)) {
      stop("Cannot infer row order from heatmap. Pass a plot built by msaHeatmap().",
        call. = FALSE)
    }
    row_order <- if (is.factor(hm_data$Sequence)) {
      levels(hm_data$Sequence)
    } else {
      unique(as.character(hm_data$Sequence))
    }
  }
  if (is.null(names(groups))) {
    stop("`groups` must be a named character or factor vector, where the ",
      "names are the sequence names from `heatmap`.", call. = FALSE)
  }
  missing_names <- setdiff(row_order, names(groups))
  if (length(missing_names)) {
    stop("`groups` is missing entries for sequence(s): ",
      paste(head(missing_names, 5), collapse = ", "),
      if (length(missing_names) > 5L) " (and more)" else "",
      call. = FALSE)
  }
  groups_vec <- groups[row_order]
  df <- data.frame(
    Sequence = factor(row_order, levels = rev(row_order)),
    x = 1,
    Group = factor(groups_vec)
  )
  p <- ggplot(df, aes(x = .data[["x"]], y = .data[["Sequence"]],
      fill = .data[["Group"]])) +
    geom_tile() +
    scale_x_continuous(expand = c(0, 0), name = NULL,
      breaks = NULL) +
    # Match msaHeatmap()'s default discrete y-axis expansion so the strip's
    # row heights are identical to the heatmap's row heights when stacked
    # side by side via patchwork.
    scale_y_discrete(name = NULL) +
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  if (!is.null(group.colours)) {
    p <- p + scale_fill_manual(values = group.colours)
  }
  p
}

#' Compose an MSA heatmap with a PID track, dendrogram, and/or group strip
#'
#' Lays out [msaHeatmap()] together with optional companion plots using
#' [patchwork] so x- and y-axes line up cleanly:
#'
#' * `top` — a ggplot to stack above the heatmap (typically [msaPID()]).
#'   Its x-axis must cover the same range as the heatmap.
#' * `left` — a ggplot to place to the left of the heatmap (typically
#'   [msaDendro()]). Its y-coordinates must match the heatmap's row order.
#'   If `left` carries an `"order"` attribute (as `msaDendro()` adds) and
#'   that order does not match `heatmap`, a warning fires.
#' * `groups` — a named character or factor vector keyed by sequence name.
#'   Drawn as a thin coloured strip just to the left of the heatmap.
#'
#' @param heatmap A ggplot returned by [msaHeatmap()].
#' @param top Optional ggplot to stack above the heatmap.
#' @param left Optional ggplot to place to the left of the heatmap.
#' @param groups Optional named character or factor vector of group labels.
#'   Names must match the sequence names in `heatmap`.
#' @param group.colours Named character vector of fill colours, one per
#'   group level. When `NULL`, ggplot's default discrete palette is used.
#' @param heights Numeric vector of row heights. Defaults to `c(1, 4)`
#'   when `top` is supplied, else `1`.
#' @param widths Numeric vector of column widths. Defaults are chosen to
#'   give the heatmap roughly 4× the width of the dendrogram and 80× the
#'   width of the group strip.
#'
#' @return A patchwork object (which prints and saves like a ggplot).
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz")
#' aln <- seqinr::read.alignment(aln_path, format = "fasta")
#' alnDF <- msa2DF(aln)
#'
#' d <- msaDendro(aln)
#' p <- msaHeatmap(alnDF, row.order = attr(d, "order"))
#' composeMSA(
#'   heatmap = p,
#'   top  = msaPID(alnDF),
#'   left = d
#' )
#' }
#'
#' @export
composeMSA <- function(heatmap, top = NULL, left = NULL,
  groups = NULL, group.colours = NULL,
  heights = NULL, widths = NULL) {

  ord <- if (!is.null(left)) attr(left, "order") else NULL
  if (!is.null(ord)) {
    hm_row_order <- attr(heatmap, "row.order")
    if (is.null(hm_row_order)) {
      hm_data <- heatmap$data
      hm_row_order <- if (is.factor(hm_data$Sequence)) {
        levels(hm_data$Sequence)
      } else {
        unique(as.character(hm_data$Sequence))
      }
    }
    if (!identical(as.character(hm_row_order), as.character(ord))) {
      warning("`left` carries an `order` attribute that does not match the ",
        "heatmap's row order. For correct vertical alignment, build the ",
        "heatmap with `msaHeatmap(..., row.order = attr(left, \"order\"))`.",
        call. = FALSE)
    }
  }

  strip <- if (!is.null(groups)) make_group_strip(heatmap, groups, group.colours) else NULL

  has_top   <- !is.null(top)
  has_left  <- !is.null(left)
  has_strip <- !is.null(strip)

  if (!has_top && !has_left && !has_strip) return(heatmap)

  plots <- list()
  letters_pool <- LETTERS
  next_idx <- 1L
  use_letter <- function() {
    L <- letters_pool[next_idx]
    next_idx <<- next_idx + 1L
    L
  }
  letter_top <- letter_left <- letter_strip <- letter_hm <- NA_character_

  if (has_top)   { plots <- c(plots, list(top));     letter_top   <- use_letter() }
  if (has_left)  { plots <- c(plots, list(left));    letter_left  <- use_letter() }
  if (has_strip) { plots <- c(plots, list(strip));   letter_strip <- use_letter() }
  plots <- c(plots, list(heatmap));                  letter_hm    <- use_letter()

  cols <- character()
  if (has_left)  cols <- c(cols, letter_left)
  if (has_strip) cols <- c(cols, letter_strip)
  cols <- c(cols, letter_hm)

  design_rows <- character()
  if (has_top) {
    top_row <- ifelse(cols == letter_hm, letter_top, "#")
    design_rows <- c(design_rows, paste(top_row, collapse = ""))
  }
  design_rows <- c(design_rows, paste(cols, collapse = ""))
  design <- paste(design_rows, collapse = "\n")

  if (is.null(heights)) heights <- if (has_top) c(1, 4) else 1
  if (is.null(widths)) {
    w <- numeric()
    if (has_left)  w <- c(w, 1)
    if (has_strip) w <- c(w, 0.1)
    w <- c(w, 4)
    widths <- w
  }

  patchwork::wrap_plots(plots, design = design,
    heights = heights, widths = widths)
}
