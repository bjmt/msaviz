#' Hierarchical-clustering dendrogram of an alignment
#'
#' Computes a pairwise distance matrix between sequences with
#' [seqinr::dist.alignment()], runs [stats::hclust()] on it, and renders the
#' resulting dendrogram as a ggplot via [ggdendro::dendro_data()]. Returns a
#' ggplot whose y-coordinates line up with [msaHeatmap()] when the heatmap
#' uses the same row order, so it can be placed to the left of a heatmap
#' via [composeMSA()] (or any patchwork composition).
#'
#' The clustered order of sequences is attached to the returned plot as the
#' `"order"` attribute, ready to pass to `msaHeatmap(row.order = ...)`.
#'
#' Requires the `seqinr` and `ggdendro` packages (declared in `Suggests`).
#'
#' @param aln Alignment, accepted in the same forms as [msa2DF()].
#' @param dist.matrix Passed through as `matrix =` to
#'   [seqinr::dist.alignment()]. One of `"identity"` (default) or
#'   `"similarity"`.
#' @param hclust.method Linkage method passed to [stats::hclust()].
#'   Defaults to `"average"`.
#' @param line.colour Colour of dendrogram segments.
#' @param flip Logical. If `TRUE` (default) the dendrogram is rotated so
#'   the root is on the left and the leaves are on the right, ready to
#'   sit alongside a heatmap. If `FALSE`, the conventional bottom-up
#'   layout is used.
#'
#' @return A ggplot object with an `"order"` attribute giving the leaf
#'   order as a character vector of sequence names.
#'
#' @examples
#' \dontrun{
#' aln_path <- system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz")
#' aln <- seqinr::read.alignment(aln_path, format = "fasta")
#' d <- msaDendro(aln)
#' head(attr(d, "order"))
#' }
#'
#' @export
msaDendro <- function(aln,
  dist.matrix = c("identity", "similarity"),
  hclust.method = "average",
  line.colour = "black",
  flip = TRUE) {

  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("Package 'seqinr' is required for msaDendro(). ",
      "Install with install.packages(\"seqinr\").", call. = FALSE)
  }
  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    stop("Package 'ggdendro' is required for msaDendro(). ",
      "Install with install.packages(\"ggdendro\").", call. = FALSE)
  }
  dist.matrix <- match.arg(dist.matrix)

  if (is(aln, "alignment")) {
    aln_obj <- aln
  } else {
    if (!is.character(aln)) aln <- as.character(aln)
    if (is.null(names(aln))) names(aln) <- as.character(seq_len(length(aln)))
    aln_obj <- seqinr::as.alignment(
      nb = length(aln),
      nam = names(aln),
      seq = unname(aln),
      com = NA
    )
  }

  d <- seqinr::dist.alignment(aln_obj, matrix = dist.matrix)
  hc <- stats::hclust(d, method = hclust.method)
  dd <- ggdendro::dendro_data(hc, type = "rectangle")
  segs <- dd$segments
  leaf_order <- hc$labels[hc$order]
  n <- length(leaf_order)

  if (flip) {
    p <- ggplot(segs,
        aes(x = -.data[["y"]], y = n + 1 - .data[["x"]],
            xend = -.data[["yend"]], yend = n + 1 - .data[["xend"]])) +
      geom_segment(colour = line.colour) +
      scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0))
  } else {
    p <- ggplot(segs,
        aes(x = .data[["x"]], y = .data[["y"]],
            xend = .data[["xend"]], yend = .data[["yend"]])) +
      geom_segment(colour = line.colour) +
      scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
  }

  p <- p +
    theme_void() +
    theme(plot.margin = margin())

  attr(p, "order") <- leaf_order
  p
}
