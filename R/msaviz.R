#' msaviz: Multiple sequence alignment visualization
#'
#' A small set of functions for turning a multiple sequence alignment into a
#' tidy data frame and then visualizing it, either as a ggplot2 heatmap or as
#' a per-sequence lollipop chart, with a few companion plots that compose
#' cleanly alongside it.
#'
#' @keywords internal
#' @aliases msaviz-package
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom png readPNG
#' @importFrom grid convertUnit unit viewport convertWidth convertHeight
#' @importFrom grid pushViewport grid.draw popViewport rasterGrob
#' @importFrom grDevices colorRampPalette palette.colors dev.off rgb
#' @importFrom grDevices dev.cur dev.set png
#' @importFrom methods is
#' @importFrom utils object.size head
#' @importFrom tibble as_tibble
"_PACKAGE"
