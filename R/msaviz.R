#' msaviz: Multiple sequence alignment visualization
#'
#' Provides a handful of useful functions for simplifying visualizations of
#' multiple sequence alignments using heatmaps and lollipop charts.
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
#' @importFrom utils object.size
#' @importFrom tibble as_tibble
"_PACKAGE"
