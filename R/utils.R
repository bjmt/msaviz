#' Integer-valued axis breaks for alignment heatmaps
#'
#' Returns a vector of integer break points suitable for the position axis
#' of an alignment heatmap. The first break is always `1` (the first
#' position of an alignment, not `0`), and remaining breaks come from
#' [base::pretty()] coerced to integers.
#'
#' @param limits Numeric length-2 vector giving the axis range. Only the
#'   upper bound is used; the lower bound is treated as `1`.
#' @param n Target number of breaks.
#'
#' @return An integer vector.
#'
#' @examples
#' integer_breaks(c(1, 9))
#' integer_breaks(c(1, 100))
#' integer_breaks(c(1, 2596))
#'
#' @export
integer_breaks <- function(limits, n = 6) {
  m <- as.integer(max(limits, na.rm = TRUE))
  if (m <= 1L) return(1L)
  if (m <= 3L) return(seq.int(1L, m))
  if (m <= 9L) return(as.integer(c(1L, ceiling(m / 2), m)))
  raw <- pretty(c(1, m), n = n)
  raw <- raw[raw >= 0 & raw <= m]
  if (length(raw) == 0L) return(c(1L, m))
  if (raw[1L] == 0) {
    raw[1L] <- 1
  } else if (raw[1L] != 1) {
    raw <- c(1, raw)
  }
  as.integer(unique(round(raw)))
}

#' Rasterize a ggplot geom layer
#'
#' Wraps a ggplot2 geom so that, at draw time, the geom is rendered into an
#' in-memory PNG via [grDevices::png()] and the result is laid back onto the
#' plot as a [grid::rasterGrob()]. This produces compact, fast-rendering
#' output for plots with many tiles or rectangles (e.g. large
#' [msaHeatmap()]s), at the cost of resolution. Antialiasing and bilinear
#' interpolation are both off by default, giving crisp blocky output well
#' suited to alignment heatmaps; turn them on if you need smooth gradients.
#'
#' @param geomObj A geom layer to rasterize, e.g. the return value of
#'   [ggplot2::geom_tile()].
#' @param png.dpi Resolution in dots per inch for the rasterized layer.
#' @param png.bg Background colour painted under the rasterized geom. Defaults
#'   to `NA` (transparent).
#' @param png.type Graphics device type passed through to [grDevices::png()].
#'   Defaults to `getOption("bitmapType")`. When that is `NULL` (e.g. on a
#'   headless Linux server), [grDevices::png()]'s own built-in default is
#'   used.
#' @param png.antialias Antialiasing mode passed through to [grDevices::png()].
#' @param raster.interpolate Logical. Bilinear interpolation when the raster
#'   is drawn back onto the panel.
#'
#' @return A ggproto layer object that can be added to a ggplot the same way
#'   the original `geomObj` would be.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + rasterGeom(geom_point())
#' }
#'
#' @export
rasterGeom <- function(geomObj, png.dpi = 150, png.bg = NA,
  png.type = getOption("bitmapType"), png.antialias = "none",
  raster.interpolate = FALSE) {

  ggproto(NULL, geomObj, geom = ggproto(NULL, geomObj$geom,
    draw_panel = function(...) {
      g <- geomObj$geom$draw_panel(...)
      vp <- if (is.null(g$vp)) viewport() else g$vp
      npc1 <- unit(1, "npc")
      width <- convertWidth(npc1, "in", valueOnly = TRUE)
      height <- convertHeight(npc1, "in", valueOnly = TRUE)
      cdev <- dev.cur()
      on.exit(dev.set(cdev), add = TRUE)
      f <- tempfile(fileext = ".png")
      on.exit(unlink(f), add = TRUE)
      png_args <- list(filename = f, width = width, height = height,
        units = "in", res = png.dpi, bg = png.bg, antialias = png.antialias)
      if (!is.null(png.type) && nzchar(png.type)) {
        png_args$type <- png.type
      }
      do.call(png, png_args)
      pushViewport(vp)
      grid.draw(g)
      popViewport()
      dev.off()
      img <- readPNG(f)
      img <- matrix(rgb(red = as.vector(img[, , 1]),
        green = as.vector(img[, , 2]), blue = as.vector(img[, , 3]),
        alpha = as.vector(img[, , 4])), dim(img)[1], dim(img)[2])
      rasterGrob(img, height = npc1, width = npc1,
        interpolate = raster.interpolate)
    }
  ))

}
