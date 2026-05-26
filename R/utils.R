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
      png(f, width = width, height = height, type = png.type,
        units = "in", res = png.dpi, bg = png.bg, antialias = png.antialias)
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
