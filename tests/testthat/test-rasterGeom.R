test_that("rasterGeom wraps a geom and produces a rasterGrob at draw time", {
  skip_if_not_installed("png")
  skip_if(identical(getOption("bitmapType"), "Xlib"),
    "Headless graphics device is unusable for rasterGeom")

  df <- msa2DF(c(a = "ACGT", b = "ACTT"), reference = "a")
  p <- ggplot2::ggplot(df, ggplot2::aes(Position, Sequence, fill = Aln)) +
    rasterGeom(ggplot2::geom_tile())
  expect_s3_class(p, "ggplot")

  out <- withr::local_tempfile(fileext = ".png")
  ggplot2::ggsave(out, p, width = 3, height = 1.5, dpi = 72)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})
