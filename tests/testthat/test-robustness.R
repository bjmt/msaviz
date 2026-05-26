test_that("msa2DF errors on duplicate sequence names with a helpful message", {
  bad <- c(s1 = "ACGT", s2 = "ACTT", s1 = "ACGA")
  expect_error(msa2DF(bad), "duplicate sequence names.*s1")
})

test_that("rasterGeom can be constructed under a NULL bitmapType option", {
  skip_if_not_installed("withr")
  # The fix in 1.1.0 lets rasterGeom() default-construct cleanly even
  # when getOption("bitmapType") is NULL (common on headless Linux). The
  # construction path runs the default-argument expression; this exercises
  # that without needing png() to actually render.
  withr::local_options(bitmapType = NULL)
  expect_no_error({
    rasterGeom(ggplot2::geom_tile())
  })
})
