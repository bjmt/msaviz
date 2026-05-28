test_that("msa2DF errors on duplicate sequence names with a helpful message", {
  bad <- c(s1 = "ACGT", s2 = "ACTT", s1 = "ACGA")
  expect_error(msa2DF(bad), "duplicate sequence names.*s1")
})

test_that("msa2DF rejects empty and zero-width alignments with a helpful message", {
  expect_error(msa2DF(character(0)), "empty")
  expect_error(msa2DF(c(a = "", b = "")), "zero width")
  expect_error(consensusSeq(character(0)), "empty")
  expect_error(consensusSeq(c(a = "", b = "")), "zero width")
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
