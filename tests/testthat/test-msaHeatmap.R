make_df <- function() {
  msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACG-ACGT"
  ), reference = "seq1")
}

test_that("msaHeatmap returns a ggplot with expected layers and scales", {
  p <- msaHeatmap(make_df(), raster = FALSE)
  expect_s3_class(p, "ggplot")
  # Exactly one tile layer
  expect_length(p$layers, 1)
  # Both axis scales are present
  scale_classes <- vapply(p$scales$scales, function(s) class(s)[1], character(1))
  expect_true(any(grepl("ScaleContinuous", scale_classes)))
  expect_true(any(grepl("ScaleDiscrete", scale_classes)))
})

test_that("msaHeatmap validates required columns", {
  df <- make_df()
  expect_error(msaHeatmap(df[, c("Sequence", "Aln")]),
    "must contain columns Sequence and Position")
  expect_error(msaHeatmap(df[, c("Sequence", "Position", "Aln")], column = "Letter"),
    "Letter column")
})

test_that("msaHeatmap can be built and printed without warnings", {
  expect_no_warning({
    p <- msaHeatmap(make_df(), raster = FALSE)
    ggplot2::ggplot_build(p)
  })
})

test_that("snapshot of a small heatmap stays stable", {
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    "small-aln-heatmap",
    msaHeatmap(make_df(), raster = FALSE)
  )
})
