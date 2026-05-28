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

test_that("x.labels.rotate and y.labels.rotate set the corresponding theme angles", {
  p <- msaHeatmap(make_df(), raster = FALSE,
    x.labels.rotate = 90, y.labels.rotate = 45)
  expect_equal(p$theme$axis.text.x$angle, 90)
  expect_equal(p$theme$axis.text.y$angle, 45)
})

test_that("emphasize.size = 1 adds no overlay; > 1 adds a second layer", {
  p1 <- msaHeatmap(make_df(), raster = FALSE)
  expect_length(p1$layers, 1L)

  p2 <- msaHeatmap(make_df(), raster = FALSE, emphasize.size = 1.5)
  expect_length(p2$layers, 2L)
})

test_that("per-axis emphasize.size.x and emphasize.size.y reach the overlay tile", {
  p <- msaHeatmap(make_df(), raster = FALSE,
    emphasize.size.x = 2, emphasize.size.y = 1.5)
  expect_length(p$layers, 2L)
  # The overlay is now a geom_rect with per-row xmin/xmax/ymin/ymax columns
  # so the size knobs are visible in the layer's data widths/heights.
  ld <- p$layers[[2]]$data
  widths <- ld$emph_xmax - ld$emph_xmin
  heights <- ld$emph_ymax - ld$emph_ymin
  # Interior cells take the full requested width/height; edge cells get
  # clipped to the panel limits (smaller).
  expect_true(any(widths == 2))
  expect_true(any(heights == 1.5))
})

test_that("emphasis overlay does not extend past the plot limits", {
  # Pick an alignment whose Alt cells include positions 1 and aln.size,
  # so the clipping path is exercised on both sides.
  df <- msa2DF(c(seq1 = "ACGTACGT", seq2 = "TCGTACGA", seq3 = "ACGTTCGT"),
    reference = "seq1", drop.gaps = FALSE)
  aln_size <- attr(df, "aln.size")
  nrows <- nlevels(df$Sequence)
  p <- msaHeatmap(df, raster = FALSE, emphasize.size = 3)
  ld <- p$layers[[2]]$data
  expect_true(all(ld$emph_xmin >= 0.5))
  expect_true(all(ld$emph_xmax <= aln_size + 0.5))
  expect_true(all(ld$emph_ymin >= 0.5))
  expect_true(all(ld$emph_ymax <= nrows + 0.5))
})

test_that("emphasize.by switches the overlay filter column", {
  df <- make_df()
  # Filter on Letter instead of Aln. Three sequences * one position with "T" =
  # the overlay data should contain exactly the rows whose Letter is "T".
  p <- msaHeatmap(df, raster = FALSE,
    emphasize = "T", emphasize.by = "Letter", emphasize.size = 1.5)
  expect_equal(
    nrow(p$layers[[2]]$data),
    sum(df$Letter == "T", na.rm = TRUE)
  )
})

test_that("emphasize.by errors clearly on a missing column", {
  expect_error(
    msaHeatmap(make_df(), emphasize.by = "NoSuch", emphasize.size = 1.5),
    "is not a column in alnDF"
  )
})

test_that("snapshot of a small heatmap stays stable", {
  # Announce the snapshot before the skip so testthat does not prune it from
  # tests/testthat/_snaps/ when vdiffr is unavailable locally.
  testthat::announce_snapshot_file(name = "small-aln-heatmap.svg")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    "small-aln-heatmap",
    msaHeatmap(make_df(), raster = FALSE)
  )
})
