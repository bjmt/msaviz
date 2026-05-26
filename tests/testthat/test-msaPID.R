make_df <- function() {
  msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACGAACGT",
    seq4 = "ACGTACGT"
  ), reference = "seq1")
}

test_that("msaPID returns a ggplot with one layer", {
  p <- msaPID(make_df())
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 1L)
})

test_that("style = 'line' uses geom_step", {
  p <- msaPID(make_df(), style = "line")
  expect_s3_class(p, "ggplot")
  geom_name <- class(p$layers[[1]]$geom)[1]
  expect_match(geom_name, "GeomStep")
})

test_that("x-limits line up with msaHeatmap", {
  df <- make_df()
  p_pid <- msaPID(df)
  p_hm <- msaHeatmap(df)
  pid_x <- ggplot2::ggplot_build(p_pid)$layout$panel_params[[1]]$x.range
  hm_x  <- ggplot2::ggplot_build(p_hm)$layout$panel_params[[1]]$x.range
  expect_equal(pid_x, hm_x, tolerance = 1e-6)
})
