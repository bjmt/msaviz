make_df <- function() {
  msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACGAACGT",
    seq4 = "ACGTACGT"
  ), reference = "seq1")
}

test_that("composeMSA with only a heatmap returns the heatmap unchanged", {
  p_hm <- msaHeatmap(make_df())
  expect_identical(composeMSA(p_hm), p_hm)
})

test_that("composeMSA returns a patchwork object when extras are passed", {
  df <- make_df()
  p_hm <- msaHeatmap(df)
  p_pid <- msaPID(df)
  out <- composeMSA(heatmap = p_hm, top = p_pid)
  expect_s3_class(out, "patchwork")
})

test_that("groups must be named and cover every sequence in the heatmap", {
  p_hm <- msaHeatmap(make_df())
  expect_error(composeMSA(heatmap = p_hm, groups = c("A","B","A","B")),
    "must be a named")
  partial <- c(seq1 = "A", seq2 = "B")
  expect_error(composeMSA(heatmap = p_hm, groups = partial),
    "missing entries for sequence")
})

test_that("groups produce a strip with one tile per sequence", {
  df <- make_df()
  p_hm <- msaHeatmap(df)
  groups <- c(seq1 = "A", seq2 = "A", seq3 = "B", seq4 = "B")
  out <- composeMSA(heatmap = p_hm, groups = groups)
  expect_s3_class(out, "patchwork")
})
