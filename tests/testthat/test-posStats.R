make_df <- function() {
  msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACGAACGT",
    seq4 = "ACG-ACGT"
  ), reference = "seq1", drop.gaps = FALSE)
}

test_that("posStats returns the documented columns and one row per position", {
  ps <- posStats(make_df())
  expect_s3_class(ps, "tbl_df")
  expect_equal(nrow(ps), 8L)
  expect_setequal(colnames(ps),
    c("Position", "n", "nGaps", "pctIdentity", "mostCommon",
      "mostCommonPct", "entropy"))
  expect_type(ps$n, "integer")
  expect_type(ps$nGaps, "integer")
})

test_that("pctIdentity is 1 on perfectly conserved positions", {
  df <- make_df()
  ps <- posStats(df)
  # Positions 1, 2, 3 are "A", "C", "G" across all four input rows
  expect_equal(ps$pctIdentity[1:3], c(1, 1, 1))
})

test_that("nGaps counts gap rows even when alnDF drops them", {
  # drop.gaps = TRUE removes the gap entry at position 4 of seq4
  df_dropped <- msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACGAACGT",
    seq4 = "ACG-ACGT"
  ), reference = "seq1", drop.gaps = TRUE)
  ps <- posStats(df_dropped)
  expect_equal(ps$nGaps[4], 1L)
  expect_equal(ps$n[4], 3L)
})

test_that("entropy is 0 on a perfectly conserved column and positive on a mixed one", {
  df <- make_df()
  ps <- posStats(df)
  expect_equal(ps$entropy[1], 0)
  expect_gt(ps$entropy[5], 0)
})
