tiny_aln <- c(
  seq1 = "ACGTACGT",
  seq2 = "ACGTTCGT",
  seq3 = "ACG-ACGT"
)

test_that("character-vector input produces a tibble with the documented columns and attrs", {
  df <- msa2DF(tiny_aln, reference = "seq1")
  expect_s3_class(df, "tbl_df")
  expect_setequal(colnames(df), c("Sequence", "Position", "Letter", "Aln"))
  expect_s3_class(df$Sequence, "factor")
  expect_type(df$Position, "integer")
  expect_s3_class(df$Letter, "factor")
  expect_s3_class(df$Aln, "factor")
  expect_setequal(levels(df$Aln), c("Ref", "Alt", "Gap"))
  expect_equal(attr(df, "aln.size"), 8L)
  expect_equal(attr(df, "drop.gaps"), TRUE)
  expect_equal(attr(df, "reference"), "seq1")
})

test_that("drop.gaps = TRUE removes gap rows; FALSE keeps them as Aln = 'Gap' with NA letter", {
  kept <- msa2DF(tiny_aln, reference = "seq1", drop.gaps = TRUE)
  expect_false(any(kept$Aln == "Gap"))

  full <- msa2DF(tiny_aln, reference = "seq1", drop.gaps = FALSE)
  gap_rows <- full[full$Aln == "Gap", ]
  expect_gt(nrow(gap_rows), 0)
  expect_true(all(is.na(gap_rows$Letter)))
})

test_that("reference = NULL builds a consensus and drops it again by default", {
  df <- msa2DF(tiny_aln)
  expect_equal(attr(df, "reference"), "consensus")
  expect_false("consensus" %in% as.character(df$Sequence))

  df_keep <- msa2DF(tiny_aln, keep.consensus = TRUE)
  expect_true("consensus" %in% as.character(df_keep$Sequence))
})

test_that("uppercase = TRUE normalises mixed-case input", {
  mixed <- c(s1 = "acgt", s2 = "ACGA")
  df <- msa2DF(mixed, reference = "s1", uppercase = TRUE)
  expect_true(all(toupper(as.character(df$Letter)) == as.character(df$Letter)))
  # The single Alt should still be detected after case folding
  expect_equal(sum(df$Aln == "Alt"), 1)
})

test_that("mismatched sequence lengths raise the documented error", {
  bad <- c(a = "AAAA", b = "AAA")
  expect_error(msa2DF(bad), "must be the same size")
})

test_that("missing reference name raises the documented error", {
  expect_error(msa2DF(tiny_aln, reference = "does-not-exist"),
    "No matching reference")
})

test_that("seqinr::alignment input works", {
  skip_if_not_installed("seqinr")
  aln <- seqinr::read.alignment(
    system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz"),
    format = "fasta"
  )
  df <- msa2DF(aln)
  expect_s3_class(df, "tbl_df")
  expect_equal(attr(df, "aln.size"), 250L)
  expect_equal(nlevels(df$Sequence), 20L)
})

test_that("Biostrings multiple-alignment input works", {
  skip_if_not_installed("Biostrings")
  aln <- Biostrings::AAMultipleAlignment(c(
    "ACDEFGH-",
    "ACDEAGHK",
    "ACDEFGHK"
  ))
  df <- msa2DF(aln, reference = NULL)
  expect_s3_class(df, "tbl_df")
  expect_equal(attr(df, "aln.size"), 8L)
})
