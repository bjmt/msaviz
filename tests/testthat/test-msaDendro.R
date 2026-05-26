test_that("msaDendro returns a ggplot with an order attribute that permutes the input", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("ggdendro")
  aln <- seqinr::read.alignment(
    system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz"),
    format = "fasta"
  )
  d <- msaDendro(aln)
  expect_s3_class(d, "ggplot")
  ord <- attr(d, "order")
  expect_type(ord, "character")
  expect_setequal(ord, aln$nam)
})

test_that("flip = FALSE returns a standard bottom-up dendrogram", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("ggdendro")
  aln <- c(s1 = "ACGT", s2 = "ACGA", s3 = "ACTT", s4 = "TCTT")
  d <- msaDendro(aln, flip = FALSE)
  expect_s3_class(d, "ggplot")
})
