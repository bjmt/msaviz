tiny_aln <- c(
  s1 = "ACGTACGT",
  s2 = "ACGTTCGT",
  s3 = "ACG-ACGT"
)

test_that("consensusSeq returns a single string of alignment width", {
  cs <- consensusSeq(tiny_aln)
  expect_type(cs, "character")
  expect_length(cs, 1L)
  expect_equal(nchar(cs), unname(nchar(tiny_aln[1])))
})

test_that("include.gaps = FALSE falls back to a letter when most are gaps", {
  aln <- c(
    s1 = "A-",
    s2 = "-G",
    s3 = "-G"
  )
  # Position 1: A vs - vs - -> include.gaps=TRUE picks "-", include.gaps=FALSE picks "A"
  expect_equal(consensusSeq(aln, include.gaps = TRUE), "-G")
  expect_equal(consensusSeq(aln, include.gaps = FALSE), "AG")
})

test_that("seqinr alignment input works", {
  skip_if_not_installed("seqinr")
  aln <- seqinr::read.alignment(
    system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz"),
    format = "fasta")
  cs <- consensusSeq(aln)
  expect_length(cs, 1L)
  expect_equal(nchar(cs), 250L)
})

test_that("ties are broken alphabetically", {
  aln <- c(s1 = "AC", s2 = "CA")
  # Each position has a tie between A and C. Alphabetical wins -> "AA".
  expect_equal(consensusSeq(aln), "AA")
})
