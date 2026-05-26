test_that("encode_aln round-trips letters", {
  m <- rbind(
    s1 = strsplit("ACGTAC", "")[[1]],
    s2 = strsplit("ACGTAA", "")[[1]],
    s3 = strsplit("ACG-AC", "")[[1]]
  )
  enc <- msaviz:::encode_aln(m)
  expect_equal(dim(enc$codes), dim(m))
  expect_equal(rownames(enc$codes), rownames(m))
  reconstructed <- matrix(enc$letters[enc$codes], nrow = nrow(m),
    dimnames = dimnames(m))
  expect_equal(reconstructed, m)
})

test_that("vectorized consensus matches the naive table-based one", {
  set.seed(1)
  m <- rbind(
    s1 = sample(c("A","C","G","T","-"), 30, replace = TRUE, prob = c(4,2,2,1,1)),
    s2 = sample(c("A","C","G","T","-"), 30, replace = TRUE, prob = c(3,3,2,1,1)),
    s3 = sample(c("A","C","G","T","-"), 30, replace = TRUE, prob = c(2,3,3,1,1)),
    s4 = sample(c("A","C","G","T","-"), 30, replace = TRUE, prob = c(1,2,2,5,0)),
    s5 = sample(c("A","C","G","T","-"), 30, replace = TRUE, prob = c(2,2,2,2,1))
  )
  naive_with_gaps <- apply(m, 2, function(x) names(sort(table(x), decreasing = TRUE))[1])
  naive_no_gaps <- apply(m, 2, function(x) {
    x <- x[x != "-"]
    if (!length(x)) "-" else names(sort(table(x), decreasing = TRUE))[1]
  })

  enc <- msaviz:::encode_aln(m)
  expect_equal(
    msaviz:::consensus_from_codes(enc$codes, enc$letters, "-", include.gaps = TRUE),
    naive_with_gaps
  )
  expect_equal(
    msaviz:::consensus_from_codes(enc$codes, enc$letters, "-", include.gaps = FALSE),
    naive_no_gaps
  )
})
