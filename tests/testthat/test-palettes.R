test_that("palettes cover their advertised alphabets in both cases", {
  expect_true(all(c("A","C","G","T","N","a","c","g","t","n")
                  %in% names(msa_palette_DNA)))
  expect_true(all(c("A","C","G","U","N","a","c","g","u","n")
                  %in% names(msa_palette_RNA)))
  upper <- c("A","R","N","D","C","Q","E","G","H","I",
             "L","K","M","F","P","S","T","W","Y","V","X")
  expect_true(all(upper %in% names(msa_palette_AA)))
  expect_true(all(tolower(upper) %in% names(msa_palette_AA)))
})

test_that("upper- and lowercase keys map to the same colour", {
  expect_equal(msa_palette_DNA[["A"]], msa_palette_DNA[["a"]])
  expect_equal(msa_palette_AA[["L"]],  msa_palette_AA[["l"]])
})

test_that("every palette value is a valid colour", {
  for (palette in list(msa_palette_DNA, msa_palette_RNA, msa_palette_AA)) {
    expect_silent(grDevices::col2rgb(palette))
  }
})
