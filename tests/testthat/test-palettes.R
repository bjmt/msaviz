test_that("palettes cover their advertised alphabets", {
  expect_setequal(names(msa_palette_DNA), c("A", "C", "G", "T", "N"))
  expect_setequal(names(msa_palette_RNA), c("A", "C", "G", "U", "N"))
  expect_true(all(c("A","R","N","D","C","Q","E","G","H","I",
                    "L","K","M","F","P","S","T","W","Y","V","X")
                  %in% names(msa_palette_AA)))
})

test_that("every palette value is a valid colour", {
  for (palette in list(msa_palette_DNA, msa_palette_RNA, msa_palette_AA)) {
    expect_silent(grDevices::col2rgb(palette))
  }
})
