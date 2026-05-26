test_that("saveHeatmap writes a non-empty file at the chosen path", {
  skip_if_not_installed("withr")
  df <- msa2DF(c(a = "ACGT", b = "ACTT"), reference = "a")
  p <- msaHeatmap(df, raster = FALSE)
  out <- withr::local_tempfile(fileext = ".pdf")
  suppressMessages(saveHeatmap(p, out))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})

test_that("explicit height/width are honoured", {
  skip_if_not_installed("withr")
  df <- msa2DF(c(a = "ACGT", b = "ACTT"), reference = "a")
  p <- msaHeatmap(df, raster = FALSE)
  out <- withr::local_tempfile(fileext = ".pdf")
  expect_invisible(suppressMessages(
    saveHeatmap(p, out, height = 30, width = 80)
  ))
  expect_true(file.exists(out))
})
