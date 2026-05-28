make_df <- function() {
  msa2DF(c(
    seq1 = "ACGTACGT",
    seq2 = "ACGTTCGT",
    seq3 = "ACGAACGT",
    seq4 = "AAGTACGT"
  ), reference = "seq1", drop.gaps = FALSE, keep.consensus = TRUE)
}

test_that("msaLollipop returns a ggplot with the expected default layers", {
  p <- msaLollipop(make_df())
  expect_s3_class(p, "ggplot")
  # baseline + stems + heads + labels
  expect_length(p$layers, 4L)
})

test_that("labels = FALSE drops the label layer; baseline = FALSE drops baseline", {
  p1 <- msaLollipop(make_df(), labels = FALSE)
  expect_length(p1$layers, 3L)
  p2 <- msaLollipop(make_df(), labels = FALSE, baseline = FALSE)
  expect_length(p2$layers, 2L)
})

test_that("drop.empty controls whether the reference row appears", {
  df <- make_df()
  p_drop <- msaLollipop(df, drop.empty = TRUE)
  expect_false("seq1" %in% attr(p_drop, "row.order"))
  p_keep <- msaLollipop(df, drop.empty = FALSE)
  expect_true("seq1" %in% attr(p_keep, "row.order"))
})

test_that("row.order accepts names not in alnDF and keeps them as empty rows", {
  df <- make_df()
  custom <- c("seq2", "spacer1", "seq3", "spacer2", "seq4")
  p <- msaLollipop(df, drop.empty = FALSE, row.order = custom)
  expect_equal(attr(p, "row.order"), custom)
})

test_that("highlight.by switches the filter column", {
  df <- make_df()
  # Filter on Letter == "T" instead of Aln == "Alt"
  p <- msaLollipop(df,
    highlight = "T", highlight.by = "Letter", labels = FALSE)
  expect_s3_class(p, "ggplot")
  # The stem layer (layer index 2 when baseline is on; 1 when off) has data
  # filtered to Letter == "T"
  stem_data <- p$layers[[2]]$data
  expect_true(all(as.character(stem_data$Letter) == "T"))
})

test_that("highlight.by, y, and head.fill.by errors are clear", {
  df <- make_df()
  expect_error(
    msaLollipop(df, highlight.by = "DoesNotExist"),
    "not a column in alnDF"
  )
  expect_error(
    msaLollipop(df, y = "DoesNotExist"),
    "not a column in alnDF"
  )
  expect_error(
    msaLollipop(df, head.fill.by = "DoesNotExist"),
    "not a column in alnDF"
  )
})

test_that("head.fill.by maps fill to a column on the head layer", {
  df <- make_df()
  p <- msaLollipop(df, head.fill.by = "Letter",
    head.fill.colours = msa_palette_DNA)
  # Layer order: baseline, stem, head, label
  head_layer <- p$layers[[3]]
  expect_true("fill" %in% names(head_layer$mapping))
  # Default (no head.fill.by) does not add a fill aesthetic to the head
  p2 <- msaLollipop(df)
  expect_false("fill" %in% names(p2$layers[[3]]$mapping))
})

test_that("custom y column threads per-cell heights into the stem", {
  df <- make_df()
  set.seed(1)
  df$score <- runif(nrow(df))
  p <- msaLollipop(df, y = "score", labels = FALSE)
  expect_s3_class(p, "ggplot")
  # The stem/head layer data inherits the parent ggplot's data, where the
  # normalised per-cell heights are stored as y_norm.
  expect_true("y_norm" %in% colnames(p$data))
  expect_lte(max(p$data$y_norm, na.rm = TRUE), 1)
})

test_that("the reference-letter fallback fires when reference is absent", {
  df <- msa2DF(c(
    seq1 = "ACGT",
    seq2 = "ACTT",
    seq3 = "AAGT"
  ), drop.gaps = FALSE)  # no keep.consensus -> consensus dropped
  expect_message(
    msaLollipop(df),
    "labels fall back to position-only"
  )
})

test_that("msaLollipop carries the attributes composeMSA needs", {
  p <- msaLollipop(make_df())
  expect_true(!is.null(attr(p, "row.order")))
  expect_equal(attr(p, "legend.pos"), "none")
  expect_equal(attr(p$data, "aln.size"), 8L)
})

test_that("composeMSA accepts msaLollipop as the heatmap argument", {
  df <- make_df()
  lol <- msaLollipop(df)
  pid <- msaPID(df)

  # Heatmap-only collapses to the lollipop itself
  expect_identical(composeMSA(lol), lol)

  # With top
  expect_s3_class(composeMSA(heatmap = lol, top = pid), "patchwork")

  # With row groups
  expect_s3_class(
    composeMSA(heatmap = lol,
      groups = c(seq1 = "A", seq2 = "A", seq3 = "B", seq4 = "B")),
    "patchwork"
  )

  # With column groups
  expect_s3_class(
    composeMSA(heatmap = lol,
      col.groups = rep(c("L", "R"), c(4, 4))),
    "patchwork"
  )
})

test_that("snapshot of a small lollipop stays stable", {
  testthat::announce_snapshot_file(name = "small-aln-lollipop.svg")
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    "small-aln-lollipop",
    msaLollipop(make_df())
  )
})
