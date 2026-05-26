## Compare the 1.0.0 msa2DF implementation against the 1.1.0 rewrite.
##
## Run with: Rscript bench/msa2DF.R
##
## This is not part of `R CMD check`; the directory is .Rbuildignore'd.

stopifnot(requireNamespace("bench", quietly = TRUE))
stopifnot(requireNamespace("seqinr", quietly = TRUE))
stopifnot(requireNamespace("ggplot2", quietly = TRUE))
stopifnot(requireNamespace("tibble", quietly = TRUE))

devtools::load_all(".", quiet = TRUE)

# ---- A faithful copy of the 1.0.0 hot paths, with no other changes ----
matrix2DF_old <- function(x, rowLevels = NULL, colNames = NULL) {
  if (is.null(rowLevels)) rowLevels <- rownames(x)
  df <- data.frame(
    row.names = NULL, check.names = FALSE,
    Var1 = rep.int(factor(rownames(x), levels = rowLevels), ncol(x)),
    Var2 = rep.int(seq_len(ncol(x)), rep.int(nrow(x), ncol(x)))
  )
  df$value <- as.vector(x)
  if (!is.null(colNames)) colnames(df) <- colNames
  df
}

msa2DF_v1 <- function(aln, reference = NULL, drop.gaps = TRUE, gap.chars = "-",
  keep.consensus = FALSE, consensus.gaps = TRUE, uppercase = FALSE) {
  calcConsensus <- is.null(reference)
  if (methods::is(aln, "alignment")) {
    aln <- structure(unlist(aln$seq), names = aln$nam)
  } else if (!is.character(aln)) {
    aln <- as.character(aln)
  }
  if (is.null(names(aln))) names(aln) <- as.character(seq_len(length(aln)))
  aln2 <- t(as.matrix(as.data.frame(strsplit(aln, split = "", fixed = TRUE),
    check.names = FALSE)))
  if (is.null(reference)) {
    if (consensus.gaps) {
      consensus <- apply(aln2, 2, function(x) names(sort(table(x), decreasing = TRUE))[1])
    } else {
      consensus <- apply(aln2, 2, function(x) {
        x <- x[!x %in% gap.chars]
        if (!length(x)) gap.chars[1] else names(sort(table(x), decreasing = TRUE))[1]
      })
    }
    aln2 <- rbind(consensus = consensus, aln2)
    reference <- "consensus"
  }
  alnDF <- matrix2DF_old(aln2, rowLevels = rownames(aln2),
    colNames = c("Sequence", "Position", "Letter"))
  ref <- as.character(aln2[reference, ])
  alnDF$Aln <- "Ref"
  alnDF$Aln[alnDF$Letter != ref[alnDF$Position]] <- "Alt"
  alnDF$Aln[alnDF$Letter %in% gap.chars] <- "Gap"
  if (drop.gaps) {
    alnDF <- alnDF[!alnDF$Letter %in% gap.chars, ]
  } else {
    alnDF$Letter[alnDF$Letter %in% gap.chars] <- NA
  }
  alnDF$Aln <- factor(alnDF$Aln, levels = c("Ref", "Alt", "Gap"))
  alnDF$Letter <- factor(alnDF$Letter)
  if (calcConsensus && !keep.consensus) {
    alnDF <- alnDF[alnDF$Sequence != "consensus", ]
    alnDF$Sequence <- factor(alnDF$Sequence, levels = levels(alnDF$Sequence)[-1])
  }
  tibble::as_tibble(alnDF, rownames = NULL)
}

# ---- Bundled slice (20 x 250) ----
small_aln <- seqinr::read.alignment(
  system.file("extdata", "hug_ribosomal.fa.gz", package = "msaviz"),
  format = "fasta"
)

cat("==== bundled slice (20 x 250) ====\n")
print(bench::mark(
  v1 = msa2DF_v1(small_aln),
  v2 = msa2DF(small_aln),
  iterations = 5,
  check = FALSE
))

# ---- Full Hug alignment (3083 x 2596) ----
full_path <- "41564_2016_BFnmicrobiol201648_MOESM206_ESM.txt"
if (file.exists(full_path)) {
  cat("\n==== full Hug alignment (3083 x 2596) ====\n")
  full_aln <- seqinr::read.alignment(full_path, format = "fasta")
  print(bench::mark(
    v1 = msa2DF_v1(full_aln),
    v2 = msa2DF(full_aln),
    iterations = 2,
    check = FALSE
  ))
} else {
  cat("\nFull alignment not available locally; skipping large benchmark.\n")
}
