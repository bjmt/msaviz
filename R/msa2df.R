#' @noRd
matrix2DF <- function(x, rowLevels = NULL, colNames = NULL) {
  if (is.null(rowLevels)) {
    rowLevels <- rownames(x)
  }
  df <- data.frame(
    row.names = NULL, check.names = FALSE,
    Var1 = rep.int(factor(rownames(x), levels = rowLevels), ncol(x)),
    Var2 = rep.int(seq_len(ncol(x)), rep.int(nrow(x), ncol(x)))
  )
  df$value <- as.vector(x)
  if (!is.null(colNames)) {
    colnames(df) <- colNames
  }
  df
}

#' @noRd
encode_aln <- function(aln_matrix) {
  # Sort the alphabet so ties in the consensus break alphabetically — this
  # matches the naive `sort(table(x), decreasing = TRUE)` semantics.
  letters_vec <- sort(unique.default(c(aln_matrix)))
  codes <- matrix(
    match(c(aln_matrix), letters_vec),
    nrow = nrow(aln_matrix),
    dimnames = dimnames(aln_matrix)
  )
  list(codes = codes, letters = letters_vec)
}

#' @noRd
consensus_from_codes <- function(codes, letters_vec, gap.chars = "-",
  include.gaps = TRUE) {
  L <- length(letters_vec)
  counts <- apply(codes, 2, tabulate, nbins = L)
  if (!is.matrix(counts)) counts <- matrix(counts, nrow = L)
  if (!include.gaps) {
    gap_codes <- which(letters_vec %in% gap.chars)
    if (length(gap_codes)) {
      counts[gap_codes, ] <- 0L
    }
    all_zero <- colSums(counts) == 0L
    if (any(all_zero)) {
      gap_idx <- if (length(gap_codes)) gap_codes[1] else {
        stop("Cannot compute gapless consensus: an entire column is gaps and ",
          "no gap character is present in the alphabet to fall back to.")
      }
      idx <- integer(ncol(counts))
      idx[!all_zero] <- max.col(t(counts[, !all_zero, drop = FALSE]),
        ties.method = "first")
      idx[all_zero] <- gap_idx
    } else {
      idx <- max.col(t(counts), ties.method = "first")
    }
  } else {
    idx <- max.col(t(counts), ties.method = "first")
  }
  letters_vec[idx]
}

#' Convert a multiple sequence alignment to a tidy data frame
#'
#' Reshapes an aligned set of equal-length sequences into a long-format tibble
#' with one row per (sequence, position), suitable for passing to
#' [msaHeatmap()] or any other ggplot2-based visualization.
#'
#' `aln` can be a named character vector (one element per sequence), an
#' `alignment` object as returned by [seqinr::read.alignment()], or any
#' object that yields a character vector through [base::as.character()]
#' (e.g. a `Biostrings` `DNAMultipleAlignment` or `AAMultipleAlignment`).
#' When `reference` is `NULL` (the default) a consensus sequence is computed
#' and used as the reference for the `Aln` column.
#'
#' Sequence names must be unique; duplicates are rejected with an error.
#'
#' @param aln Alignment to convert. See Details.
#' @param reference Name of the sequence in `aln` to use as reference. When
#'   `NULL`, a consensus sequence is computed and used.
#' @param drop.gaps Logical. If `TRUE` (default), rows where the reference or
#'   the sequence holds a gap character are dropped from the returned tibble.
#'   If `FALSE` gaps are kept and marked with `Aln = "Gap"` and
#'   `Letter = NA`.
#' @param gap.chars Character vector of one or more single-character strings
#'   treated as gaps.
#' @param keep.consensus Logical. When `reference` is `NULL`, controls whether
#'   the computed consensus row is kept in the output. Defaults to `FALSE`.
#' @param consensus.gaps Logical. When computing the consensus, count gaps
#'   when picking the most-common character (`TRUE`) or ignore them
#'   (`FALSE`).
#' @param uppercase Logical. If `TRUE`, force every letter to upper case
#'   before comparing to the reference. Useful when input sequences mix cases.
#' @param verbose Logical. Print progress and size diagnostics to the console.
#'
#' @return A [tibble][tibble::tibble-package] with columns `Sequence`
#'   (factor), `Position` (integer), `Letter` (factor), and `Aln` (factor with
#'   levels `Ref`, `Alt`, `Gap`). The result carries the attributes
#'   `aln.size`, `drop.gaps` and `reference`.
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACG-ACGT"
#' )
#' msa2DF(aln, reference = "seq1")
#'
#' @export
msa2DF <- function(aln, reference = NULL, drop.gaps = TRUE, gap.chars = "-",
  keep.consensus = FALSE, consensus.gaps = TRUE, uppercase = FALSE, verbose = FALSE) {
  calcConsensus <- is.null(reference)
  if (verbose) {
    message("Input aln object size: ", format(object.size(aln), units = "Mb"))
  }
  if (is(aln, "alignment")) {
    aln <- structure(unlist(aln$seq), names = aln$nam)
  } else if (!is.character(aln)) {
    aln <- as.character(aln)
  }
  if (is.null(names(aln))) {
    names(aln) <- as.character(seq_len(length(aln)))
  }
  dup_idx <- which(duplicated(names(aln)))
  if (length(dup_idx)) {
    dupes <- unique(names(aln)[dup_idx])
    stop("Alignment contains duplicate sequence names: ",
      paste(head(dupes, 5), collapse = ", "),
      if (length(dupes) > 5) paste0(" (and ", length(dupes) - 5L, " more)") else "",
      call. = FALSE)
  }
  if (!is.null(reference) && !reference %in% names(aln)) {
    stop("No matching reference in alignment", call. = FALSE)
  }
  if (length(unique(nchar(aln))) != 1L) {
    stop("All sequences in the alignment must be the same size", call. = FALSE)
  }
  if (uppercase) aln <- toupper(aln)
  if (verbose) {
    message("Parsing alignment containing ", format(length(aln), big.mark = ","),
      " sequences and ", format(unique(nchar(aln)), big.mark = ","), " positions...")
  }
  aln2 <- do.call(rbind, strsplit(aln, "", fixed = TRUE))
  rownames(aln2) <- names(aln)
  if (verbose) {
    message("Matrix aln object size: ", format(object.size(aln2), units = "Mb"))
  }

  if (is.null(reference)) {
    if ("consensus" %in% rownames(aln2)) {
      stop("Alignment already contains a sequence named 'consensus'", call. = FALSE)
    }
    enc <- encode_aln(aln2)
    consensus <- consensus_from_codes(enc$codes, enc$letters,
      gap.chars = gap.chars, include.gaps = consensus.gaps)
    aln2 <- rbind(consensus = consensus, aln2)
    reference <- "consensus"
  }

  alnDF <- matrix2DF(aln2, rowLevels = rownames(aln2),
    colNames = c("Sequence", "Position", "Letter"))
  if (verbose) {
    message("Initial alnDF object size: ", format(object.size(alnDF), units = "Mb"))
  }
  attributes(alnDF)$aln.size <- ncol(aln2)
  attributes(alnDF)$drop.gaps <- drop.gaps

  ref <- as.character(aln2[reference, ])
  alnDF$Aln <- "Ref"
  alnDF$Aln[alnDF$Letter != ref[alnDF$Position]] <- "Alt"
  alnDF$Aln[alnDF$Letter %in% gap.chars] <- "Gap"
  attributes(alnDF)$reference <- reference
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
  if (verbose) {
    message("Alphabet contains ", length(levels(alnDF$Letter)), " letters.")
    if (!is.null(reference)) {
      nAlt <- as.vector(table(alnDF$Aln)["Alt"])
      message("Alignment contains ", format(nAlt, big.mark = ","),
        " Alt letters (", round(100 * (nAlt / sum(!is.na(alnDF$Letter))), 2),
        "% of all non-gap letters).")
    }
    message("Final alnDF object size: ", format(object.size(alnDF), units = "Mb"))
  }
  as_tibble(alnDF, rownames = NULL)
}
