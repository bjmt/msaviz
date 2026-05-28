#' Compute the consensus sequence of a multiple sequence alignment
#'
#' Returns the most common character at each position. By default, gap
#' characters compete with letters when picking the most common; set
#' `include.gaps = FALSE` to ignore gaps when the alphabet supports a clear
#' winner at a position.
#'
#' Ties between equally common characters are broken by occurrence order in
#' the input (the first letter to appear wins).
#'
#' @param aln Alignment, accepted in the same forms as [msa2DF()].
#' @param gap.chars Character vector of one or more single-character strings
#'   treated as gaps.
#' @param include.gaps Logical. If `TRUE` (default), gap characters count
#'   toward the per-position consensus. If `FALSE`, gaps are dropped from
#'   the count and the most common non-gap letter wins; positions that are
#'   entirely gaps fall back to `gap.chars[1]`.
#'
#' @return A single character string of length equal to the alignment width.
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACG-ACGT"
#' )
#' consensusSeq(aln)
#' consensusSeq(aln, include.gaps = FALSE)
#'
#' @export
consensusSeq <- function(aln, gap.chars = "-", include.gaps = TRUE) {
  if (is(aln, "alignment")) {
    aln <- structure(unlist(aln$seq), names = aln$nam)
  } else if (!is.character(aln)) {
    aln <- as.character(aln)
  }
  if (length(aln) == 0L) {
    stop("Alignment is empty (no sequences)", call. = FALSE)
  }
  widths <- nchar(aln)
  if (length(unique(widths)) != 1L) {
    stop("All sequences in the alignment must be the same size", call. = FALSE)
  }
  if (widths[1L] == 0L) {
    stop("Alignment has zero width (sequences are empty)", call. = FALSE)
  }
  aln2 <- do.call(rbind, strsplit(aln, "", fixed = TRUE))
  enc <- encode_aln(aln2)
  paste(
    consensus_from_codes(enc$codes, enc$letters,
      gap.chars = gap.chars, include.gaps = include.gaps),
    collapse = ""
  )
}
