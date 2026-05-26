#' Per-position summary statistics for an alignment
#'
#' Computes a small set of per-position summaries from a tidy alignment data
#' frame as produced by [msa2DF()]: how many non-gap entries cover the
#' position, how many gaps, the fraction matching the reference (`Aln ==
#' "Ref"`), the most common letter, and the Shannon entropy of the column.
#'
#' The result has one row per alignment position, including positions that
#' are entirely gaps (which receive `NA` for letter-based stats).
#'
#' @param alnDF A tibble produced by [msa2DF()] containing the columns
#'   `Sequence`, `Position`, `Letter`, and `Aln`.
#' @param include.gaps Logical. If `TRUE`, gap characters are included when
#'   computing `mostCommon` and `entropy`. If `FALSE` (default) gaps are
#'   excluded from those statistics.
#'
#' @return A tibble with one row per position and columns `Position`, `n`
#'   (non-gap count), `nGaps`, `pctIdentity` (`n_Ref / n`), `mostCommon`,
#'   `mostCommonPct`, and `entropy` (Shannon, natural log).
#'
#' @examples
#' aln <- c(
#'   seq1 = "ACGTACGT",
#'   seq2 = "ACGTTCGT",
#'   seq3 = "ACGAACGT",
#'   seq4 = "ACGTACGT"
#' )
#' posStats(msa2DF(aln, reference = "seq1"))
#'
#' @export
posStats <- function(alnDF, include.gaps = FALSE) {
  if (any(!c("Sequence", "Position", "Letter", "Aln") %in% colnames(alnDF))) {
    stop("alnDF must contain columns Sequence, Position, Letter, Aln",
      call. = FALSE)
  }

  aln_size <- attributes(alnDF)$aln.size
  if (is.null(aln_size)) aln_size <- max(alnDF$Position)
  all_positions <- seq_len(aln_size)

  nseq <- nlevels(alnDF$Sequence)
  if (!nseq) nseq <- length(unique(alnDF$Sequence))

  is_gap <- is.na(alnDF$Letter) | as.character(alnDF$Aln) == "Gap"
  nongap <- alnDF[!is_gap, , drop = FALSE]

  pos_nongap <- factor(nongap$Position, levels = all_positions)
  n_nongap <- as.integer(tabulate(pos_nongap, nbins = aln_size))
  n_gap <- as.integer(nseq) - n_nongap

  ref_pos <- factor(nongap$Position[nongap$Aln == "Ref"], levels = all_positions)
  n_ref <- as.integer(tabulate(ref_pos, nbins = aln_size))
  pctIdentity <- ifelse(n_nongap > 0L, n_ref / n_nongap, NA_real_)

  if (include.gaps) {
    letter_str <- ifelse(is_gap, "-", as.character(alnDF$Letter))
    pos_combined <- as.integer(alnDF$Position)
  } else {
    letter_str <- as.character(nongap$Letter)
    pos_combined <- as.integer(nongap$Position)
  }

  by_pos <- split(letter_str, factor(pos_combined, levels = all_positions))
  mostCommon <- character(aln_size)
  mostCommonPct <- numeric(aln_size)
  entropy <- numeric(aln_size)
  for (i in seq_along(by_pos)) {
    letters_at <- by_pos[[i]]
    if (length(letters_at) == 0L) {
      mostCommon[i] <- NA_character_
      mostCommonPct[i] <- NA_real_
      entropy[i] <- NA_real_
      next
    }
    tab <- tabulate(factor(letters_at))
    labels <- levels(factor(letters_at))
    max_i <- which.max(tab)
    n_here <- length(letters_at)
    mostCommon[i] <- labels[max_i]
    mostCommonPct[i] <- tab[max_i] / n_here
    p_i <- tab / n_here
    p_i <- p_i[p_i > 0]
    entropy[i] <- -sum(p_i * log(p_i))
  }

  tibble::tibble(
    Position = all_positions,
    n = n_nongap,
    nGaps = n_gap,
    pctIdentity = pctIdentity,
    mostCommon = mostCommon,
    mostCommonPct = mostCommonPct,
    entropy = entropy
  )
}
