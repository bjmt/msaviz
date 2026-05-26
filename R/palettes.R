#' Default colour palettes for nucleic acid and amino acid alignments
#'
#' Named character vectors that map each letter of an alphabet to a hex
#' colour, suitable for passing as `letter.colours =` to
#' [msaHeatmap()] when `column = "Letter"`. Each palette ships with both
#' upper- and lower-case keys so the same palette works whether the
#' upstream alignment reader (e.g. [seqinr::read.alignment()], which
#' lowercases by default) preserved case or not.
#'
#' `msa_palette_DNA` covers `A C G T N`. `msa_palette_RNA` is the same
#' palette with `T` replaced by `U`. `msa_palette_AA` covers the 20
#' standard amino acids plus `X` (unknown) and follows a Clustal-X-style
#' physico-chemical colour scheme.
#'
#' @return A named character vector of hex colour strings.
#'
#' @examples
#' msa_palette_DNA
#' msa_palette_AA
#'
#' aln <- c(seq1 = "ACGTACGT", seq2 = "ACGTTCGT")
#' msaHeatmap(msa2DF(aln, reference = "seq1"),
#'   column = "Letter", letter.colours = msa_palette_DNA)
#'
#' @name msa_palettes
NULL

#' @noRd
both_cases <- function(x) {
  out <- c(x, x)
  names(out) <- c(toupper(names(x)), tolower(names(x)))
  out[!duplicated(names(out))]
}

#' @rdname msa_palettes
#' @export
msa_palette_DNA <- both_cases(c(
  A = "#2ca02c",
  C = "#1f77b4",
  G = "#ff7f0e",
  T = "#d62728",
  N = "grey60"
))

#' @rdname msa_palettes
#' @export
msa_palette_RNA <- both_cases(c(
  A = "#2ca02c",
  C = "#1f77b4",
  G = "#ff7f0e",
  U = "#d62728",
  N = "grey60"
))

#' @rdname msa_palettes
#' @export
msa_palette_AA <- both_cases(c(
  # Hydrophobic (Clustal blue)
  A = "#80a0f0", I = "#80a0f0", L = "#80a0f0", M = "#80a0f0",
  F = "#80a0f0", W = "#80a0f0", V = "#80a0f0",
  # Positive charge (red)
  K = "#f01505", R = "#f01505",
  # Negative charge (magenta)
  D = "#c048c0", E = "#c048c0",
  # Polar (green)
  N = "#15c015", Q = "#15c015", S = "#15c015", T = "#15c015",
  # Cysteine (pink)
  C = "#f08080",
  # Glycine (orange)
  G = "#f09048",
  # Proline (yellow)
  P = "#c0c000",
  # Aromatic (cyan)
  H = "#15a4a4", Y = "#15a4a4",
  # Unknown
  X = "grey60"
))
