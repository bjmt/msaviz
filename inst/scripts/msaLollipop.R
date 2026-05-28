#!/usr/bin/env Rscript --no-save --no-restore

args <- commandArgs(TRUE)

helpFun <- function(quit = !interactive()) {
  message(appendLF = FALSE, paste0(collapse = "\n", c(
"msaLollipop v1.3.0  Copyright (C) 2023-2026 Benjamin Jean-Marie Tremblay",
"",
"Usage:  msaLollipop [options] --input file.fa",
"",
" --input       <str>  Input file.",
" --output      <str>  Output file. Default: lollipop.pdf",
" --read-fun    <str>  Alignment read function. Possible functions:",
"                      - seqinr::read.alignment (default)",
"                      - Biostrings::readDNAMultipleAlignment",
"                      - Biostrings::readAAMultipleAlignment",
" --aln-fmt     <str>  Alignment file format. One of: fasta (default), clustal,",
"                      phylip, stockholm (Biostrings only), mase (seqinr only),",
"                      msf (seqinr only).",
" --reference   <str>  Sequence name to use as the reference. Default: a",
"                      consensus sequence is computed (RefPosAlt labels then",
"                      fall back to position-only).",
" --highlight   <str>  Comma-separated list of values to draw lollipops for.",
"                      Default: Alt",
" --highlight-by <str> Column to match --highlight against. Default: Aln",
" --no-labels          Suppress the RefPosAlt labels above each head.",
" --no-baseline        Suppress the per-position baseline segments.",
" --no-drop-empty      Keep sequences with no highlighted cells as empty rows.",
" --width       <num>  Output width in inches. Default: 11.5",
" --height      <num>  Output height in inches. Default: 7",
" --verbose            Print progress information.",
" --help               Show this help message.",
""
  )))
  if (quit) quit(save = "no")
}

stop2 <- function(...) {
  msg <- paste0(c(...), collapse = "")
  msg <- paste0(msg, "\nRun msaLollipop --help to see usage.")
  stop(msg, call. = FALSE)
}

getArg <- function(x, default = NULL, inputs = NULL, isFlag = FALSE,
  deleteArg = TRUE, optional = TRUE) {
  x <- paste0("--", x)
  i <- which(args == x)
  if (!length(i)) {
    if (!optional) {
      stop2("Missing argument ", x)
    }
    if (isFlag) return(FALSE) else return(default)
  }
  if (length(i) > 1) {
    stop2("Found more than one argument for ", x)
  }
  if (isFlag) {
    if (deleteArg) args <<- args[-i]
    return(TRUE)
  }
  if (i == length(args)) {
    stop2("Missing value for argument ", x)
  }
  arg <- args[i + 1]
  if (!nchar(arg) || substr(arg, 1, 2) == "--") {
    stop2("Missing value for argument ", x)
  }
  if (length(args) > (i + 1)) {
    if (substr(args[i + 2], 1, 2) != "--") {
      stop2("Found more than one value for argument ", x)
    }
  }
  if (deleteArg) args <<- args[-c(i, i + 1)]
  if (!is.null(inputs) && !arg %in% inputs) {
    stop2("Incorrect value for argument ", x,
      " [", paste0(inputs, collapse = ", "), "]")
  }
  arg
}

if (getArg("help", isFlag = TRUE)) helpFun()

readFun <- getArg("read-fun", default = "seqinr::read.alignment",
  inputs = c("seqinr::read.alignment", "Biostrings::readDNAMultipleAlignment",
    "Biostrings::readAAMultipleAlignment"))
if (readFun == "seqinr::read.alignment") {
  alnFmt <- getArg("aln-fmt", default = "fasta",
    inputs = c("mase", "clustal", "phylip", "fasta", "msf"))
} else {
  alnFmt <- getArg("aln-fmt", default = "fasta",
    inputs = c("clustal", "phylip", "fasta", "stockholm"))
}
inFile <- getArg("input", optional = FALSE)
outFile <- getArg("output", default = "lollipop.pdf")
reference <- getArg("reference", default = NULL)
highlightStr <- getArg("highlight", default = "Alt")
highlightBy <- getArg("highlight-by", default = "Aln")
noLabels <- getArg("no-labels", isFlag = TRUE)
noBaseline <- getArg("no-baseline", isFlag = TRUE)
noDropEmpty <- getArg("no-drop-empty", isFlag = TRUE)
widthArg <- getArg("width", default = "11.5")
heightArg <- getArg("height", default = "7")
verbose <- getArg("verbose", isFlag = TRUE)

if (length(args)) {
  stop2("Found ", length(args), " unused argument(s) [",
    paste0(args, collapse = " "), "]")
}

widthArg <- suppressWarnings(as.numeric(widthArg))
if (is.na(widthArg)) stop2("--width must be numeric")
heightArg <- suppressWarnings(as.numeric(heightArg))
if (is.na(heightArg)) stop2("--height must be numeric")

highlight <- strsplit(highlightStr, ",", fixed = TRUE)[[1]]

if (verbose) {
  message("msaLollipop v1.3.0")
  message("------------------")
  message("Using the following library paths to look for packages:\n",
    paste0(paste0("  ", .libPaths()), collapse = "\n"))
}

if (!requireNamespace("msaviz", quietly = TRUE)) {
  stop2("Could not find the msaviz package. Install it with ",
    "remotes::install_github(\"bjmt/msaviz\").")
}

readFunParts <- strsplit(readFun, "::", fixed = TRUE)[[1]]
if (!requireNamespace(readFunParts[1], quietly = TRUE)) {
  stop2("Could not find package ", readFunParts[1])
}
readFn <- get(readFunParts[2], envir = asNamespace(readFunParts[1]))

if (verbose) message("Reading alignment from ", inFile, " (", alnFmt, ")...")
aln <- readFn(inFile, format = alnFmt)

alnDF <- msaviz::msa2DF(aln, reference = reference, drop.gaps = FALSE,
  keep.consensus = TRUE, verbose = verbose)

if (verbose) message("Building lollipop chart...")
p <- msaviz::msaLollipop(alnDF,
  highlight = highlight,
  highlight.by = highlightBy,
  labels = !noLabels,
  baseline = !noBaseline,
  drop.empty = !noDropEmpty)

if (verbose) message("Saving lollipop chart to ", outFile, "...")
ggplot2::ggsave(outFile, p, width = widthArg, height = heightArg)

if (verbose) message("Done.")
