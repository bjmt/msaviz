## Build the bundled example alignment shipped in `inst/extdata/`.
##
## The source is the 3,083-sequence ribosomal-protein alignment from:
##
##   Hug, L. A. et al. (2016). A new view of the tree of life.
##   Nature Microbiology 1, 16048. doi:10.1038/nmicrobiol.2016.48
##   (Supplementary Data file `41564_2016_BFnmicrobiol201648_MOESM206_ESM.txt`)
##
## The full alignment is ~8 MB, which is far too big to ship inside an R
## package. We slice 20 sequences spread evenly across the input and a 250-
## column window from the middle of the alignment, then gzip the result. The
## slice is diverse enough to demonstrate `msa2DF()` / `msaHeatmap()` without
## ballooning the tarball.
##
## Re-run with `Rscript data-raw/prepare_example.R` from the package root.

src <- "41564_2016_BFnmicrobiol201648_MOESM206_ESM.txt"
dst <- "inst/extdata/hug_ribosomal.fa.gz"
n_seqs <- 20
win_start <- 800L
win_end <- 1049L

stopifnot(file.exists(src))

lines <- readLines(src)
hdr_idx <- grep("^>", lines)
n_total <- length(hdr_idx)
stopifnot(n_total > n_seqs)

pick <- as.integer(round(seq(1, n_total, length.out = n_seqs)))
pick <- unique(pick)

names_only <- sub("^>", "", lines[hdr_idx[pick]])

starts <- hdr_idx[pick] + 1L
ends <- c(hdr_idx[-1] - 1L, length(lines))[pick]
seqs <- vapply(seq_along(pick), function(i) {
  paste0(lines[starts[i]:ends[i]], collapse = "")
}, character(1))

stopifnot(length(unique(nchar(seqs))) == 1L)
seqs <- substr(seqs, win_start, win_end)

out_lines <- character(2L * length(seqs))
out_lines[seq.int(1L, length(out_lines), by = 2L)] <- paste0(">", names_only)
out_lines[seq.int(2L, length(out_lines), by = 2L)] <- seqs

con <- gzfile(dst, open = "wt")
writeLines(out_lines, con)
close(con)

cat("Wrote", dst, "with", length(seqs), "sequences x", unique(nchar(seqs)),
    "positions (", file.info(dst)$size, "bytes ).\n")
