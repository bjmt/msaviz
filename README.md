
# msaviz

<!-- badges: start -->

[![R-CMD-check](https://github.com/bjmt/msaviz/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bjmt/msaviz/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`msaviz` is a small R package for visualizing multiple sequence
alignments (MSAs). It reshapes an alignment into a tidy `tibble` and
turns it into a ggplot2 heatmap, with a `rasterGeom()` wrapper that
keeps file sizes and render times sane even for alignments with millions
of cells.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("bjmt/msaviz")
```

`msaviz` does not depend on any particular alignment reader. Use either
[seqinr](https://cran.r-project.org/package=seqinr) or
[Biostrings](https://bioconductor.org/packages/Biostrings/) to load your
sequences.

## Quick example

``` r
library(msaviz)

aln_path <- system.file(
  "extdata", "hug_ribosomal.fa.gz", package = "msaviz"
)
aln <- seqinr::read.alignment(aln_path, format = "fasta")

alnDF <- msa2DF(aln)
msaHeatmap(alnDF)
```

<img src="man/figures/README-example-1.png" alt="" width="100%" />

By default `msaHeatmap()` colours every cell by whether it matches the
reference sequence (a computed consensus if you don’t pass
`reference =`). Switch to per-letter colouring with `column = "Letter"`.

## Saving plots

`saveHeatmap()` is a `ggsave()` wrapper that picks a reasonable size
based on the number of sequences and axis text size:

``` r
saveHeatmap(msaHeatmap(alnDF), "aln.pdf")
```

## Shell script

A standalone CLI is shipped in `inst/scripts/msaHeatmap.R`. After
installing the package:

``` r
system.file("scripts", "msaHeatmap.R", package = "msaviz")
```

Use it as:

``` sh
Rscript msaHeatmap.R --input aln.fa --output aln.pdf --column Aln
```

See `--help` for all available options.

## Composing with PID, groups, and a dendrogram

`composeMSA()` lines the heatmap up with optional companion plots — a
percent-identity track, a hierarchical-clustering dendrogram, and a
left-side group annotation strip — using `patchwork`:

``` r
d <- msaDendro(aln)
groups <- setNames(rep(c("Bacteria", "Archaea", "Eukarya"),
                       length.out = length(aln$nam)),
                   aln$nam)
group_cols <- c(Bacteria = "#1f77b4", Archaea = "#ff7f0e",
                Eukarya  = "#2ca02c")

composeMSA(
  heatmap = msaHeatmap(alnDF, row.order = attr(d, "order"),
                       legend.pos = "none"),
  top  = msaPID(alnDF),
  left = d,
  groups = groups[attr(d, "order")],
  group.colours = group_cols
)
```

<img src="man/figures/README-compose-1.png" alt="" width="100%" />

## More

The package vignette walks through `Aln` vs `Letter` colouring, custom
palettes, the `raster = TRUE/FALSE` trade-off, `posStats()`, and the
companion plotters:

``` r
vignette("msaviz")
```
