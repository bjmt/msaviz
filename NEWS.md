# msaviz 1.0.0

* First proper release. Public API:
  * `msa2DF()` reshapes a character vector, `seqinr::alignment`, or
    `Biostrings::*MultipleAlignment` into a tidy tibble.
  * `msaHeatmap()` draws the alignment as a ggplot2 heatmap, coloured by
    similarity to a reference sequence (`column = "Aln"`) or by per-letter
    identity (`column = "Letter"`).
  * `saveHeatmap()` is a `ggsave()` wrapper that auto-sizes the output based
    on the number of sequences.
  * `rasterGeom()` rasterizes a ggplot geom layer at draw time to keep large
    heatmaps compact and fast.
* Bundled example slice of the Hug *et al.* (2016) ribosomal-protein
  alignment in `inst/extdata/hug_ribosomal.fa.gz`.
* Standalone CLI in `inst/scripts/msaHeatmap.R`.
* Vignette `vignette("msaviz")` walking through `Aln` vs `Letter` colouring,
  custom palettes, and the rasterization trade-off.
