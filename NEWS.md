# msaviz 1.1.0

## New features

* `msaPID()` draws a per-position percent-identity track designed to stack
  above [`msaHeatmap()`] with a matching x-axis.
* `msaDendro()` runs `seqinr::dist.alignment()` + `hclust()` on the
  alignment and returns a ggplot dendrogram whose y-axis lines up with
  the heatmap. The leaf order is attached as the `"order"` attribute,
  ready to pass through to `msaHeatmap(row.order = ...)`.
* `composeMSA()` lays out a heatmap with optional `top` (e.g. `msaPID()`)
  and `left` (e.g. `msaDendro()`) companion plots, plus an optional
  `groups` row-annotation strip, using `patchwork`. This finishes the
  group-annotation feature that 1.0.0 advertised but did not implement.
* `consensusSeq()` exposes the per-position consensus calculation as a
  standalone function. `include.gaps = FALSE` ignores gap characters when
  picking the most-common letter.
* `posStats()` returns one row per position with non-gap count, gap count,
  fraction matching the reference, most-common letter, that letter's
  frequency, and Shannon entropy.
* `msa_palette_DNA`, `msa_palette_RNA`, `msa_palette_AA` — drop-in
  `letter.colours =` palettes for the standard alphabets.
* `integer_breaks()` is exported as a tick-position helper for the
  position axis.

## Performance

* `msa2DF()` no longer round-trips through `data.frame` to build the
  per-letter matrix, and the consensus is computed via a vectorized
  `tabulate()` + `max.col()` pass instead of one `table()` per column.
  Result: roughly 10× faster on small/medium alignments. On very large
  alignments the long-format reshape dominates and the speedup is more
  modest (~10%).

## Robustness

* `msa2DF()` now errors with a helpful message on duplicate sequence
  names instead of silently deduping them inside the data.frame
  conversion.
* `rasterGeom()` and `msaHeatmap(raster = TRUE)` no longer fail when
  `getOption("bitmapType")` is `NULL` (common on headless Linux). The
  `png()` device's compiled-in default is used instead.
* `msaHeatmap()` x-axis breaks are now produced by the new
  `integer_breaks()` helper, replacing a private-API call into
  `ggplot2::scale_x_continuous(...)$get_breaks()` and guaranteeing
  integer-valued ticks.

## API changes

* `msaHeatmap()` lost the never-implemented `groups`, `group.box.width`,
  `group.colours`, `group.labels`, and `group.pos` arguments. Pass
  `groups =` and `group.colours =` to `composeMSA()` instead. No code
  that worked under 1.0.0 stops working — those arguments were silently
  ignored.

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
