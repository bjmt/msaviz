## Submission

This is a new submission of msaviz, version 1.3.0. The package reshapes
multiple sequence alignments into tidy data frames and renders them as
ggplot2 heatmaps or per-sequence lollipop charts.

## R CMD check results

0 errors | 0 warnings | 1 note

The single note is the expected "New submission" message, as this is the
first release of the package on CRAN.

During local checking one further note may appear, "unable to verify current
time", which is an artefact of the checking machine being unable to reach a
time server and does not reflect anything in the package itself. It does not
arise on a machine with a verifiable clock.

## Test environments

* local macOS (darwin), R 4.4.1
* win-builder, R-release and R-devel
* GitHub Actions: ubuntu-latest, windows-latest and macos-latest, across R
  release, R-devel and R oldrel-1

## Notes for the reviewer

One of the suggested packages, Biostrings, comes from Bioconductor rather than
CRAN. It is used only optionally, as a convenience for reading alignment files
and for accepting Biostrings alignment objects as input, and every reference to
it in the tests, the vignette and the command-line scripts is guarded with
requireNamespace() or skip_if_not_installed(), so the package builds, checks
and runs fully without it.

## Downstream dependencies

There are no downstream dependencies, as this is a new package.
