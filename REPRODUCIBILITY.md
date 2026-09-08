# Reproducibility

This document records how the results in `results/` were produced, what was
verified, and how the analysis differs from the previously published version
of this repository.

## Environment

| Component | Value |
|---|---|
| R | 4.5.2 (2025-10-31), `x86_64-apple-darwin20` |
| BLAS | R reference BLAS (`libRblas.0.dylib`), single-threaded |
| LipidTrend | 1.3.1, from `renv/cellar/LipidTrend_1.3.1.tar.gz` |
| Bioconductor | 3.22 |
| Package versions | pinned in `renv.lock` (122 packages) |
| Random seed | `1234`, set before every `analyzeLipidRegion()` call |

`run_all.R` sets `OMP_NUM_THREADS=1` (and the OpenBLAS/MKL/Accelerate
equivalents) before launching any script, so matrix products are always summed
in the same order rather than one chosen by thread scheduling.

## Changes from the previously published scripts

The scripts in `scripts/` are derived from the previously published version of
this repository, which was generated for LipidTrend v1.0.0. Three changes were
made, and nothing else. Analysis parameters — `radius`, `test`,
`abund_weight`, `permute_time`, `p_cutoff`, `ref_group`, `split_chain` — are
unchanged.

1. **Fixed random seed.** The published scripts set no seed. Because
   `analyzeLipidRegion()` draws its permutations from R's global RNG, every
   run of those scripts was an independent random sample, and borderline
   features could change significance between runs. Each script now pins
   `RNGkind()` and calls `set.seed(1234)` immediately before every
   `analyzeLipidRegion()` call — inside loops as well, so each analysis is
   reproducible independently of loop order.

2. **Project-relative paths.** The published scripts carried a `~/data/...`
   path template that had to be edited by hand. Paths are now relative to the
   repository root, and each script refuses to run from anywhere else.

3. **`regionalTestFC()` output.** New in LipidTrend v1.3.1: a region-level
   paired t-test and fold-change summary for each Increase/Decrease region,
   written as `*_regionalTestFC_result.csv` alongside each
   `*_smoothing_result.csv`. Its columns are `direction`, `n.features`,
   `regional.test.pval`, `regional.FC` and `log2.regional.FC`, where
   `regional.FC` is the **median** fold-change across the significant
   features in that direction — a median rather than a mean, so that a single
   extreme feature cannot dominate the summary of the region.

One small fix was also required. In
`scripts/HCC/tsne1_tsne0/TG/coordinate_correlation.R`, four references read
`cor.res$estimate` / `cor.res$p.value`, where the objects in scope are
`x.cor.res` and `y.cor.res`. They were corrected to the matching `x.`/`y.`
objects. These values are summary statistics that the script does not pass to
either plot, so the figures are unaffected.

## Verification

### 1. Determinism under a fixed seed

`run_all.R` was executed twice from scratch, each script in its own R process:

| Output | Files | Result |
|---|---|---|
| `*_smoothing_result.csv` | 47 | byte-for-byte identical |
| `*_regionalTestFC_result.csv` | 47 | byte-for-byte identical |
| `*.pdf` | 51 | identical content |

No file differed, and no file was missing from either run. The PDFs are not
byte-identical because R's PDF device stamps `/CreationDate` and `/ModDate` at
write time; with those two fields excluded, every PDF is identical, at the
same file size.

### 2. Environment restores from the lockfile

`renv::restore()` was tested into a clean library with an empty renv cache, to
simulate a machine that has never seen these packages. LipidTrend 1.3.1 was
built from the bundled `renv/cellar/` tarball and all 41 dependencies were
installed at their locked versions.

Two corrections were made to the generated lockfile:

* **LipidTrend was recorded as `Source: Bioconductor`.** `renv` infers this
  from the `biocViews` field in the package DESCRIPTION. When these results
  were generated, Bioconductor did not yet carry version 1.3.1 (release
  carried 1.0.0 and devel 1.3.0), so restoring on another machine would have
  silently installed a different version of the package rather than the one
  these results were produced with. The record is now `Source: Cellar`, which
  resolves to the tarball shipped in `renv/cellar/`.

  This remains the correct record even once 1.3.1 is published to
  Bioconductor. Bioconductor advances its release branch twice a year and
  moves superseded versions out of the main repository, so a lockfile that
  points at a Bioconductor version has a limited shelf life. The bundled
  tarball does not: it pins the exact bits this analysis was run against, for
  as long as the repository exists.

* **Six Bioconductor packages were recorded against
  `https://bioc-release.r-universe.dev`.** r-universe serves only current
  builds and is not a version archive. Each of those six versions
  (DelayedArray 0.36.1, GenomicRanges 1.62.1, IRanges 2.44.0, S4Arrays
  1.10.1, SparseArray 1.10.10, XVector 0.50.0) was confirmed present in the
  Bioconductor 3.22 release repository, and the records now point there.

> **Expected `renv::status()` output.** Because of the first correction,
> `renv::status()` always reports one discrepancy:
>
> ```
> The following package(s) are out of sync [lockfile != library]:
> - LipidTrend   [1.3.1: Cellar != Bioconductor]
> ```
>
> This is expected and correct. `renv` infers the installed package's source
> from its `biocViews` field, while the lockfile deliberately records the
> bundled tarball. The other 121 packages are in sync. Do not "fix" it by
> running `renv::snapshot()` — see below.

> **Maintenance note.** `renv::snapshot()` re-derives both of these fields
> from the installed packages and will silently revert them: LipidTrend goes
> back to `Source: Bioconductor`, and those six packages back to r-universe.
> After any future `renv::snapshot()`, re-apply both corrections and confirm
> the lockfile still reads `"Source": "Cellar"` for LipidTrend and
> `"Repository": "Bioconductor 3.22"` for the six packages above.
>
> Installing LipidTrend needs the same care. `renv::install("LipidTrend")`
> resolves the name against the configured repositories and will install the
> Bioconductor version, not the bundled one — it is the version number in the
> lockfile that binds `renv::restore()` to the tarball. To install it
> directly, give the path explicitly:
> `renv::install("./renv/cellar/LipidTrend_1.3.1.tar.gz")`.

## Run outcome

All 27 scripts complete successfully, producing 47 `*_smoothing_result.csv`
files, 47 `*_regionalTestFC_result.csv` files and 51 PDFs in `results/`, with
per-script output and timings in `logs/`.

## Reproducing this from scratch

Clone the repository, then start R with its root as the working directory —
`.Rprofile` activates `renv` on startup — and run:

```r
install.packages("renv")
renv::restore()
source("run_all.R")
```

Data provenance: the 26 `.rds` inputs are unmodified copies of the previously
published `data/` directory. Their SHA-256 checksums are listed in
`docs/data_checksums.txt`.
