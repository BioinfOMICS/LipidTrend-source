# LipidTrend Manuscript – Source Code and Data

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18765955.svg)](https://doi.org/10.5281/zenodo.18765955)
<!-- badges: end -->

This repository contains all source code and processed data required to
reproduce the figures and supplementary results presented in:

**LipidTrend: A Structure-Aware Framework for Detecting Continuous Trends in
Lipidomic Remodeling**

## Overview

LipidTrend is a structure-aware statistical framework for detecting continuous
lipidomic remodeling in chain-length × double-bond structural space.

This repository provides:

* All R scripts used to generate main and supplementary figures
* Processed input data (as described in the manuscript)
* A fully locked computational environment via `renv`
* A fixed random seed, so every figure and table is bit-for-bit reproducible

No new datasets were generated in this study. All lipidomics datasets were
obtained from previously published studies as detailed in the manuscript.

## Software Environment

All analyses were performed under:

* R 4.5.2
* LipidTrend 1.3.1 (shipped as a source tarball in `renv/cellar/`)
* Bioconductor 3.22
* All other package versions pinned in `renv.lock`

To reproduce the computational environment:

```r
install.packages("renv")
renv::restore()
```

`renv::restore()` installs every package at the exact version recorded in
`renv.lock`, including LipidTrend 1.3.1 from `renv/cellar/`. Start R from the
repository root so that `.Rprofile` activates `renv` automatically.

### Recommended BLAS

`renv` pins R package versions, but not the BLAS library that R uses for
matrix arithmetic. Different BLAS implementations sum the same matrix product
in different orders and disagree in the last few bits of the result.

**These results were produced with R's own reference (netlib) BLAS, which is
the recommended configuration for reproducing them.** Check which BLAS your R
is using:

```r
extSoftVersion()[["BLAS"]]
```

A path ending in `libRblas` is the reference BLAS, and is the default for the
official R builds on macOS and Windows — nothing further to do. Many Linux
distributions instead route R through an optimised library. If the value above
names FlexiBLAS, switch to the netlib backend at the start of your session:

```r
flexiblas::flexiblas_switch(flexiblas::flexiblas_load_backend("NETLIB"))
```

Otherwise, select the reference build system-wide before starting R — on
Debian and Ubuntu, `update-alternatives --config libblas.so.3` offers it.

### Notes for Linux and other platforms

`renv.lock` is platform-independent; the `renv/library/` directory is not, and
is deliberately excluded from version control. Each machine builds its own
library from the lockfile.

On Linux, most CRAN packages will be compiled from source unless you point
`renv` at a binary repository. Posit Package Manager serves Linux binaries and
is much faster:

```r
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
renv::restore()
```

Substitute your distribution's codename for `jammy`. LipidTrend itself is a
pure-R package with no compiled code, so it installs from the bundled tarball
on any platform.

## Reproducing the Figures

To reproduce every figure and table:

```bash
Rscript run_all.R
```

This runs all 27 analysis scripts, each in its own R process, writing results
to `results/`, per-script output to `logs/`, and a status table to
`logs/run_summary.tsv`.

To reproduce a single figure, run its script from the repository root:

```r
source("scripts/<dataset>/<comparison>/<lipid class>/<script>.R")
```

All scripts:

* Load data from `data/`
* Set a fixed random seed before each permutation test
* Perform structure-aware smoothing and permutation testing
* Apply BH-FDR correction
* Save figures and tables to `results/`

## Reproducibility

Every script sets `RNGkind()` explicitly and calls `set.seed(1234)`
immediately before each `analyzeLipidRegion()` call. Re-running any script
reproduces its CSV output byte-for-byte. PDF outputs differ only in the
`/CreationDate` and `/ModDate` metadata that R's PDF device embeds at write
time; their rendered content is identical.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the exact environment, the
verification that was carried out, and the changes made relative to the
previously published scripts.

## Methodological Notes

Each analysis follows the workflow described in the manuscript:

1. Species-level statistical testing.
2. Signed log-transformed regional scoring.
3. Gaussian kernel smoothing in structural space.
4. Permutation-based null estimation.
5. Benjamini–Hochberg FDR control.
6. Significant region identification (FDR < 0.05).
7. Regional testing and regional effect size (`regionalTestFC()`).

Supported analysis modes demonstrated in this repository:

* One-dimensional (chain length or unsaturation)
* Two-dimensional (chain length × double bonds)
* Abundance-weighted and unweighted smoothing
* Odd–even chain stratification

## Repository Structure

```
data/        processed input data, one .rds per analysis (26 files)
scripts/     analysis scripts, mirroring the data layout (27 files)
results/     figures and tables produced by run_all.R
logs/        per-script output and run_summary.tsv
renv/        renv infrastructure; renv/cellar/ holds the LipidTrend tarball
docs/        repository diagram and input-data checksums
run_all.R    run every analysis script
```

The diagram below summarizes the relationship between datasets, analysis
modules, and manuscript figures.

![](docs/repository_structure.png)

## Data Sources

All datasets analyzed in this study were obtained from previously published
work and are described in the manuscript and Supplementary Table 1.

No raw data redistribution beyond published material is included.

## Code Availability

The LipidTrend R package is publicly available:

* Bioconductor: [https://bioconductor.org/packages/LipidTrend](https://bioconductor.org/packages/LipidTrend)
* GitHub: [https://github.com/BioinfOMICS/LipidTrend](https://github.com/BioinfOMICS/LipidTrend)

This repository specifically documents the exact scripts and processed inputs
used to generate the figures in the manuscript.

The analyses use LipidTrend 1.3.1, released through both channels above. The
same version is also bundled as a source tarball in `renv/cellar/`, so that
`renv::restore()` always installs precisely the version these results were
produced with — Bioconductor advances its release and devel branches on a
six-month cycle, and pinning the tarball keeps the environment reproducible
well beyond the lifetime of any one branch.

Releases of this repository are numbered after the LipidTrend version whose
results they contain, so that the two can be matched at a glance.
