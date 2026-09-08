# ============================================================================
# Analysis Script for supp.Fig4AC
# ============================================================================
# Dataset: Dietary_Intervention
# Comparison: ALDR vs. DR
# Class: TG
# Dimension: 1D
# Characteristic: chainLength
# Chain: total chain
# Split chain: TRUE
# Weighted: unweighted
# ============================================================================
# Generated: 2026-01-16
# ============================================================================

# ============================================================================
# Re-run with LipidTrend v1.3.1
# ============================================================================
# This script is derived from the upstream repository
# (github.com/BioinfOMICS/LipidTrend-source, generated for LipidTrend v1.0.0).
# Changes from upstream, and nothing else:
#   1. Fixed RNG seed. The upstream script set no seed, so its permutation
#      p-values were not reproducible between runs. RNGkind() is pinned and
#      set.seed(SEED) is called immediately before every analyzeLipidRegion().
#   2. Project-relative paths, replacing the upstream path template.
#   3. regionalTestFC() tables written alongside each smoothing result
#      (region-level paired t-test / fold-change, new in v1.3.1).
# Analysis parameters are unchanged.
# ============================================================================


# Load packages
library(LipidTrend)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)

# Pin the RNG algorithm so permutation draws are identical across R versions
# and platforms, rather than inherited from the session default.
RNGkind(kind="Mersenne-Twister", normal.kind="Inversion",
        sample.kind="Rejection")


# ============================================================================
# Configuration
# ============================================================================

# Paths are relative to the project root. Run this script from the root of
# the LipidTrend-source project, so that renv activates and these paths
# resolve:
#     Rscript scripts/<dataset>/<comparison>/<class>/<script>.R
# or, from an R session started at the root:
#     source("scripts/<dataset>/<comparison>/<class>/<script>.R")
if (!dir.exists("data") || !dir.exists("scripts")) {
    stop("Working directory must be the LipidTrend-source project root; got: ",
         getwd())
}
dataPATH <- "data/Dietary_Intervention/ALDR_DR/TG"
outPATH <- "results/Dietary_Intervention/ALDR_DR/TG"
dir.create(outPATH, recursive=TRUE, showWarnings=FALSE)

# Random seed for the permutation test. analyzeLipidRegion() draws its
# permutations from R's global RNG, so a fixed seed -- set immediately before
# every call below -- is what makes this script reproducible.
SEED <- 1234

# Dataset parameters
split <- TRUE
weight <- FALSE
# Reference group
ref <- "DR"
# Column name
colName <- "Total.C"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_1D_chainLength.rds"))

# ============================================================================
# Run Analysis
# ============================================================================

set.seed(SEED)
res <- LipidTrend::analyzeLipidRegion(
    se, ref_group = ref, split_chain = split, chain_col = colName, 
    test = "t.test", abund_weight = weight, permute_time = 100000)
plot1 <- LipidTrend::plotRegion1D(res, p_cutoff=0.05)
plot2 <- LipidTrend::plotRegion1D(res, p_cutoff=0.05, y_scale='log10')

# ============================================================================
# Save Results
# ============================================================================

# Smoothing results
even.tab <- LipidTrend::even_chain_result(res)

# Save tables
data.table::fwrite(even.tab, file.path(outPATH, "totalChain_1D_evenChain_smoothing_result.csv"))

# Region-level paired t-test and fold-change (new in LipidTrend v1.3.1)
reg.tab <- LipidTrend::regionalTestFC(res, p_cutoff=0.05)
data.table::fwrite(reg.tab$even_result, file.path(outPATH, "totalChain_1D_evenChain_regionalTestFC_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot.pdf"), plot1$even_result, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot_logScale.pdf"), plot2$even_result, width = 6.5, height = 4.5)

# End of script

