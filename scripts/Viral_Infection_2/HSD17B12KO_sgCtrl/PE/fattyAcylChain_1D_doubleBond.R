# ============================================================================
# Analysis Script for supp.Fig3D_1D_db
# ============================================================================
# Dataset: Viral_Infection_2
# Comparison: HSD17B12KO vs. sgCtrl
# Class: PE
# Dimension: 1D
# Characteristic: doubleBond
# Chain: fatty acyl chain
# Split chain: FALSE
# Weighted: weighted
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
dataPATH <- "data/Viral_Infection_2/HSD17B12KO_sgCtrl/PE"
outPATH <- "results/Viral_Infection_2/HSD17B12KO_sgCtrl/PE"
dir.create(outPATH, recursive=TRUE, showWarnings=FALSE)

# Random seed for the permutation test. analyzeLipidRegion() draws its
# permutations from R's global RNG, so a fixed seed -- set immediately before
# every call below -- is what makes this script reproducible.
SEED <- 1234

# Dataset parameters
split <- FALSE
weight <- TRUE
# Reference group
ref <- "sgCtrl"
# Column name
colName <- "DB"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "fattyAcylChain_1D_doubleBond.rds"))

# ============================================================================
# Run Analysis
# ============================================================================

set.seed(SEED)
res <- LipidTrend::analyzeLipidRegion(
    se, ref_group = ref, split_chain = split, chain_col = colName, 
    test = "t.test", abund_weight = weight, permute_time = 100000)
plot <- LipidTrend::plotRegion1D(res, p_cutoff=0.05)

# ============================================================================
# Save Results
# ============================================================================

# Smoothing results
res.tab <- LipidTrend::result(res)

# Save tables
data.table::fwrite(res.tab, file.path(outPATH, "fattyAcylChain_1D_doubleBond_smoothing_result.csv"))

# Region-level paired t-test and fold-change (new in LipidTrend v1.3.1)
reg.tab <- LipidTrend::regionalTestFC(res, p_cutoff=0.05)
data.table::fwrite(reg.tab, file.path(outPATH, "fattyAcylChain_1D_doubleBond_regionalTestFC_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "fattyAcylChain_1D_doubleBond_plot.pdf"), plot, width = 6.5, height = 4.5)


# End of script

