# ============================================================================
# Analysis Script for supp.Fig3B_1D_db
# ============================================================================
# Dataset: Viral_Infection_2
# Comparison: HSD17B12KO vs. sgCtrl
# Class: TG
# Dimension: 1D
# Characteristic: doubleBond
# Chain: fatty acyl chain
# Split chain: FALSE
# Weighted: weighted
# ============================================================================
# Generated: 2026-01-16
# ============================================================================


# Load packages
library(LipidTrend)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)


# ============================================================================
# Configuration
# ============================================================================

# IMPORTANT:
# This is a path template for illustration.
# Please modify `~` and the placeholders according to your local environment.
dataPATH <- "~/data/Viral_Infection_2/HSD17B12KO_sgCtrl/TG"
outPATH <- "~/results/Viral_Infection_2/HSD17B12KO_sgCtrl/TG"

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

# Save plots
ggplot2::ggsave(file.path(outPATH, "fattyAcylChain_1D_doubleBond_plot.pdf"), plot, width = 6.5, height = 4.5)


# End of script

