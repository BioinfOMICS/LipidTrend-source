# ============================================================================
# Analysis Script for Fig2AB
# ============================================================================
# Dataset: Mice_Development_2
# Comparison: PM-P21 vs. PM-P0
# Class: TG
# Dimension: 1D
# Characteristic: doubleBond
# Chain: total chain
# Split chain: FALSE
# Weighted: unweighted
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
dataPATH <- "~/data/Mice_Development_2/PM-P21_PM-P0/TG"
outPATH <- "~/results/Mice_Development_2/PM-P21_PM-P0/TG"

# Dataset parameters
split <- FALSE
weight <- FALSE
# Reference group
ref <- "PM-P0"
# Column name
colName <- "Total.DB"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_1D_doubleBond.rds"))

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
data.table::fwrite(res.tab, file.path(outPATH, "totalChain_1D_doubleBond_smoothing_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_doubleBond_plot.pdf"), plot, width = 6.5, height = 4.5)


# End of script

