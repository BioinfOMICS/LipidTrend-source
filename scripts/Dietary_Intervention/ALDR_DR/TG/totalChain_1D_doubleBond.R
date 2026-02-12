# ============================================================================
# Analysis Script for supp.Fig4BD
# ============================================================================
# Dataset: Dietary_Intervention
# Comparison: ALDR vs. DR
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
dataPATH <- "~/data/Dietary_Intervention/ALDR_DR/TG"
outPATH <- "~/results/Dietary_Intervention/ALDR_DR/TG"

# Dataset parameters
split <- FALSE
weight <- FALSE
# Reference group
ref <- "DR"
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
plot1 <- LipidTrend::plotRegion1D(res, p_cutoff=0.05)
plot2 <- LipidTrend::plotRegion1D(res, p_cutoff=0.05, y_scale='log10')

# ============================================================================
# Save Results
# ============================================================================

# Smoothing results
res.tab <- LipidTrend::result(res)

# Save tables
data.table::fwrite(res.tab, file.path(outPATH, "totalChain_1D_doubleBond_smoothing_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_doubleBond_plot.pdf"), plot1, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_doubleBond_plot_logScale.pdf"), plot2, width = 6.5, height = 4.5)

# End of script

