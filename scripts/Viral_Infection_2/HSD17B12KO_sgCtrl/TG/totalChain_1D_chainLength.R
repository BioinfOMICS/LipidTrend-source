# ============================================================================
# Analysis Script for Fig1CD
# ============================================================================
# Dataset: Viral_Infection_2
# Comparison: HSD17B12KO vs. sgCtrl
# Class: TG
# Dimension: 1D
# Characteristic: chainLength
# Chain: total chain
# Split chain: TRUE
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
split <- TRUE
weight <- TRUE
# Reference group
ref <- "sgCtrl"
# Column name
colName <- "Total.C"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_1D_chainLength.rds"))

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
even.tab <- LipidTrend::even_chain_result(res)
odd.tab <- LipidTrend::odd_chain_result(res)

# Save tables
data.table::fwrite(even.tab, file.path(outPATH, "totalChain_1D_evenChain_smoothing_result.csv"))
data.table::fwrite(odd.tab, file.path(outPATH, "totalChain_1D_oddChain_smoothing_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot.pdf"), plot$even_result, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_oddChain_plot.pdf"), plot$odd_result, width = 6.5, height = 4.5)


# End of script

