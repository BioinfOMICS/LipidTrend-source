# ============================================================================
# Analysis Script for Fig1AB
# ============================================================================
# Dataset: HCC
# Comparison: tsne1 vs. tsne0
# Class: TG
# Dimension: 2D
# Characteristic: chainLength and doubleBond
# Chain: total chain
# Split chain: TRUE
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
dataPATH <- "~/data/HCC/tsne1_tsne0/TG"
outPATH <- "~/results/HCC/tsne1_tsne0/TG"

# Dataset parameters
split <- TRUE
weight <- FALSE
# Reference group
ref <- "tsne0"
# Column name
colName <- "Total.C"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_2D.rds"))

# ============================================================================
# Run Analysis
# ============================================================================

res <- LipidTrend::analyzeLipidRegion(
    se, ref_group = ref, split_chain = split, chain_col = colName, 
    test = "t.test", abund_weight = weight, permute_time = 100000)
plot <- LipidTrend::plotRegion2D(res, p_cutoff=0.05)

# ============================================================================
# Save Results
# ============================================================================

# Smoothing results
even.tab <- LipidTrend::even_chain_result(res)
odd.tab <- LipidTrend::odd_chain_result(res)

# Save tables
data.table::fwrite(even.tab, file.path(outPATH, "totalChain_2D_evenChain_smoothing_result.csv"))
data.table::fwrite(odd.tab, file.path(outPATH, "totalChain_2D_oddChain_smoothing_result.csv"))

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_2D_evenChain_plot.pdf"), plot$even_result, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(outPATH, "totalChain_2D_oddChain_plot.pdf"), plot$odd_result, width = 6.5, height = 4.5)


# End of script

