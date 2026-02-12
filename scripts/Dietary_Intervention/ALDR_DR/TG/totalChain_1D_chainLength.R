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

# Save plots
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot.pdf"), plot1$even_result, width = 6.5, height = 4.5)
ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot_logScale.pdf"), plot2$even_result, width = 6.5, height = 4.5)

# End of script

