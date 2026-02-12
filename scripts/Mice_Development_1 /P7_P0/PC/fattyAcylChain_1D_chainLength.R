# ============================================================================
# Analysis Script for Fig3A
# ============================================================================
# Dataset: Mice_Development_1 
# Comparison: P7 vs. P0
# Class: PC
# Dimension: 1D
# Characteristic: chainLength
# Chain: fatty acyl chain
# Split chain: TRUE
# Weighted: both
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
dataPATH <- "~/data/Mice_Development_1 /P7_P0/PC"
outPATH <- "~/results/Mice_Development_1 /P7_P0/PC"

# Dataset parameters
split <- TRUE
# Reference group
ref <- "P0"
# Column name
colName <- "C"

# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "fattyAcylChain_1D_chainLength.rds"))

for(weight in c('weighted', 'unweighted')){
  
  abund_weight <- switch(weight, 
                         'weighted'=TRUE, 
                         'unweighted'=FALSE)
  
  # ============================================================================
  # Run Analysis
  # ============================================================================
  
  res <- LipidTrend::analyzeLipidRegion(
    se, ref_group = ref, split_chain = split, chain_col = colName, 
    test = "t.test", abund_weight = abund_weight, permute_time = 100000)
  plot <- LipidTrend::plotRegion1D(res, p_cutoff=0.05)
  
  # ============================================================================
  # Save Results
  # ============================================================================
  
  # Smoothing results
  even.tab <- LipidTrend::even_chain_result(res)
  
  # Save tables
  data.table::fwrite(even.tab, file.path(outPATH, paste0(weight, "_fattyAcylChain_1D_evenChain_smoothing_result.csv")))
  
  # Save plots
  ggplot2::ggsave(file.path(outPATH, paste0(weight, "_fattyAcylChain_1D_evenChain_plot.pdf")), plot$even_result, width = 6.5, height = 4.5)
  
}

# End of script

