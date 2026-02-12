# ============================================================================
# Analysis Script for Fig4AB_1D
# ============================================================================
# Dataset: Viral_Infection_1
# Comparison: Day14-p16 vs. Day14-Mock
# Class: DG
# Dimension: 1D
# Characteristic: chainLength
# Chain: total chain
# Split chain: TRUE/FALSE
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
dataPATH <- "~/data/Viral_Infection_1/Day14-p16_Day14-Mock/DG"
outPATH <- "~/results/Viral_Infection_1/Day14-p16_Day14-Mock/DG"

# Dataset parameters
weight <- FALSE
# Reference group
ref <- "Day14-Mock"


# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_1D_chainLength.rds"))

for(split in c(TRUE, FALSE)){
  
  # Column name
  colName <- if(split) "Total.C" else NULL
  
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
  
  if(split){
    # Smoothing results
    even.tab <- LipidTrend::even_chain_result(res)
    odd.tab <- LipidTrend::odd_chain_result(res)
    
    # Save tables
    data.table::fwrite(even.tab, file.path(outPATH, "totalChain_1D_evenChain_smoothing_result.csv"))
    data.table::fwrite(odd.tab, file.path(outPATH, "totalChain_1D_oddChain_smoothing_result.csv"))
    
    # Save plots
    ggplot2::ggsave(file.path(outPATH, "totalChain_1D_evenChain_plot.pdf"), plot$even_result, width = 6.5, height = 4.5)
    ggplot2::ggsave(file.path(outPATH, "totalChain_1D_oddChain_plot.pdf"), plot$odd_result, width = 6.5, height = 4.5)
  }else{
    # Smoothing results
    res.tab <- LipidTrend::result(res)
    
    # Save tables
    data.table::fwrite(res.tab, file.path(outPATH, "totalChain_1D_smoothing_result.csv"))

    # Save plots
    ggplot2::ggsave(file.path(outPATH, "totalChain_1D_plot.pdf"), plot, width = 8, height = 4.5)
  }
  
}


# End of script

