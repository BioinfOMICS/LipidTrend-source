# ============================================================================
# Analysis Script for supp.Fig1
# ============================================================================
# Dataset: HCC
# Comparison: tsne1 vs. tsne0
# Class: TG
# Characteristic: chainLength and doubleBond
# ============================================================================
# Generated: 2026-02-13
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
library(tidyverse)
library(data.table)
library(SummarizedExperiment)


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
dataPATH <- "data/HCC/tsne1_tsne0/TG"
outPATH <- "results/HCC/tsne1_tsne0/TG"
dir.create(outPATH, recursive=TRUE, showWarnings=FALSE)


# ============================================================================
# Load Data
# ============================================================================

se <- readRDS(file.path(dataPATH, "totalChain_2D.rds"))
abundance <- SummarizedExperiment::assay(se)
group.info <- SummarizedExperiment::colData(se) %>% as.data.frame()
lipid.char <- SummarizedExperiment::rowData(se) %>% as.data.frame()

# ============================================================================
# Run Analysis
# ============================================================================

t.stat <- apply(abundance, 1, function(x) t.test(x ~ group.info$group)$statistic)
stat.df <- data.frame(x = lipid.char$Total.C,
                      y = lipid.char$Total.DB,
                      t.stat)

x.cor <- NULL
for(ii in unique(stat.df$y)){
    tmp <- stat.df[stat.df$y == ii,]
    if(nrow(tmp) > 1){
        tmp.out <- cbind(tmp$t.stat[-nrow(tmp)],
                         tmp$t.stat[-1])
        x.cor <- rbind(x.cor, tmp.out)
    }
}
x.cor.df <- data.frame(self = x.cor[,1],
                       neighbor = x.cor[,2])
x.cor.res <- cor.test(x.cor.df$self, x.cor.df$neighbor, method = "pearson")
x.r.val <- x.cor.res$estimate
x.p.val <- x.cor.res$p.value


y.cor <- NULL
for(ii in unique(stat.df$x)){
    tmp <- stat.df[stat.df$x == ii,]
    if(nrow(tmp) > 1){
        tmp.out <- cbind(tmp$t.stat[-nrow(tmp)],
                         tmp$t.stat[-1])
        y.cor <- rbind(y.cor, tmp.out)
    }
}
y.cor.df <- data.frame(self = y.cor[,1],
                       neighbor = y.cor[,2])
y.cor.res <- cor.test(y.cor.df$self, y.cor.df$neighbor, method = "pearson")
y.r.val <- y.cor.res$estimate
y.p.val <- y.cor.res$p.value


x.cor.plot <- ggplot(x.cor.df) +
    geom_point(aes(x = self, y = neighbor)) +
    geom_abline(intercept = 0, slope = 1, col = "brown") + 
    scale_x_continuous(limits = c(min(x.cor.df), max(x.cor.df))) +
    scale_y_continuous(limits = c(min(x.cor.df), max(x.cor.df))) +
    theme_bw() +
    labs(title = "x-coordinate") + 
    theme(plot.title = element_text(size=18, hjust = 0.5), 
          axis.title = element_text(size=16), 
          axis.text = element_text(size=14))

y.cor.plot <- ggplot(y.cor.df) +
    geom_point(aes(x = self, y = neighbor)) +
    geom_abline(intercept = 0, slope = 1, col = "brown") + 
    scale_x_continuous(limits = c(min(y.cor.df), max(y.cor.df))) +
    scale_y_continuous(limits = c(min(y.cor.df), max(y.cor.df))) +
    theme_bw() +
    labs(title = "y-coordinate") + 
    theme(plot.title = element_text(size=18, hjust = 0.5), 
          axis.title = element_text(size=16), 
          axis.text = element_text(size=14))


# ============================================================================
# Save Results
# ============================================================================

# Save plots
ggplot2::ggsave(file.path(outPATH, "x_coordinate_correlation.pdf"), x.cor.plot, width = 5.8, height = 5.5)
ggplot2::ggsave(file.path(outPATH, "y_coordinate_correlation.pdf"), y.cor.plot, width = 5.8, height = 5.5)


# End of script

