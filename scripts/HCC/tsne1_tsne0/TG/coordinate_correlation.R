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


# Load packages
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
x.r.val <- cor.res$estimate
x.p.val <- cor.res$p.value


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
y.r.val <- cor.res$estimate
y.p.val <- cor.res$p.value


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

