library(ggrepel)
library(SummarizedExperiment)
library(preprocessCore)
library(tidyverse)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(corrplot)
library(dplyr)
rm(list = ls())
date = '20250626'

dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_prenatal.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

cpm_qnorm_list <- list()

cell_types <- setdiff(unique(colData(dt)$CellType), "None")

# Raw counts
counts <- assay(dt)

zero_count_report <- list()

for (ct in cell_types) {
  message("Processing CellType: ", ct)
  
  idx <- which(colData(dt)$CellType == ct)
  counts_ct <- counts[, idx]
  
  zero_peak_count <- sum(rowSums(counts_ct) == 0)
  message(" - Number of all-zero peaks: ", zero_peak_count)
  
  zero_count_report[[ct]] <- zero_peak_count
  
  nonzero_peaks <- rowSums(counts_ct) > 0
  counts_ct <- counts_ct[nonzero_peaks, ]
  
  cpm_ct <- t(t(counts_ct) / colSums(counts_ct)) * 1e6
  
  # Quantile normalization
  cpm_qnorm_ct <- normalize.quantiles(as.matrix(cpm_ct))
  
  colnames(cpm_qnorm_ct) <- colnames(cpm_ct)
  rownames(cpm_qnorm_ct) <- rownames(cpm_ct)
  
  cpm_qnorm_list[[ct]] <- cpm_qnorm_ct
}

