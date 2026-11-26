rm(list = ls())
library(rtracklayer)
library(ggrepel)
library(SummarizedExperiment)
library(preprocessCore)
library(tidyverse)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(corrplot)
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)

dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_prenatal.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

peaks <- rowData(dt) %>% as.data.frame()

peaks_gr <- GRanges(
  seqnames = peaks$seqnames,
  ranges = IRanges(start = peaks$start, end = peaks$end)
)

blacklist <- import.bed("~/Dropbox/Resources/hg38-blacklist.v2.bed.gz")

# Remove all overlapping regions
filtered_peaks_gr <- peaks_gr[!peaks_gr %over% blacklist]

cat("Before filtering:", length(peaks_gr), "\n")
cat("After filtering:", length(filtered_peaks_gr), "\n")

keep_idx <- which(!peaks_gr %over% blacklist)

dt_filtered <- dt[keep_idx, ]

saveRDS(dt_filtered, "~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_prenatal.blacklist_filtered.rds")












dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_postnatal.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

peaks <- rowData(dt) %>% as.data.frame()

peaks_gr <- GRanges(
  seqnames = peaks$seqnames,
  ranges = IRanges(start = peaks$start, end = peaks$end)
)

blacklist <- import.bed("~/Dropbox/Resources/hg38-blacklist.v2.bed.gz")

# Remove all overlapping regions
filtered_peaks_gr <- peaks_gr[!peaks_gr %over% blacklist]

cat("Before filtering:", length(peaks_gr), "\n")
cat("After filtering:", length(filtered_peaks_gr), "\n")

keep_idx <- which(!peaks_gr %over% blacklist)

dt_filtered <- dt[keep_idx, ]

saveRDS(dt_filtered, "~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_postnatal.blacklist_filtered.rds")
