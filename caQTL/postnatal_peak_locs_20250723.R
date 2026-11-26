rm(list = ls())
library(SummarizedExperiment)
library(tidyverse)
library(edgeR)
library(data.table)
library(stringr)

# ----------- READ IN DATA -----------
exprMat = readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_postnatal.blacklist_filtered.rds")

peak_df <- rowData(exprMat) %>%
  as.data.frame() %>%
  mutate(
    ATAC_ID = rownames(.),  # f1, f2, ...
    chr = gsub("chr", "", as.character(seqnames))  # chr1 → 1
  ) %>%
  dplyr::select(ATAC_ID, chr, start, end) %>%
  dplyr::rename(`#chr` = chr) %>%
  dplyr::select(`#chr`, start, end, ATAC_ID)

write.table(
  peak_df,
  file = "~/Dropbox/BrainVar2_WD/Tables/postnatal.peak_locs.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
