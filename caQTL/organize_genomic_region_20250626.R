rm(list = ls())
library(tidyverse)
date = '20250626'

library(GenomicRanges)
library(tidyverse)

cell_types <- c("ExN", "InN", "Glia", "NNC", "Progenitor")

base_path <- "~/Dropbox/BrainVar2_WD/Data/prenatal_peaks/"

annotation_summary <- list()

for (ct in cell_types) {
  gr <- readRDS(file.path(base_path, paste0(ct, "-reproduciblePeaks.gr.rds")))
  
  freq <- as.data.frame(table(gr$peakType)) %>%
    rename(Region = Var1, Count = Freq) %>%
    mutate(CellType = ct,
           Percent = Count / sum(Count) * 100)
  
  annotation_summary[[ct]] <- freq
}

annotation_long <- bind_rows(annotation_summary)

annotation_long$CellType <- factor(annotation_long$CellType, levels = cell_types)


library(ggplot2)
library(RColorBrewer)

pdf(paste('~/Dropbox/BrainVar2_WD/Figures/plot.prenatal.genomic_annotation_prop_celltype', date, 'pdf',
          sep = '.'),
    height = 3, width = 4.5)

ggplot(annotation_long, aes(x = CellType, y = Percent, fill = Region)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Prenatal: Genomic Annotation by Cell Type",
       x = "Cell Type",
       y = "Proportion of Peaks",
       fill = "Region") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(color = 'black', size = 8),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'))

dev.off()










library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)
library(GenomicFeatures)
library(Biostrings)
library(org.Hs.eg.db)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_prenatal.blacklist_filtered.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

valid_cells <- !(colData(dt)$CellType %in% c("None"))
dt_filtered <- dt[, valid_cells]

# Extract peak info (rowData is the same across cells)
peaks <- rowData(dt_filtered) %>% as.data.frame()

peaks_gr <- GRanges(
  seqnames = peaks$seqnames,
  ranges = IRanges(start = peaks$start, end = peaks$end)
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

genes <- genes(txdb)
exons <- exons(txdb)
tss <- promoters(genes, upstream = 0, downstream = 1)
promoterRegion <- c(2000, 100)
promotersGR <- promoters(genes, upstream = promoterRegion[1], downstream = promoterRegion[2])

peakSummits <- resize(peaks_gr, width = 1, fix = "center")

inPromoter <- overlapsAny(peakSummits, promotersGR)
inGene <- overlapsAny(peakSummits, genes)
inExon <- overlapsAny(peakSummits, exons)

peakType <- rep("Distal", length(peaks_gr))
peakType[inGene & !inExon] <- "Intronic"
peakType[inGene & inExon] <- "Exonic"
peakType[inPromoter] <- "Promoter"
mcols(peaks_gr)$peakType <- peakType

peak_df = peaks_gr %>% as.data.frame()

anno_simple_table <- table(peak_df$peakType)
anno_simple_percent <- round(100 * anno_simple_table / sum(anno_simple_table), 2)
anno_simple_percent

anno_df <- as.data.frame(anno_simple_percent)
colnames(anno_df) <- c("Region", "Percent")

pdf(paste('~/Dropbox/BrainVar2_WD/Figures/plot.prenatal.genomic_annotation_prop', date, 'pdf',
          sep = '.'),
    height = 4, width = 5)
ggplot(anno_df, aes(x = reorder(Region, -Percent), y = Percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Prenatal: Genomic Annotation of Peaks",
       x = "Genomic Region",
       y = "Percentage (%)") +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'))
dev.off()












library(GenomicRanges)
library(tidyverse)

cell_types <- c("ExN", "InN", "Glia", "NNC")

base_path <- "~/Dropbox/BrainVar2_WD/Data/postnatal_peaks/"

annotation_summary <- list()

for (ct in cell_types) {
  gr <- readRDS(file.path(base_path, paste0(ct, "-reproduciblePeaks.gr.rds")))
  
  freq <- as.data.frame(table(gr$peakType)) %>%
    rename(Region = Var1, Count = Freq) %>%
    mutate(CellType = ct,
           Percent = Count / sum(Count) * 100)
  
  annotation_summary[[ct]] <- freq
}

annotation_long <- bind_rows(annotation_summary)

annotation_long$CellType <- factor(annotation_long$CellType, levels = cell_types)


library(ggplot2)
library(RColorBrewer)

pdf(paste('~/Dropbox/BrainVar2_WD/Figures/plot.postnatal.genomic_annotation_prop_celltype', date, 'pdf',
          sep = '.'),
    height = 3, width = 4.5)

ggplot(annotation_long, aes(x = CellType, y = Percent, fill = Region)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Postnatal: Genomic Annotation by Cell Type",
       x = "Cell Type",
       y = "Proportion of Peaks",
       fill = "Region") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(color = 'black', size = 8),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 8, color = 'black'))

dev.off()













library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationHub)
library(GenomicFeatures)
library(Biostrings)
library(org.Hs.eg.db)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_postnatal.blacklist_filtered.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

# Filter out 'Progenitor' and 'None' samples
valid_cells <- !(colData(dt)$CellType %in% c("Progenitor", "None"))
dt_filtered <- dt[, valid_cells]

# Extract peak info (rowData is the same across cells)
peaks <- rowData(dt_filtered) %>% as.data.frame()

peaks_gr <- GRanges(
  seqnames = peaks$seqnames,
  ranges = IRanges(start = peaks$start, end = peaks$end)
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

genes <- genes(txdb)
exons <- exons(txdb)
tss <- promoters(genes, upstream = 0, downstream = 1)
promoterRegion <- c(2000, 100)
promotersGR <- promoters(genes, upstream = promoterRegion[1], downstream = promoterRegion[2])

peakSummits <- resize(peaks_gr, width = 1, fix = "center")

inPromoter <- overlapsAny(peakSummits, promotersGR)
inGene <- overlapsAny(peakSummits, genes)
inExon <- overlapsAny(peakSummits, exons)

peakType <- rep("Distal", length(peaks_gr))
peakType[inGene & !inExon] <- "Intronic"
peakType[inGene & inExon] <- "Exonic"
peakType[inPromoter] <- "Promoter"
mcols(peaks_gr)$peakType <- peakType

peak_df = peaks_gr %>% as.data.frame()



anno_simple_table <- table(peak_df$peakType)
nrow(peak_df)
anno_simple_percent <- round(100 * anno_simple_table / sum(anno_simple_table), 2)
anno_simple_percent

anno_df <- as.data.frame(anno_simple_percent)
colnames(anno_df) <- c("Region", "Percent")

pdf(paste('~/Dropbox/BrainVar2_WD/Figures/plot.postnatal.genomic_annotation_prop', date, 'pdf',
          sep = '.'),
    height = 4, width = 5)
ggplot(anno_df, aes(x = reorder(Region, -Percent), y = Percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Postnatal: Genomic Annotation of Peaks",
       x = "Genomic Region",
       y = "Percentage (%)") +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 8, color = 'black'),
        plot.title = element_text(size = 10, color = 'black'))
dev.off()


