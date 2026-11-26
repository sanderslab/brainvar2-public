rm(list = ls())
library(tidyverse)
library(edgeR)
library("optparse")
library(data.table)
library(stringr)
date = '20250709'

# ##read in data
# exprMat = readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_postnatal.blacklist_filtered.rds")
# 
# mat = assay(exprMat)
# 
# meta <- tibble(
#   sample = colnames(mat),
#   donor = sub("^(.*)_.*\\.Rep.*$", "\\1", colnames(mat)),
#   celltype = sub("^.*_(.*)\\.Rep.*$", "\\1", colnames(mat))
# )
# 
# mat_long <- as.data.frame(mat) %>%
#   rownames_to_column("feature") %>%
#   pivot_longer(-feature, names_to = "sample", values_to = "count") %>%
#   left_join(meta, by = "sample")
# 
# donor_counts <- mat_long %>%
#   group_by(donor, celltype) %>%
#   summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = celltype, values_from = total_count, values_fill = 0)
# 
# print(donor_counts)
# donor_counts = donor_counts %>%
#   dplyr::select(-None, -Progenitor)
# 
# donor_counts <- donor_counts %>%
#   filter(rowSums(across(where(is.numeric))) > 0)
# 
# write.table(donor_counts,
#             '~/Dropbox/BrainVar2_WD/Tables/table.postnatal.donor_count.txt',
#             quote = F, row.names = F, col.names = T, sep = '\t')
# 
# 
# mat_long <- as.data.frame(mat) %>%
#   rownames_to_column("feature") %>%
#   pivot_longer(-feature, names_to = "sample", values_to = "count") %>%
#   mutate(celltype = sub(".*_(.*)\\.Rep.*", "\\1", sample))
# 
# gene_counts <- mat_long %>%
#   group_by(celltype, feature) %>%
#   summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
#   pivot_wider(names_from = celltype, values_from = total, values_fill = 0) %>%
#   column_to_rownames("feature")
# 
# gene_counts = gene_counts %>%
#   dplyr::select(-None, -Progenitor)
# 
# gene_counts <- gene_counts %>%
#   filter(rowSums(across(where(is.numeric))) > 0)
# 
# gene_counts <- gene_counts[rownames(mat), ]
# 
# write.table(gene_counts,
#             '~/Dropbox/BrainVar2_WD/Tables/table.postnatal.peak_count.txt',
#             quote = F, row.names = F, col.names = T, sep = '\t')
# 
# 
# cts = unique(colData(exprMat)$CellType[!(colData(exprMat)$CellType %in% c('None', 'Progenitor'))])
# 
# for(celltype in cts){
#   print(paste('Cell type:', celltype))
#   ## filter cell type
#   filter_exp <- colData(exprMat)$CellType == celltype
#   
#   exprMat_ft <- exprMat[, filter_exp]
#   
#   geneMat = gene_counts
#   donorMat = donor_counts %>% as.data.frame()
#   rownames(donorMat) = donorMat$donor
#   
#   #filter donors with too few cells
#   # ncell=50
#   # any(colData(exprMat_ft)$nCells < ncell)
#   
#   geneMat_ft = geneMat %>%
#     dplyr::select(all_of(celltype))
#   donorMat_ft = donorMat %>%
#     dplyr::select(all_of(celltype))
#   mat_ft = assay(exprMat_ft)
#   
#   colnames(mat_ft) <- sub("^(.*?)_.*$", "\\1", colnames(mat_ft))
#   donorMat_ft <- donorMat_ft[colnames(mat_ft), , drop = FALSE]
#   table(colnames(mat_ft) == rownames(donorMat_ft))
#   x=DGEList(counts=mat_ft,lib.size=as.numeric(donorMat_ft[,1]),samples=rownames(donorMat_ft),genes=geneMat_ft)
#   
#   ##filter genes without expression
#   print(table(rowSums(x$counts==0)==length(rownames(donorMat))))
#   #keep.exprs <- filterByExpr(x)
#   
#   ###MINIMUM UMI COUNTS HARDCODED HERE
#   # keep.exprs <- rowSums(x$counts) >= 15
#   # x <- x[keep.exprs, keep.lib.sizes=FALSE]
#   x <- calcNormFactors(x, method = "TMM")
#   
#   print(x$samples$lib.size)
#   print(x$samples$norm.factors)
#   
#   lcpm <- edgeR::cpm(x, log=TRUE)
#   s_lcpm <- t(scale(t(lcpm)))
#   n <- ncol(lcpm)
#   zvalues <- qnorm(ppoints(n))
#   z <- lcpm
#   for (i in 1:nrow(z)) z[i,] <- zvalues[order(order(z[i,]))]
#   
#   outdir = '~/Dropbox/BrainVar2_WD/Tables/'
#   write.table(lcpm,file=paste0(outdir,"/postnatal.",celltype,"_atac.lcpm"),quote=FALSE)
#   write.table(s_lcpm,file=paste0(outdir,"/postnatal.",celltype,"_scaled_atac.lcpm"),quote=FALSE)
#   write.table(z,file=paste0(outdir,"/postnatal.",celltype,"_atac.lcpm.qnorm"),quote=FALSE)
#   
# }









pca_results <- list()

cell_types = c('ExN', 'InN', 'Glia', 'NNC')

exclude_samples <- list(
  ExN = c('H463DFCm'), # H463DFCm
  Glia = character(0),
  InN = c('H470DFCm', 'H481DFCm'), # H470DFCm
  NNC = character(0),
  Progenitor = character(0)
)


for (ct in cell_types) {
  message("Processing CellType: ", ct)
  
  mat = data.table::fread(paste('~/Dropbox/BrainVar2_WD/Tables/postnatal.', ct, '_atac.lcpm.qnorm',
                                sep = ''))
  mat = mat %>%
    dplyr::select(-V1)
  
  mat_ct <- mat %>% as.data.frame()
  
  if (ct %in% names(exclude_samples) && length(exclude_samples[[ct]]) > 0) {
    exclude_list <- exclude_samples[[ct]]
    mat_ct <- mat_ct[, !colnames(mat_ct) %in% exclude_list]
    message(sprintf("  Excluded %d samples for %s", length(exclude_list), ct))
  }
  
  zero_var <- apply(mat_ct, 1, function(x) var(x) == 0)
  num_zero_var <- sum(zero_var)
  num_total <- nrow(mat_ct)
  num_remaining <- num_total - num_zero_var
  
  message(sprintf("  Total features: %d | Removed 0-variance: %d | Remaining: %d",
                  num_total, num_zero_var, num_remaining))
  
  mat_ct_filtered <- mat_ct[!zero_var, ]
  
  if (nrow(mat_ct_filtered) < 2) {
    message("  (Too few variable features, skipped)")
    next
  }
  
  pca <- prcomp(t(mat_ct_filtered), center = TRUE, scale. = TRUE)
  
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  cum_exp <- cumsum(var_exp)
  
  pca_results[[ct]] <- list(
    pca = pca,
    var_explained = var_exp,
    cum_var_explained = cum_exp
  )
}


pdf(paste('~/Dropbox/BrainVar2_WD/Figures/plot.postnatal.elbow.by_celltype', date, 'pdf',
          sep = '.'), height = 5)
par(mfrow = c(2, 3))

for (ct in names(pca_results)) {
  var_exp <- pca_results[[ct]]$var_explained
  
  plot(var_exp[1:30], type = "b", pch = 19,
       main = paste("Elbow Plot:", ct),
       xlab = "Principal Component",
       ylab = "Variance Explained")
  
  abline(v = which.max(diff(diff(var_exp[1:50]))), col = "darkgreen", lty = 2)
  abline(h = 0.05, col = "gray", lty = 3)
}
dev.off()






pc_use <- list(
  ExN = 5,
  Glia = 5,
  InN = 5,
  NNC = 5,
  Progenitor = 5
)

mg = read.delim('~/Dropbox/BrainVar2_WD/Data/brainvar2.postnatal.metadata.20250625.txt')

library(ggrepel)

plot_list <- list()

for (ct in names(pca_results)) {
  pca_obj <- pca_results[[ct]]$pca
  n_pc <- pc_use[[ct]]
  
  # Variance explained
  var_exp <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2))[1:2] * 100
  x_lab <- sprintf("PC1 (%.1f%%)", var_exp[1])
  y_lab <- sprintf("PC2 (%.1f%%)", var_exp[2])
  
  df.pc <- as.data.frame(pca_obj$x[, 1:min(n_pc, 6)])
  colnames(df.pc) <- paste0("PC", 1:ncol(df.pc))
  df.pc$BioSampleID <- sapply(strsplit(rownames(df.pc), "_"), `[`, 1)
  
  df.pc <- left_join(df.pc, mg %>%
                       distinct(BioSampleID, .keep_all = TRUE),
                     by = "BioSampleID")
  
  # Outlier
  iqr <- IQR(df.pc$PC2, na.rm = TRUE)
  q1 <- quantile(df.pc$PC2, 0.25, na.rm = TRUE)
  q3 <- quantile(df.pc$PC2, 0.75, na.rm = TRUE)
  upper <- q3 + 1.5 * iqr
  lower <- q1 - 1.5 * iqr
  df.pc$outlier <- df.pc$PC2 > upper | df.pc$PC2 < lower
  
  # PC1 vs PC2 by Period
  gg1 <- ggplot(df.pc, aes(PC1, PC2, color = as.factor(Period))) +
    geom_point(size = 2) +
    geom_text_repel(data = df.pc[df.pc$outlier, ],
                    aes(label = paste0(matching_WGS_library, " (", BioSampleID, ")")),
                    size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_minimal() +
    ggtitle(paste0(ct, " | PC1 vs PC2 by period")) +
    theme(axis.title = element_text(size = 8, color = 'black'),
          axis.text = element_text(size = 8, color = 'black'),
          plot.title = element_text(size = 8, color = 'black'),
          legend.title = element_text(size = 8, color = 'black'),
          legend.text = element_text(size = 8, color = 'black'))
  
  # PC1 vs PC2 by Sex
  gg2 <- ggplot(df.pc, aes(PC1, PC2, color = SexForAnalysis)) +
    geom_point(size = 2) +
    geom_text_repel(data = df.pc[df.pc$outlier, ],
                    aes(label = paste0(matching_WGS_library, " (", BioSampleID, ")")),
                    size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_minimal() +
    ggtitle(paste0(ct, " | PC1 vs PC2 by sex")) +
    theme(axis.title = element_text(size = 8, color = 'black'),
          axis.text = element_text(size = 8, color = 'black'),
          plot.title = element_text(size = 8, color = 'black'),
          legend.title = element_text(size = 8, color = 'black'),
          legend.text = element_text(size = 8, color = 'black'))
  
  # PC1 vs PC2 by ancestry
  gg3 <- ggplot(df.pc, aes(PC1, PC2, color = AncestryPCA)) +
    geom_point(size = 2) +
    geom_text_repel(data = df.pc[df.pc$outlier, ],
                    aes(label = paste0(matching_WGS_library, " (", BioSampleID, ")")),
                    size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_minimal() +
    ggtitle(paste0(ct, " | PC1 vs PC2 by ancestry")) +
    theme(axis.title = element_text(size = 8, color = 'black'),
          axis.text = element_text(size = 8, color = 'black'),
          plot.title = element_text(size = 8, color = 'black'),
          legend.title = element_text(size = 8, color = 'black'),
          legend.text = element_text(size = 8, color = 'black'))
  
  plot_list[[paste0(ct, "_PC1vs2_Period")]] <- gg1
  plot_list[[paste0(ct, "_PC1vs2_Sex")]] <- gg2
  plot_list[[paste0(ct, "_PC1vs2_Ancestry")]] <- gg3
  
  write.table(
    df.pc,
    file = paste('~/Dropbox/BrainVar2_WD/Tables/postnatal.PCA', ct, date, 'txt',
                 sep = '.'),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}



out_dir <- "~/Dropbox/BrainVar2_WD/Figures/post_pca_by_celltype"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in names(pca_results)) {
  plot1 <- plot_list[[paste0(ct, "_PC1vs2_Period")]]
  plot2 <- plot_list[[paste0(ct, "_PC1vs2_Sex")]]
  plot3 <- plot_list[[paste0(ct, "_PC1vs2_Ancestry")]]
  
  combined_plot <- plot_grid(plot1, plot2, plot3, ncol = 3)
  
  out_file <- file.path(out_dir, paste0("PCA_", ct, ".", date, ".round2.pdf"))
  
  ggsave(filename = out_file, plot = combined_plot, width = 14, height = 4, units = "in")
}



