#########################################################################
# Bulk RNA QC                                                           #
# Edited by Narjes, based on Emilie's code for Bulk RNA Quality Control #
# Version 1.0 - 27 March 2025                                           #
#########################################################################

# Set the random seed for reproducibility
set.seed(1)

# Install necessary libraries (if not already installed)
required_libraries <- c("data.table", "RNAseqQC", "DESeq2", "ensembldb", "dplyr", "ggplot2", 
                        "purrr", "tidyr", "tibble", "magrittr", "tidyverse", "limma", "edgeR", 
                        "pcaMethods", "pheatmap", "sva")

# Install missing libraries
for (pkg in required_libraries) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Define paths for the input and output files
input_count_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/bulkRNA_qc/250401_htseq_counts_BrainVar_v2.tsv"  # Update to your actual path
library_metadata_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/Meta-data-analysis/brainvar2_library_metadata_v1.3.tsv"  # Update to your actual path
input_gtf_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/bulkRNA_qc/gencode.v41.annotation.gtf"  # Update to your actual path
cpm_output_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/bulkRNA_qc/brainvar_bulk_RNA_CPM_gene_count_matrix.txt"  # Update to your actual path
tpm_output_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/bulkRNA_qc/brainvar_bulk_RNA_TPM_gene_count_matrix.txt"  # Update to your actual path
output_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/bulkRNA_qc/"  # Update to your actual output path
meta_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/Meta-data-analysis/brainvar2_subject_metadata_v1.3.tsv"  # Subject metadata path

# Number of samples: 262 - 5 samples without subject id
n_samples <- 262 - 5

# Read count file and library metadata
raw_gene_mapped_matrix <- fread(input_count_file, data.table = FALSE)
library_metadata <- fread(library_metadata_file, data.table = FALSE)

# Filter library metadata to only include bulk RNA libraries
library_metadata <- library_metadata[library_metadata$LibraryModality == 'bulkRNA', ]
library_metadata <- library_metadata[!(library_metadata$SubjectID == ''), ]

# Modify column names (RNA library ids) according to the pattern in the metadata file
modify_column_name <- function(col_name) {
  if (grepl("^HSB", col_name)) {
    # If it starts with 'HSB', remove 'HSB' and add '_R1'
    return(paste0(sub("^HSB", "", col_name), "_R1"))
  } else if (grepl("^[0-9]", col_name)) {
    # If it starts with a number, add '_R2'
    return(paste0(col_name, "_R2"))
  } else if (grepl("^H", col_name)) {
    # If it starts with 'H' but not 'HSB', remove 'H' and add '_R2'
    return(paste0(sub("^H", "", col_name), "_R2"))
  } else {
    # For other column names, leave unchanged
    return(col_name)
  }
}

# Apply the modify_column_name function to convert columns into library id format
colnames(raw_gene_mapped_matrix) <- sapply(colnames(raw_gene_mapped_matrix), modify_column_name)

# Ensure the first column is "gene" (this column contains Ensembl gene IDs)
colnames(raw_gene_mapped_matrix)[1] <- "gene"

# Read GTF file to get gene annotations (Ensembl to Gene Symbol mapping)
gtf_data <- fread(input_gtf_file, sep = "\t", header = FALSE, quote = "")

# Extract relevant gene information from GTF (Ensembl to gene symbol)
gtf_data <- gtf_data[V3 == "gene", .(
  chr = V1,
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_name = sub('.*gene_name "([^"]+)".*', '\\1', V9),
  gene_type = sub('.*gene_type "([^"]+)".*', '\\1', V9),
  start = V4,
  end = V5,
  gene_length = V5 - V4 + 1
)]

# Remove rows with NA gene_id
gtf_data <- gtf_data[gene_id != "NA"]

# Map Ensembl gene IDs to gene symbols using GTF data
# We will merge the raw counts with gene symbols from the GTF file
genes <- as.data.frame(raw_gene_mapped_matrix$gene)
colnames(genes)[1] <- "gene_id"
gene_symbols <- merge(genes, gtf_data, by = "gene_id", all.x = TRUE)

# Rename the new column appropriately
raw_gene_mapped_matrix$gene <- gene_symbols$gene_name

# Define CPM calculation function
cpm <- function(df) {
  # Calculate counts per million (CPM)
  total_counts <- colSums(df[, -1], na.rm = TRUE)
  cp <- as.data.frame(t(t(df[, -1]) / total_counts))
  cpm <- cp * 1e6
  final_cpm <- cbind(df[, 1], cpm)
  colnames(final_cpm)[1] <- "gene"
  return(final_cpm)
}

# Calculate CPM
CPM <- cpm(raw_gene_mapped_matrix)

# Save CPM results to a file
fwrite(CPM, cpm_output_path, col.names = TRUE, row.names = FALSE, sep = "\t")

# Add gene lengths to the raw counts data
gene_lengths <- merge(genes, gtf_data, by = "gene_id", all.x = TRUE)

# Define TPM calculation function
tpm <- function(df, lengths) {
  # Calculate transcripts per million (TPM)
  rpk <- (df[, -1]) / lengths
  m_scaling <- colSums(rpk, na.rm = TRUE) / 1000000
  print(any(is.na(rpk)))
  
  tpm <- as.data.frame(t(t(rpk) / m_scaling))
  final_tpm <- cbind(df[, 1], tpm)
  colnames(final_tpm)[1] <- "gene"
  return(final_tpm)
}

# Calculate TPM
TPM <- tpm(raw_gene_mapped_matrix, gene_lengths$gene_length)

# Save TPM results to a file
fwrite(TPM, tpm_output_path, col.names = TRUE, row.names = FALSE, sep = "\t")

# Perform QC Analysis
raw_gene_mapped_matrix$gene <- gene_symbols$gene_id

# Add gene annotations to raw counts
qc_table <- merge(raw_gene_mapped_matrix, gtf_data, by.y = "gene_id", by.x = "gene")

# Calculate total counts for each sample
all_counts <- colSums(qc_table[, 2:n_samples])

# Subset to lncRNA and protein coding genes
lp <- qc_table[which(qc_table$gene_type == "lncRNA" | qc_table$gene_type == "protein_coding"), ]
lp_counts <- colSums(lp[, 2:n_samples])
lp_to_all <- lp_counts / all_counts

# Create histogram of protein-coding/lncRNA to total gene mapped reads ratio
pdf(paste0(output_path, "bulk_brainvar_ratio_of_protein_coding_and_lncRNA_gene_mapped_reads_to_total_reads.pdf"), height = 8, width = 8)
hist(t(lp_to_all), col = "cornflowerblue", breaks = 10, main = "Protein-coding/lncRNA to \n Total Gene Mapped Reads Ratio", xlab = "Ratio of Protein-coding/lncRNA to Total Gene Mapped Reads", ylab = "Frequency")
dev.off()

# Set expression cutoff for TPM (e.g., 0.5)
expression_cutoff <- 0.5
TPM$gene <- gene_symbols$gene_id
tpm_table <- merge(TPM, gtf_data, by.y = "gene_id", by.x = "gene")
TPM$gene <- gene_symbols$gene_name

# Calculate fraction of expressed genes
expressed <- lapply(tpm_table[, 2:n_samples], function(x) { sum(x > expression_cutoff) })
fraction_expressed <- t(as.data.frame(expressed)) / nrow(qc_table)

# Create a data frame with lp_to_all and fraction_expressed
df <- as.data.frame(cbind(lp_to_all, fraction_expressed))
colnames(df)[2] <- "fraction_expressed"

# Create a scatter plot of lp_to_all vs fraction_expressed
pdf(paste0(output_path, "bulk_brainvar_ratio_of_protein_coding_and_lncRNA_gene_mapped_reads_to_total_reads_vs_fraction_of_expressed_genes.pdf"), height = 8, width = 8)
plot(df$lp_to_all, df$fraction_expressed, pch = 19, col = "cornflowerblue", main = "Protein-coding/lncRNA to \n Total Gene Mapped Reads Ratio vs Fraction of Expressed Genes", xlab = "Ratio of Protein-coding/lncRNA to Total Gene Mapped Reads", ylab = "Fraction of Expressed Genes")
dev.off()

# Analyze haemoglobin genes (if listed on the HGNC website)
haemoglobin_genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

# Get rRNA genes
rRNA_genes <- gtf_data$gene_name[which(gtf_data$gene_type == "rRNA")]

# Calculate haemoglobin to protein-coding/lncRNA ratio
haemoglobin <- qc_table[which(qc_table$gene_name %in% haemoglobin_genes), ]
haem_counts <- colSums(haemoglobin[, 2:n_samples])
haem_to_lp <- haem_counts / lp_counts

# Create a histogram of haemoglobin to protein-coding/lncRNA ratio
df <- as.data.frame(haem_to_lp)
colnames(df) <- "haem_to_lp"
pdf(paste0(output_path, "bulk_brainvar_ratio_of_rRNA_gene_mapped_reads_to_protein_coding_and_lncRNA_mapped_reads_hemo.pdf"), height = 5, width = 5)
hist(df$haem_to_lp, col = "cornflowerblue", breaks = 10, main = "Haemoglobin to Protein-coding/lncRNA Ratio", xlab = "Ratio of Haemoglobin to Protein-coding/lncRNA", ylab = "Frequency")
dev.off()

# Calculate rRNA to protein-coding/lncRNA ratio
rRNA <- qc_table[which(qc_table$gene_name %in% rRNA_genes), ]
rRNA_counts <- colSums(rRNA[, 2:n_samples])
rRNA_to_lp <- rRNA_counts / lp_counts

# Create a histogram of rRNA to protein-coding/lncRNA ratio
df <- as.data.frame(rRNA_to_lp)
colnames(df) <- "rRNA_to_lp"
pdf(paste0(output_path, "bulk_brainvar_ratio_of_rRNA_gene_mapped_reads_to_protein_coding_and_lncRNA_mapped_reads.pdf"), height = 5, width = 5)
hist(df$rRNA_to_lp, col = "cornflowerblue", breaks = 10, main = "rRNA to Protein-coding/lncRNA Ratio", xlab = "Ratio of rRNA to Protein-coding/lncRNA", ylab = "Frequency")
dev.off()

# Analyze Y chromosome genes; do not use all Y chrom genes but the below selected ones
Y_chromosome <- tpm_table[which(tpm_table$chr == "chrY" & (tpm_table$gene_type == "lncRNA" | tpm_table$gene_type == "protein_coding")), ]
Y_mean_exp <- colMeans(Y_chromosome[, 2:n_samples])

# Define Y chromosome genes
Y_only <- c("DAZ1", "DAZ2", "DAZ3", "DAZ4", "SRY", "TSPY1", "RBMY1A1", "RBMY1B", "RBMY1D", "RBMY1E", "RBMY1F", "RBMY1J", "CDY1", "CDY2A", "CDY2B", "USP9Y", "VCY", "BPY2", "PRY", "ZFY", "EIF1AY", "RPS4Y1", "RPS4Y2", "AMELY", "TMSB4Y", "XKRY")

# Calculate XIST expression
XIST <- tpm_table[which(tpm_table$gene_name == "XIST"), ]

# Calculate Y chromosome gene expression
Y_genes <- tpm_table[which(tpm_table$gene_name %in% Y_only), ]
Y_genes_exp <- colMeans(Y_genes[, 2:n_samples])

# Create a data frame with XIST and Y chromosome gene expression
df <- as.data.frame(cbind(t(XIST[2:n_samples]), Y_genes_exp))
df$LibraryID <- rownames(df)
colnames(df)[1] <- "XIST"

# Merge with library metadata
merged_df <- merge(df, library_metadata[, c('LibraryID', 'SexForAnalysis')], by.x = "LibraryID", by.y = "LibraryID", all.x = TRUE)

# Create a scatter plot of XIST vs Y chromosome gene expression
pdf(paste0(output_path, "bulk_brainvar_mean_Y_vs_XIST_TPM.pdf"), height = 5, width = 5)
ggplot(merged_df, aes(x = XIST, y = Y_genes_exp, color = SexForAnalysis)) +
  geom_point() +
  scale_color_manual(values = c("XY" = "blue", "XX" = "pink")) +
  theme_minimal() +
  labs(title = "XIST vs Y Chromosome Gene Expression", x = "XIST Expression (TPM)", y = "Y Chromosome Gene Expression (TPM)")
dev.off()

# Calculate QC metrics
qc_metrics <- cbind(df, all_counts, lp_counts, lp_to_all, rRNA_counts, rRNA_to_lp, haem_counts, haem_to_lp, fraction_expressed)
qc_metrics$sample <- rownames(qc_metrics)
qc_metrics <- qc_metrics[c(ncol(qc_metrics), 1:(ncol(qc_metrics) - 1))]

# Write the result to a file
fwrite(qc_metrics, file = paste0(output_path, "brainvar_bulk_RNA_per_sample_qc_metrics.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Further QC plots using RNAseq QC package

count_mat <- raw_gene_mapped_matrix  # Input raw count matrix
meta <- read_delim(meta_path)  # Read metadata
meta <- meta[!is.na(meta$matching_bulkRNA_library), ]  # Remove rows with NA in matching_bulkRNA_library column
rownames(count_mat) <- count_mat$gene
count_mat <- count_mat[, !colnames(count_mat) %in% "gene"]
sample_to_select<-intersect(meta$matching_bulkRNA_library, colnames(count_mat))
meta_selected <- meta[meta$matching_bulkRNA_library %in% sample_to_select, ]
count_mat<-count_mat[,sample_to_select]
meta_selected<-data.frame(meta_selected)

# Create a DESeq2 object for analysis (need for some functions)
dds <- make_dds(counts = count_mat, metadata = meta_selected, ah_record = "AH89426")

# Plot total counts across samples (raw counts)
pdf(paste0(output_path, "bulkRNA_brainvar_total_counts.pdf"), height = 5, width = 5)
plot_total_counts(dds)  # Title the plot as Total Counts
dev.off()

# Plot library complexity (sequencing depth across samples)
pdf(paste0(output_path, "bulkRNA_brainvar_library_complexity.pdf"), height = 5, width = 5)
plot_library_complexity(dds)  # Title the plot as Library Complexity
dev.off()

# Plot gene detection rate (proportion of genes detected across samples)
pdf(paste0(output_path, "bulkRNA_brainvar_gene_detection.pdf"), height = 5, width = 5)
plot_gene_detection(dds)  # Title the plot as Gene Detection Rate
dev.off()

# Plot gene biotypes (categorize genes by their type, e.g., protein-coding)
pdf(paste0(output_path, "bulkRNA_brainvar_biotypes.pdf"), height = 10, width = 10)
plot_biotypes(dds)  # Title the plot as Gene Biotypes Distribution
dev.off()

# Filter genes based on minimum count threshold 
dds <- filter_genes(dds, min_count = 5, min_rep = 20)
vsd <- vst(dds)

# Plot Mean vs Standard Deviation (to assess variance across genes)
pdf(paste0(output_path, "bulkRNA_brainvar_mean_sd_plot.pdf"), height = 5, width = 5)
mean_sd_plot(vsd)  # Title the plot accordingly
dev.off()

# Plot sample clustering (Euclidean distance) based on gene expression
pdf(paste0(output_path, "bulkRNA_brainvar_sample_clustering.pdf"), height = 10, width = 10)
plot_sample_clustering(vsd, anno_vars = c("Period"), distance = "euclidean")  # Annotate by sex for clustering
dev.off()

normalized_counts<-CPM[,2:n_samples]

# PCA for normalized counts (already using CPM)
pca_result <- prcomp(t(normalized_counts))  # Transpose so samples are in rows

# Plot PCA for normalized counts (using CPM)
pca_df <- data.frame(PC1 = pca_result$x[, 1], 
                     PC2 = pca_result$x[, 2], 
                     Sample = colnames(normalized_counts))

pdf(paste0(output_path, "bulkRNA_brainvar_pca_normalized_counts_CPM.pdf"), height = 10, width = 10)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) + 
  geom_point() + 
  geom_text(vjust = 1.5, hjust = 1.5) + 
  theme_minimal() + 
  labs(title = "PCA Plot of Normalized Counts (CPM)", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)")  # PCA of CPM normalized counts
dev.off()

#------------------
# Run PCA
pca_result <- prcomp(t(normalized_counts))

# Variance explained
explained_var <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Build PCA dataframe
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = colnames(normalized_counts)
)

metadata <- as.data.frame(meta_path)

# Optional: add grouping info if you have a metadata table
pca_df <- left_join(
  pca_df,
  meta,
  by = c("Sample" = "matching_bulkRNA_library")
)

# Save PDF
pdf(paste0(output_path, "bulkRNA_brainvar_pca_normalized_counts_CPM.pdf"), height = 6, width = 7)

library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 4, aes(color = Sample)) +  # replace Sample with Group if available
  geom_text(vjust = -1, size = 3.5, show.legend = FALSE) +
  theme_classic(base_size = 16) +
  labs(
    title = "PCA of Bulk RNA-seq (CPM normalized counts)",
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  ) +
  scale_color_brewer(palette = "Set2")  # nice color scheme for publication

dev.off()



















# Step 4: Heatmap of correlations between samples using CPM
cor_matrix <- cor(normalized_counts)

# Generate a heatmap of sample correlations
pdf(paste0(output_path, "bulkRNA_brainvar_correlation_heatmap.pdf"), height = 40, width = 40)
pheatmap(cor_matrix, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete", 
         main = "Sample Correlation Heatmap")  # Title for heatmap
dev.off()

# Batch correction using ComBat (if needed, based on sample batch information)
batch_variable <- ifelse(grepl("R1", colnames(normalized_counts)), "brainvar1", "brainvar2")
batch_variable <- factor(batch_variable)

# Apply ComBat batch correction (optional, depending on if batch effect is present)
corrected_counts <- ComBat(dat = normalized_counts, batch = batch_variable, mod = NULL, par.prior = TRUE, prior.plots = FALSE)

# Before correction PCA (using corrected CPM data)
pca_before <- prcomp(t(normalized_counts))

# Calculate percentage variance explained
explained_var_before <- round(100 * summary(pca_before)$importance[2, 1:2], 1)

# Plot PCA before batch correction
pdf(paste0(output_path, "bulkRNA_brainvar_pca_before_correction.pdf"), height = 5, width = 5)
par(cex.lab = 1, cex.axis = 1, cex.main = 1.2)
plot(pca_before$x[, 1:2], col = batch_variable, pch = 16,
     main = "PCA Before Batch Correction",
     xlab = paste0("PC1 (", explained_var_before[1], "%)"),
     ylab = paste0("PC2 (", explained_var_before[2], "%)"))
text(pca_before$x[, 1], pca_before$x[, 2], labels = colnames(normalized_counts), cex = 0.5, pos = 2)
dev.off()

# After correction PCA (using corrected CPM data)
pca_after <- prcomp(t(corrected_counts))

# Calculate percentage variance explained
explained_var_after <- round(100 * summary(pca_after)$importance[2, 1:2], 1)

# Plot PCA after batch correction
pdf(paste0(output_path, "bulkRNA_brainvar_pca_after_correction.pdf"), height = 5, width = 5)
par(cex.lab = 1, cex.axis = 1, cex.main = 1.2)
plot(pca_after$x[, 1:2], col = batch_variable, pch = 16,
     main = "PCA After Batch Correction",
     xlab = paste0("PC1 (", explained_var_after[1], "%)"),
     ylab = paste0("PC2 (", explained_var_after[2], "%)"))
text(pca_after$x[, 1], pca_after$x[, 2], labels = colnames(corrected_counts), cex = 0.8, pos = 4)
dev.off()

cor_matrix$Period<-library_metadata[library_metadata$LibraryID %in% rownames(cor_matrix),c('Period')]
