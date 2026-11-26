#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Build a gene-level counts matrix from salmon quant.sf files.
# Uses tximport so transcript-length differences across samples are
# handled correctly (this is what DESeq2 expects downstream, but
# this script only writes the matrix - no DESeq2 code here).
#
# Requires: tximport (BiocManager::install("tximport"))
# ------------------------------------------------------------------

if (!requireNamespace("tximport", quietly = TRUE)) {
    stop("Install tximport first: BiocManager::install('tximport')")
}
library(tximport)

# ------------------------------------------------------------------
# EDIT THESE
# ------------------------------------------------------------------

salmon_dir   <- path.expand("~/adar_rinaldi/salmon_fastq_quant")
tx2gene_file <- path.expand("~/adar_rinaldi/tx2gene.tsv")
gene2symbol_file <- path.expand("~/adar_rinaldi/gene2symbol.tsv")   # <-- from make_gene2symbol.sh
out_counts_file <- "salmon_fastq_gene_counts_matrix.tsv"
out_tpm_file    <- "salmon_fastq_gene_TPM_matrix.tsv"

# ------------------------------------------------------------------
# Find every sample subfolder that has a quant.sf
# ------------------------------------------------------------------

sample_dirs <- list.dirs(salmon_dir, recursive = FALSE)
samples     <- basename(sample_dirs)
files       <- file.path(sample_dirs, "quant.sf")
names(files) <- samples

files <- files[file.exists(files)]
cat("Found", length(files), "samples:\n")
print(names(files))

# ------------------------------------------------------------------
# tx2gene mapping (from make_tx2gene.sh)
# ------------------------------------------------------------------

tx2gene <- read.table(tx2gene_file, sep = "\t", header = FALSE,
                       col.names = c("TXNAME", "GENEID"))

gene2symbol <- read.table(gene2symbol_file, sep = "\t", header = FALSE,
                           col.names = c("GeneID", "gene_symbol"))

# ------------------------------------------------------------------
# Import and roll up to gene level
# ------------------------------------------------------------------

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# ------------------------------------------------------------------
# Raw (estimated) counts matrix
# ------------------------------------------------------------------

counts_matrix <- as.data.frame(txi$counts)
counts_matrix <- cbind(GeneID = rownames(counts_matrix), counts_matrix)
counts_matrix <- merge(gene2symbol, counts_matrix, by = "GeneID", all.y = TRUE)
write.table(counts_matrix, out_counts_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote counts matrix:", nrow(counts_matrix), "genes x", ncol(counts_matrix) - 1,
    "samples to", out_counts_file, "\n")

# ------------------------------------------------------------------
# TPM matrix (gene-level, transcript TPMs summed by tximport's abundance)
# ------------------------------------------------------------------

tpm_matrix <- as.data.frame(txi$abundance)
tpm_matrix <- cbind(GeneID = rownames(tpm_matrix), tpm_matrix)
tpm_matrix <- merge(gene2symbol, tpm_matrix, by = "GeneID", all.y = TRUE)
write.table(tpm_matrix, out_tpm_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote TPM matrix:   ", nrow(tpm_matrix), "genes x", ncol(tpm_matrix) - 1,
    "samples to", out_tpm_file, "\n")
