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

salmon_dir   <- path.expand("~/adar_rinaldi/gene_quantification/salmon")
tx2gene_file <- path.expand("~/adar_rinaldi/tx2gene.tsv")
out_file     <- "salmon_gene_counts_matrix.tsv"

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

# ------------------------------------------------------------------
# Import and roll up to gene level
# ------------------------------------------------------------------

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

counts_matrix <- as.data.frame(txi$counts)
counts_matrix <- cbind(GeneID = rownames(counts_matrix), counts_matrix)

write.table(counts_matrix, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote", nrow(counts_matrix), "genes x", ncol(counts_matrix) - 1,
    "samples to", out_file, "\n")
