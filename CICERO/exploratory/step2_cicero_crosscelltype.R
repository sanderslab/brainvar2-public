#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 2: Cross-cell-type comparison of Cicero results.
#
# Usage:
#   Rscript step2_cicero_crosscelltype.R --indir /path/to/v3_coacc_cicero --outdir cicero_crosscelltype
#
# What this does:
#   A) Peak-level CCAN overlap across cell types
#        - binary peak-in-any-CCAN membership matrix
#        - pairwise Jaccard similarity -> clustered heatmap
#        - UpSet plot of peak set overlaps
#   B) Connection-level (edge) overlap across cell types
#        - each connection = sorted(Peak1,Peak2) pair, thresholded
#        - pairwise Jaccard similarity -> clustered heatmap
#        - UpSet plot of edge set overlaps
#
# By default uses ALL 14 groups (so ExN/InN/Micro prenatal + postnatal
# are separate "cell types" here). Pass --postnatal_only to restrict
# to the 11 postnatal groups for a same-timepoint cell-type comparison
# (recommended as the primary read; keep the full version as a sanity
# check that prenatal/postnatal versions of the same lineage cluster
# near each other rather than randomly).
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(UpSetR)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--indir", type = "character", default = "."),
  make_option("--outdir", type = "character", default = "cicero_crosscelltype"),
  make_option("--threshold", type = "double", default = 0.2,
              help = "abs(coaccess) threshold for defining an edge [default %default]"),
  make_option("--postnatal_only", action = "store_true", default = FALSE,
              help = "Restrict to postnatal groups only for a same-timepoint comparison"),
  make_option("--upset_nsets", type = "integer", default = 14,
              help = "Max number of sets to show in UpSet plots [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

get_groups <- function(indir) {
  files <- list.files(indir, pattern = "^cicero_.*_connections\\.tsv$", full.names = FALSE)
  gsub("^cicero_|_connections\\.tsv$", "", files) |> sort()
}

groups <- get_groups(opt$indir)
if (opt$postnatal_only) {
  groups <- groups[grepl("postnatal$", groups)]
}
cat(sprintf("Using %d groups: %s\n", length(groups), paste(groups, collapse = ", ")))

# ==================================================================
# A) PEAK-LEVEL CCAN OVERLAP
# ==================================================================
cat("\n--- Loading CCAN peak sets ---\n")
ccan_peak_sets <- lapply(groups, function(g) {
  f <- file.path(opt$indir, paste0("cicero_", g, "_ccans.tsv"))
  df <- fread(f, sep = "\t")
  unique(df$Peak)
})
names(ccan_peak_sets) <- groups

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

n <- length(groups)
jmat_peaks <- matrix(1, n, n, dimnames = list(groups, groups))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i < j) {
      j_val <- jaccard(ccan_peak_sets[[i]], ccan_peak_sets[[j]])
      jmat_peaks[i, j] <- j_val
      jmat_peaks[j, i] <- j_val
    }
  }
}
write.csv(jmat_peaks, file.path(opt$outdir, "A_ccan_peak_jaccard_matrix.csv"))
cat("Wrote A_ccan_peak_jaccard_matrix.csv\n")

plot_jaccard_heatmap <- function(jmat, outpath, main, low_color) {
  m <- jmat
  diag(m) <- NA  # diagonal (self-similarity = 1) is uninformative; render as black
  max_val <- max(m, na.rm = TRUE)
  breaks <- seq(0, max_val, length.out = 51)
  png(outpath, width = 1400, height = 1200, res = 150)
  pheatmap(m,
           main = main,
           color = colorRampPalette(c("white", low_color))(50),
           breaks = breaks,
           na_col = "black",
           display_numbers = TRUE, number_format = "%.2f", fontsize_number = 6)
  dev.off()
  cat("Wrote", outpath, "\n")
}

plot_jaccard_heatmap(jmat_peaks, file.path(opt$outdir, "A_ccan_peak_jaccard_heatmap.png"),
                      "Jaccard similarity: peaks in CCANs (cell type vs cell type)", "firebrick")

# UpSet plot (limit to a manageable number of sets if requested)
upset_groups <- groups[seq_len(min(opt$upset_nsets, length(groups)))]
png(file.path(opt$outdir, "A_ccan_peak_upset.png"), width = 1800, height = 900, res = 150)
print(upset(fromList(ccan_peak_sets[upset_groups]),
            nsets = length(upset_groups), order.by = "freq",
            mainbar.y.label = "Peaks in intersection", sets.x.label = "Total peaks in CCANs"))
dev.off()
cat("Wrote A_ccan_peak_upset.png\n")

# ==================================================================
# B) CONNECTION (EDGE)-LEVEL OVERLAP
# ==================================================================
cat("\n--- Loading connection edge sets (this can take a while for large groups) ---\n")

load_edges <- function(indir, group, threshold) {
  f <- file.path(indir, paste0("cicero_", group, "_connections.tsv"))
  df <- fread(f, sep = "\t", select = c("Peak1", "Peak2", "coaccess"))
  df <- df[abs(coaccess) >= threshold]
  # canonicalize edge direction so Peak1-Peak2 == Peak2-Peak1
  p1 <- pmin(df$Peak1, df$Peak2)
  p2 <- pmax(df$Peak1, df$Peak2)
  unique(paste(p1, p2, sep = "||"))
}

edge_sets <- lapply(groups, function(g) {
  cat("  edges:", g, "\n")
  load_edges(opt$indir, g, opt$threshold)
})
names(edge_sets) <- groups

jmat_edges <- matrix(1, n, n, dimnames = list(groups, groups))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i < j) {
      j_val <- jaccard(edge_sets[[i]], edge_sets[[j]])
      jmat_edges[i, j] <- j_val
      jmat_edges[j, i] <- j_val
    }
  }
}
write.csv(jmat_edges, file.path(opt$outdir, "B_connection_edge_jaccard_matrix.csv"))
cat("Wrote B_connection_edge_jaccard_matrix.csv\n")

plot_jaccard_heatmap(jmat_edges, file.path(opt$outdir, "B_connection_edge_jaccard_heatmap.png"),
                      "Jaccard similarity: individual connections (cell type vs cell type)", "steelblue")

png(file.path(opt$outdir, "B_connection_edge_upset.png"), width = 1800, height = 900, res = 150)
print(upset(fromList(edge_sets[upset_groups]),
            nsets = length(upset_groups), order.by = "freq",
            mainbar.y.label = "Connections in intersection", sets.x.label = "Total connections"))
dev.off()
cat("Wrote B_connection_edge_upset.png\n")

# ==================================================================
# C) Simple summary: how much is "core"/shared vs cell-type-specific
# ==================================================================
peak_membership_count <- table(unlist(lapply(ccan_peak_sets, unique)))
peak_specificity <- data.table(
  n_groups_present = as.integer(names(table(peak_membership_count))),
  n_peaks = as.integer(table(peak_membership_count))
)
fwrite(peak_specificity, file.path(opt$outdir, "C_peak_sharing_across_groups.tsv"), sep = "\t")
cat("Wrote C_peak_sharing_across_groups.tsv (how many peaks are private vs shared across N groups)\n")

p_spec <- ggplot(peak_specificity, aes(x = n_groups_present, y = n_peaks)) +
  geom_col(fill = "darkorchid") +
  labs(x = "Number of cell types a peak's CCAN membership is seen in",
       y = "Number of peaks",
       title = "Peak sharing across cell types (CCAN membership)") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "C_peak_sharing_barplot.png"), p_spec, width = 8, height = 5)
cat("Wrote C_peak_sharing_barplot.png\n")

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
