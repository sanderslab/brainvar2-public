#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3: Prenatal vs postnatal comparison of Cicero results.
# Runs for every lineage that has both a "<lineage>_prenatal" and
# "<lineage>_postnatal" group present (in your data: ExN, InN, Micro).
#
# Usage:
#   Rscript step3_cicero_prenatal_postnatal.R --indir /path/to/v3_coacc_cicero --outdir cicero_dev_comparison
#
# For each lineage, produces:
#   <lineage>_edge_overlap_summary.tsv        - n shared / prenatal-only / postnatal-only edges
#   <lineage>_edge_venn_bar.png                - bar chart of the above
#   <lineage>_score_distribution_overlay.png   - coaccess score histograms, prenatal vs postnatal
#   <lineage>_distance_distribution_overlay.png- log10(distance) histograms, prenatal vs postnatal
#   <lineage>_ccan_size_overlay.png            - CCAN size histograms, prenatal vs postnatal
#   <lineage>_pairtype_comparison.png          - pair_type fractions, prenatal vs postnatal
#   <lineage>_shared_edge_score_scatter.png    - for edges present at both timepoints,
#                                                 coaccess score prenatal vs postnatal (concordance)
# Plus a combined summary across all lineages:
#   ALL_lineages_dev_summary.tsv
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--indir", type = "character", default = "."),
  make_option("--outdir", type = "character", default = "cicero_dev_comparison"),
  make_option("--threshold", type = "double", default = 0.2,
              help = "abs(coaccess) threshold for defining an edge [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- discover lineages with both timepoints ----------------------
get_groups <- function(indir) {
  files <- list.files(indir, pattern = "^cicero_.*_connections\\.tsv$", full.names = FALSE)
  gsub("^cicero_|_connections\\.tsv$", "", files)
}
all_groups <- get_groups(opt$indir)
prenatal_groups <- all_groups[grepl("_prenatal$", all_groups)]
lineages <- gsub("_prenatal$", "", prenatal_groups)
lineages <- lineages[paste0(lineages, "_postnatal") %in% all_groups]
cat(sprintf("Found %d lineages with both timepoints: %s\n", length(lineages), paste(lineages, collapse = ", ")))

# ---------------------- helpers ----------------------
peak_midpoint <- function(peak) {
  parts <- str_split_fixed(peak, "_", 3)
  (as.numeric(parts[, 2]) + as.numeric(parts[, 3])) / 2
}
peak_chrom <- function(peak) str_split_fixed(peak, "_", 3)[, 1]

load_connections <- function(indir, group, threshold) {
  f <- file.path(indir, paste0("cicero_", group, "_connections.tsv"))
  df <- fread(f, sep = "\t")
  df <- df[!is.na(coaccess) & abs(coaccess) >= threshold]
  chrom1 <- peak_chrom(df$Peak1); chrom2 <- peak_chrom(df$Peak2)
  mid1 <- peak_midpoint(df$Peak1); mid2 <- peak_midpoint(df$Peak2)
  df$distance <- ifelse(chrom1 == chrom2, abs(mid1 - mid2), NA_real_)
  p1 <- pmin(df$Peak1, df$Peak2); p2 <- pmax(df$Peak1, df$Peak2)
  df$edge_id <- paste(p1, p2, sep = "||")
  df
}

load_ccans <- function(indir, group) {
  f <- file.path(indir, paste0("cicero_", group, "_ccans.tsv"))
  fread(f, sep = "\t")
}

overall_summary <- list()

for (lin in lineages) {
  cat(sprintf("\n=== Lineage: %s ===\n", lin))
  pre_group  <- paste0(lin, "_prenatal")
  post_group <- paste0(lin, "_postnatal")

  conn_pre  <- load_connections(opt$indir, pre_group, opt$threshold)
  conn_post <- load_connections(opt$indir, post_group, opt$threshold)
  ccan_pre  <- load_ccans(opt$indir, pre_group)
  ccan_post <- load_ccans(opt$indir, post_group)

  # ---------------- A) edge overlap (shared / gained / lost) ----------------
  edges_pre  <- unique(conn_pre$edge_id)
  edges_post <- unique(conn_post$edge_id)
  shared     <- intersect(edges_pre, edges_post)
  pre_only   <- setdiff(edges_pre, edges_post)
  post_only  <- setdiff(edges_post, edges_pre)

  edge_summary <- data.table(
    lineage = lin,
    category = c("shared", "prenatal_only", "postnatal_only"),
    n_edges = c(length(shared), length(pre_only), length(post_only))
  )
  fwrite(edge_summary, file.path(opt$outdir, paste0(lin, "_edge_overlap_summary.tsv")), sep = "\t")
  cat(sprintf("  Edges: %d shared, %d prenatal-only, %d postnatal-only\n",
              length(shared), length(pre_only), length(post_only)))

  p_venn <- ggplot(edge_summary, aes(x = category, y = n_edges, fill = category)) +
    geom_col() +
    scale_fill_manual(values = c(shared = "gray50", prenatal_only = "purple", postnatal_only = "darkgreen")) +
    labs(title = paste0(lin, ": connection overlap, prenatal vs postnatal"),
         x = NULL, y = "Number of connections") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(opt$outdir, paste0(lin, "_edge_venn_bar.png")), p_venn, width = 6, height = 5)

  # ---------------- B) score distribution overlay ----------------
  conn_pre$timepoint  <- "prenatal"
  conn_post$timepoint <- "postnatal"
  combined_conn <- rbind(conn_pre[, .(coaccess, distance, pair_type, timepoint)],
                          conn_post[, .(coaccess, distance, pair_type, timepoint)])

  p_score <- ggplot(combined_conn, aes(x = coaccess, fill = timepoint)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, position = "identity", alpha = 0.5) +
    scale_fill_manual(values = c(prenatal = "purple", postnatal = "darkgreen")) +
    labs(title = paste0(lin, ": co-accessibility score, prenatal vs postnatal"),
         x = "Co-accessibility score", y = "Density") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$outdir, paste0(lin, "_score_distribution_overlay.png")), p_score, width = 7, height = 5)

  # ---------------- C) distance distribution overlay ----------------
  dist_df <- combined_conn[!is.na(distance) & distance > 0]
  dist_df[, log10_distance := log10(distance)]
  p_dist <- ggplot(dist_df, aes(x = log10_distance, fill = timepoint)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, position = "identity", alpha = 0.5) +
    scale_fill_manual(values = c(prenatal = "purple", postnatal = "darkgreen")) +
    labs(title = paste0(lin, ": connection distance, prenatal vs postnatal"),
         x = "log10(distance, bp)", y = "Density") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$outdir, paste0(lin, "_distance_distribution_overlay.png")), p_dist, width = 7, height = 5)

  # ---------------- D) CCAN size overlay ----------------
  ccan_pre$timepoint  <- "prenatal"
  ccan_post$timepoint <- "postnatal"
  size_pre  <- ccan_pre[, .N, by = .(CCAN, timepoint)]
  size_post <- ccan_post[, .N, by = .(CCAN, timepoint)]
  size_combined <- rbind(size_pre, size_post)

  p_ccan <- ggplot(size_combined, aes(x = N, fill = timepoint)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, position = "identity", alpha = 0.5) +
    scale_fill_manual(values = c(prenatal = "purple", postnatal = "darkgreen")) +
    labs(title = paste0(lin, ": CCAN size, prenatal vs postnatal"),
         x = "Peaks per CCAN", y = "Density") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$outdir, paste0(lin, "_ccan_size_overlay.png")), p_ccan, width = 7, height = 5)

  # ---------------- E) pair_type composition comparison ----------------
  if ("pair_type" %in% names(combined_conn)) {
    pt <- combined_conn[, .N, by = .(timepoint, pair_type)]
    pt[, frac := N / sum(N), by = timepoint]
    p_pt <- ggplot(pt, aes(x = timepoint, y = frac, fill = pair_type)) +
      geom_col() +
      labs(title = paste0(lin, ": pair type composition"), x = NULL, y = "Fraction of connections") +
      theme_bw(base_size = 11)
    ggsave(file.path(opt$outdir, paste0(lin, "_pairtype_comparison.png")), p_pt, width = 6, height = 5)
  }

  # ---------------- F) shared-edge score concordance ----------------
  if (length(shared) > 0) {
    pre_scores  <- conn_pre[edge_id %in% shared, .(edge_id, coaccess_prenatal = coaccess)]
    post_scores <- conn_post[edge_id %in% shared, .(edge_id, coaccess_postnatal = coaccess)]
    merged <- merge(pre_scores, post_scores, by = "edge_id")
    corr_val <- cor(merged$coaccess_prenatal, merged$coaccess_postnatal, use = "complete.obs")

    p_scatter <- ggplot(merged, aes(x = coaccess_prenatal, y = coaccess_postnatal)) +
      geom_point(alpha = 0.15, size = 0.5, color = "steelblue") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste0(lin, ": shared-connection score concordance (r = ", round(corr_val, 2), ")"),
           x = "Co-accessibility score (prenatal)", y = "Co-accessibility score (postnatal)") +
      theme_bw(base_size = 11)
    ggsave(file.path(opt$outdir, paste0(lin, "_shared_edge_score_scatter.png")), p_scatter, width = 6, height = 6)
  } else {
    corr_val <- NA
  }

  overall_summary[[lin]] <- data.table(
    lineage = lin,
    n_edges_prenatal = length(edges_pre),
    n_edges_postnatal = length(edges_post),
    n_shared = length(shared),
    n_prenatal_only = length(pre_only),
    n_postnatal_only = length(post_only),
    frac_prenatal_retained = length(shared) / length(edges_pre),
    frac_postnatal_novel = length(post_only) / length(edges_post),
    shared_edge_score_corr = corr_val,
    n_ccans_prenatal = ccan_pre[, uniqueN(CCAN)],
    n_ccans_postnatal = ccan_post[, uniqueN(CCAN)],
    mean_ccan_size_prenatal = mean(size_pre$N),
    mean_ccan_size_postnatal = mean(size_post$N)
  )
}

summary_df <- rbindlist(overall_summary)
fwrite(summary_df, file.path(opt$outdir, "ALL_lineages_dev_summary.tsv"), sep = "\t")
cat("\nWrote ALL_lineages_dev_summary.tsv\n")
print(summary_df)

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
