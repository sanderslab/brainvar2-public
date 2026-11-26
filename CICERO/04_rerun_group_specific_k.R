#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 5: Rerun Cicero for under-capped groups with k rescaled to
# maintain a consistent metacell-overlap ratio across all groups,
# WITHOUT downsampling (uses each group's full native cell count).
#
# Default targets (n_cells < 4000): Micro_prenatal, SMC_postnatal,
# Immune_postnatal, VLMC_postnatal. k is set to round(ratio * n_cells),
# matching the ~1% overlap ratio your 5000-cell groups have at k=50.
#
# Usage:
#   Rscript 04_rerun_group_specific_k.R \
#     --export_dir ./cicero_export/ \
#     --outdir ./cicero_corrected \
#     --ratio 0.01 --min_k 4
#
# Outputs (in --outdir):
#   cicero_<group>_connections.tsv, cicero_<group>_ccans.tsv
#     -- same format as your original per-group files, drop-in
#        replacements for the 4 corrected groups
#   step5_rescale_k_summary.tsv
#     -- before/after comparison: original (k=50) vs corrected (rescaled
#        k) values for n_connections, n_ccans, mean_ccan_size
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export/"),
  make_option("--outdir", type = "character", default = "./cicero_corrected"),
  make_option("--groups", type = "character",
              default = "Micro_prenatal,SMC_postnatal,Immune_postnatal,VLMC_postnatal",
              help = "comma-separated groups to rerun with rescaled k [default %default]"),
  make_option("--ratio", type = "double", default = 0.01,
              help = "target fraction of a group's cell pool per metacell (i.e. k/n_cells) [default %default]"),
  make_option("--min_k", type = "integer", default = 4,
              help = "floor on k, since Cicero becomes unstable at very low k [default %default]"),
  make_option("--k_override", type = "character", default = NA_character_,
              help = "optional comma-separated group:k pairs to use exact values instead of the ratio formula, e.g. 'Micro_prenatal:10,SMC_postnatal:25'"),
  make_option("--distance_constraint", type = "double", default = 5e5),
  make_option("--window", type = "double", default = 1e6),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--seed", type = "integer", default = 123),
  make_option("--original_summary", type = "character", default = NA_character_,
              help = "optional path to your ORIGINAL step1_summary_stats.tsv, for a before/after comparison table")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

source("03_cicero_core_pipeline.R")

genomic_info <- load_genomic_and_peaktype(opt$export_dir)
genomic_coords <- genomic_info$genomic_coords
peak_type_lookup <- genomic_info$peak_type_lookup

groups <- trimws(strsplit(opt$groups, ",")[[1]])
cat(sprintf("Rescaling k for %d groups: %s\n", length(groups), paste(groups, collapse = ", ")))

k_override_map <- list()
if (!is.na(opt$k_override)) {
  pairs <- trimws(strsplit(opt$k_override, ",")[[1]])
  for (p in pairs) {
    kv <- trimws(strsplit(p, ":")[[1]])
    k_override_map[[kv[1]]] <- as.integer(kv[2])
  }
  cat("Using k_override for:", paste(names(k_override_map), unlist(k_override_map), sep = "=", collapse = ", "), "\n")
} else {
  cat(sprintf("Using ratio-based k (ratio = %.4f, min_k = %d) for all groups\n", opt$ratio, opt$min_k))
}
cat("\n")

results <- list()

for (grp in groups) {
  grp_dir <- file.path(opt$export_dir, grp)
  cell_ids_all <- readLines(file.path(grp_dir, "cell_barcodes.txt"))
  cell_ids_all <- cell_ids_all[nzchar(cell_ids_all)]
  n_cells <- length(cell_ids_all)

  k_new <- if (grp %in% names(k_override_map)) {
    k_override_map[[grp]]
  } else {
    max(opt$min_k, round(opt$ratio * n_cells))
  }
  cat(sprintf("=== %s: n_cells = %d, k = %d (was 50) ===\n", grp, n_cells, k_new))

  res <- run_cicero_group(
    grp_dir = grp_dir, cell_subset = NULL,
    genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup,
    k = k_new, distance_constraint = opt$distance_constraint, window = opt$window,
    threshold = opt$threshold, umap_seed = opt$seed,
    out_prefix = file.path(opt$outdir, paste0("cicero_", grp)),
    verbose = TRUE
  )

  if (!res$ok) {
    cat(sprintf("  !! FAILED: %s\n", res$error))
  } else {
    cat(sprintf("  OK: n_connections=%d, n_ccans=%d, mean_ccan_size=%.2f\n",
                res$summary$n_connections, res$summary$n_ccans, res$summary$mean_ccan_size))
  }
  res$summary$k_used <- k_new
  results[[grp]] <- res$summary
}

summary_df <- rbindlist(results, fill = TRUE)

if (!is.na(opt$original_summary) && file.exists(opt$original_summary)) {
  orig <- fread(opt$original_summary)[group %in% groups,
    .(group, n_connections_original_k50 = n_connections, n_ccans_original_k50 = n_ccans,
      mean_ccan_size_original_k50 = mean_ccan_size)]
  summary_df <- merge(summary_df, orig, by = "group", all.x = TRUE)
  summary_df <- summary_df[, .(group, n_cells, k_used,
                                n_connections_original_k50, n_connections_corrected = n_connections,
                                n_ccans_original_k50, n_ccans_corrected = n_ccans,
                                mean_ccan_size_original_k50, mean_ccan_size_corrected = mean_ccan_size)]
}

fwrite(summary_df, file.path(opt$outdir, "step5_rescale_k_summary.tsv"), sep = "\t")
cat("\nWrote step5_rescale_k_summary.tsv\n")
print(summary_df)

cat(sprintf("\nDone. Corrected connection/ccan files for %d groups written to: %s/\n",
            length(groups), opt$outdir))
cat("Next: assemble a merged dataset (original 10 groups + these 4 corrected ones)\n")
cat("with 05_assemble_corrected_dataset.sh before rerunning downstream analyses.\n")
