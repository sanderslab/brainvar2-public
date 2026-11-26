#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 5c: k-sweep for groups where the ratio-based k rescaling
# overcorrected (too-small absolute k -> noisy/fragmented CCANs).
#
# Runs Cicero at several candidate k values for one or more groups
# (at their FULL native cell count, no downsampling) and reports
# n_connections / n_ccans / mean_ccan_size for each, so you can pick
# the k where all three land closest to your normal range (roughly
# what your 10 cell-count-matched n=5000 groups show).
#
# Usage:
#   Rscript step5c_k_sweep.R \
#     --export_dir ./cicero_export/ \
#     --groups Micro_prenatal,SMC_postnatal \
#     --k_values 6,10,15,20,25,30 \
#     --outdir cicero_k_sweep
#
# Outputs:
#   step5c_k_sweep_results.tsv   -- one row per group x k
#   step5c_k_sweep_<metric>.png  -- metric vs k, one line per group,
#                                   with a shaded band showing the
#                                   normal range you supply
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export/"),
  make_option("--groups", type = "character", default = "Micro_prenatal,SMC_postnatal"),
  make_option("--k_values", type = "character", default = "6,10,15,20,25,30,40,50",
              help = "comma-separated k values to try [default %default]"),
  make_option("--outdir", type = "character", default = "cicero_k_sweep"),
  make_option("--distance_constraint", type = "double", default = 5e5),
  make_option("--window", type = "double", default = 1e6),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--seed", type = "integer", default = 123),
  # normal range from your 10 matched n=5000 groups -- adjust if needed
  make_option("--normal_n_connections_lo", type = "double", default = 60000),
  make_option("--normal_n_connections_hi", type = "double", default = 330000),
  make_option("--normal_n_ccans_lo", type = "double", default = 2600),
  make_option("--normal_n_ccans_hi", type = "double", default = 5000),
  make_option("--normal_ccan_size_lo", type = "double", default = 8),
  make_option("--normal_ccan_size_hi", type = "double", default = 33)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

source("03_cicero_core_pipeline.R")

genomic_info <- load_genomic_and_peaktype(opt$export_dir)
genomic_coords <- genomic_info$genomic_coords
peak_type_lookup <- genomic_info$peak_type_lookup

groups <- trimws(strsplit(opt$groups, ",")[[1]])
k_values <- as.integer(trimws(strsplit(opt$k_values, ",")[[1]]))

cat(sprintf("Sweeping k = %s for groups: %s\n\n", paste(k_values, collapse = ", "), paste(groups, collapse = ", ")))

results <- list()

for (grp in groups) {
  grp_dir <- file.path(opt$export_dir, grp)
  for (k in k_values) {
    cat(sprintf("=== %s, k = %d ===\n", grp, k))
    res <- run_cicero_group(
      grp_dir = grp_dir, cell_subset = NULL,
      genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup,
      k = k, distance_constraint = opt$distance_constraint, window = opt$window,
      threshold = opt$threshold, umap_seed = opt$seed,
      out_prefix = NULL, verbose = FALSE
    )
    if (!res$ok) {
      cat(sprintf("  FAILED: %s\n", res$error))
      next
    }
    cat(sprintf("  n_connections=%d, n_ccans=%d, mean_ccan_size=%.2f\n",
                res$summary$n_connections, res$summary$n_ccans, res$summary$mean_ccan_size))
    row <- res$summary
    row$k <- k
    results[[paste(grp, k)]] <- row
  }
}

sweep_df <- rbindlist(results, fill = TRUE)

# flag whether each metric falls inside the "normal" band you supplied
sweep_df[, in_range_connections := n_connections >= opt$normal_n_connections_lo & n_connections <= opt$normal_n_connections_hi]
sweep_df[, in_range_ccans := n_ccans >= opt$normal_n_ccans_lo & n_ccans <= opt$normal_n_ccans_hi]
sweep_df[, in_range_ccan_size := mean_ccan_size >= opt$normal_ccan_size_lo & mean_ccan_size <= opt$normal_ccan_size_hi]
sweep_df[, n_metrics_in_range := in_range_connections + in_range_ccans + in_range_ccan_size]

fwrite(sweep_df, file.path(opt$outdir, "step5c_k_sweep_results.tsv"), sep = "\t")
cat("\nWrote step5c_k_sweep_results.tsv\n\n")
print(sweep_df[, .(group, k, n_connections, n_ccans, mean_ccan_size, n_metrics_in_range)])

cat("\nBest k per group (most metrics simultaneously in normal range, ties broken by n_connections closest to range midpoint):\n")
mid_conn <- (opt$normal_n_connections_lo + opt$normal_n_connections_hi) / 2
best <- sweep_df[, .SD[order(-n_metrics_in_range, abs(n_connections - mid_conn))][1], by = group]
print(best[, .(group, k, n_connections, n_ccans, mean_ccan_size, n_metrics_in_range)])
fwrite(best, file.path(opt$outdir, "step5c_recommended_k.tsv"), sep = "\t")
cat("Wrote step5c_recommended_k.tsv\n")

# ---------------- plots ----------------
metrics_plot <- list(
  n_connections = c(opt$normal_n_connections_lo, opt$normal_n_connections_hi),
  n_ccans = c(opt$normal_n_ccans_lo, opt$normal_n_ccans_hi),
  mean_ccan_size = c(opt$normal_ccan_size_lo, opt$normal_ccan_size_hi)
)
for (m in names(metrics_plot)) {
  band <- metrics_plot[[m]]
  p <- ggplot(sweep_df, aes(x = k, y = .data[[m]], color = group)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
             alpha = 0.1, fill = "forestgreen") +
    geom_line() + geom_point(size = 2) +
    labs(title = paste0(m, " vs k"), subtitle = "shaded band = normal range from matched n=5000 groups") +
    theme_bw(base_size = 11)
  ggsave(file.path(opt$outdir, paste0("step5c_k_sweep_", m, ".png")), p, width = 8, height = 6)
}
cat("\nDone. Plots + tables in:", opt$outdir, "\n")
