#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3g: Downsample the larger groups and rerun Cicero end-to-end
# (Approach #3).
#
# For any group with more cells than --target_n, subsamples its cells down
# to --target_n (seeded, without replacement) and reprocesses it through
# the SAME Cicero pipeline as 02_run_cicero_from_export.R (same k /
# distance_constraint / window / threshold -- see 03_cicero_core_pipeline.R,
# which this script sources). Groups already at or below --target_n are
# passed through unchanged. Compares the resulting metrics and cross-group
# rankings against the full-N step1 summary, so you can see directly
# whether conclusions from steps 1-3 survive at matched N or were driven
# by cell count.
#
# Run this from the directory containing 03_cicero_core_pipeline.R, in the
# same environment as 02_run_cicero_from_export.R (monocle 2.x, monocle3,
# cicero 1.29.2, Matrix, VGAM, stringr, uwot; plus optparse, data.table,
# ggplot2, ggrepel).
#
# Usage:
#   Rscript step3g_downsample_reprocess.R \
#     --export_dir ./cicero_export/ \
#     --target_n 571 \
#     --full_summary cicero_qc_plots/step1_summary_stats.tsv \
#     --outdir cicero_downsampled_571
#
# (571 = Micro_prenatal's own N in the running example -- pass whatever
# your smallest included group's N is, or any N you want everyone matched
# to. Omit --target_n and pass --cell_counts instead to auto-pick the min.)
#
# Outputs (in --outdir):
#   step3g_downsample_summary.tsv    - per-group metrics at target_n
#   step3g_full_vs_downsampled.tsv   - full-N vs downsampled side by side
#   step3g_rank_correlation.tsv      - Spearman rho per metric
#   step3g_rank_scatter_<metric>.png - one scatter plot per key metric
# Plus the actual reprocessed Cicero connections/CCANs in
# --out_data_dir, ready to feed back into step1/step2/step3.
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export/",
              help = "cicero_export/ directory used by 02_run_cicero_from_export.R [default %default]"),
  make_option("--target_n", type = "integer", default = NA_integer_,
              help = "downsample every larger group to this many cells (required unless --cell_counts is given)"),
  make_option("--cell_counts", type = "character", default = NA_character_,
              help = "optional tsv (group, n_cells); if given and --target_n is unset, target_n = min(n_cells) among --groups"),
  make_option("--groups", type = "character", default = NA_character_,
              help = "optional comma-separated subset of groups; default = every group in --export_dir"),
  make_option("--full_summary", type = "character", default = NA_character_,
              help = "optional step1_summary_stats.tsv (full-N) to compare against"),
  make_option("--outdir", type = "character", default = "cicero_downsampled",
              help = "output directory for comparison tables/plots [default %default]"),
  make_option("--out_data_dir", type = "character", default = NA_character_,
              help = "where to write reprocessed connections/ccans tsvs [default: <outdir>/v3_coacc_cicero_downsampled]"),
  make_option("--k", type = "integer", default = 50),
  make_option("--distance_constraint", type = "double", default = 5e5),
  make_option("--window", type = "double", default = 1e6),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--seed", type = "integer", default = 123)
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
if (is.na(opt$out_data_dir)) opt$out_data_dir <- file.path(opt$outdir, "v3_coacc_cicero_downsampled")
dir.create(opt$out_data_dir, showWarnings = FALSE, recursive = TRUE)

source("03_cicero_core_pipeline.R")

genomic_info <- load_genomic_and_peaktype(opt$export_dir)
genomic_coords <- genomic_info$genomic_coords
peak_type_lookup <- genomic_info$peak_type_lookup

groups_dirs <- list.dirs(opt$export_dir, full.names = FALSE, recursive = FALSE)
if (!is.na(opt$groups)) {
  wanted <- trimws(strsplit(opt$groups, ",")[[1]])
  groups_dirs <- intersect(groups_dirs, wanted)
}
cat(sprintf("[CHECK] Groups to process: %d (%s)\n", length(groups_dirs), paste(groups_dirs, collapse = ", ")))
if (length(groups_dirs) == 0) stop("No groups found in ", opt$export_dir)

# ---------------------- resolve target_n ----------------------
target_n <- opt$target_n
if (is.na(target_n)) {
  if (is.na(opt$cell_counts)) stop("Provide either --target_n or --cell_counts")
  cc <- fread(opt$cell_counts)
  if (!all(c("group", "n_cells") %in% names(cc))) stop("--cell_counts must contain columns: group, n_cells")
  cc <- cc[group %in% groups_dirs]
  if (nrow(cc) == 0) stop("--cell_counts has no rows matching the selected groups")
  target_n <- min(cc$n_cells)
  cat(sprintf("[CHECK] --target_n not given; using min(n_cells) across selected groups = %d\n", target_n))
}
cat(sprintf("[CHECK] Downsampling target: %d cells\n", target_n))

# ---------------------- process each group ----------------------
summary_rows <- list()

for (grp in groups_dirs) {
  grp_dir <- file.path(opt$export_dir, grp)
  cell_ids_all <- readLines(file.path(grp_dir, "cell_barcodes.txt"))
  cell_ids_all <- cell_ids_all[nzchar(cell_ids_all)]
  n_full <- length(cell_ids_all)

  outfile <- file.path(opt$out_data_dir, paste0("cicero_", grp, "_connections.tsv"))
  if (file.exists(outfile)) {
    cat(sprintf("[SKIP] %s already processed.\n", grp))
    next
  }

  if (n_full <= target_n) {
    cat(sprintf("[NOTE] %s has only %d cells (<= target_n = %d); using all cells, not subsampling.\n",
                grp, n_full, target_n))
    cell_subset <- cell_ids_all
  } else {
    set.seed(opt$seed + sum(utf8ToInt(grp)))  # per-group but reproducible seed
    cell_subset <- sample(cell_ids_all, target_n)
  }

  cat(sprintf("\n==== %s: %d -> %d cells ====\n", grp, n_full, length(cell_subset)))
  res <- run_cicero_group(
    grp_dir = grp_dir, cell_subset = cell_subset,
    genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup,
    k = opt$k, distance_constraint = opt$distance_constraint, window = opt$window,
    threshold = opt$threshold, umap_seed = opt$seed,
    out_prefix = file.path(opt$out_data_dir, paste0("cicero_", grp))
  )
  if (!res$ok) cat(sprintf("[FAILED] %s: %s\n", grp, res$error))
  summary_rows[[grp]] <- res$summary
}

downsample_summary <- rbindlist(summary_rows, fill = TRUE)
fwrite(downsample_summary, file.path(opt$outdir, "step3g_downsample_summary.tsv"), sep = "\t")
cat("\nWrote step3g_downsample_summary.tsv\n")
print(downsample_summary)

# ---------------------- compare to full-N summary, if given ----------------------
if (!is.na(opt$full_summary)) {

  full_summary <- fread(opt$full_summary)
  metrics <- intersect(
    c("n_connections", "n_ccans", "mean_ccan_size", "median_ccan_size", "median_distance_bp"),
    intersect(names(full_summary), names(downsample_summary))
  )

  merged <- merge(
    full_summary[, c("group", metrics), with = FALSE],
    downsample_summary[, c("group", metrics), with = FALSE],
    by = "group", suffixes = c("_fullN", "_downsampled")
  )
  merged <- merged[group %in% downsample_summary$group]
  fwrite(merged, file.path(opt$outdir, "step3g_full_vs_downsampled.tsv"), sep = "\t")
  cat("Wrote step3g_full_vs_downsampled.tsv\n")

  rho_rows <- list()
  for (m in metrics) {
    col_full <- paste0(m, "_fullN")
    col_down <- paste0(m, "_downsampled")
    sub <- merged[is.finite(get(col_full)) & is.finite(get(col_down))]
    if (nrow(sub) < 3) next

    rho <- suppressWarnings(cor(sub[[col_full]], sub[[col_down]], method = "spearman"))
    rho_rows[[m]] <- data.table(metric = m, n_groups = nrow(sub), spearman_rho_fullN_vs_downsampled = rho)

    p <- ggplot(sub, aes(x = get(col_full), y = get(col_down))) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      geom_point(color = "steelblue", size = 2.5) +
      geom_text_repel(aes(label = group), size = 2.8, max.overlaps = 20) +
      labs(x = paste0(m, " (full N)"), y = paste0(m, " (downsampled to ", target_n, ")"),
           title = sprintf("%s: full N vs downsampled (Spearman rho = %.2f)", m, rho)) +
      theme_bw(base_size = 11)
    ggsave(file.path(opt$outdir, paste0("step3g_rank_scatter_", m, ".png")), p, width = 7, height = 6)
    cat(sprintf("Wrote step3g_rank_scatter_%s.png (rho = %.2f)\n", m, rho))
  }

  fwrite(rbindlist(rho_rows), file.path(opt$outdir, "step3g_rank_correlation.tsv"), sep = "\t")
  cat("Wrote step3g_rank_correlation.tsv\n")
  cat("\nA Spearman rho close to 1 for a metric means the cross-group ranking for that\n",
      "metric survives at matched N -- i.e. it's likely a real cell-type effect, not a\n",
      "cell-count artifact. A low or negative rho means treat comparisons on that\n",
      "metric with caution.\n", sep = "")

} else {
  cat("\n[NOTE] --full_summary not given; skipping full-N vs downsampled comparison.\n",
      "Rerun with --full_summary cicero_qc_plots/step1_summary_stats.tsv to get it.\n", sep = "")
}

cat(sprintf("\nDone. Comparison outputs in: %s/  |  reprocessed Cicero outputs in: %s/\n",
            opt$outdir, opt$out_data_dir))
