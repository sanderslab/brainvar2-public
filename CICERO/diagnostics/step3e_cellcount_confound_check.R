#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3e: Cell-count confound check (Approach #1).
#
# Uses data you already have -- no reruns needed. For each key metric,
# fits metric ~ n_cells across your existing groups and reports the
# residuals: a group far off the fitted line is a real outlier after
# accounting for cell count; a group that sits on the line is just
# following the expected cell-count relationship, not showing distinct
# biology.
#
# Usage:
#   Rscript step3e_cellcount_confound_check.R \
#     --summary_stats cicero_qc_plots/step1_summary_stats.tsv \
#     --cell_counts cell_counts_vs_k.tsv \
#     --outdir cicero_cellcount_confound_check
#
# --cell_counts must be a tsv with at least columns: group, n_cells
#
# Outputs (in --outdir):
#   residuals_summary.tsv   - group x metric: observed, fitted, residual,
#                              standardized residual, flagged_outlier
#   fit_summary.tsv         - per metric: R^2, Spearman rho, n groups, n flagged
#   residuals_<metric>.png  - metric vs n_cells, fitted line, outliers in red
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--summary_stats", type = "character", default = "cicero_qc_plots/step1_summary_stats.tsv",
              help = "step1_summary_stats.tsv [default %default]"),
  make_option("--cell_counts", type = "character", default = "cell_counts_vs_k.tsv",
              help = "tsv with at least columns: group, n_cells [default %default]"),
  make_option("--outdir", type = "character", default = "cicero_cellcount_confound_check",
              help = "output directory [default %default]"),
  make_option("--sd_flag", type = "double", default = 2,
              help = "flag a group if |standardized residual| exceeds this many SD [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- load + merge ----------------------
summary_df <- fread(opt$summary_stats)
cell_counts <- fread(opt$cell_counts)
if (!all(c("group", "n_cells") %in% names(cell_counts))) {
  stop("--cell_counts must contain columns: group, n_cells")
}
cell_counts <- unique(cell_counts[, .(group, n_cells)])

df <- merge(summary_df, cell_counts, by = "group")
missing_groups <- setdiff(summary_df$group, df$group)
if (length(missing_groups) > 0) {
  warning("No cell count found for: ", paste(missing_groups, collapse = ", "), " -- these are dropped.")
}
cat(sprintf("Fitting cell-count trends across %d groups.\n", nrow(df)))

# ---------------------- fit metric ~ n_cells for each metric ----------------------
metrics <- intersect(c("n_connections", "n_ccans", "mean_ccan_size", "median_distance_bp"), names(df))
if (length(metrics) == 0) stop("None of the expected metrics (n_connections, n_ccans, mean_ccan_size, median_distance_bp) were found in --summary_stats")

residual_rows <- list()
fit_rows <- list()

for (m in metrics) {
  sub <- df[is.finite(get(m)) & is.finite(n_cells)]
  if (nrow(sub) < 3) {
    cat("Skipping", m, "- fewer than 3 groups with finite values\n")
    next
  }

  fit <- lm(as.formula(paste(m, "~ n_cells")), data = sub)
  r2 <- summary(fit)$r.squared
  rho <- suppressWarnings(cor(sub$n_cells, sub[[m]], method = "spearman"))
  resid_std <- rstandard(fit)

  sub[, `:=`(
    fitted_value = fitted(fit),
    residual = get(m) - fitted(fit),
    residual_sd = resid_std,
    flagged_outlier = abs(resid_std) > opt$sd_flag
  )]

  residual_rows[[m]] <- sub[, .(
    group, metric = m, n_cells,
    observed = get(m), fitted_value, residual, residual_sd, flagged_outlier
  )]

  fit_rows[[m]] <- data.table(
    metric = m, r_squared = r2, spearman_rho = rho,
    n_groups = nrow(sub), n_flagged = sum(sub$flagged_outlier)
  )
}

residuals_summary <- rbindlist(residual_rows)
fit_summary <- rbindlist(fit_rows)

fwrite(residuals_summary, file.path(opt$outdir, "residuals_summary.tsv"), sep = "\t")
fwrite(fit_summary, file.path(opt$outdir, "fit_summary.tsv"), sep = "\t")
cat("Wrote residuals_summary.tsv and fit_summary.tsv\n")

cat("\n=== Fit summary (high R^2 = mostly a cell-count effect, not biology) ===\n")
print(fit_summary)

flagged <- residuals_summary[flagged_outlier == TRUE]
cat(sprintf("\n=== Groups flagged as outliers beyond cell count (|residual| > %.1f SD) ===\n", opt$sd_flag))
if (nrow(flagged) > 0) {
  print(flagged[, .(group, metric, n_cells, observed, fitted_value, residual_sd)])
} else {
  cat("None.\n")
}

# ---------------------- plots ----------------------
for (m in metrics) {
  sub <- residuals_summary[metric == m]
  if (nrow(sub) == 0) next

  p <- ggplot(sub, aes(x = n_cells, y = observed)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.6) +
    geom_point(aes(color = flagged_outlier), size = 2.5) +
    geom_text_repel(aes(label = group), size = 2.8, max.overlaps = 20) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue"), guide = "none") +
    labs(x = "n_cells", y = m,
         title = sprintf("%s vs n_cells (red = outlier beyond %.1f SD from trend)", m, opt$sd_flag)) +
    theme_bw(base_size = 11)

  ggsave(file.path(opt$outdir, paste0("residuals_", m, ".png")), p, width = 7, height = 5)
  cat("Wrote residuals_", m, ".png\n", sep = "")
}

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
