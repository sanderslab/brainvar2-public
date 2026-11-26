#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# APPROACH 1: Does cell count explain cross-group differences?
# (Same script as before -- included here for completeness so all 4
#  approaches are in one set.)
#
# Usage:
#   Rscript approach1_cellcount_regression.R \
#     --summary_stats cicero_qc_plots/step1_summary_stats.tsv \
#     --cell_counts cell_counts_vs_k.tsv \
#     --outdir approach1_regression_check
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

option_list <- list(
  make_option("--summary_stats", type = "character",
              default = "cicero_qc_plots/step1_summary_stats.tsv"),
  make_option("--cell_counts", type = "character", default = "cell_counts_vs_k.tsv"),
  make_option("--outdir", type = "character", default = "approach1_regression_check")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

stats <- fread(opt$summary_stats)
cells <- fread(opt$cell_counts)[, .(group, n_cells)]
merged <- merge(stats, cells, by = "group")
fwrite(merged, file.path(opt$outdir, "merged_metrics_with_cellcounts.tsv"), sep = "\t")

metrics_to_check <- intersect(
  c("n_connections", "n_ccans", "mean_ccan_size", "median_distance_bp"),
  names(merged)
)

residual_list <- list()
for (m in metrics_to_check) {
  d <- merged[, .(group, n_cells, value = get(m))]
  d <- d[is.finite(value) & is.finite(n_cells)]
  fit <- lm(value ~ n_cells, data = d)
  d[, fitted := predict(fit)]
  d[, resid_raw := value - fitted]
  d[, resid_sd := resid_raw / sd(resid_raw)]
  d[, is_outlier := abs(resid_sd) > 2]
  r2 <- summary(fit)$r.squared
  spearman_r <- cor(d$n_cells, d$value, method = "spearman")
  cat(sprintf("\n[%s] R^2 = %.3f | Spearman rho = %.3f\n", m, r2, spearman_r))
  outliers <- d[is_outlier == TRUE][order(-abs(resid_sd))]
  if (nrow(outliers) > 0) print(outliers[, .(group, n_cells, value, resid_sd = round(resid_sd, 2))])
  else cat("  (no outliers -- metric fully explained by cell count)\n")
  d[, metric := m]
  residual_list[[m]] <- d

  p <- ggplot(d, aes(x = n_cells, y = value)) +
    geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
    geom_point(aes(color = is_outlier), size = 3) +
    ggrepel::geom_text_repel(aes(label = group), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "steelblue")) +
    labs(title = paste0(m, " vs. n_cells"),
         subtitle = sprintf("R^2=%.2f, rho=%.2f | red = >2SD outlier from trend", r2, spearman_r),
         x = "Number of cells", y = m) +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(opt$outdir, paste0(m, "_vs_ncells.png")), p, width = 8, height = 6)
}

fwrite(rbindlist(residual_list), file.path(opt$outdir, "residuals_summary.tsv"), sep = "\t")
cat("\nDone. Outputs in:", opt$outdir, "\n")
