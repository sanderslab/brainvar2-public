#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3c: Check for a near-threshold co-accessibility score pile-up.
#
# Rationale: an unstable graphical-lasso fit (e.g. too little independent
# information per genomic window) tends to retain many weak/borderline
# edges rather than fewer strong ones. If Micro_prenatal's excess edges
# are disproportionately clustered just above the 0.2 threshold compared
# to other groups, that's evidence of an unstable fit rather than a
# genuinely denser regulatory landscape.
#
# Usage:
#   Rscript step3c_score_pileup_check.R --indir /path/to/v3_coacc_cicero \
#       --groups Micro_prenatal,Micro_postnatal,ExN_prenatal,ExN_postnatal \
#       --outdir cicero_pileup_check
#
# Outputs:
#   pileup_check_overlay_density.png   - density-normalized score histograms,
#                                         shared (non-free) x/y axes, all groups
#                                         on one panel for direct comparison
#   pileup_check_zoomed_0.15_0.4.png   - same, zoomed to the threshold region
#   pileup_check_ecdf.png              - empirical CDF of |coaccess| per group;
#                                         a steep rise right at the threshold
#                                         with a long flat run after = pile-up
#   pileup_check_fraction_near_threshold.tsv
#                                       - % of each group's edges falling in
#                                         [threshold, threshold+band] vs above
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option("--indir", type = "character", default = "."),
  make_option("--outdir", type = "character", default = "cicero_pileup_check"),
  make_option("--groups", type = "character",
              default = "Micro_prenatal,Micro_postnatal,ExN_prenatal,ExN_postnatal",
              help = "Comma-separated list of groups to compare [default %default]"),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--band", type = "double", default = 0.1,
              help = "Width above threshold counted as 'near-threshold' [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

groups <- trimws(strsplit(opt$groups, ",")[[1]])
cat(sprintf("Comparing %d groups: %s\n", length(groups), paste(groups, collapse = ", ")))

load_scores <- function(indir, group) {
  f <- file.path(indir, paste0("cicero_", group, "_connections.tsv"))
  df <- fread(f, sep = "\t", select = c("coaccess"))
  df <- df[!is.na(coaccess)]
  df[, abs_coaccess := abs(coaccess)]
  df[, group := group]
  df
}

all_scores <- rbindlist(lapply(groups, function(g) {
  cat("Loading:", g, "\n")
  load_scores(opt$indir, g)
}))

# ---------------- overlay density histogram (shared axes) ----------------
p_overlay <- ggplot(all_scores, aes(x = abs_coaccess, color = group, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, position = "identity", alpha = 0.15) +
  geom_vline(xintercept = opt$threshold, linetype = "dashed", color = "black") +
  labs(title = "abs(co-accessibility score) distribution, shared axes",
       subtitle = paste0("Dashed line = threshold (", opt$threshold, ")"),
       x = "abs(coaccess)", y = "Density") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "pileup_check_overlay_density.png"), p_overlay, width = 8, height = 6)
cat("Wrote pileup_check_overlay_density.png\n")

# ---------------- zoomed view near threshold ----------------
p_zoom <- ggplot(all_scores, aes(x = abs_coaccess, color = group, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, position = "identity", alpha = 0.15) +
  geom_vline(xintercept = opt$threshold, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(opt$threshold - 0.05, opt$threshold + 0.2)) +
  labs(title = "Zoomed: abs(co-accessibility score) near threshold",
       x = "abs(coaccess)", y = "Density") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "pileup_check_zoomed_0.15_0.4.png"), p_zoom, width = 8, height = 6)
cat("Wrote pileup_check_zoomed_0.15_0.4.png\n")

# ---------------- ECDF ----------------
p_ecdf <- ggplot(all_scores, aes(x = abs_coaccess, color = group)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  geom_vline(xintercept = opt$threshold, linetype = "dashed", color = "black") +
  labs(title = "Empirical CDF of abs(co-accessibility score)",
       subtitle = "A steep rise right at threshold + long flat run after = pile-up of weak edges",
       x = "abs(coaccess)", y = "Cumulative fraction of edges") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "pileup_check_ecdf.png"), p_ecdf, width = 8, height = 6)
cat("Wrote pileup_check_ecdf.png\n")

# ---------------- quantitative pile-up summary ----------------
pileup_summary <- all_scores[, .(
  n_total = .N,
  n_near_threshold = sum(abs_coaccess >= opt$threshold & abs_coaccess < (opt$threshold + opt$band)),
  n_above_band = sum(abs_coaccess >= (opt$threshold + opt$band)),
  median_abs_score = median(abs_coaccess)
), by = group]
pileup_summary[, frac_near_threshold := n_near_threshold / n_total]
pileup_summary[, frac_above_band := n_above_band / n_total]
pileup_summary <- pileup_summary[order(-frac_near_threshold)]

fwrite(pileup_summary, file.path(opt$outdir, "pileup_check_fraction_near_threshold.tsv"), sep = "\t")
cat("\nWrote pileup_check_fraction_near_threshold.tsv\n")
print(pileup_summary)

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
cat("\nHow to read this:\n")
cat(sprintf("  - frac_near_threshold = fraction of edges with abs(coaccess) in [%.2f, %.2f)\n",
            opt$threshold, opt$threshold + opt$band))
cat("  - If Micro_prenatal's frac_near_threshold is much higher than the other groups,\n")
cat("    its excess edges are disproportionately weak/borderline -> supports an unstable fit.\n")
cat("  - If frac_near_threshold is similar across groups but Micro_prenatal simply has\n")
cat("    proportionally more edges everywhere (including strong ones), that argues against\n")
cat("    a pure threshold-instability explanation and the difference is more structural.\n")
