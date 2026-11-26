#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3f: Restrict primary conclusions to the cell-count-matched
# subset (Approach #2).
#
# Does NOT rerun Cicero -- relabels your existing step1 summary table so
# reporting draws primary conclusions from the groups that share the same
# cell number, and treats the unequal-N groups as an explicit secondary /
# exploratory tier rather than pretending all groups are on equal footing.
#
# Usage:
#   Rscript step3f_matched_subset_tier.R \
#     --summary_stats cicero_qc_plots/step1_summary_stats.tsv \
#     --cell_counts cell_counts_vs_k.tsv \
#     --outdir cicero_matched_tier
#
# --cell_counts must be a tsv with at least columns: group, n_cells
#
# Outputs (in --outdir):
#   step3f_summary_with_tier.tsv     - step1 summary + n_cells + tier column
#   step3f_matched_tier_rankings.tsv - per-metric group ranking, matched tier only
#   step3f_summary_by_tier.png       - all groups, colored by tier
#   step3f_matched_tier_only.png     - matched tier only (cite this one)
#   step3f_caveat.txt                - plain-text caveat for methods/notes
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--summary_stats", type = "character", default = "cicero_qc_plots/step1_summary_stats.tsv",
              help = "step1_summary_stats.tsv [default %default]"),
  make_option("--cell_counts", type = "character", default = "cell_counts_vs_k.tsv",
              help = "tsv with at least columns: group, n_cells [default %default]"),
  make_option("--outdir", type = "character", default = "cicero_matched_tier",
              help = "output directory [default %default]"),
  make_option("--cap_n", type = "integer", default = 5000,
              help = "the cell-number cap used at export; groups at/near this N are 'matched' [default %default]"),
  make_option("--tol", type = "integer", default = 0,
              help = "a group counts as 'matched' if n_cells >= cap_n - tol [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- load + tier ----------------------
summary_df <- fread(opt$summary_stats)
cell_counts <- fread(opt$cell_counts)
if (!all(c("group", "n_cells") %in% names(cell_counts))) {
  stop("--cell_counts must contain columns: group, n_cells")
}
cell_counts <- unique(cell_counts[, .(group, n_cells)])

df <- merge(summary_df, cell_counts, by = "group", all.x = TRUE)
if (any(is.na(df$n_cells))) {
  missing_groups <- df$group[is.na(df$n_cells)]
  warning("No cell count found for: ", paste(missing_groups, collapse = ", "),
          " -- excluded from tiering.")
  df <- df[!is.na(n_cells)]
}

df[, tier := ifelse(n_cells >= (opt$cap_n - opt$tol), "matched", "exploratory")]
setorder(df, tier, -n_connections)

fwrite(df, file.path(opt$outdir, "step3f_summary_with_tier.tsv"), sep = "\t")
cat("Wrote step3f_summary_with_tier.tsv\n")
cat(sprintf("\nMatched tier (n=%d): %s\n", sum(df$tier == "matched"),
            paste(df$group[df$tier == "matched"], collapse = ", ")))
cat(sprintf("Exploratory tier (n=%d): %s\n", sum(df$tier == "exploratory"),
            paste(df$group[df$tier == "exploratory"], collapse = ", ")))

# ---------------------- matched-tier-only rankings ----------------------
metrics <- intersect(c("n_connections", "n_ccans", "mean_ccan_size", "median_ccan_size",
                        "median_distance_bp", "frac_pass_thresh"), names(df))

matched <- df[tier == "matched"]
rankings <- rbindlist(lapply(metrics, function(m) {
  sub <- matched[, .(group, value = get(m))]
  setorder(sub, -value)
  sub[, `:=`(metric = m, rank = seq_len(.N))]
  sub[, .(metric, rank, group, value)]
}))
fwrite(rankings, file.path(opt$outdir, "step3f_matched_tier_rankings.tsv"), sep = "\t")
cat("Wrote step3f_matched_tier_rankings.tsv",
    "(these rankings are safe to read as cell-type effects: n_cells is identical",
    "across every group being compared)\n", fill = TRUE)

# ---------------------- plots ----------------------
summary_long <- df |>
  select(group, tier, all_of(metrics)) |>
  pivot_longer(-c(group, tier), names_to = "metric", values_to = "value")

p_all <- ggplot(summary_long, aes(x = reorder(group, value), y = value, fill = tier)) +
  geom_col() +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  coord_flip() +
  scale_fill_manual(values = c(matched = "steelblue", exploratory = "firebrick")) +
  labs(x = NULL, y = NULL, fill = "Tier", title = "All groups, colored by cell-count tier") +
  theme_bw(base_size = 9)
ggsave(file.path(opt$outdir, "step3f_summary_by_tier.png"), p_all, width = 18, height = 9)
cat("Wrote step3f_summary_by_tier.png\n")

matched_long <- summary_long |> filter(tier == "matched")
p_matched <- ggplot(matched_long, aes(x = reorder(group, value), y = value)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  coord_flip() +
  labs(x = NULL, y = NULL,
       title = sprintf("Cell-count-matched tier only (n_cells >= %d)", opt$cap_n - opt$tol)) +
  theme_bw(base_size = 9)
ggsave(file.path(opt$outdir, "step3f_matched_tier_only.png"), p_matched,
       width = 18, height = max(6, 3 * ceiling(length(metrics) / 3)))
cat("Wrote step3f_matched_tier_only.png\n")

# ---------------------- caveat text ----------------------
caveat <- c(
  sprintf("Cell-count tiering (generated %s)", Sys.time()),
  "",
  sprintf("Matched tier (n_cells >= %d, n=%d groups): %s",
          opt$cap_n - opt$tol, sum(df$tier == "matched"),
          paste(df$group[df$tier == "matched"], collapse = ", ")),
  "  -> Primary comparisons should be drawn from this tier. Cell number is",
  "     identical (or effectively identical) across these groups, so any",
  "     differences in Cicero output reflect biology, not sampling depth.",
  "",
  sprintf("Exploratory tier (n_cells < %d, n=%d groups): %s",
          opt$cap_n - opt$tol, sum(df$tier == "exploratory"),
          paste(df$group[df$tier == "exploratory"], collapse = ", ")),
  "  -> Treat conclusions involving these groups as exploratory/secondary",
  "     until cross-checked at matched N. See:",
  "       step3e_cellcount_confound_check.R  - quantitative residual test",
  "       step3g_downsample_reprocess.R      - do the ranking / conclusions",
  "                                             survive when everyone is",
  "                                             downsampled to match?",
  "       step3h_bootstrap_ci.R              - is an apparent difference",
  "                                             bigger than sampling noise",
  "                                             at that cell count?"
)
writeLines(caveat, file.path(opt$outdir, "step3f_caveat.txt"))
cat("Wrote step3f_caveat.txt\n")

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
