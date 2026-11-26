#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3h: Within-group subsampling confidence intervals for Cicero
# metrics (Approach #4).
#
# For each group, repeatedly subsamples its cells and reruns the full
# Cicero pipeline (via 03_cicero_core_pipeline.R), so you get a distribution
# -- not just a single point estimate -- for each metric. This tells you
# how much a metric would plausibly vary from sampling noise alone at
# that cell count, so you can judge whether an observed difference
# between two groups is bigger than that.
#
# Note on terminology: this resamples WITHOUT replacement (subsampling),
# not classic with-replacement bootstrap. Cicero's KNN-based metacell
# aggregation (make_cicero_cds) needs unique cell barcodes as matrix
# column names, so duplicating cells the way a real bootstrap would isn't
# compatible with the pipeline. Repeated subsampling at a fixed size is
# the standard substitute and answers the same question: how much would
# this metric move around just from which cells happened to be sampled?
#
# Run this from the directory containing 03_cicero_core_pipeline.R, in the
# same environment as 02_run_cicero_from_export.R.
#
# Usage (bootstrap each group at its own N, 20 reps each):
#   Rscript step3h_bootstrap_ci.R \
#     --export_dir ./cicero_export/ \
#     --n_reps 20 --frac 0.8 \
#     --outdir cicero_bootstrap_ci
#
# Usage (bootstrap specific groups at one SHARED matched N -- directly
# comparable to step3g's downsampled point estimates):
#   Rscript step3h_bootstrap_ci.R \
#     --export_dir ./cicero_export/ \
#     --groups Micro_prenatal,SMC_postnatal,Immune_postnatal,VLMC_postnatal,ExN_postnatal \
#     --n_reps 20 --target_n 571 \
#     --outdir cicero_bootstrap_ci_matched571
#
# Outputs (in --outdir):
#   step3h_bootstrap_raw_reps.tsv    - one row per group x rep
#   step3h_bootstrap_ci_summary.tsv  - per group x metric: mean, sd, CI
#   step3h_ci_<metric>.png           - CI errorbar plot per metric
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export/",
              help = "cicero_export/ directory used by 02_run_cicero_from_export.R [default %default]"),
  make_option("--groups", type = "character", default = NA_character_,
              help = "optional comma-separated subset of groups; default = every group in --export_dir"),
  make_option("--n_reps", type = "integer", default = 20,
              help = "number of subsampling repeats per group [default %default]"),
  make_option("--frac", type = "double", default = 0.8,
              help = "fraction of each group's own cells sampled per rep, used when --target_n is not set [default %default]"),
  make_option("--target_n", type = "integer", default = NA_integer_,
              help = "if set, every rep for every group samples exactly this many cells; groups smaller than this are skipped"),
  make_option("--outdir", type = "character", default = "cicero_bootstrap_ci",
              help = "output directory [default %default]"),
  make_option("--k", type = "integer", default = 50),
  make_option("--distance_constraint", type = "double", default = 5e5),
  make_option("--window", type = "double", default = 1e6),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--seed", type = "integer", default = 123),
  make_option("--ci", type = "double", default = 0.95,
              help = "confidence interval width, e.g. 0.95 for 2.5/97.5 percentiles [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

source("03_cicero_core_pipeline.R")

genomic_info <- load_genomic_and_peaktype(opt$export_dir)
genomic_coords <- genomic_info$genomic_coords
peak_type_lookup <- genomic_info$peak_type_lookup

groups_dirs <- list.dirs(opt$export_dir, full.names = FALSE, recursive = FALSE)
if (!is.na(opt$groups)) {
  wanted <- trimws(strsplit(opt$groups, ",")[[1]])
  groups_dirs <- intersect(groups_dirs, wanted)
}
if (length(groups_dirs) == 0) stop("No groups found in ", opt$export_dir)

cat(sprintf("[CHECK] Groups to bootstrap: %d (%s)\n", length(groups_dirs), paste(groups_dirs, collapse = ", ")))
cat(sprintf(
  "[CHECK] %d reps per group, %s\n", opt$n_reps,
  if (is.na(opt$target_n)) sprintf("sampling %.0f%% of each group's own cells per rep", 100 * opt$frac)
  else sprintf("sampling exactly %d cells per rep", opt$target_n)
))

# ---------------------- run reps ----------------------
raw_rows <- list()

for (grp in groups_dirs) {
  grp_dir <- file.path(opt$export_dir, grp)
  cell_ids_all <- readLines(file.path(grp_dir, "cell_barcodes.txt"))
  cell_ids_all <- cell_ids_all[nzchar(cell_ids_all)]
  n_full <- length(cell_ids_all)

  if (!is.na(opt$target_n)) {
    if (n_full < opt$target_n) {
      cat(sprintf("[SKIP] %s has only %d cells (< target_n = %d).\n", grp, n_full, opt$target_n))
      next
    }
    n_rep <- opt$target_n
  } else {
    n_rep <- max(2, round(opt$frac * n_full))
  }

  cat(sprintf("\n==== %s: %d reps of %d / %d cells ====\n", grp, opt$n_reps, n_rep, n_full))

  for (rep_i in seq_len(opt$n_reps)) {
    set.seed(opt$seed + sum(utf8ToInt(grp)) + rep_i * 997)
    cell_subset <- sample(cell_ids_all, n_rep)

    res <- run_cicero_group(
      grp_dir = grp_dir, cell_subset = cell_subset,
      genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup,
      k = opt$k, distance_constraint = opt$distance_constraint, window = opt$window,
      threshold = opt$threshold, umap_seed = opt$seed + rep_i,
      out_prefix = NULL, verbose = FALSE
    )
    row <- res$summary
    row$rep <- rep_i
    row$ok <- res$ok
    raw_rows[[paste(grp, rep_i)]] <- row
    cat(sprintf("  rep %d/%d: %s\n", rep_i, opt$n_reps,
                if (res$ok) "ok" else paste0("FAILED (", res$error, ")")))
  }
}

raw_reps <- rbindlist(raw_rows, fill = TRUE)
fwrite(raw_reps, file.path(opt$outdir, "step3h_bootstrap_raw_reps.tsv"), sep = "\t")
cat("\nWrote step3h_bootstrap_raw_reps.tsv\n")

# ---------------------- CI summary ----------------------
metrics <- intersect(
  c("n_connections", "n_ccans", "mean_ccan_size", "median_ccan_size", "median_distance_bp"),
  names(raw_reps)
)
alpha <- 1 - opt$ci

ci_rows <- list()
for (grp in unique(raw_reps$group)) {
  sub_grp <- raw_reps[group == grp & ok == TRUE]
  if (nrow(sub_grp) == 0) next
  for (m in metrics) {
    vals <- sub_grp[[m]]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2) next
    ci_rows[[paste(grp, m)]] <- data.table(
      group = grp, metric = m, n_reps_ok = length(vals),
      n_cells_per_rep = sub_grp$n_cells[1],
      mean = mean(vals), sd = sd(vals),
      ci_lo = as.numeric(quantile(vals, alpha / 2)),
      ci_hi = as.numeric(quantile(vals, 1 - alpha / 2))
    )
  }
}
ci_summary <- rbindlist(ci_rows)
fwrite(ci_summary, file.path(opt$outdir, "step3h_bootstrap_ci_summary.tsv"), sep = "\t")
cat("Wrote step3h_bootstrap_ci_summary.tsv\n")
print(ci_summary)

# ---------------------- plots: one errorbar plot per metric ----------------------
for (m in metrics) {
  sub <- ci_summary[metric == m]
  if (nrow(sub) == 0) next

  p <- ggplot(sub, aes(x = reorder(group, mean), y = mean)) +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), color = "steelblue") +
    coord_flip() +
    labs(x = NULL, y = m,
         title = sprintf("%s: subsampling %.0f%% CI (up to %d reps)", m, opt$ci * 100, opt$n_reps)) +
    theme_bw(base_size = 10)
  ggsave(file.path(opt$outdir, paste0("step3h_ci_", m, ".png")), p,
         width = 8, height = max(3, 0.4 * nrow(sub)))
  cat(sprintf("Wrote step3h_ci_%s.png\n", m))
}

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
cat("If two groups' CIs for a metric don't overlap, the difference between them is\n",
    "unlikely to be explained by sampling noise at that cell count. If they overlap\n",
    "substantially, be cautious calling that a real difference.\n", sep = "")
