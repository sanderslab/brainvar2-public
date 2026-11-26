#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3h (fast): Parallelized within-group subsampling CIs.
#
# Speed changes vs. the original step3h_bootstrap_ci.R:
#   1. Each group's matrix/coords are loaded from disk ONCE
#      (load_group_assets), not once per rep.
#   2. Reps run in parallel via parallel::mclapply (fork-based -- workers
#      share the preloaded matrix via copy-on-write memory, no per-rep
#      re-reading or re-serializing).
#   3. By default, each group's distance_parameter is estimated ONCE from
#      its full cell set and reused across all reps for that group
#      (--reuse_distance_param, on by default). Pass
#      --no_reuse_distance_param for the exact original (slower) behavior
#      of re-estimating it fresh on every subsample.
#
# Usage:
#   Rscript step3h_bootstrap_ci_fast.R \
#     --export_dir ./cicero_export/ \
#     --n_reps 20 --frac 0.8 --n_cores 8 \
#     --outdir cicero_bootstrap_ci_fast
#
#   Rscript step3h_bootstrap_ci_fast.R \
#     --export_dir ./cicero_export/ \
#     --groups Micro_prenatal,SMC_postnatal,Immune_postnatal,VLMC_postnatal,ExN_postnatal \
#     --n_reps 20 --target_n 571 --n_cores 8 \
#     --outdir cicero_bootstrap_ci_matched571_fast
#
# Outputs: identical set/format to step3h_bootstrap_ci.R
#   step3h_bootstrap_raw_reps.tsv, step3h_bootstrap_ci_summary.tsv,
#   step3h_ci_<metric>.png
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export/"),
  make_option("--groups", type = "character", default = NA_character_),
  make_option("--n_reps", type = "integer", default = 20),
  make_option("--frac", type = "double", default = 0.8),
  make_option("--target_n", type = "integer", default = NA_integer_),
  make_option("--outdir", type = "character", default = "cicero_bootstrap_ci_fast"),
  make_option("--k", type = "integer", default = 50),
  make_option("--distance_constraint", type = "double", default = 5e5),
  make_option("--window", type = "double", default = 1e6),
  make_option("--threshold", type = "double", default = 0.2),
  make_option("--seed", type = "integer", default = 123),
  make_option("--ci", type = "double", default = 0.95),
  make_option("--n_cores", type = "integer", default = max(1, parallel::detectCores() - 1),
              help = "parallel workers (fork-based; Linux/Mac only) [default: detectCores()-1]"),
  make_option("--reuse_distance_param", action = "store_true", default = TRUE,
              help = "estimate distance_parameter once per group (full cell set) and reuse across reps [default TRUE]"),
  make_option("--no_reuse_distance_param", action = "store_false", dest = "reuse_distance_param",
              help = "re-estimate distance_parameter fresh on every rep (slower, matches original exactly)")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

if (.Platform$OS.type == "windows" && opt$n_cores > 1) {
  warning("mclapply forking is not supported on Windows; forcing n_cores = 1.")
  opt$n_cores <- 1
}

source("exploratory/cicero_core_pipeline_fast.R")

genomic_info <- load_genomic_and_peaktype(opt$export_dir)
genomic_coords <- genomic_info$genomic_coords
peak_type_lookup <- genomic_info$peak_type_lookup

groups_dirs <- list.dirs(opt$export_dir, full.names = FALSE, recursive = FALSE)
if (!is.na(opt$groups)) {
  wanted <- trimws(strsplit(opt$groups, ",")[[1]])
  groups_dirs <- intersect(groups_dirs, wanted)
}
if (length(groups_dirs) == 0) stop("No groups found in ", opt$export_dir)

cat(sprintf("[CHECK] Groups: %d (%s)\n", length(groups_dirs), paste(groups_dirs, collapse = ", ")))
cat(sprintf("[CHECK] n_cores = %d | reuse_distance_param = %s\n", opt$n_cores, opt$reuse_distance_param))

# ---------------------- preload assets + build task list ----------------------
tasks <- list()
fixed_dp_by_group <- list()

for (grp in groups_dirs) {
  grp_dir <- file.path(opt$export_dir, grp)
  t0 <- Sys.time()
  assets <- load_group_assets(grp_dir)
  cat(sprintf("[LOAD] %s: %d cells, %d peaks (%.1fs)\n",
              grp, length(assets$cell_ids_all), nrow(assets$mat),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  n_full <- length(assets$cell_ids_all)
  if (!is.na(opt$target_n)) {
    if (n_full < opt$target_n) {
      cat(sprintf("[SKIP] %s has only %d cells (< target_n = %d).\n", grp, n_full, opt$target_n))
      next
    }
    n_rep <- opt$target_n
  } else {
    n_rep <- max(2, round(opt$frac * n_full))
  }

  if (isTRUE(opt$reuse_distance_param)) {
    t0 <- Sys.time()
    fixed_dp <- tryCatch(
      estimate_group_distance_parameter(
        assets, genomic_coords, k = opt$k, distance_constraint = opt$distance_constraint,
        window = opt$window, umap_seed = opt$seed, verbose = TRUE
      ),
      error = function(e) {
        cat(sprintf("[WARN] %s: full-group distance_parameter estimation failed (%s); ",
                     grp, conditionMessage(e)))
        cat("falling back to per-rep estimation for this group.\n")
        NULL
      }
    )
    cat(sprintf("[DIST_PARAM] %s: %s (%.1fs)\n", grp,
                if (is.null(fixed_dp)) "per-rep fallback" else signif(fixed_dp, 6),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  } else {
    fixed_dp <- NULL
  }
  fixed_dp_by_group[[grp]] <- fixed_dp

  for (rep_i in seq_len(opt$n_reps)) {
    tasks[[length(tasks) + 1]] <- list(grp = grp, rep_i = rep_i, n_rep = n_rep, assets = assets, fixed_dp = fixed_dp)
  }
}

cat(sprintf("\n[CHECK] Total tasks (group x rep) to run: %d across %d cores\n", length(tasks), opt$n_cores))

# ---------------------- run all tasks in parallel ----------------------
run_one_task <- function(task) {
  set.seed(opt$seed + sum(utf8ToInt(task$grp)) + task$rep_i * 997)
  cell_subset <- sample(task$assets$cell_ids_all, task$n_rep)

  res <- run_cicero_group_fast(
    assets = task$assets, cell_subset = cell_subset,
    genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup,
    k = opt$k, distance_constraint = opt$distance_constraint, window = opt$window,
    threshold = opt$threshold, umap_seed = opt$seed + task$rep_i,
    fixed_distance_parameter = task$fixed_dp, out_prefix = NULL, verbose = FALSE
  )
  row <- res$summary
  row$rep <- task$rep_i
  row$ok <- res$ok
  row
}

t0 <- Sys.time()
raw_rows <- parallel::mclapply(tasks, run_one_task, mc.cores = opt$n_cores, mc.preschedule = TRUE)
cat(sprintf("\n[TIMING] All %d tasks finished in %.1f minutes\n",
            length(tasks), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

raw_reps <- rbindlist(raw_rows, fill = TRUE)
fwrite(raw_reps, file.path(opt$outdir, "step3h_bootstrap_raw_reps.tsv"), sep = "\t")
cat("Wrote step3h_bootstrap_raw_reps.tsv\n")

n_failed <- sum(!raw_reps$ok)
if (n_failed > 0) cat(sprintf("[WARN] %d / %d reps failed -- see 'ok'/error rows for details.\n", n_failed, nrow(raw_reps)))

# ---------------------- CI summary (identical logic to the original) ----------------------
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

for (m in metrics) {
  sub <- ci_summary[metric == m]
  if (nrow(sub) == 0) next
  p <- ggplot(sub, aes(x = reorder(group, mean), y = mean)) +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), color = "steelblue") +
    coord_flip() +
    labs(x = NULL, y = m,
         title = sprintf("%s: subsampling %.0f%% CI (up to %d reps)", m, opt$ci * 100, opt$n_reps)) +
    theme_bw(base_size = 10)
  ggsave(file.path(opt$outdir, paste0("step3h_ci_", m, ".png")), p, width = 8, height = max(3, 0.4 * nrow(sub)))
}

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
