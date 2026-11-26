#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 6: Assemble one master metadata table for reporting -- pulls
# together everything a methods section / supplementary table needs:
# per-group cell count, k used (and whether it was corrected from the
# default 50), coaccessibility threshold, final summary metrics, and
# whether that group was flagged as a cell-count-driven outlier in the
# regression check.
#
# Usage:
#   Rscript step6_build_metadata_table.R \
#     --cell_counts cell_counts_vs_k.tsv \
#     --summary_stats cicero_qc_plots_corrected/step1_summary_stats.tsv \
#     --residuals approach1_regression_check_corrected/residuals_summary.tsv \
#     --k_corrections cicero_corrected/step5_rescale_k_summary.tsv \
#     --threshold 0.2 \
#     --outdir cicero_reporting
#
# Only --summary_stats is required; the others are optional and are
# simply skipped (with a note) if not supplied or not found.
#
# Outputs (in --outdir):
#   metadata_table_full.tsv   -- full precision, every column, one row
#                                 per group -- use this as your
#                                 supplementary table
#   metadata_table_report.md  -- rounded, publication-style markdown
#                                 table with the columns a methods
#                                 section / figure legend typically needs
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--cell_counts", type = "character", default = NA_character_,
              help = "cell_counts_vs_k.tsv (group, n_cells, k, pct_cells_per_metacell)"),
  make_option("--summary_stats", type = "character",
              default = "cicero_qc_plots_corrected/step1_summary_stats.tsv",
              help = "the (corrected) step1_summary_stats.tsv -- required"),
  make_option("--residuals", type = "character", default = NA_character_,
              help = "residuals_summary.tsv from approach1_cellcount_regression.R"),
  make_option("--k_corrections", type = "character", default = NA_character_,
              help = "step5_rescale_k_summary.tsv, to note which groups had k corrected and from/to what"),
  make_option("--threshold", type = "double", default = 0.2,
              help = "coaccessibility threshold used (constant across groups) [default %default]"),
  make_option("--default_k", type = "integer", default = 50,
              help = "the default/uncorrected k value, for flagging which groups deviate from it [default %default]"),
  make_option("--outdir", type = "character", default = "cicero_reporting")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(opt$summary_stats)) stop("--summary_stats file not found: ", opt$summary_stats)
meta <- fread(opt$summary_stats)

# ---------------- derive lineage / timepoint from group name ----------------
# group names are "<lineage>_prenatal" or "<lineage>_postnatal"
meta[, timepoint := ifelse(grepl("_prenatal$", group), "prenatal",
                     ifelse(grepl("_postnatal$", group), "postnatal", NA_character_))]
meta[, lineage := sub("_(prenatal|postnatal)$", "", group)]

# ---------------- attach cell counts ----------------
if (!is.na(opt$cell_counts) && file.exists(opt$cell_counts)) {
  cells <- fread(opt$cell_counts)[, .(group, n_cells)]
  meta <- merge(meta, cells, by = "group", all.x = TRUE)
  cat("Attached cell counts from", opt$cell_counts, "\n")
} else {
  meta[, n_cells := NA_integer_]
  cat("NOTE: --cell_counts not supplied/found -- n_cells will be NA. Reporting table\n",
      "      will still work but is much more useful with this attached.\n", sep = "")
}

# ---------------- attach k used / correction status ----------------
if (!is.na(opt$k_corrections) && file.exists(opt$k_corrections)) {
  kcorr <- fread(opt$k_corrections)
  k_col <- if ("k_used" %in% names(kcorr)) "k_used" else if ("k" %in% names(kcorr)) "k" else NA
  if (!is.na(k_col)) {
    kcorr_small <- kcorr[, .(group, k_used = get(k_col))]
    meta <- merge(meta, kcorr_small, by = "group", all.x = TRUE)
  } else {
    meta[, k_used := NA_integer_]
  }
  cat("Attached corrected k values from", opt$k_corrections, "\n")
} else {
  meta[, k_used := NA_integer_]
  cat("NOTE: --k_corrections not supplied/found -- assuming default k for all groups.\n")
}
meta[, k_used := ifelse(is.na(k_used), opt$default_k, k_used)]
meta[, k_corrected := k_used != opt$default_k]
meta[, coaccess_threshold := opt$threshold]

# ---------------- attach regression outlier flags ----------------
if (!is.na(opt$residuals) && file.exists(opt$residuals)) {
  resid <- fread(opt$residuals)
  resid_wide <- dcast(resid, group ~ metric, value.var = "resid_sd")
  setnames(resid_wide, setdiff(names(resid_wide), "group"), paste0("resid_sd_", setdiff(names(resid_wide), "group")))
  meta <- merge(meta, resid_wide, by = "group", all.x = TRUE)

  outlier_wide <- dcast(resid, group ~ metric, value.var = "is_outlier")
  setnames(outlier_wide, setdiff(names(outlier_wide), "group"), paste0("outlier_", setdiff(names(outlier_wide), "group")))
  meta <- merge(meta, outlier_wide, by = "group", all.x = TRUE)

  outlier_cols <- grep("^outlier_", names(meta), value = TRUE)
  meta[, any_cellcount_outlier := Reduce(`|`, lapply(.SD, function(x) isTRUE(x) | (is.logical(x) & x)), accumulate = FALSE), .SDcols = outlier_cols]
  meta[, any_cellcount_outlier := apply(meta[, ..outlier_cols], 1, function(r) any(as.logical(r), na.rm = TRUE))]
  cat("Attached regression residuals/outlier flags from", opt$residuals, "\n")
} else {
  cat("NOTE: --residuals not supplied/found -- outlier flag columns omitted.\n")
}

# ---------------- reorder columns sensibly ----------------
front_cols <- c("group", "lineage", "timepoint", "n_cells", "k_used", "k_corrected", "coaccess_threshold")
front_cols <- intersect(front_cols, names(meta))
setcolorder(meta, c(front_cols, setdiff(names(meta), front_cols)))
meta <- meta[order(lineage, timepoint)]

fwrite(meta, file.path(opt$outdir, "metadata_table_full.tsv"), sep = "\t")
cat("\nWrote metadata_table_full.tsv (", nrow(meta), "groups,", ncol(meta), "columns)\n")

# ---------------- rounded publication-style markdown table ----------------
report_cols <- intersect(c(
  "group", "n_cells", "k_used", "k_corrected", "coaccess_threshold",
  "n_connections", "n_ccans", "mean_ccan_size", "median_ccan_size",
  "median_distance_bp", "any_cellcount_outlier"
), names(meta))
report <- meta[, ..report_cols]

round_if_numeric <- function(x, digits = 2) if (is.numeric(x)) round(x, digits) else x
report[, (names(report)) := lapply(.SD, round_if_numeric)]

md_lines <- c(
  paste0("| ", paste(names(report), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(report)), collapse = "|"), "|"),
  apply(report, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
)
writeLines(md_lines, file.path(opt$outdir, "metadata_table_report.md"))
cat("Wrote metadata_table_report.md (rounded, publication-style)\n")

cat("\n--- Preview ---\n")
print(report)

cat(sprintf("\nDone. Outputs in: %s/\n", opt$outdir))
