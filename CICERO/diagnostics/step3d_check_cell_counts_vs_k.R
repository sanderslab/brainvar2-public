#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3d: Check cell counts per group vs. the fixed Cicero k parameter.
#
# Rationale: 02_run_cicero_from_export.R uses a FIXED k=50 for every group's
# make_cicero_cds() aggregation step, regardless of how many cells that
# group has (01_export_for_cicero.R allows anywhere from 500 to 5000 cells
# per group). A group near the 500-cell floor has each k=50 metacell
# drawn from 10% of its cell pool -> heavy overlap between metacells ->
# artificially inflated co-accessibility. A group near the 5000-cell cap
# has each metacell drawn from only 1% of its pool -> much less overlap.
#
# This script just counts cell_barcodes.txt lines per group directory
# in the cicero_export/ folder and reports k/n_cells ratio, sorted so
# you can immediately see if Micro_prenatal is a floor-count outlier.
#
# Usage:
#   Rscript step3d_check_cell_counts_vs_k.R --export_dir ./cicero_export --k 50
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--export_dir", type = "character", default = "./cicero_export",
              help = "Path to the cicero_export directory (per-group subfolders) [default %default]"),
  make_option("--k", type = "integer", default = 50,
              help = "k used in make_cicero_cds() [default %default]"),
  make_option("--out", type = "character", default = "cell_counts_vs_k.tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!dir.exists(opt$export_dir)) {
  stop("export_dir not found: ", opt$export_dir,
       "\nPoint --export_dir at the folder containing per-group subfolders ",
       "(each with cell_barcodes.txt), i.e. wherever 01_export_for_cicero.R wrote its output.")
}

group_dirs <- list.dirs(opt$export_dir, full.names = TRUE, recursive = FALSE)
cat(sprintf("Found %d group directories in %s\n", length(group_dirs), opt$export_dir))

results <- rbindlist(lapply(group_dirs, function(d) {
  grp <- basename(d)
  barcode_file <- file.path(d, "cell_barcodes.txt")
  if (!file.exists(barcode_file)) {
    cat(sprintf("  [SKIP] %s -- no cell_barcodes.txt found\n", grp))
    return(NULL)
  }
  n_cells <- length(readLines(barcode_file, warn = FALSE))
  data.table(group = grp, n_cells = n_cells)
}))

results[, k := opt$k]
results[, pct_cells_per_metacell := 100 * k / n_cells]
results <- results[order(-pct_cells_per_metacell)]

fwrite(results, opt$out, sep = "\t")
cat(sprintf("\nWrote %s\n\n", opt$out))
print(results)

cat("\nHow to read this:\n")
cat("  - pct_cells_per_metacell = what fraction of that group's ENTIRE cell pool\n")
cat("    goes into each single k-nearest-neighbor metacell.\n")
cat("  - Higher values (near the 500-cell floor) => metacells overlap heavily with\n")
cat("    each other => inflated, less-independent co-accessibility estimates.\n")
cat("  - Lower values (near the 5000-cell cap) => metacells are much more distinct\n")
cat("    from each other => more reliable estimates.\n")
cat("  - If Micro_prenatal sits at/near the top of this table (fewest cells,\n")
cat("    highest overlap fraction), that supports the metacell-overlap explanation\n")
cat("    for its inflated edge count and CCAN sizes.\n")
