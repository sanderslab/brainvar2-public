#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 3b: CCAN genomic span / contiguity check.
#
# For each requested group, computes per-CCAN:
#   - n_peaks
#   - chromosome(s) involved (flags any CCAN spanning >1 chromosome,
#     which should never happen and would indicate a real bug)
#   - genomic span (max peak end - min peak start)
#   - density = n_peaks / (span_kb)  -> low density + huge span suggests
#     a "collapsed mega-region" rather than a real tight regulatory module
#
# Usage:
#   Rscript step3b_ccan_span_check.R --indir /path/to/v3_coacc_cicero \
#       --groups Micro_prenatal,Micro_postnatal --outdir cicero_ccan_span_check
#
# Outputs (in --outdir):
#   <group>_ccan_span_table.tsv   - full per-CCAN table, sorted by n_peaks desc
#   ALL_top20_ccans_combined.tsv  - top 20 largest CCANs per group, combined
#   span_vs_npeaks_scatter.png    - n_peaks vs span_bp, colored by group,
#                                   log-log scale (flat/shallow trend for a
#                                   group = big CCANs are mostly just wide,
#                                   sparse regions, not densely-linked modules)
#   density_comparison_boxplot.png- peaks-per-kb density distribution per group
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(dplyr)
})

option_list <- list(
  make_option("--indir", type = "character", default = "."),
  make_option("--outdir", type = "character", default = "cicero_ccan_span_check"),
  make_option("--groups", type = "character", default = "Micro_prenatal,Micro_postnatal",
              help = "Comma-separated list of groups to check [default %default]"),
  make_option("--top_n", type = "integer", default = 20,
              help = "How many largest CCANs per group to report in detail [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

groups <- str_split(opt$groups, ",")[[1]] |> str_trim()
cat(sprintf("Checking %d groups: %s\n", length(groups), paste(groups, collapse = ", ")))

parse_peak <- function(peak) {
  parts <- str_split_fixed(peak, "_", 3)
  data.table(chrom = parts[, 1], start = as.numeric(parts[, 2]), end = as.numeric(parts[, 3]))
}

all_span_tables <- list()
all_top20 <- list()

for (g in groups) {
  cat(sprintf("\n=== %s ===\n", g))
  f <- file.path(opt$indir, paste0("cicero_", g, "_ccans.tsv"))
  ccan <- fread(f, sep = "\t")

  parsed <- parse_peak(ccan$Peak)
  ccan <- cbind(ccan, parsed)

  span_table <- ccan[, .(
    n_peaks = .N,
    n_chroms = uniqueN(chrom),
    chroms_involved = paste(sort(unique(chrom)), collapse = ","),
    min_start = min(start),
    max_end = max(end)
  ), by = CCAN]

  span_table[, span_bp := max_end - min_start]
  span_table[, span_kb := span_bp / 1000]
  span_table[, density_peaks_per_kb := n_peaks / span_kb]
  span_table[, group := g]
  span_table <- span_table[order(-n_peaks)]

  # flag any CCAN spanning more than one chromosome (should not happen; real bug if so)
  n_multichrom <- sum(span_table$n_chroms > 1)
  if (n_multichrom > 0) {
    cat(sprintf("  !! WARNING: %d CCAN(s) span multiple chromosomes -- investigate separately.\n", n_multichrom))
  } else {
    cat("  OK: every CCAN is confined to a single chromosome.\n")
  }

  fwrite(span_table, file.path(opt$outdir, paste0(g, "_ccan_span_table.tsv")), sep = "\t")
  cat(sprintf("  Wrote %s_ccan_span_table.tsv (%d CCANs total)\n", g, nrow(span_table)))

  top20 <- head(span_table, opt$top_n)
  cat(sprintf("  Top %d CCANs by peak count:\n", opt$top_n))
  print(top20[, .(CCAN, n_peaks, chroms_involved, span_kb = round(span_kb, 1),
                   density_peaks_per_kb = round(density_peaks_per_kb, 3))])

  all_span_tables[[g]] <- span_table
  all_top20[[g]] <- top20
}

combined_span <- rbindlist(all_span_tables)
combined_top20 <- rbindlist(all_top20)
fwrite(combined_top20, file.path(opt$outdir, "ALL_top20_ccans_combined.tsv"), sep = "\t")
cat("\nWrote ALL_top20_ccans_combined.tsv\n")

# ---------------- plot: span vs n_peaks (log-log), colored by group ----------------
p_scatter <- ggplot(combined_span, aes(x = span_bp, y = n_peaks, color = group)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "CCAN size vs genomic span",
       x = "Genomic span (bp, log scale)", y = "Peaks per CCAN (log scale)") +
  theme_bw(base_size = 11)
ggsave(file.path(opt$outdir, "span_vs_npeaks_scatter.png"), p_scatter, width = 7, height = 6)
cat("Wrote span_vs_npeaks_scatter.png\n")

# ---------------- plot: density comparison ----------------
p_density <- ggplot(combined_span, aes(x = group, y = density_peaks_per_kb, fill = group)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10() +
  labs(title = "CCAN density (peaks per kb of genomic span)",
       x = NULL, y = "Peaks / kb (log scale)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(opt$outdir, "density_comparison_boxplot.png"), p_density, width = 6, height = 5)
cat("Wrote density_comparison_boxplot.png\n")

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
cat("\nHow to read this:\n")
cat("  - If the largest CCANs in one group have MUCH larger span but LOWER density\n")
cat("    than equivalently-sized CCANs in other groups, that group's big modules are\n")
cat("    likely collapsed/sparse mega-regions rather than tight regulatory modules.\n")
cat("  - Any n_chroms > 1 would indicate a parsing/pipeline bug and should be examined directly.\n")
