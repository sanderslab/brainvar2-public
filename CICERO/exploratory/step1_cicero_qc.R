#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Step 1: Descriptive characterization of Cicero co-accessibility results.
#
# Usage:
#   Rscript step1_cicero_qc.R --indir /path/to/v3_coacc_cicero --outdir cicero_qc_plots
#
# Expects (per group, e.g. "ExN_postnatal", "Micro_prenatal", ...):
#   cicero_<group>_connections.tsv   columns: Peak1  Peak2  coaccess  pair_type
#   cicero_<group>_ccans.tsv         columns: Peak   CCAN
#
# Peak format assumed: chr_start_end  (e.g. chr10_100009708_100010208)
#
# Outputs (in --outdir):
#   step1_summary_stats.tsv
#   01_coaccess_score_distributions.png
#   02_distance_distributions.png
#   03_ccan_size_distributions.png
#   04_pairtype_composition.png
#   05_summary_metrics_across_groups.png
# ------------------------------------------------------------------

suppressMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(tidyr)
})

# ---------------------- args ----------------------
option_list <- list(
  make_option("--indir", type = "character", default = ".",
              help = "Directory containing cicero_*_connections.tsv / _ccans.tsv"),
  make_option("--outdir", type = "character", default = "cicero_qc_plots",
              help = "Output directory for plots and summary table"),
  make_option("--threshold", type = "double", default = 0.2,
              help = "coaccess threshold used to compute 'pass' fraction [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- helpers ----------------------
peak_midpoint <- function(peak) {
  parts <- str_split_fixed(peak, "_", 3)
  (as.numeric(parts[, 2]) + as.numeric(parts[, 3])) / 2
}

peak_chrom <- function(peak) {
  str_split_fixed(peak, "_", 3)[, 1]
}

get_groups <- function(indir) {
  files <- list.files(indir, pattern = "^cicero_.*_connections\\.tsv$", full.names = FALSE)
  gsub("^cicero_|_connections\\.tsv$", "", files) |> sort()
}

load_connections <- function(indir, group) {
  f <- file.path(indir, paste0("cicero_", group, "_connections.tsv"))
  df <- fread(f, sep = "\t")
  df <- df[!is.na(coaccess)]
  chrom1 <- peak_chrom(df$Peak1)
  chrom2 <- peak_chrom(df$Peak2)
  mid1 <- peak_midpoint(df$Peak1)
  mid2 <- peak_midpoint(df$Peak2)
  df$distance <- ifelse(chrom1 == chrom2, abs(mid1 - mid2), NA_real_)
  df$group <- group
  df
}

load_ccans <- function(indir, group) {
  f <- file.path(indir, paste0("cicero_", group, "_ccans.tsv"))
  df <- fread(f, sep = "\t")
  df$group <- group
  df
}

# ---------------------- load everything ----------------------
groups <- get_groups(opt$indir)
cat(sprintf("Found %d groups: %s\n", length(groups), paste(groups, collapse = ", ")))

all_conn <- rbindlist(lapply(groups, function(g) {
  cat("Loading connections:", g, "\n")
  load_connections(opt$indir, g)
}), fill = TRUE)

all_ccan <- rbindlist(lapply(groups, function(g) {
  cat("Loading ccans:", g, "\n")
  load_ccans(opt$indir, g)
}), fill = TRUE)

# ---------------------- summary table ----------------------
ccan_sizes <- all_ccan[, .N, by = .(group, CCAN)]

summary_df <- all_conn[, .(
  n_connections = as.numeric(.N),
  n_connections_pass_thresh = as.numeric(sum(abs(coaccess) >= opt$threshold)),
  median_coaccess = as.numeric(median(coaccess)),
  median_distance_bp = as.numeric(median(distance, na.rm = TRUE))
), by = group]
summary_df[, frac_pass_thresh := n_connections_pass_thresh / n_connections]

ccan_summary <- ccan_sizes[, .(
  n_ccans = as.numeric(uniqueN(CCAN)),
  mean_ccan_size = as.numeric(mean(N)),
  median_ccan_size = as.numeric(median(N)),
  max_ccan_size = as.numeric(max(N))
), by = group]

peaks_in_ccans <- all_ccan[, .(n_peaks_in_ccans = uniqueN(Peak)), by = group]

summary_df <- summary_df |>
  left_join(ccan_summary, by = "group") |>
  left_join(peaks_in_ccans, by = "group") |>
  arrange(desc(n_connections))

fwrite(summary_df, file.path(opt$outdir, "step1_summary_stats.tsv"), sep = "\t")
cat("Wrote step1_summary_stats.tsv\n")
print(summary_df)

# ---------------------- plot theme ----------------------
theme_qc <- theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 7),
        axis.text = element_text(size = 6))

# ---------------------- 01: coaccess score distributions ----------------------
p1 <- ggplot(all_conn, aes(x = coaccess)) +
  geom_histogram(bins = 50, fill = "steelblue", color = NA) +
  facet_wrap(~group, ncol = 5, scales = "free_y") +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Co-accessibility score", y = "Count") +
  theme_qc

ggsave(file.path(opt$outdir, "01_coaccess_score_distributions.png"), p1,
       width = 20, height = 4 * ceiling(length(groups) / 5), limitsize = FALSE)
cat("Wrote 01_coaccess_score_distributions.png\n")

# ---------------------- 02: distance distributions (log10) ----------------------
dist_df <- all_conn[!is.na(distance) & distance > 0]
dist_df[, log10_distance := log10(distance)]

p2 <- ggplot(dist_df, aes(x = log10_distance)) +
  geom_histogram(bins = 50, fill = "seagreen", color = NA) +
  facet_wrap(~group, ncol = 5, scales = "free_y") +
  labs(x = "log10(distance, bp)", y = "Count") +
  theme_qc

ggsave(file.path(opt$outdir, "02_distance_distributions.png"), p2,
       width = 20, height = 4 * ceiling(length(groups) / 5), limitsize = FALSE)
cat("Wrote 02_distance_distributions.png\n")

# ---------------------- 03: CCAN size distributions ----------------------
p3 <- ggplot(ccan_sizes, aes(x = N)) +
  geom_histogram(bins = 30, fill = "darkorange", color = NA) +
  facet_wrap(~group, ncol = 5, scales = "free_y") +
  labs(x = "Peaks per CCAN", y = "Count") +
  theme_qc

ggsave(file.path(opt$outdir, "03_ccan_size_distributions.png"), p3,
       width = 20, height = 4 * ceiling(length(groups) / 5), limitsize = FALSE)
cat("Wrote 03_ccan_size_distributions.png\n")

# ---------------------- 04: pair type composition ----------------------
if ("pair_type" %in% names(all_conn)) {
  pt <- all_conn[, .N, by = .(group, pair_type)]
  pt[, frac := N / sum(N), by = group]

  p4 <- ggplot(pt, aes(x = group, y = frac, fill = pair_type)) +
    geom_col() +
    labs(x = NULL, y = "Fraction of connections", fill = "Pair type") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1))

  ggsave(file.path(opt$outdir, "04_pairtype_composition.png"), p4,
         width = max(8, 0.4 * length(groups)), height = 6)
  cat("Wrote 04_pairtype_composition.png\n")
} else {
  cat("No pair_type column found — skipping plot 04.\n")
}

# ---------------------- 05: summary metrics across groups ----------------------
metrics <- c("n_connections", "n_connections_pass_thresh", "median_coaccess",
             "median_distance_bp", "n_ccans", "mean_ccan_size")

summary_long <- summary_df |>
  select(group, all_of(metrics)) |>
  pivot_longer(-group, names_to = "metric", values_to = "value")

p5 <- ggplot(summary_long, aes(x = reorder(group, value), y = value)) +
  geom_col(fill = "darkorange") +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9)

ggsave(file.path(opt$outdir, "05_summary_metrics_across_groups.png"), p5,
       width = 18, height = 9)
cat("Wrote 05_summary_metrics_across_groups.png\n")

cat(sprintf("\nDone. All outputs in: %s/\n", opt$outdir))
