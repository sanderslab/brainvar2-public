# =============================================================================
# violin_plot.R  –  sQTL PSI violin plots with real allele labels
#
# Input TSV columns (written by plot_sqtl.py):
#   sample_id | PSI | genotype | allele_label | stage
#
# Usage (command line):
#   Rscript violin_plot.R <input_dir_or_file> [output_dir]
#
# Usage (interactive):
#   source("violin_plot.R")
#   plot_violin_batch("~/violin_export", output_dir = "~/violin_export/R_plots")
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

use_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

# -----------------------------------------------------------------
# Colour palette: one colour per unique allele label.
# Labels are ordered by genotype (ref/ref → het → alt/alt).
# -----------------------------------------------------------------
GENO_PALETTE <- c(
  "0" = "#4DAF4A",   # homozygous ref  → green
  "1" = "#FF7F00",   # heterozygous    → orange
  "2" = "#E41A1C"    # homozygous alt  → red
)

MIN_N_TEST <- 3

# =============================================================================
# SIGNIFICANCE HELPERS
# =============================================================================

sig_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

sig_label <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return(sprintf("p=%.1e", p))
  if (p < 0.05)  return(sprintf("p=%.3f", p))
  return(sprintf("p=%.2f", p))
}

# =============================================================================
# build_allele_levels()
#   Returns allele labels ordered by median genotype (ref→het→alt)
# =============================================================================

build_allele_levels <- function(df) {
  df %>%
    group_by(allele_label) %>%
    summarise(median_num = median(genotype, na.rm = TRUE), .groups = "drop") %>%
    arrange(median_num) %>%
    pull(allele_label)
}

# =============================================================================
# compute_pairwise()
# =============================================================================

compute_pairwise <- function(sub_df, allele_levels, y_max, y_range) {
  
  present <- allele_levels[allele_levels %in% unique(as.character(sub_df$allele_label))]
  if (length(present) < 2) return(data.frame())
  
  pairs <- combn(present, 2, simplify = FALSE)
  widths <- sapply(pairs, function(p)
    abs(which(present == p[2]) - which(present == p[1])))
  pairs <- pairs[order(widths)]
  
  bracket_base <- y_max * 1.02
  bracket_step <- y_range * 0.14
  tick_drop    <- y_range * 0.015
  
  rows <- list()
  
  for (k in seq_along(pairs)) {
    g1 <- pairs[[k]][1]
    g2 <- pairs[[k]][2]
    
    v1 <- sub_df$PSI[sub_df$allele_label == g1]
    v2 <- sub_df$PSI[sub_df$allele_label == g2]
    
    if (length(v1) < MIN_N_TEST || length(v2) < MIN_N_TEST) next
    
    pval <- tryCatch(wilcox.test(v1, v2)$p.value, error = function(e) NA_real_)
    
    x1      <- which(allele_levels == g1)
    x2      <- which(allele_levels == g2)
    y_brack <- bracket_base + (length(rows)) * bracket_step
    
    rows[[length(rows) + 1]] <- data.frame(
      group1    = g1,
      group2    = g2,
      x1        = x1,
      x2        = x2,
      x_mid     = (x1 + x2) / 2,
      pval      = pval,
      stars     = sig_stars(pval),
      p_label   = sig_label(pval),
      y_bracket = y_brack,
      y_tick_lo = y_brack - tick_drop,
      y_stars   = y_brack + y_range * 0.03,
      y_pval    = y_brack + y_range * 0.07,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

# =============================================================================
# plot_violin()
#   df columns: sample_id, PSI, genotype, allele_label, stage
# =============================================================================

plot_violin <- function(df, title = NULL, output_path = NULL) {
  
  allele_levels <- build_allele_levels(df)
  
  df <- df %>%
    mutate(allele_label = factor(allele_label, levels = allele_levels))
  
  if (nrow(df) == 0) return(NULL)
  
  # map each allele_label to a palette colour via genotype
  num_map <- df %>%
    group_by(allele_label) %>%
    summarise(num = as.character(round(median(genotype, na.rm = TRUE))),
              .groups = "drop")
  
  geno_colors <- setNames(
    GENO_PALETTE[num_map$num],
    num_map$allele_label
  )
  
  y_min_global <- floor(min(df$PSI, na.rm = TRUE))
  y_max_global <- ceiling(max(df$PSI, na.rm = TRUE))
  y_range      <- max(y_max_global - y_min_global, 0.01)
  y_buffer     <- y_range * 0.10
  
  # Per-stage significance brackets
  stages <- unique(df$stage)
  
  pval_list <- lapply(stages, function(s) {
    sub <- filter(df, stage == s)
    res <- compute_pairwise(sub, allele_levels, y_max_global, y_range)
    if (nrow(res) > 0) res$stage <- s
    res
  })
  pval_df <- do.call(rbind, Filter(function(x) nrow(x) > 0, pval_list))
  
  y_top <- if (!is.null(pval_df) && nrow(pval_df) > 0)
    max(pval_df$y_pval, na.rm = TRUE) + y_range * 0.08
  else
    y_max_global * 1.10
  
  n_df <- df %>%
    group_by(stage, allele_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(y_pos = y_min_global - y_buffer * 0.5,
           label = paste0("n=", n))
  
  p <- ggplot(df, aes(x = allele_label, y = PSI,
                      fill = allele_label, colour = allele_label)) +
    
    geom_violin(trim = FALSE, alpha = 0.40, scale = "width", colour = NA) +
    
    stat_summary(fun = median, geom = "crossbar",
                 width = 0.4, linewidth = 0.55,
                 color = "grey20", fatten = 0) +
    
    {
      if (use_beeswarm)
        ggbeeswarm::geom_beeswarm(size = 0.9, alpha = 0.55,
                                  colour = "black", show.legend = FALSE)
      else
        geom_jitter(width = 0.12, size = 0.9, alpha = 0.55,
                    colour = "black", show.legend = FALSE)
    } +
    
    geom_text(data = n_df,
              aes(x = allele_label, y = y_pos, label = label),
              inherit.aes = FALSE, size = 2.4, colour = "grey50") +
    
    scale_fill_manual(values = geno_colors, name = "Genotype") +
    scale_colour_manual(values = geno_colors, guide = "none") +
    
    coord_cartesian(ylim = c(y_min_global - y_buffer, y_top)) +
    
    facet_wrap(~stage, nrow = 1) +
    
    labs(
      title = title,
      x     = "Genotype",
      y     = "Normalized PSI"
    ) +
    
    theme_bw(base_size = 10) +
    theme(
      plot.title      = element_text(size = 8, hjust = 0.5),
      legend.position = "right"
    )
  
  # Add significance brackets if available
  if (!is.null(pval_df) && nrow(pval_df) > 0 &&
      all(c("x1", "x2") %in% colnames(pval_df))) {
    
    p <- p +
      geom_segment(data = pval_df,
                   aes(x = x1, xend = x1, y = y_tick_lo, yend = y_bracket),
                   inherit.aes = FALSE, linewidth = 0.35) +
      geom_segment(data = pval_df,
                   aes(x = x2, xend = x2, y = y_tick_lo, yend = y_bracket),
                   inherit.aes = FALSE, linewidth = 0.35) +
      geom_segment(data = pval_df,
                   aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
                   inherit.aes = FALSE, linewidth = 0.35) +
      geom_text(data = pval_df,
                aes(x = x_mid, y = y_stars, label = stars),
                inherit.aes = FALSE, size = 3) +
      geom_text(data = pval_df,
                aes(x = x_mid, y = y_pval, label = p_label),
                inherit.aes = FALSE, size = 2.5)
  }
  
  if (!is.null(output_path)) {
    pdf(output_path, width = 8, height = 5)
    print(p)
    dev.off()
    message("Saved: ", output_path)
  }
  
  p
}

# =============================================================================
# FILE WRAPPERS
# =============================================================================

plot_violin_from_file <- function(tsv_path, output_dir = NULL) {
  
  df <- read.delim(tsv_path, stringsAsFactors = FALSE)
  
  # derive a plot title from filename
  title <- sub("_violin_data\\.tsv$", "", basename(tsv_path))
  title <- gsub("__", "\n", title)
  
  out_pdf <- NULL
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    base    <- sub("_violin_data\\.tsv$", "", basename(tsv_path))
    out_pdf <- file.path(output_dir, paste0(base, "_R.pdf"))
  }
  
  if (!is.null(out_pdf) && file.exists(out_pdf)) {
    message("Skipping (exists): ", out_pdf)
    return(invisible(NULL))
  }
  
  plot_violin(df, title = title, output_path = out_pdf)
}

plot_violin_batch <- function(input_dir, output_dir = NULL) {
  
  files <- list.files(input_dir,
                      pattern = "_violin_data\\.tsv$",
                      full.names = TRUE)
  
  if (length(files) == 0) {
    message("No *_violin_data.tsv files found in: ", input_dir)
    return(invisible(NULL))
  }
  
  if (is.null(output_dir))
    output_dir <- file.path(input_dir, "R_plots")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf("Processing %d file(s) → %s", length(files), output_dir))
  
  for (f in files) {
    message("  ", basename(f))
    try(plot_violin_from_file(f, output_dir))
  }
}

# =============================================================================
# CLI
# =============================================================================

if (!interactive()) {
  args  <- commandArgs(trailingOnly = TRUE)
  input <- args[1]
  out   <- if (length(args) >= 2) args[2] else NULL
  
  if (dir.exists(input)) {
    plot_violin_batch(input, out)
  } else {
    plot_violin_from_file(input, out)
  }
}

# =============================================================================
# Interactive examples (uncomment to run):
plot_violin_batch("/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/sQTL/violin_export",
                  output_dir = "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/sQTL/violin_export/R_plots")
# =============================================================================