#!/usr/bin/env Rscript
# =============================================================================
# Leafcutter sashimi plot — genotype-stratified per stage
# Usage: Rscript leafcutter_sashimi.R <TSV>
#
# TSV filename format: ...rsXXXXXX__CHROM_START_END_clu_XXXXX...
# TSV columns:         sample_id, genotype, stage
#
# Uses make_differential_splicing_plot() which supports any number of groups.
# make_cluster_plot() only supports exactly 2 groups so is not used here.
#
# Inputs from leafviz RData:
#   counts  — junctions x samples integer matrix
#   introns — data frame with clusterID, start, end, verdict, deltapsi cols
#   (exons_table loaded separately from gencode file)
#
# Output: one PDF per stage, stacked panels one per genotype (HomRef/Het/HomAlt)
# =============================================================================

suppressPackageStartupMessages({
  library(foreach)
  library(reshape2)
  library(Hmisc)
  library(gridExtra)
  library(ggplot2)
})

# Source only the function we need directly from the cloned repo
# This avoids needing rstan/TailRank which are only used by other leafcutter functions
LEAFCUTTER_R_DIR <- path.expand("~/leafcutter/leafcutter/R")
source(file.path(LEAFCUTTER_R_DIR, "junction_plot.R"))   # make_differential_splicing_plot
source(file.path(LEAFCUTTER_R_DIR, "utils.R"))           # get_intron_meta

# Patch make_differential_splicing_plot to fix "invalid value 0 for 'digits'"
# This happens when a junction has zero mean counts across a group.
# The bug is in: format(intron_meta$counts[i], digits=2, scientific=FALSE)
# when counts[i] == 0, digits resolves to 0 which format() rejects.
# Fix: wrap summary_func to add a tiny pseudocount preventing exact zeros.
.make_differential_splicing_plot_orig <- make_differential_splicing_plot
make_differential_splicing_plot <- function(y, x = numeric(nrow(y)) + 1,
                                            exons_table = NULL, len = 500,
                                            length_transform = function(g) log(g + 1),
                                            main_title = NA, snp_pos = NA,
                                            summary_func = colMeans,
                                            legend_title = "Mean count") {
  # Wrap summary_func to ensure no exact zeros reach format()
  safe_summary <- function(mat) {
    vals <- summary_func(mat)
    pmax(vals, 0.01)   # floor at 0.01 — has 1 sig digit, format(digits=2) works
  }
  .make_differential_splicing_plot_orig(
    y = y, x = x, exons_table = exons_table, len = len,
    length_transform = length_transform, main_title = main_title,
    snp_pos = snp_pos, summary_func = safe_summary,
    legend_title = legend_title
  )
}

# =============================================================================
# CONFIG — update these paths
# =============================================================================

LEAFVIZ_RDATA <- "/mnt/s3/FXN/SoowheeKim/leafviz_filtered.RData"
COUNTS_FILE   <- "/mnt/s3/FXN/SoowheeKim/cluster/leafcutter_perind.counts.nochr.gz"
EXONS_FILE    <- "/mnt/s3/FXN/SoowheeKim/gencode_v41_all_exons.txt.gz"

# =============================================================================
# STEP 1: Parse TSV filename for rsid, chrom, coords, cluster_id
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript leafcutter_sashimi.R <TSV>")
tsv_file <- args[1]

m <- regmatches(basename(tsv_file),
                regexec("(rs\\d+)__(\\d+)_(\\d+)_(\\d+).*(clu_\\d+)",
                        basename(tsv_file)))[[1]]
if (length(m) == 0) stop("Filename format not recognised: ", basename(tsv_file))

rsid       <- m[2]
chrom      <- m[3]
sqtl_start <- as.integer(m[4])
sqtl_end   <- as.integer(m[5])
cluster_id <- m[6]

cat(sprintf("[parse] rsid=%s chrom=%s start=%d end=%d cluster=%s\n",
            rsid, chrom, sqtl_start, sqtl_end, cluster_id))

# Output directory named after sQTL: rsid_chr:start-end_clusterid
OUT_DIR <- file.path(".", sprintf("%s_chr%s_%d-%d_%s",
                                  rsid, chrom, sqtl_start, sqtl_end, cluster_id))
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
cat(sprintf("[parse] Output directory: %s\n", OUT_DIR))

# =============================================================================
# STEP 2: Load genotype TSV
# =============================================================================

cat("[load] Genotype TSV:", tsv_file, "\n")
geno_df <- read.table(tsv_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
geno_df$genotype <- as.integer(geno_df$genotype)
cat(sprintf("[load] %d samples | stages: %s\n",
            nrow(geno_df), paste(unique(geno_df$stage), collapse = ", ")))

# =============================================================================
# STEP 3: Load leafviz RData
# =============================================================================

cat("[load] Loading RData:", LEAFVIZ_RDATA, "\n")
e <- new.env()
load(LEAFVIZ_RDATA, envir = e)
cat("[load] Objects:", paste(ls(e), collapse = ", "), "\n\n")

# Print summary of each object so we know exactly what's available
for (nm in ls(e)) {
  obj <- get(nm, envir = e)
  cat(sprintf("  [%s] class=%s", nm, paste(class(obj), collapse = "/")))
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat(sprintf(" dim=%dx%d cols: %s",
                nrow(obj), ncol(obj),
                paste(head(colnames(obj), 8), collapse = ", ")))
    if (is.data.frame(obj) && nrow(obj) > 0)
      cat(sprintf("\n         row1: %s",
                  paste(head(as.character(unlist(obj[1,])), 6), collapse = " | ")))
  } else if (is.list(obj)) {
    cat(sprintf(" length=%d names: %s",
                length(obj), paste(head(names(obj), 6), collapse = ", ")))
  }
  cat("\n")
}
cat("\n")

# =============================================================================
# STEP 4: Extract counts matrix and introns data frame
#
# leafviz typically stores:
#   counts  — matrix, rows = junctions (chr:start:end:clu_ID), cols = samples
#   introns — data frame with clusterID, start, end, verdict, deltapsi
# =============================================================================

# =============================================================================
# STEP 4: Extract introns from RData + load counts from perind file
#
# This RData contains:
#   all.introns  — data frame: chr, start, end, clusterID (285k junctions)
#   verdict.list — list per cluster of junction verdicts
#   coord.list   — list per cluster of coordinates
# Counts are in the separate perind file (not in RData).
# =============================================================================

# Introns come from all.introns
introns_df   <- e$all.introns
clusterID_col <- "clusterID"

cat(sprintf("\n[extract] all.introns: %d rows | cols: %s\n",
            nrow(introns_df), paste(colnames(introns_df), collapse = ", ")))

# Filter to our cluster
cluster_introns <- introns_df[introns_df[[clusterID_col]] == cluster_id, ]
if (nrow(cluster_introns) == 0)
  cluster_introns <- introns_df[grepl(cluster_id, introns_df[[clusterID_col]]), ]

cat(sprintf("[filter] %d introns for %s\n", nrow(cluster_introns), cluster_id))
if (nrow(cluster_introns) == 0)
  stop(sprintf("No introns for %s. Example clusterIDs: %s",
               cluster_id,
               paste(head(introns_df[[clusterID_col]], 5), collapse = ", ")))
print(cluster_introns)

# Load counts from perind file
# Format: space-separated, row name = junction ID (e.g. 11:3039931:3040448:clu_20443_NA)
# Values are "numerator/denominator" — numerator = junction read count
cat("\n[load] Reading perind counts:", COUNTS_FILE, "\n")

# Read header preview
con <- gzfile(COUNTS_FILE, "r")
header_line <- readLines(con, n = 1)
close(con)
cat("[load] Header preview:", substr(header_line, 1, 200), "\n")

counts_raw <- read.table(gzfile(COUNTS_FILE), header = TRUE, sep = " ",
                         check.names = FALSE, stringsAsFactors = FALSE,
                         row.names = 1)

cat(sprintf("[load] Raw dims: %d x %d\n", nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("[load] Example value: %s\n", as.character(counts_raw[1, 1])))

# Extract numerator from "num/denom" fraction strings -> integer read counts
cat("[load] Extracting numerator read counts ...\n")
counts_mat <- apply(counts_raw, 2, function(col) {
  as.integer(sub("/.*", "", col))
})
rownames(counts_mat) <- rownames(counts_raw)

cat(sprintf("[load] counts_mat: %d junctions x %d samples\n",
            nrow(counts_mat), ncol(counts_mat)))
cat(sprintf("[load] rowname example: %s\n", rownames(counts_mat)[1]))
cat(sprintf("[load] colname examples: %s\n",
            paste(head(colnames(counts_mat), 3), collapse = " | ")))
cat(sprintf("[load] value example (should be integer): %d\n", counts_mat[1, 1]))

# Filter to cluster rows
cluster_count_rows <- grep(cluster_id, rownames(counts_mat), value = TRUE)
cat(sprintf("[filter] %d count rows for %s\n", length(cluster_count_rows), cluster_id))

if (length(cluster_count_rows) == 0)
  stop(sprintf("No rows in counts_mat match %s.\nExample rownames:\n  %s",
               cluster_id,
               paste(head(rownames(counts_mat), 5), collapse = "\n  ")))

cat(paste(" ", cluster_count_rows), sep = "\n")
counts_cluster <- counts_mat[cluster_count_rows, , drop = FALSE]

# =============================================================================
# STEP 5: Load exons table
# =============================================================================

cat("[load] Loading exons from file:", EXONS_FILE, "\n")
exons_table <- read.table(gzfile(EXONS_FILE), header = TRUE,
                          sep = "\t", stringsAsFactors = FALSE)
colnames(exons_table)[1] <- "chr"

# Junction IDs use no chr prefix (e.g. "11:3039931:...") so strip "chr" from
# exons_table to match — make_differential_splicing_plot filters on
# exons_table[exons_table$chr == intron_meta$chr[1], ] where intron_meta$chr
# is parsed directly from colnames(y)
exons_table$chr <- sub("^chr", "", exons_table$chr)

cat(sprintf("[load] Exons: %d rows | chr examples: %s\n",
            nrow(exons_table),
            paste(head(unique(exons_table$chr), 3), collapse = ", ")))

# =============================================================================
# STEP 6: Plot per stage using make_differential_splicing_plot
#
# make_differential_splicing_plot(y, x, exons_table, snp_pos, ...)
#   y   — [samples x junctions] matrix, colnames must be junction IDs
#   x   — [samples] grouping vector, one panel per unique value
#   snp_pos — single numeric position to mark with a vertical line
#
# Note: make_cluster_plot only supports exactly 2 groups (hardcodes plots[[1]]
# and plots[[2]]), so make_differential_splicing_plot is used here for 3
# genotype groups.
# =============================================================================

gt_labels <- c("0" = "HomRef (0/0)", "1" = "Het (0/1)", "2" = "HomAlt (1/1)")
snp_pos <- as.numeric((sqtl_start + sqtl_end) / 2)[1]

# Note: snp_pos marks the sQTL position with a vertical line.
# Set to NA to disable if causing errors with multi-junction clusters.
# snp_pos <- NA

for (stage in c("Prenatal", "Postnatal")) {

  cat(sprintf("\n[plot] ========== Stage: %s ==========\n", stage))

  df_stage <- geno_df[geno_df$stage == stage, ]
  if (nrow(df_stage) == 0) { cat("[plot] No samples, skipping.\n"); next }

  # ---- Match sample IDs to count matrix columns ----
  # Count matrix columns may be full paths or have BAM suffixes — normalise both
  strip_id <- function(x) {
    x <- sub("\\.STARAligned.*", "", x)   # strip STAR suffix first
    x <- sub("\\.bam$", "", x, ignore.case = TRUE)
    x <- basename(x)                      # strip any path
    x <- gsub("[-. ]", "_", x)           # normalise separators
    x
  }
  count_cols   <- colnames(counts_cluster)
  count_cols_n <- strip_id(count_cols)
  stage_ids_n  <- strip_id(df_stage$sample_id)

  cat("[plot] Count col examples (normalised):",
      paste(head(count_cols_n, 3), collapse = ", "), "\n")
  cat("[plot] TSV sample examples (normalised):",
      paste(head(stage_ids_n, 3), collapse = ", "), "\n")

  matched <- match(stage_ids_n, count_cols_n)
  found   <- !is.na(matched)

  if (sum(found) == 0) {
    cat("[plot] WARNING: no sample IDs matched count matrix columns.\n")
    cat("[plot] Count col examples: ", paste(head(count_cols, 3), collapse = ", "), "\n")
    cat("[plot] TSV sample examples:", paste(head(df_stage$sample_id, 3), collapse = ", "), "\n")
    next
  }
  cat(sprintf("[plot] Matched %d / %d samples.\n", sum(found), nrow(df_stage)))

  # ---- Build y [samples x junctions] and x [samples grouping] ----
  y <- t(counts_cluster[, matched[found], drop = FALSE])

  x <- factor(
    gt_labels[as.character(df_stage$genotype[found])],
    levels = c("HomRef (0/0)", "Het (0/1)", "HomAlt (1/1)")
  )
  x <- droplevels(x)

  cat("[plot] Group sizes:\n"); print(table(x))

  if (nlevels(x) < 2) {
    cat("[plot] Fewer than 2 genotype groups — skipping.\n"); next
  }

  cat(sprintf("[plot] y matrix: %d samples x %d junctions\n", nrow(y), ncol(y)))

  # ---- Compute PSI matrix [samples x junctions] ----
  # PSI = junction_count / cluster_total per sample
  # This is what shows the sQTL effect — genotype shifts the proportion, not raw counts
  row_totals <- rowSums(y)
  row_totals[row_totals == 0] <- 1
  psi_mat <- y / row_totals

  psi_summary <- do.call(rbind, lapply(levels(x), function(gt) {
    idx <- x == gt
    if (sum(idx) == 0) return(NULL)
    colMeans(psi_mat[idx, , drop = FALSE])
  }))
  rownames(psi_summary) <- levels(x)
  cat("[plot] Mean PSI per genotype:\n"); print(round(psi_summary, 3))

  # Delta PSI summary (HomAlt - HomRef)
  if (all(c("HomRef (0/0)", "HomAlt (1/1)") %in% rownames(psi_summary))) {
    delta_psi <- psi_summary["HomAlt (1/1)", ] - psi_summary["HomRef (0/0)", ]
    cat("[plot] Delta PSI (HomAlt - HomRef):\n")
    print(round(delta_psi, 3))
    cat(sprintf("[plot] sQTL junction delta PSI: %.3f\n",
                delta_psi[grep(paste0(sqtl_start, ".*", sqtl_end, "|",
                                      sqtl_end, ".*", sqtl_start),
                               names(delta_psi))]))
  }

  # ---- Scale PSI to counts-like values for plotting ----
  # Multiply PSI by a fixed constant (100) so arc sizes reflect
  # proportion not absolute counts — comparable across genotypes and stages
  # This means arc thickness = % of reads at this junction (0-100 scale)
  PSI_SCALE  <- 100
  y_plot     <- psi_mat * PSI_SCALE

  # ---- Junction labels: show which is the sQTL junction ----
  junc_names <- colnames(y)
  sqtl_junc  <- grep(paste0(sqtl_start, ".*", sqtl_end, "|",
                             sqtl_end, ".*", sqtl_start), junc_names)
  cat(sprintf("[plot] sQTL junction index: %s\n",
              if (length(sqtl_junc)) junc_names[sqtl_junc] else "not found"))

  out_pdf <- file.path(OUT_DIR, sprintf("%s_%s_sashimi.pdf", rsid, stage))
  cat("[plot] Writing:", out_pdf, "\n")

  tryCatch({
    pdf(out_pdf, width = 14, height = 4 * nlevels(x) + 3)

    # summary_func on PSI*100 values — colMeans gives mean % usage per junction
    # All panels now on same 0-100 scale so thickness is directly comparable
    norm_summary <- function(mat) pmax(colMeans(mat), 0.01)

    make_differential_splicing_plot(
      y              = y_plot,
      x              = x,
      exons_table    = exons_table,
      snp_pos        = snp_pos,
      summary_func   = norm_summary,
      legend_title   = "Mean PSI (%)",
      main_title     = sprintf("%s | %s | chr%s:%d-%d | %s",
                               rsid, stage, chrom, sqtl_start, sqtl_end, cluster_id)
    )

    dev.off()
    cat("[plot] SUCCESS ->", out_pdf, "\n")

    # ---- PSI annotation table using ggplot2 (no pdftools needed) ----
    psi_pdf <- file.path(OUT_DIR, sprintf("%s_%s_sashimi_annotated.pdf", rsid, stage))

    tryCatch({
      junc_short <- sub(":clu_.*", "", colnames(psi_summary))
      junc_short <- sub("^[^:]+:", "", junc_short)   # strip chr

      # Mark sQTL junction
      is_sqtl <- grepl(paste0(sqtl_start, ".*", sqtl_end, "|",
                               sqtl_end, ".*", sqtl_start), colnames(psi_summary))
      junc_short[is_sqtl] <- paste0(junc_short[is_sqtl], "\n(sQTL)")

      # Build long-form data frame for ggplot table
      psi_df <- as.data.frame(psi_summary)
      psi_df$genotype <- rownames(psi_df)
      psi_long <- reshape2::melt(psi_df, id.vars = "genotype",
                                 variable.name = "junction",
                                 value.name = "PSI")
      psi_long$junction_short <- junc_short[match(psi_long$junction,
                                                   colnames(psi_summary))]
      psi_long$label   <- sprintf("%.3f", psi_long$PSI)
      psi_long$is_sqtl <- is_sqtl[match(psi_long$junction, colnames(psi_summary))]
      psi_long$genotype <- factor(psi_long$genotype,
                                  levels = c("HomRef (0/0)", "Het (0/1)", "HomAlt (1/1)"))

      gt_cols <- c("HomRef (0/0)" = "#337fc2",
                   "Het (0/1)"    = "#c6a52d",
                   "HomAlt (1/1)" = "#bc005e")

      tbl_plot <- ggplot(psi_long, aes(x = junction_short, y = genotype,
                                        fill = genotype, label = label)) +
        geom_tile(colour = "white", linewidth = 0.5, alpha = 0.15) +
        geom_text(aes(colour = genotype,
                      fontface = ifelse(is_sqtl, "bold", "plain")),
                  size = 3.2) +
        scale_fill_manual(values = gt_cols, guide = "none") +
        scale_colour_manual(values = gt_cols, guide = "none") +
        theme_minimal(base_size = 10) +
        theme(
          axis.text.x  = element_text(size = 7, angle = 30, hjust = 1,
                                      face = ifelse(is_sqtl, "bold", "plain")),
          axis.text.y  = element_text(colour = gt_cols[levels(psi_long$genotype)],
                                      face = "bold"),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = 9, face = "bold")
        ) +
        labs(x = "Junction", y = NULL,
             title = "Mean PSI per genotype  (sQTL junction in bold)")

      # Delta PSI row
      if (all(c("HomRef (0/0)", "HomAlt (1/1)") %in% rownames(psi_summary))) {
        dpsi <- psi_summary["HomAlt (1/1)", ] - psi_summary["HomRef (0/0)", ]
        dpsi_df <- data.frame(
          junction_short = junc_short,
          label          = sprintf("%+.3f", dpsi),
          is_sqtl        = is_sqtl,
          y              = 0.5
        )
        dpsi_plot <- ggplot(dpsi_df, aes(x = junction_short, y = 1, label = label)) +
          geom_tile(aes(fill = is_sqtl), colour = "white", linewidth = 0.5, alpha = 0.1) +
          geom_text(aes(colour = is_sqtl,
                        fontface = ifelse(is_sqtl, "bold", "plain")),
                    size = 3.2) +
          scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "#e63946"),
                            guide = "none") +
          scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "#e63946"),
                              guide = "none") +
          theme_minimal(base_size = 10) +
          theme(axis.text   = element_blank(),
                panel.grid  = element_blank(),
                plot.title  = element_text(size = 9, face = "bold")) +
          labs(x = NULL, y = NULL, title = "ΔΨ (HomAlt − HomRef)")
      } else {
        dpsi_plot <- NULL
      }

      # Combine and save
      plots_list <- if (!is.null(dpsi_plot)) {
        list(tbl_plot, dpsi_plot)
      } else {
        list(tbl_plot)
      }
      heights <- if (!is.null(dpsi_plot)) c(3, 1) else c(3)

      pdf(psi_pdf, width = max(8, length(junc_short) * 1.5), height = 5)
      do.call(gridExtra::grid.arrange, c(plots_list, list(ncol = 1, heights = heights)))
      dev.off()
      cat("[plot] PSI table ->", psi_pdf, "\n")

    }, error = function(e2) {
      try(dev.off(), silent = TRUE)
      cat("[plot] PSI table skipped:", conditionMessage(e2), "\n")
    })

  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("[plot] ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n[done] Finished.\n")
