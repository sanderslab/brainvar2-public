#!/usr/bin/env Rscript

library(data.table)
library(Matrix)
library(stringr)

###===========================================================
### 0. DIRECTORIES & SETUP
###===========================================================

ld_input <- "/mnt/s3/eqtl/wgs_bv1_and_bv2/upload/LD.SELECT"
xref_file <- file.path(ld_input, "rs-2-chr.pos.ref.alt.xref")
ld_output <- "data/LD"

cat("\n=========================================================\n")
cat("   BrainVar LD Construction Script - DEBUG MODE\n")
cat("=========================================================\n\n")

if (!dir.exists(ld_output)) dir.create(ld_output, recursive = TRUE)

###===========================================================
### 1. LOAD ALT.ID → RSID MAPPING
###===========================================================

cat("Step 1: Loading SNP mapping (alt.id → rsID)...\n")
stopifnot(file.exists(xref_file))

xref <- fread(xref_file, data.table = FALSE)
colnames(xref)[colnames(xref) == "id"] <- "rsid"
colnames(xref)[colnames(xref) == "alt.id"] <- "altid"

cat(" - Number of mapping rows:", nrow(xref), "\n")
cat(" - Column names:", paste(colnames(xref), collapse=", "), "\n")
cat(" - First few rows:\n")
print(head(xref))

alt_to_rsid <- setNames(xref$rsid, xref$altid)

###===========================================================
### 2. PARSE BLOCK FILENAMES
###===========================================================

parse_block <- function(filename) {
  base <- basename(filename)
  core <- sub("\\.vars\\.gz$|\\.correlation\\.gz$", "", base)
  m <- str_match(core, "^chr([0-9]+):([0-9]+)-([0-9]+)$")
  if (any(is.na(m))) return(NULL)
  list(chr = as.integer(m[2]), start = as.integer(m[3]), end = as.integer(m[4]), key = core)
}

###===========================================================
### 3. REBUILD LD BLOCK
###===========================================================

build_block_matrix <- function(vars_file, cor_file) {
  cat("\n  Loading block vars:", vars_file, "\n")
  vars <- fread(vars_file, data.table = FALSE)
  colnames(vars)[1] <- "altid"

  cat("   - Variants in block:", nrow(vars), "\n")
  print(head(vars))

  cat("  Loading block correlations:", cor_file, "\n")
  cor_vec <- fread(cor_file, header = FALSE, data.table = FALSE)[[1]]

  n <- nrow(vars)
  expected <- n*(n+1)/2

  cat("   - Expecting", expected, "values; received:", length(cor_vec), "\n")
  if (length(cor_vec) != expected) {
    stop("Correlation vector length mismatch for block ", vars_file)
  }

  cat("  Reconstructing LD matrix (", n, "x", n, ") ...\n")
  R <- matrix(0, n, n)
  R[lower.tri(R, diag = TRUE)] <- cor_vec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]

  cat("  - First 5x5 of LD block:\n")
  print(R[1:5,1:5])

  rownames(R) <- vars$altid
  colnames(R) <- vars$altid

  return(R)
}

###===========================================================
### 4. SCAN LD.SELECT FOR BLOCK FILES
###===========================================================

cat("\n\nStep 2: Scanning LD.SELECT directory for block files...\n")

vars_files <- list.files(ld_input, pattern="\\.vars\\.gz$", full.names = TRUE)
cor_files  <- list.files(ld_input, pattern="\\.correlation\\.gz$", full.names = TRUE)

cat(" - Found", length(vars_files), "vars files and", length(cor_files), "correlation files\n")

blocks <- list()

# Register vars files
for (vf in vars_files) {
  parsed <- parse_block(vf)
  if (!is.null(parsed)) {
    blocks[[parsed$key]] <- list(vars=vf, chr=parsed$chr, start=parsed$start, end=parsed$end)
  }
}

# Register correlation files
for (cf in cor_files) {
  parsed <- parse_block(cf)
  if (!is.null(parsed) && parsed$key %in% names(blocks)) {
    blocks[[parsed$key]]$cor <- cf
  }
}

cat(" - Total usable LD blocks:", sum(sapply(blocks, function(x) !is.null(x$cor))), "\n")

###===========================================================
### 5. BUILD LD PER CHROMOSOME
###===========================================================

for (chr in 1:22) {

  cat("\n===========================================================\n")
  cat("Building LD for chromosome", chr, "\n")
  cat("===========================================================\n")

  blk_list <- Filter(function(b) !is.null(b$cor) && b$chr == chr, blocks)

  if (length(blk_list) == 0) {
    cat(" - No blocks for chromosome", chr, "\n")
    next
  }

  cat(" - Number of blocks:", length(blk_list), "\n")

  # Sort blocks by coordinate
  blk_list <- blk_list[order(sapply(blk_list, `[[`, "start"))]

  block_mats <- list()

  for (i in seq_along(blk_list)) {
    b <- blk_list[[i]]
    cat("\n>>> BLOCK", i, ":", b$start, "-", b$end, "\n")
    block_mats[[i]] <- build_block_matrix(b$vars, b$cor)
  }

  cat("\nMerging blocks into one LD matrix...\n")
  LD_chr <- as.matrix(bdiag(block_mats))

  cat(" - Combined LD dimension:", dim(LD_chr)[1], "x", dim(LD_chr)[2], "\n")

  ###-----------------------------------------------------------
  ###  CONVERT alt.id → rsID
  ###-----------------------------------------------------------

  old_ids <- rownames(LD_chr)
  rs_ids <- alt_to_rsid[old_ids]

  cat(" - Variants before rsID mapping:", length(old_ids), "\n")
  cat(" - Variants with valid rsID:", sum(!is.na(rs_ids)), "\n")
  cat(" - Variants dropped (no rsID):", sum(is.na(rs_ids)), "\n")

  keep <- which(!is.na(rs_ids))

  LD_chr <- LD_chr[keep, keep]
  rs_ids <- rs_ids[keep]

  rownames(LD_chr) <- rs_ids
  colnames(LD_chr) <- rs_ids

  cat(" - Final LD dimension (after filtering):",
      dim(LD_chr)[1], "x", dim(LD_chr)[2], "\n")

  cat(" - First few rsIDs:\n")
  print(head(rs_ids))

  ###-----------------------------------------------------------
  ### SAVE RESULT
  ###-----------------------------------------------------------

  out_file <- file.path(ld_output, paste0("LD_with_blocks_chr", chr, ".rds"))
  cat("\nSaving LD matrix to:", out_file, "\n")
  saveRDS(LD_chr, out_file)

  cat("Chromosome", chr, "completed.\n\n")
}

cat("\n===========================================================\n")
cat("ALL CHROMOSOMES PROCESSED SUCCESSFULLY\n")
cat("LD MATRICES SAVED IN:", ld_output, "\n")
cat("===========================================================\n")
