#!/usr/bin/env Rscript
library(data.table)

## ============================ CONFIG ============================
ld_dir    <- "/mnt/s3/FXN/wgs_bv1_and_bv2/upload/LD.SELECT"
xref_file <- file.path(ld_dir, "rs-2-chr.pos.ref.alt.xref")
out_dir   <- "/mnt/s3/FXN/wgs_bv1_and_bv2/upload/LD_exact_blocks"

# Storage order of the flattened lower triangle in *.correlation.gz
#   "row"  -> r21, r31, r32, r41, ...  (most common)
#   "col"  -> column-major lower.tri() order
LOWER_ORDER <- "row"

# rsID(s) to check for after building
CHECK_RSIDS <- c("rs13197384")
## ================================================================

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\n==================================================\n")
cat("   Building LD matrices\n")
cat("==================================================\n\n")

### 1. Load mapping alt.id -> rsID
cat("Loading mapping file...\n")
xref <- fread(xref_file, data.table = FALSE)
colnames(xref)[colnames(xref) == "id"]     <- "rsid"
colnames(xref)[colnames(xref) == "alt.id"] <- "altid"
alt_to_rsid <- setNames(xref$rsid, xref$altid)

# Quick check: is the target rsID even in the mapping file?
for (target in CHECK_RSIDS) {
  in_xref <- target %in% xref$rsid
  cat("  ", target, "present in xref mapping:", in_xref, "\n")
}

### 2. Get all block files
vars_files <- list.files(ld_dir, pattern = "\\.vars\\.gz$", full.names = TRUE)
cat("Found", length(vars_files), "block vars files\n")

### 3. Process each block
for (vf in vars_files) {
  block <- sub("\\.vars\\.gz$", "", basename(vf))
  cf    <- file.path(ld_dir, paste0(block, ".correlation.gz"))

  cat("\n--------------------------------------------------\n")
  cat("Processing block:", block, "\n")

  if (!file.exists(cf)) {
    cat("  Missing correlation file, skipping\n")
    next
  }

  ## ---- load variants ----
  vars <- fread(vf, header = TRUE, data.table = FALSE)
  colnames(vars)[1] <- "ID"      # BrainVar format "chr:pos:ref:alt"
  n <- nrow(vars)

  ## ---- load flattened lower triangle ----
  r    <- fread(cf, header = FALSE, data.table = FALSE)
  vals <- r[[1]]                 # one value per line

  if (ncol(r) > 1) {
    cat("  WARNING: correlation file has", ncol(r),
        "columns; using column 1 only\n")
  }

  expected <- n * (n + 1) / 2
  if (length(vals) != expected) {
    cat("  ERROR: expected", expected, "values (n*(n+1)/2) but got",
        length(vals), "- skipping block\n")
    next
  }

  ## ---- build full symmetric LD matrix ----
  R <- matrix(0, nrow = n, ncol = n)
  rownames(R) <- vars$ID
  colnames(R) <- vars$ID

  if (LOWER_ORDER == "row") {
    R[upper.tri(R, diag = TRUE)] <- vals
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
  } else {
    R[lower.tri(R, diag = TRUE)] <- vals
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
  }

  cat("  Raw LD matrix size:", nrow(R), "x", ncol(R), "\n")

  ## ---- map to rsID ----
  rs <- alt_to_rsid[rownames(R)]
  cat("  Mapped", sum(!is.na(rs)), "variants to rsID\n")
  cat("  Unmapped:", sum(is.na(rs)), "\n")

  keep <- !is.na(rs)
  R2   <- R[keep, keep, drop = FALSE]
  rs   <- rs[keep]

  if (any(duplicated(rs))) {
    cat("  WARNING:", sum(duplicated(rs)),
        "duplicate rsID(s) after mapping - keeping first occurrence\n")
    dup <- !duplicated(rs)
    R2  <- R2[dup, dup, drop = FALSE]
    rs  <- rs[dup]
  }

  rownames(R2) <- rs
  colnames(R2) <- rs
  cat("  Final LD size:", nrow(R2), "x", ncol(R2), "\n")

  ## ---- report if any target rsID is in this block ----
  hit <- CHECK_RSIDS[CHECK_RSIDS %in% rs]
  if (length(hit) > 0) {
    cat("  *** FOUND target rsID(s) in this block:",
        paste(hit, collapse = ", "), "***\n")
  }

  ## ---- sanity check symmetry ----
  if (!isTRUE(all.equal(R2, t(R2), check.attributes = FALSE))) {
    cat("  WARNING: final matrix is not symmetric\n")
  }

  ## ---- save RDS ----
  out_file <- file.path(out_dir, paste0(block, ".rds"))
  saveRDS(R2, out_file)
  cat("  Saved:", out_file, "\n")
}

cat("\n==================================================\n")
cat("     ALL BLOCKS PROCESSED\n")
cat("==================================================\n")

### 4. Final scan: which saved block(s) contain each target rsID?
cat("\n==================================================\n")
cat("   Checking saved blocks for target rsID(s)\n")
cat("==================================================\n")

rds_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

for (target in CHECK_RSIDS) {
  found_in <- character(0)
  for (f in rds_files) {
    m <- readRDS(f)
    if (target %in% rownames(m)) {
      found_in <- c(found_in, sub("\\.rds$", "", basename(f)))
    }
  }
  if (length(found_in) > 0) {
    cat("\n", target, "FOUND in", length(found_in), "block(s):\n")
    cat("   ", paste(found_in, collapse = "\n    "), "\n")
  } else {
    cat("\n", target, "NOT found in any saved block\n")
  }
}
