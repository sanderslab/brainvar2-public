#!/usr/bin/env Rscript
library(data.table)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: make_ld_block.R <block>")
}
block = args[1]

cat("\n==================================================\n")
cat("  STARTING BLOCK:", block, "\n")
cat("  Time:", timestamp(), "\n")
cat("==================================================\n\n")

ld_dir <- "/mnt/s3/eqtl/wgs_bv1_and_bv2/upload/LD.SELECT"
xref_file <- file.path(ld_dir, "rs-2-chr.pos.ref.alt.xref")
out_dir <- "data/LD_exact_blocks_v2"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, paste0(block, ".rds"))

### Skip if already exists
if (file.exists(out_file)) {
  cat("[", timestamp(), "]  Skipping", block, "→ RDS already exists.\n")
  quit(save = "no")
}

cat("[", timestamp(), "]  Loading mapping file...\n")
xref <- fread(xref_file, data.table = FALSE)
colnames(xref)[colnames(xref) == "id"] <- "rsid"
colnames(xref)[colnames(xref) == "alt.id"] <- "altid"
alt_to_rsid <- setNames(xref$rsid, xref$altid)

### Input files
vf <- file.path(ld_dir, paste0(block, ".vars.gz"))
cf <- file.path(ld_dir, paste0(block, ".correlation.gz"))

cat("[", timestamp(), "]  Checking input files...\n")
if (!file.exists(vf)) stop("Missing vars.gz file for block: ", block)
if (!file.exists(cf)) stop("Missing correlation.gz file for block: ", block)

cat("[", timestamp(), "]  Loading variant list...\n")
vars <- fread(vf, header = TRUE, data.table = FALSE)
colnames(vars)[1] <- "ID"
cat("   • Loaded", nrow(vars), "variants\n")

cat("[", timestamp(), "]  Loading correlation matrix (lower triangle)...\n")
r <- fread(cf, header = FALSE, data.table = FALSE)
cat("   • Loaded", nrow(r), "correlation entries\n")

### Build LD matrix
cat("[", timestamp(), "]  Building LD matrix...\n")
R <- matrix(0, nrow = nrow(vars), ncol = nrow(vars))
rownames(R) <- vars$ID
colnames(R) <- vars$ID

R[lower.tri(R, diag = TRUE)] <- r[,1]
ii <- which(upper.tri(R, diag = TRUE))
R[ii] <- t(R)[ii]

cat("   • LD matrix:", nrow(R), "x", ncol(R), "\n")

### Map to rsID
cat("[", timestamp(), "]  Mapping variant IDs to rsIDs...\n")
rs <- alt_to_rsid[rownames(R)]

mapped <- sum(!is.na(rs))
unmapped <- sum(is.na(rs))
cat("   • Mapped:", mapped, "\n")
cat("   • Unmapped:", unmapped, "\n")

keep <- !is.na(rs)
R2 <- R[keep, keep, drop = FALSE]
rs <- rs[keep]
rownames(R2) <- rs
colnames(R2) <- rs

cat("   • Final LD size:", nrow(R2), "x", ncol(R2), "\n")

### Save
cat("[", timestamp(), "]  Saving RDS file...\n")
saveRDS(R2, out_file)
cat("[", timestamp(), "]  Saved:", out_file, "\n")

cat("\n==================================================\n")
cat("  BLOCK COMPLETED:", block, "\n")
cat("  Time:", timestamp(), "\n")
cat("==================================================\n\n")
