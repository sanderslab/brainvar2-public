# ============================================================
# make_homer_input_caQTL.R
#
# Converts significant caQTL tables (one file per cell type x stage)
# into HOMER-format peak files for findMotifsGenome.pl.
#
# HOMER peak file format (tab-separated, no header):
#   PeakID  chr  start  end  strand
# ATAC peaks have no strand, so strand is set to "+" for every peak
# (findMotifsGenome.pl scans both strands regardless).
#
# Usage:
#   Rscript make_homer_input_caQTL.R
#   (or source in RStudio after editing the CONFIG paths below)
#
# Then run HOMER per file, e.g.:
#   findMotifsGenome.pl ExN_Prenatal_homer.txt hg38 ExN_Prenatal_motifs/ -size given
# ============================================================

library(data.table)

# ---- CONFIG: edit these paths ----
caqtl_dir <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/caQTL/sig_caQTL"
out_dir   <- file.path(caqtl_dir, "homer_input_caqtl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

CELL_TYPES <- c("ExN", "InN", "Glia", "NNC")
STAGES     <- c("Prenatal", "Postnatal")

# same two possible filename patterns as the motif enrichment script
FILE_TEMPLATE_ANNOTATED <- file.path(caqtl_dir, "{stage}_{celltype}_all_chr_1e-5_peakanno.tsv")
FILE_TEMPLATE_PLAIN     <- file.path(caqtl_dir, "{stage}_{celltype}_all_chr_1e-5_annotated.tsv")

resolve_path <- function(celltype, stage) {
  fill <- function(tpl) gsub("\\{celltype\\}", celltype, gsub("\\{stage\\}", stage, tpl))
  p1 <- fill(FILE_TEMPLATE_ANNOTATED)
  if (file.exists(p1)) return(p1)
  p2 <- fill(FILE_TEMPLATE_PLAIN)
  if (file.exists(p2)) return(p2)
  NA_character_
}

# ---- conversion: one caQTL file -> one HOMER peak file ----
write_homer_input <- function(path, out_file) {
  df <- fread(path, sep = "\t")
  stopifnot(all(c("chr", "start", "end", "peak_coordinate") %in% names(df)))
  
  # collapse multi-variant rows down to unique peaks
  peaks <- unique(df[, .(chr, start, end, peak_coordinate)])
  
  homer_df <- data.table(
    PeakID = peaks$peak_coordinate,
    chr    = peaks$chr,
    start  = peaks$start,
    end    = peaks$end,
    strand = "+"
  )
  
  fwrite(homer_df, out_file, sep = "\t", col.names = FALSE)
  nrow(homer_df)
}

# ---- run over all cell type x stage combinations ----
for (ct in CELL_TYPES) {
  for (st in STAGES) {
    path <- resolve_path(ct, st)
    if (is.na(path)) {
      warning("No caQTL file found for ", st, " / ", ct)
      next
    }
    
    out_file <- file.path(out_dir, paste0(ct, "_", st, "_homer.txt"))
    n <- write_homer_input(path, out_file)
    cat(sprintf("%s / %s: wrote %d unique peaks -> %s\n", st, ct, n, out_file))
  }
}

cat("\nHOMER input files written to: ", out_dir, "\n", sep = "")