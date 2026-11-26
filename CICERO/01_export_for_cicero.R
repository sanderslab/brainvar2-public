# ==============================================================================
# 01_export_for_cicero.R
# RUN THIS IN: macs2_env (has ArchR)
#
# Rebuilds the same analysis_group definitions as the ArchR runs, then for
# each group exports: binarized peak matrix, peak coords, cell barcodes,
# reduced-dim coords (ArchR UMAP if present, else IterativeLSI), plus one
# shared genomic_coords.tsv and peak_type_lookup.tsv.
# No monocle3/cicero needed here -- those only get installed in cicero_env.
#
# Checkpointed: a group is skipped if its DONE marker already exists.
#
# Run with: Rscript 01_export_for_cicero.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Matrix)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

addArchRThreads(threads = 3)

SUBTYPE_COL <- "rna_celltype_v1"
MIN_CELLS_PER_GROUP <- 500
MAX_CELLS_PER_GROUP <- 5000
export_dir <- "./cicero_export/"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
message("[DEBUG] Export dir: ", export_dir)

# ------------------------------------------------------------------------------
# STEP 1: Load project + rebuild the same analysis_group assignment
# ------------------------------------------------------------------------------
proj <- loadArchRProject("./ArchR_all_with_coaccessibility_complete")
ccd <- getCellColData(proj)

stopifnot("celltype_group" %in% colnames(ccd))
stopifnot(SUBTYPE_COL %in% colnames(ccd))

period_suffix <- ifelse(grepl("_prenatal$",  ccd$celltype_group), "prenatal",
                  ifelse(grepl("_postnatal$", ccd$celltype_group), "postnatal", NA_character_))
is_exn_inn  <- grepl("^(ExN|InN)_", ccd$celltype_group)
is_glia_nnc <- grepl("^(Glia|NNC)_", ccd$celltype_group)

analysis_group <- rep(NA_character_, nrow(ccd))
analysis_group[is_exn_inn] <- as.character(ccd$celltype_group[is_exn_inn])

subtype_raw   <- as.character(ccd[[SUBTYPE_COL]])
subtype_clean <- gsub("[^A-Za-z0-9]+", "_", trimws(subtype_raw))
glia_nnc_valid <- is_glia_nnc & !is.na(subtype_clean) & subtype_clean != "" & !is.na(period_suffix)
analysis_group[glia_nnc_valid] <- paste0(subtype_clean[glia_nnc_valid], "_", period_suffix[glia_nnc_valid])

proj$analysis_group <- analysis_group

group_counts <- sort(table(analysis_group, useNA = "no"), decreasing = TRUE)
groups_kept  <- names(group_counts)[group_counts >= MIN_CELLS_PER_GROUP]
message("[CHECK] Groups to export: ", length(groups_kept))
print(group_counts[groups_kept])

# ------------------------------------------------------------------------------
# STEP 2: Shared files -- genomic coords + peak type lookup (written once)
# ------------------------------------------------------------------------------
std_chroms <- paste0("chr", c(1:22, "X"))
seq_lens <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genomic_coords <- data.frame(chr = std_chroms, length = as.numeric(seq_lens[std_chroms]))
write.table(genomic_coords, file.path(export_dir, "genomic_coords.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("[DEBUG] Wrote genomic_coords.tsv (", nrow(genomic_coords), " chroms)")

peakSet_full <- getPeakSet(proj)
stopifnot("peakType" %in% colnames(mcols(peakSet_full)))
peak_ids_full <- paste0(as.character(seqnames(peakSet_full)), "_", start(peakSet_full), "_", end(peakSet_full))
peak_type_df <- data.frame(peak_id = peak_ids_full,
                            peakType = as.character(mcols(peakSet_full)$peakType))
write.table(peak_type_df, file.path(export_dir, "peak_type_lookup.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("[DEBUG] Wrote peak_type_lookup.tsv (", nrow(peak_type_df), " peaks)")

# ------------------------------------------------------------------------------
# STEP 3: Per-group export
# ------------------------------------------------------------------------------
for (grp in groups_kept) {

  grp_safe <- gsub("[^A-Za-z0-9]", "_", grp)
  grp_dir  <- file.path(export_dir, grp_safe)

  if (file.exists(file.path(grp_dir, "DONE"))) {
    message("[SKIP] ", grp, " already exported.")
    next
  }
  dir.create(grp_dir, showWarnings = FALSE, recursive = TRUE)
  message("\n---- Exporting group: ", grp, " ----")

  cellsOfGroup <- getCellNames(proj)[!is.na(proj$analysis_group) & proj$analysis_group == grp]

  if (length(cellsOfGroup) > MAX_CELLS_PER_GROUP) {
    set.seed(123)
    cells_use <- sample(cellsOfGroup, MAX_CELLS_PER_GROUP)
    message("  [DEBUG] Subsampled ", grp, " from ", length(cellsOfGroup), " to ", MAX_CELLS_PER_GROUP)
  } else {
    cells_use <- cellsOfGroup
  }

  proj_temp <- proj[cells_use, ]

  pm <- tryCatch(
    getMatrixFromProject(ArchRProj = proj_temp, useMatrix = "PeakMatrix", binarize = TRUE),
    error = function(e) { message("  [ERROR] getMatrixFromProject: ", conditionMessage(e)); NULL }
  )
  if (is.null(pm)) next

  mat <- assay(pm)
  peaks_gr <- rowRanges(pm)
  peak_ids <- paste0(as.character(seqnames(peaks_gr)), "_", start(peaks_gr), "_", end(peaks_gr))
  rownames(mat) <- peak_ids
  message("  [DEBUG] Peak matrix: ", nrow(mat), " peaks x ", ncol(mat), " cells")

  Matrix::writeMM(as(mat, "CsparseMatrix"), file.path(grp_dir, "peak_matrix.mtx"))
  writeLines(peak_ids, file.path(grp_dir, "peak_ids.txt"))
  writeLines(colnames(mat), file.path(grp_dir, "cell_barcodes.txt"))

  peak_df <- data.frame(chr = as.character(seqnames(peaks_gr)),
                         bp1 = start(peaks_gr), bp2 = end(peaks_gr))
  write.table(peak_df, file.path(grp_dir, "peak_coords.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Reduced dims: reuse ArchR's own UMAP if present (no monocle3 needed),
  # else export the LSI coords so cicero_env can compute UMAP itself.
  has_archr_umap <- tryCatch("UMAP" %in% names(proj_temp@embeddings), error = function(e) FALSE)
  if (has_archr_umap) {
    emb <- getEmbedding(proj_temp, embedding = "UMAP", returnDF = TRUE)
    emb <- emb[colnames(mat), , drop = FALSE]
    colnames(emb) <- c("UMAP1", "UMAP2")
    write.table(emb, file.path(grp_dir, "umap_coords.tsv"), sep = "\t", quote = FALSE)
    writeLines("umap", file.path(grp_dir, "reduced_dim_type.txt"))
    message("  [DEBUG] Exported existing ArchR UMAP")
  } else {
    lsi <- getReducedDims(proj_temp, reducedDims = "IterativeLSI")
    lsi <- lsi[colnames(mat), , drop = FALSE]
    write.table(lsi, file.path(grp_dir, "lsi_coords.tsv"), sep = "\t", quote = FALSE)
    writeLines("lsi", file.path(grp_dir, "reduced_dim_type.txt"))
    message("  [DEBUG] No ArchR UMAP -- exported LSI for cicero_env to reduce itself")
  }

  writeLines("done", file.path(grp_dir, "DONE"))
  message("  [DEBUG] Done -> ", grp_dir)
}

message("\n[DONE] Export complete -> ", export_dir)
message("Now copy/rsync '", export_dir, "' to wherever cicero_env will run, then run 02_run_cicero_from_export.R there.")
