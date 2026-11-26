# ==============================================================================
# DAR Pipeline: Differential Accessibility (Cell Type × Development Stage)
# Purpose: Find DARs comparing prenatal vs postnatal within each broad cell type
# Date: 2025-10-12
# R Version: >= 4.0.0
# Notes: Assumes combined_proj already has peaks (PeakMatrix) and metadata columns:
#        - broad.type   (ExN, InN, Glia, NNC, Progenitor, etc.)
#        - group        (prenatal, transition, postnatal)
#        - celltype_group or celltype_dev (broad.type_group)
# ==============================================================================

# ---- Library loading ---------------------------------------------------------
library(ArchR)   # Single-cell ATAC-seq analysis toolkit

# ---- Global config -----------------------------------------------------------
set.seed(123)
addArchRThreads(threads = 4)
options(future.globals.maxSize = 110595407872)
addArchRGenome("hg38")

# ---- Input: load project & metadata ------------------------------------------
# Replace path if necessary; if combined_proj already exists in R env, skip load
combined_proj <- loadArchRProject("/mnt/s3/Narjes/BrainVar2/ArchR_all_with_peaks/ArchR_all_with_peaks")

# If metadata not present in project, you can add like in your ATAC script:
# metadata <- read.delim('/path/to/meta.csv', sep = ',', header = TRUE, row.names = 1)
# combined_proj <- add_metadata_to_archr(combined_proj, metadata, "all")

# ---- Prepare groups: broad cell type × development ---------------------------
# Ensure 'broad.type' and 'group' exist:
if (!"broad.type" %in% colnames(getCellColData(combined_proj))) {
  stop("combined_proj lacks 'broad.type' column. Add it before running this pipeline.")
}
if (!"group" %in% colnames(getCellColData(combined_proj))) {
  stop("combined_proj lacks 'group' column. Add it before running this pipeline.")
}

# Remove transition; strict prenatal vs postnatal:
cells_keep <- getCellNames(combined_proj)[getCellColData(combined_proj)$group %in% c("prenatal", "postnatal")]
combined_proj <- combined_proj[cells_keep, ]

# Extract metadata from project
broad_type <- getCellColData(combined_proj, select = "broad.type")$broad.type
rna_celltype <- getCellColData(combined_proj, select = "rna_celltype_v1")$rna_celltype_v1
group <- getCellColData(combined_proj, select = "group")$group

# Initialize unified cell type
unified_celltype <- broad_type

# Replace with rna_celltype_v1 for OPC and Micro
unified_celltype[rna_celltype %in% c("OPC", "Micro")] <- rna_celltype[rna_celltype %in% c("OPC", "Micro")]

# Add unified cell type column to project
combined_proj <- addCellColData(
  ArchRProj = combined_proj,
  data = unified_celltype,
  cells = getCellNames(combined_proj),
  name = "unified_celltype",
  force = TRUE
)

# Filter to only cell types of interest
cell_types_to_analyze <- c("ExN", "InN", "Micro", "OPC")

cells_keep <- getCellNames(combined_proj)[
  getCellColData(combined_proj)$unified_celltype %in% cell_types_to_analyze &
    getCellColData(combined_proj)$group %in% c("prenatal", "postnatal")
]

combined_proj <- combined_proj[cells_keep, ]

# Create combined identity: celltype + developmental stage
celltype_dev <- paste(
  getCellColData(combined_proj)$unified_celltype,
  getCellColData(combined_proj)$group,
  sep = "_"
)

combined_proj <- addCellColData(
  ArchRProj = combined_proj,
  data = celltype_dev,
  cells = getCellNames(combined_proj),
  name = "celltype_dev",
  force = TRUE
)

# Print counts
message("Counts by celltype_dev:")
print(table(getCellColData(combined_proj)$celltype_dev))

# ---- Ensure PeakMatrix exists ------------------------------------------------
# If you previously created peaks & peak matrix, skip. Otherwise:
# combined_proj <- addGroupCoverages(ArchRProj = combined_proj, groupBy = "celltype_dev", force = TRUE)
# combined_proj <- addReproduciblePeakSet(ArchRProj = combined_proj, groupBy = "celltype_dev",
#                                        pathToMacs2 = findMacs2(), force = TRUE)
# combined_proj <- addPeakMatrix(combined_proj, force = TRUE)

if (!"PeakMatrix" %in% getAvailableMatrices(combined_proj)) {
  stop("PeakMatrix not found. Run addReproduciblePeakSet + addPeakMatrix before continuing.")
}

# ---- DARs per cell type: prenatal vs postnatal --------------------------------
# For each broad cell type, compare prenatal vs postnatal using getMarkerFeatures.

output_dir <- "./ATAC_DARs_celltype_dev_20251012/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

markers_global <- getMarkerFeatures(
  ArchRProj = combined_proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype_dev",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 10000
)

markerList <- getMarkers(markers_global, cutOff = "FDR <= 1 & Log2FC >= -Inf")

# Convert and write full table
markers_df_global <- data.frame()
for (stage in names(markerList)) {
  if (!is.null(markerList[[stage]]) && nrow(markerList[[stage]]) > 0) {
    temp_df <- as.data.frame(markerList[[stage]])
    temp_df$stage <- stage
    temp_df$peak <- rownames(temp_df)
    markers_df_global <- rbind(markers_df_global, temp_df)
  }
}

global_outfile <- file.path(output_dir, "DARs_celltype_dev_major.csv")
write.csv(markers_df_global, global_outfile, row.names = FALSE)
message("Global DARs saved to: ", global_outfile)
write.csv(markers_df_global, '/mnt/s3/Narjes/BrainVar2/DARs_celltype_dev_major.csv', row.names = FALSE)

# ---- Summary output ----------------------------------------------------------
message("DAR pipeline complete.")
message("Output directory: ", normalizePath(output_dir))
message("Global celltype_dev markers: ", global_outfile)

# Save R objects for future inspection
saveRDS(
  list(globalMarkers = markers_global, globalMarkerList = markerList),
  file = file.path(output_dir, "DARs_results_objects.rds")
)
message("Saved R objects to DARs_results_objects.rds")