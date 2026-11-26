# ==============================================================================
# Differential Accessibility Analysis Pipeline for BrainVar ATAC Data
# ==============================================================================
# Purpose: Perform differential accessibility analysis on single-nuclei ATAC-seq
#          data using ArchR. Generates marker peaks for broad, major, and fine
#          cell type classifications separately in prenatal and postnatal
#          datasets.
# Author: Narjes Rohani
# Date: Aug 17, 2025
# R Version: >= 4.0.0
# ==============================================================================

# ------------------------------------------------------------------------------
# Package Installation and Loading
# ------------------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  } else {
    message("Package already installed: ", pkg)
  }
}

# Install BiocManager if needed (required for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}

# Install devtools if needed (required to install ArchR from GitHub)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", dependencies = TRUE)
}

# Bioconductor packages
bioc_packages <- c(
  "BSgenome.Hsapiens.UCSC.hg38"   # genome reference required by addArchRGenome("hg38")
)

# CRAN packages
cran_packages <- c(
  "future"                        # parallel backend for ArchR multithreading
)

# Install ArchR from GitHub (not available on CRAN/Bioconductor)
if (!requireNamespace("ArchR", quietly = TRUE)) {
  devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())
  ArchR::installExtraPackages()
}

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load libraries
print("Loading libraries...")

suppressPackageStartupMessages({
  library(ArchR)                       # core single-cell ATAC analysis
  library(future)                      # parallel processing backend
  library(BSgenome.Hsapiens.UCSC.hg38) # hg38 genome reference
})
# ------------------------------------------------------------------------------
# Global Settings Configuration
# ------------------------------------------------------------------------------

print("Setting seed, future plan, and memory options...")

set.seed(123)

plan("multisession", workers = 4)

addArchRThreads(threads = 1)

# options(future.globals.maxSize = 110595407872)

addArchRGenome("hg38")

# ------------------------------------------------------------------------------
# Main Processing Function
# ------------------------------------------------------------------------------

# Function to process ArchR project and perform differential accessibility analysis
# @param archr_proj: ArchR project containing single-cell ATAC data
# @param group_label: String identifier for the dataset (e.g., "prenatal", "postnatal")
process_archr <- function(archr_proj, group_label) {
  print(paste("Processing", group_label, "dataset..."))
  # ------------------------------------------------------------------------------
  # Differential Accessibility Analysis - Fine Cell Types
  # ------------------------------------------------------------------------------
  
  print("  Running DA analysis by fine cell types (celltype_v1)...")
  # Process only the most abundant cell types
  celltype_counts <- table(archr_proj$rna_celltype_v1)
  # Filter out cell types with too few cells (minimum 50 cells)
  min_cells <- 50
  valid_celltypes <- names(celltype_counts)[celltype_counts >= min_cells]
  small_celltypes <- names(celltype_counts)[celltype_counts < min_cells]
  
  top_celltypes <- valid_celltypes
  print(paste("Using only top", length(top_celltypes), "cell types:", paste(top_celltypes, collapse = ", ")))
  
  cells_to_keep_alt <- which(archr_proj$rna_celltype_v1 %in% top_celltypes)
  archr_subset_alt <- archr_proj[cells_to_keep_alt, ]
  
  # Get marker peaks for fine cell types
  markersPeaks_fine <- getMarkerFeatures(
    ArchRProj = archr_subset_alt,
    useMatrix = "PeakMatrix",
    groupBy = "rna_celltype_v1",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 10000
  )
  
  # Extract marker peaks as a data frame
  markerList_fine <- getMarkers(markersPeaks_fine, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
  
  # Convert to a single data frame for all cell types
  markers_fine_df <- data.frame()
  for (cell_type in names(markerList_fine)) {
    if (nrow(markerList_fine[[cell_type]]) > 0) {
      temp_df <- markerList_fine[[cell_type]]
      temp_df <- as.data.frame(temp_df)   # convert from S4Vectors::DataFrame
      
      temp_df$cluster <- cell_type
      temp_df$gene <- rownames(temp_df)
      markers_fine_df <- rbind(markers_fine_df, temp_df)
    }
  }
  
  # Save results to CSV file
  output_file_fine <- paste0('./ATACmarkers1792025/', group_label, "_markers_minor_celltypes_all.csv")
  write.csv(markers_fine_df, output_file_fine, row.names = FALSE)
  print(paste("  Saved:", output_file_fine))

  # ------------------------------------------------------------------------------
  # Differential Accessibility Analysis - Broad Cell Types
  # ------------------------------------------------------------------------------
  
  print("  Running DA analysis by broad cell types...")
  
  # Get marker peaks for broad cell types
  markersPeaks_broad <- getMarkerFeatures(
    ArchRProj = archr_proj,
    useMatrix = "PeakMatrix",
    groupBy = "cell_broad",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 10000  # Subsample for computational efficiency
  )
  
  # Extract marker peaks as a data frame
  markerList_broad <-getMarkers(markersPeaks_broad, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")

  
  markers_broad_df <- data.frame()
  
  for (cell_type in names(markerList_broad)) {
    temp_df <- markerList_broad[[cell_type]]
    
    if (!is.null(temp_df) && nrow(temp_df) > 0) {
      temp_df <- as.data.frame(temp_df)   # convert from S4Vectors::DataFrame
      temp_df$cluster <- cell_type
      temp_df$gene <- rownames(temp_df)   # will be NULL unless rownames exist
      markers_broad_df <- rbind(markers_broad_df, temp_df)
    }
  }
  
  # Save results to CSV file
  output_file_broad <- paste0('./ATACmarkers1792025/', group_label, "_markers_broad_celltypes_all.csv")
  dir.create('./ATACmarkers1792025/', showWarnings = FALSE, recursive = TRUE)
  write.csv(markers_broad_df, output_file_broad, row.names = FALSE)
  print(paste("  Saved:", output_file_broad))
  
  # ------------------------------------------------------------------------------
  # Differential Accessibility Analysis - Major Cell Classes
  # ------------------------------------------------------------------------------
  
  print("  Running DA analysis by majorclass_v1...")
  print("  Checking if enough cells for broad cell type analysis...")

  # Get marker peaks for major cell classes
  markersPeaks_major <- getMarkerFeatures(
    ArchRProj = archr_proj,
    useMatrix = "PeakMatrix",
    groupBy = "rna_majorclass_v1",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    maxCells = 10000
  )
  
  # Extract marker peaks as a data frame
  markerList_major <- getMarkers(markersPeaks_major, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
  
  # Convert to a single data frame for all cell types
  markers_major_df <- data.frame()
  for (cell_type in names(markerList_major)) {
    if (nrow(markerList_major[[cell_type]]) > 0) {
      
      temp_df <- markerList_major[[cell_type]]
      temp_df <- as.data.frame(temp_df)   # convert from S4Vectors::DataFrame
      
      temp_df$cluster <- cell_type
      temp_df$gene <- rownames(temp_df)
      markers_major_df <- rbind(markers_major_df, temp_df)
    }
  }
  
  # Save results to CSV file
  output_file_major <- paste0('./ATACmarkers1792025/', group_label, "_markers_majorclass_v1_celltypes_all.csv")
  write.csv(markers_major_df, output_file_major, row.names = FALSE)
  print(paste("  Saved:", output_file_major))
  
  print(paste("  Completed processing of", group_label, "dataset"))
  
  # Return the processed ArchR project for further analysis if needed
  return(archr_proj)
  

}

# ------------------------------------------------------------------------------
# Main Analysis Execution
# ------------------------------------------------------------------------------

print("Loading prenatal and postnatal ArchR projects...")

# Load preprocessed ArchR projects
# These should be ArchR projects with cell type annotations already created
postnatal_proj <- loadArchRProject("./Save-proj_postnatal")
prenatal_proj <- loadArchRProject("./Save-proj_prenatal")

# Load meta data
metadata <- read.delim('/mnt/s3/Sanders_lab/BrainVar2/ArchR_brainvar/ATAC_P3-P15_meta_D240727.tsv', 
                       sep = ',', header = TRUE, row.names = 1)
# Find common cells
cellNames <- getCellNames(prenatal_proj)
common_cells <- intersect(rownames(metadata), cellNames)

print(paste("Original ArchR project cells:", length(cellNames)))
print(paste("Cells with metadata:", length(common_cells)))
print(paste("Cells to remove:", length(cellNames) - length(common_cells)))

# Now add metadata to the subsetted project
data = metadata[common_cells, ]

metadata_cols <- c(
  "rna_sdf_score", "rna_sdf_class",
  "rna_celltype_v1", "rna_majorclass_v1", "atac_celltype_v1",
  "Age", "AgeUnits", "AgeDays", "Period", "Epoch",
  "Sex", "SexForAnalysis"
)

for (colname in metadata_cols) {
  prenatal_proj <- addCellColData(
    ArchRProj = prenatal_proj,
    data = data[[colname]],
    cells = common_cells,
    name = colname,
    force = FALSE
  )
}

# Your grouping list
list.broad.type <- list(
  ExN        = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
  InN        = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
  Glia       = c('Glia'), 
  NNC        = c('NNC'), 
  Progenitor = c('Progenitor')
)

# Flatten the list to a lookup named vector:
lookup <- unlist(
  lapply(names(list.broad.type), function(broad) {
    setNames(rep(broad, length(list.broad.type[[broad]])), list.broad.type[[broad]])
  })
)
rna_classes <- ArchR::getCellColData(prenatal_proj, select = "rna_majorclass_v1")
# Create the new vector by mapping
cell_broad <- lookup[ as.character(rna_classes$rna_majorclass_v1) ]

# For cells with no match, assign "Other"
cell_broad[is.na(cell_broad)] <- "Other"

prenatal_proj <- ArchR::addCellColData(
  ArchRProj = prenatal_proj,
  data = cell_broad,
  cells = rownames(rna_classes),
  name = "cell_broad",
  force = FALSE
)


# Find common cells
cellNames <- getCellNames(postnatal_proj)
common_cells <- intersect(rownames(metadata), cellNames)

print(paste("Original ArchR project cells:", length(cellNames)))
print(paste("Cells with metadata:", length(common_cells)))
print(paste("Cells to remove:", length(cellNames) - length(common_cells)))

# Now add metadata to the subsetted project
data = metadata[common_cells, ]

metadata_cols <- c(
  "rna_sdf_score", "rna_sdf_class",
  "rna_celltype_v1", "rna_majorclass_v1", "atac_celltype_v1",
  "Age", "AgeUnits", "AgeDays", "Period", "Epoch",
  "Sex", "SexForAnalysis"
)

for (colname in metadata_cols) {
  postnatal_proj <- addCellColData(
    ArchRProj = postnatal_proj,
    data = data[[colname]],
    cells = common_cells,
    name = colname,
    force = FALSE
  )
}

# Your grouping list
list.broad.type <- list(
  ExN        = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
  InN        = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
  Glia       = c('Glia'), 
  NNC        = c('NNC'), 
  Progenitor = c('Progenitor')
)

# Flatten the list to a lookup named vector:
lookup <- unlist(
  lapply(names(list.broad.type), function(broad) {
    setNames(rep(broad, length(list.broad.type[[broad]])), list.broad.type[[broad]])
  })
)
rna_classes <- ArchR::getCellColData(postnatal_proj, select = "rna_majorclass_v1")
# Create the new vector by mapping
cell_broad <- lookup[ as.character(rna_classes$rna_majorclass_v1) ]

# For cells with no match, assign "Other" 
cell_broad[is.na(cell_broad)] <- "Other"

postnatal_proj <- ArchR::addCellColData(
  ArchRProj = postnatal_proj,
  data = cell_broad,
  cells = rownames(rna_classes),
  name = "cell_broad",
  force = TRUE
)



print(paste("Loaded prenatal dataset with", nCells(prenatal_proj), "cells"))
print(paste("Loaded postnatal dataset with", nCells(postnatal_proj), "cells"))

# ------------------------------------------------------------------------------
# Execute Differential Accessibility Analysis
# ------------------------------------------------------------------------------

# Process prenatal dataset
print("STARTING PRENATAL ANALYSIS")
prenatal_processed <- process_archr(prenatal_proj, "prenatal")

# Process postnatal dataset  
print("STARTING POSTNATAL ANALYSIS")
postnatal_processed <- process_archr(postnatal_proj, "postnatal")

print("ANALYSIS COMPLETE")
print("All differential accessibility results have been saved to CSV files.")
