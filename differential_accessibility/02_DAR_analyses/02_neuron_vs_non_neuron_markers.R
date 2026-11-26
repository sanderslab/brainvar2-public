# ==============================================================================
# ATAC Neuron vs Non-Neuron Differential Accessibility Analysis
# ==============================================================================
# Purpose: Compare chromatin accessibility between neurons and non-neuronal
#          cell types in prenatal and postnatal brain development datasets
# Author: Narjes Rohani
# Date: Aug 17, 2025
# R Version: >= 4.0.0
# Analysis: Neuron vs Non-Neuron differential accessibility using ArchR
# ==============================================================================

# Library Loading ==============================================================
print("Loading required libraries...")

library(ArchR)        # Single-cell ATAC-seq analysis toolkit
library(dplyr)        # Data manipulation and transformation
library(data.table)   # Fast data manipulation

# Global Configuration =========================================================
print("Configuring global analysis settings...")

# Set random seed for reproducible results across all analyses
set.seed(123)

# Set ArchR threading
addArchRThreads(threads = 3)

# Set maximum object size for parallel operations (~103 GB)
options(future.globals.maxSize = 110595407872)

# Configure ArchR genome annotation
addArchRGenome("hg38")

# Main Analysis Function =======================================================

#' Process ArchR object and perform differential accessibility analysis
#'
#' This function performs comprehensive cell type classification and differential
#' accessibility analysis comparing neurons vs non-neuronal cells
#'
#' @param archr_obj ArchR object containing single-cell ATAC-seq data
#' @param group_label Character string identifying the dataset (e.g., "prenatal")
#' @return ArchR object with cell type classification added
process_archr_neuron_analysis <- function(archr_obj, group_label) {
  
  print(paste("Processing", group_label, "ATAC dataset for Neuron vs Non-Neuron analysis..."))
  print(paste("Initial cell count:", nCells(archr_obj)))
  print(paste("Initial peak count:", nrow(getPeakSet(archr_obj))))
  
  # Cell Type Classification =================================================
  print("Step 1: Assigning neuron vs non-neuron categories...")
  
  # Define specific non-neuronal cell types (matching RNA analysis)
  # These are well-characterized, functionally distinct cell populations
  non_neuron_types <- c(
    "Astrocyte",   # Glial cells providing metabolic and structural support
    "Micro",       # Microglia (brain-resident immune cells)
    "Endo",        # Endothelial cells (blood-brain barrier)
    "OPC",         # Oligodendrocyte precursor cells (myelination)
    "VLMC",        # Vascular leptomeningeal cells
    "SMC"          # Smooth muscle cells
  )
  
  # Cell types to exclude from the comparison (ambiguous or low-confidence)
  exclude_types <- c("Unknown", "Oligo", "Immune")
  
  # Get cell type annotations
  cell_types <- getCellColData(archr_obj, select = "rna_celltype_v1")$rna_celltype_v1
  
  # Create neuron vs non-neuron classification
  cell_type_group <- ifelse(
    cell_types %in% exclude_types,
    NA,
    ifelse(cell_types %in% non_neuron_types,
           as.character(cell_types),  # Keep the specific non-neuronal type name
           "Neuron"                   # All other non-excluded types become "Neuron"
    )
  )
  
  print("Cell type groups identified:")
  print(table(cell_type_group, useNA = "always"))
  
  # Add classification to ArchR project
  archr_obj <- addCellColData(
    ArchRProj = archr_obj,
    data = cell_type_group,
    cells = getCellNames(archr_obj),
    name = "cell_type_group",
    force = FALSE
  )
  
  # Quality Control Filtering ================================================
  print("Step 2: Filtering cells for analysis...")
  
  initial_count <- nCells(archr_obj)
  
  # Remove cells with unknown identity for clean comparison
  valid_cells <- !is.na(cell_type_group)
  cells_to_keep <- getCellNames(archr_obj)[valid_cells]
  archr_obj <- archr_obj[cells_to_keep, ]
  
  final_count <- nCells(archr_obj)
  cells_removed <- initial_count - final_count
  
  print("Quality control applied:")
  print(paste("  - Cells removed (unknown types):", cells_removed))
  print(paste("  - Cells retained:", final_count))
  print(paste("  - Retention rate:", round(final_count / initial_count * 100, 1), "%"))
  
  # Display final cell type distribution
  cell_type_counts <- table(getCellColData(archr_obj, select = "cell_type_group")$cell_type_group)
  print("Final cell type distribution:")
  print(cell_type_counts)
  
  # Validate sufficient cells for comparison
  unique_groups <- unique(getCellColData(archr_obj, select = "cell_type_group")$cell_type_group)
  if (length(unique_groups) < 2) {
    stop(paste("Insufficient cell type diversity in", group_label,
               "dataset. Need at least 2 cell types for comparison."))
  }
  
  # Check minimum cells per group
  min_cells_threshold <- 50
  small_groups <- names(cell_type_counts)[cell_type_counts < min_cells_threshold]
  if (length(small_groups) > 0) {
    print(paste("Warning: Some cell types have fewer than", min_cells_threshold, "cells:"))
    print(small_groups)
    print("Consider grouping or filtering these cell types for robust analysis.")
  }
  
  # Differential Accessibility Analysis ======================================
  print("Step 3: Performing differential accessibility analysis...")
  print("  Analysis type: Neuron vs Non-Neuron (all pairwise comparisons)")
  
  # getMarkerFeatures compares each cell type against all others (like FindAllMarkers)
  # - Uses Wilcoxon rank-sum test for statistical testing
  # - Corrects for technical biases (TSS enrichment, fragment count)
  # - Subsamples cells for computational efficiency while maintaining power
  min_cells <- 30
  celltype_counts <- table(archr_obj$cell_type_group)
  
  valid_celltypes <- names(celltype_counts)[celltype_counts >= min_cells]
  top_celltypes <- valid_celltypes
  print(paste("Using only top", length(top_celltypes), "cell types:",
              paste(top_celltypes, collapse = ", ")))
  
  cells_to_keep_alt <- which(archr_obj$cell_type_group %in% top_celltypes)
  archr_subset_alt <- archr_obj[cells_to_keep_alt, ]
  
  markersPeaks_neuron <- getMarkerFeatures(
    ArchRProj = archr_subset_alt,
    useMatrix = "PeakMatrix",
    groupBy = "cell_type_group",
    bias = c("TSSEnrichment", "log10(nFrags)"),  # Correct for technical biases
    testMethod = "wilcoxon",                     # Non-parametric statistical test
    maxCells = 10000,                            # Subsample for efficiency
    verbose = TRUE                               # Show progress during computation
  )
  
  # Extract marker peaks as data frame
  print("  Extracting and formatting results...")
  markerList_neuron <- getMarkers(markersPeaks_neuron, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
  
  # Convert to a comprehensive data frame for all cell types
  markers_neuron_df <- data.frame()
  for (cell_type in names(markerList_neuron)) {
    if (!is.null(markerList_neuron[[cell_type]]) && nrow(markerList_neuron[[cell_type]]) > 0) {
      temp_df <- as.data.frame(markerList_neuron[[cell_type]])  # Convert from S4Vectors::DataFrame
      
      # Add metadata columns for consistency with RNA analysis
      temp_df$cluster <- cell_type
      temp_df$gene <- rownames(temp_df)
      
      markers_neuron_df <- rbind(markers_neuron_df, temp_df)
    }
  }
  
  # Generate output filename (matching RNA format)
  output_file <- paste0('./ATACmarkers/', group_label, "_markers_neuronvsnn_celltypes_all_v3.csv")
  
  # Ensure output directory exists
  dir.create('./ATACmarkers/', showWarnings = FALSE, recursive = TRUE)
  
  # Save comprehensive results
  write.csv(markers_neuron_df, output_file, row.names = FALSE)
  
  return(archr_obj)
}

# Data Loading and Analysis Execution =========================================
print("Loading preprocessed ATAC-seq datasets...")

print("Loading postnatal ArchR project...")
postnatal_proj <- loadArchRProject("./Save-proj_postnatal")

print("Loading prenatal ArchR project...")
prenatal_proj <- loadArchRProject("./Save-proj_prenatal")

# Load and add metadata
print("Loading and adding metadata...")

metadata <- read.delim('/mnt/s3/Sanders_lab/BrainVar2/ArchR_brainvar/ATAC_P3-P15_meta_D240727.tsv',
                       sep = ',', header = TRUE, row.names = 1)

# Function to add metadata to ArchR project
add_metadata_to_archr <- function(archr_proj, metadata, project_name) {
  print(paste("Adding metadata to", project_name, "project..."))
  
  # Find common cells
  cellNames <- getCellNames(archr_proj)
  common_cells <- intersect(rownames(metadata), cellNames)
  
  print(paste("Original ArchR project cells:", length(cellNames)))
  print(paste("Cells with metadata:", length(common_cells)))
  print(paste("Cells to remove:", length(cellNames) - length(common_cells)))
  
  # Subset metadata to common cells
  data <- metadata[common_cells, ]
  
  metadata_cols <- c(
    "rna_sdf_score", "rna_sdf_class",
    "rna_celltype_v1", "rna_majorclass_v1", "atac_celltype_v1",
    "Age", "AgeUnits", "AgeDays", "Period", "Epoch",
    "Sex", "SexForAnalysis"
  )
  
  for (colname in metadata_cols) {
    if (colname %in% colnames(data)) {
      archr_proj <- addCellColData(
        ArchRProj = archr_proj,
        data = data[[colname]],
        cells = common_cells,
        name = colname,
        force = FALSE
      )
    }
  }
  
  # Add broad cell type classification
  list.broad.type <- list(
    ExN        = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'),
    InN        = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'),
    Glia       = c('Glia'),
    NNC        = c('NNC'),
    Progenitor = c('Progenitor')
  )
  
  # Create lookup for broad types
  lookup <- unlist(
    lapply(names(list.broad.type), function(broad) {
      setNames(rep(broad, length(list.broad.type[[broad]])), list.broad.type[[broad]])
    })
  )
  
  rna_classes <- getCellColData(archr_proj, select = "rna_majorclass_v1")
  cell_broad <- lookup[as.character(rna_classes$rna_majorclass_v1)]
  cell_broad[is.na(cell_broad)] <- "Other"
  
  archr_proj <- addCellColData(
    ArchRProj = archr_proj,
    data = cell_broad,
    cells = rownames(rna_classes),
    name = "cell_broad",
    force = FALSE
  )
  
  return(archr_proj)
}

# Add metadata to both projects
prenatal_proj <- add_metadata_to_archr(prenatal_proj, metadata, "prenatal")
postnatal_proj <- add_metadata_to_archr(postnatal_proj, metadata, "postnatal")

print(paste("Loaded prenatal dataset with", nCells(prenatal_proj), "cells"))
print(paste("Loaded postnatal dataset with", nCells(postnatal_proj), "cells"))

# Execute Neuron vs Non-Neuron Differential Accessibility Analysis ============
print("STARTING NEURON VS NON-NEURON DIFFERENTIAL ACCESSIBILITY ANALYSIS")

print("=== PRENATAL ANALYSIS ===")
prenatal_processed <- process_archr_neuron_analysis(prenatal_proj, "prenatal")

print("=== POSTNATAL ANALYSIS ===")
postnatal_processed <- process_archr_neuron_analysis(postnatal_proj, "postnatal")

print("ATAC NEURON VS NON-NEURON ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!")