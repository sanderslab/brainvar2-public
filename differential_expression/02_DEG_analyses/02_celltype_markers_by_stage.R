# ==============================================================================
# DEG Analysis Pipeline for BrainVar Data
# ==============================================================================
# Purpose: Differential expression analysis on snRNA-seq data using Seurat.
#          Processes prenatal and postnatal BrainVar datasets and generates
#          marker genes for different cell type classifications.
# Author: Narjes Rohani
# Date: Aug 11, 2025
# R Version: >= 4.0.0
# ==============================================================================

# ------------------------------------------------------------------------------
# Package Installation and Loading
# ------------------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  } else {
    message("Package already installed: ", pkg)
  }
}

required_packages <- c(
  "Seurat",   # Single-cell analysis toolkit
  "dplyr",    # Data manipulation
  "future"    # Parallel processing
)

for (pkg in required_packages) {
  install_if_missing(pkg)
}

print("Loading libraries...")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

# ------------------------------------------------------------------------------
# Global Settings Configuration
# ------------------------------------------------------------------------------

setwd("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Enrichment_analysis")

print("Setting seed, future plan, and memory options...")

set.seed(123)
plan("multisession", workers = 3)

# options(future.globals.maxSize = 110595407872)

# ------------------------------------------------------------------------------
# Main Processing Function
# ------------------------------------------------------------------------------

# Function to process Seurat object and perform differential expression analysis
# @param seurat_obj: Seurat object containing single-nuclei data
# @param group_label: String identifier for the dataset (e.g., "prenatal", "postnatal")
process_seurat <- function(seurat_obj, group_label) {
  print(paste("Processing", group_label, "dataset..."))
  
  # Define broad cell type categories for classification
  # This mapping groups detailed cell types into broader functional categories
  list.broad.type <- list(
    ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'),
    InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'),
    Glia = c('Glia'),
    NNC = c('NNC'),
    Progenitor = c('Progenitor')
  )
  
  # Create broad cell type classification based on majorclass_v1
  seurat_obj$broad.type <- seurat_obj$majorclass_v1
  seurat_obj <- subset(seurat_obj, subset = celltype_v1 != "Unknown")
  
  # Map detailed cell types to broad categories
  for (cell_type in names(list.broad.type)) {
    seurat_obj$broad.type[seurat_obj$broad.type %in% list.broad.type[[cell_type]]] <- cell_type
  }
  
  # ------------------------------------------------------------------------------
  # Quality Control Filtering
  # ------------------------------------------------------------------------------
  
  # Filter out low-quality cells based on mitochondrial and ribosomal gene content
  # Criteria: ≤15% mitochondrial genes AND ≥0.25% ribosomal genes
  seurat_obj <- seurat_obj[, seurat_obj$percent_mt <= 15 & seurat_obj$percent_ribo >= 0.25]
  
  print(paste("  Filtered to", ncol(seurat_obj), "high-quality cells"))
  
  # ------------------------------------------------------------------------------
  # Differential Expression Analysis - Broad Cell Types
  # ------------------------------------------------------------------------------
  
  print("  Running DE analysis by broad cell types...")
  
  # Set cell identities to broad cell types for DE analysis
  Idents(seurat_obj) <- seurat_obj$broad.type
  
  # Perform differential expression analysis across all broad cell types
  # This finds markers that distinguish each broad cell type from all others
  markers_broad_all <- FindAllMarkers(
    seurat_obj, 
    min.cells.group = 10,
    min.cells.feature = 10,
    group.by = "broad.type", 
    only.pos = FALSE  # Include both up- and down-regulated genes
  )
  
  # Save results to CSV file
  output_file_broad <- paste0('./RNAmarkers1692025/',group_label, "_markers_broad_celltypes_all.csv")
  write.csv(markers_broad_all, output_file_broad, row.names = FALSE)
  print(paste("  Saved:", output_file_broad))
  
  # ------------------------------------------------------------------------------
  # Differential Expression Analysis - Major Cell Classes
  # ------------------------------------------------------------------------------
  
  print("  Running DE analysis by majorclass_v1...")
  
  # Set cell identities to major cell classes for DE analysis
  Idents(seurat_obj) <- seurat_obj$majorclass_v1
  
  # Perform differential expression analysis across all major cell classes
  markers_majorclass_all <- FindAllMarkers(
    seurat_obj, 
    min.cells.group = 10,
    min.cells.feature = 10,
    group.by = "majorclass_v1", 
    only.pos = FALSE  # Include both up- and down-regulated genes
  )
  
  # Save results to CSV file
  output_file_major <- paste0('./RNAmarkers1692025/',group_label, "_markers_majorclass_v1_celltypes_all.csv")
  write.csv(markers_majorclass_all, output_file_major, row.names = FALSE)
  print(paste("  Saved:", output_file_major))
  
  # ------------------------------------------------------------------------------
  # Differential Expression Analysis - Fine Cell Types
  # ------------------------------------------------------------------------------
  
  print("  Running DE analysis by fine cell types (celltype_v1)...")
  
  # Set cell identities to fine-grained cell types for DE analysis
  Idents(seurat_obj) <- seurat_obj$celltype_v1
  
  # Perform differential expression analysis across all fine cell types
  # This provides the most detailed view of cell type-specific gene expression
  markers_fine_all <- FindAllMarkers(
    seurat_obj, 
    min.cells.group = 10,
    min.cells.feature = 10,
    group.by = "celltype_v1", 
    only.pos = FALSE  # Include both up- and down-regulated genes
  )
  
  # Save results to CSV file
  output_file_fine <- paste0('./RNAmarkers1692025/',group_label, "_markers_minor_celltypes_all.csv")
  write.csv(markers_fine_all, output_file_fine, row.names = FALSE)
  print(paste("  Saved:", output_file_fine))
  
  print(paste("  Completed processing of", group_label, "dataset"))
}

# ------------------------------------------------------------------------------
# Main Analysis Execution
# ------------------------------------------------------------------------------

print("Loading prenatal and postnatal datasets...")

# Load preprocessed and downsampled brain single-nuclei datasets
# These should be Seurat objects with cell type annotations already assigned. Change to your paths
postnatal <- readRDS("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Data/combined_downsample_brainvar_rna_postnatal.rds")
prenatal <- readRDS("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Data/combined_downsample_brainvar_rna_prenatal.rds")

print(paste("Loaded prenatal dataset with", ncol(prenatal), "cells"))
print(paste("Loaded postnatal dataset with", ncol(postnatal), "cells"))


# ------------------------------------------------------------------------------
# Execute Differential Expression Analysis
# ------------------------------------------------------------------------------

# Process prenatal dataset
print("STARTING PRENATAL ANALYSIS")
process_seurat(prenatal, "prenatal")

# Process postnatal dataset  
print("STARTING POSTNATAL ANALYSIS")
process_seurat(postnatal, "postnatal")

print("ANALYSIS COMPLETE")
print("All differential expression results have been saved to CSV files.")