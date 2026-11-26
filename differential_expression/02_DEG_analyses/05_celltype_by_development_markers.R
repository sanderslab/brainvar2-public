# ==============================================================================
# DEG Analysis Pipeline for BrainVar Data - Combined Cell Type &
# Development
# ==============================================================================
# Purpose: Differential expression analysis comparing prenatal vs postnatal
#          across cell types (ExN, InN, NNC, Glia) in a single analysis.
# Author: Narjes Rohani
# Date: Aug 11, 2025
# R Version: >= 4.0.0
# ==============================================================================

# Package Management Functions ================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  } else {
    message("Package already installed: ", pkg)
  }
}

# Required Packages ============================================================

required_packages <- c(
  "Seurat",
  "future"
)

for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Load Libraries ===============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

# Global Configuration =========================================================

set.seed(123)
plan("multisession", workers = 3)
# Set maximum object size for parallel processing (~103 GB)
options(future.globals.maxSize = 110595407872)

# Data Loading =================================================================
print("Loading preprocessed Seurat object...")

# Load the combined brain development dataset
combined <- readRDS("/Users/nrohani/Library/Mobile Documents/com~apple~CloudDocs/combined_all_brainvar_rna.rds")

print(paste("Loaded dataset with", ncol(combined), "cells and", nrow(combined), "genes"))

# Cell Type Classification =====================================================
print("Assigning cell type categories from multiple sources...")

# Create a unified cell type column that pulls from both majorclass_v1 and celltype_v1
combined$unified_celltype <- combined$majorclass_v1


list.broad.type <- list(
  ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'),
  InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'),
  Glia = c('Glia'),
  NNC = c('NNC'),
  Progenitor = c('Progenitor')
)

# Map broad cell types (ExN, InN) from majorclass_v1
for (cell_type in names(list.broad.type)) {
  matching_cells <- combined$unified_celltype %in% list.broad.type[[cell_type]]
  combined$unified_celltype[matching_cells] <- cell_type
}

print("Cell type assignment complete.")
print("Unified cell types:")
print(table(combined$unified_celltype))

# Quality Control Filtering ===================================================
print("Applying quality control filters...")

# Apply stringent quality control filters
initial_cell_count <- ncol(combined)

combined <- combined[, 
                     combined$percent_mt <= 15 &      # Remove high mitochondrial content
                       combined$percent_ribo >= 0.25    # Remove low ribosomal content
]

filtered_cell_count <- ncol(combined)
print(paste("Quality control complete:", 
            initial_cell_count - filtered_cell_count, 
            "cells removed,", 
            filtered_cell_count, 
            "cells retained"))

# Developmental Period Grouping ================================================
print("Assigning developmental stage groups...")

# Create broader developmental stage groupings
combined$group <- ifelse(
  combined$period >= 3 & combined$period <= 6, "prenatal",      # Periods 3-6
  ifelse(combined$period == 7, "transition",                    # Period 7 (birth)
         "postnatal")                                           # Periods 8+
)

# Display group distribution
group_counts <- table(combined$group)
print("Developmental group distribution:")
print(group_counts)

# Filter to Cell Types of Interest ============================================
print("Filtering to cell types of interest...")

# Define cell types for analysis
cell_types_to_analyze <- c("ExN", "InN", "NNC", "Glia")

# Subset to only these cell types
combined <- combined[, combined$unified_celltype %in% cell_types_to_analyze]

print(paste("Cells after filtering to target cell types:", ncol(combined)))
print("Cell type distribution:")
print(table(combined$unified_celltype))

# Create Combined Identity for Cell Type + Development Stage ==================
print("Creating combined cell type and developmental stage identities...")

# Create a new metadata column combining cell type and developmental stage
# This allows FindMarkers to compare specific combinations
# Format: "CellType_DevelopmentStage" (e.g., "ExN_prenatal", "InN_postnatal")
combined$celltype_dev <- paste(combined$unified_celltype, 
                               combined$group, 
                               sep = "_")

# Display the combinations created
print("Cell type × Development stage combinations:")
combination_counts <- table(combined$celltype_dev)
print(combination_counts)

# Set this combined identity as the active identity
Idents(combined) <- combined$celltype_dev

# Differential Expression Analysis ============================================
print("\n========================================")
print("Starting combined differential expression analysis")
print("========================================")

# Create output directory if it doesn't exist
output_dir <- './RNAmarkers332026/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print(paste("Created output directory:", output_dir))
}

all_markers <- FindAllMarkers(
  combined,
  min.cells.group = 10,
  min.cells.feature = 10,
  only.pos = FALSE,             # Include both up- and down-regulated genes
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"           # Wilcoxon rank sum test
)

print(paste("FindAllMarkers complete. Found", nrow(all_markers), "markers"))

# Generate timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Save complete results
complete_file <- paste0(output_dir, "markers_celltype_development_", timestamp, ".csv")
write.csv(all_markers, complete_file, row.names = FALSE)
