# ==============================================================================
# DEG Analysis Pipeline for BrainVar Data
# ==============================================================================
# Purpose: Merge BrainVar data and downsample (random 3000 cells per cell type)
#          to make it ready for DEG analysis.
# Author: Narjes Rohani
# Date: Aug 11, 2025
# R Version: >= 4.0.0
# ==============================================================================

# Load Required Libraries =====================================================

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  } else {
    message("Package already installed: ", pkg)
  }
}

required_packages <- c(
  "Seurat",
  "dplyr",
  "future"
)

for (pkg in required_packages) {
  install_if_missing(pkg)
}

if (!requireNamespace("schard", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", dependencies = TRUE)
  }
  devtools::install_github("cellgeni/schard", upgrade = "never")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(schard)
})

# Global Configuration ========================================================

print("Setting global configuration...")
set.seed(123)
plan("multisession", workers = 8)
options(future.globals.maxSize = 110595407872)

# Data Loading and Integration ================================================
print("Loading Seurat objects for each developmental period...")

# Define developmental periods to analyze (periods 3-6 for prenatal and 8-15 for postnatal)
periods <- 3:6

# Load individual h5ad files and convert to Seurat objects
seu_list <- lapply(periods, function(p) {
  # Construct file path for each period
  path <- paste0('./v1_annotation_by_period/BrainVar_v1_anno_period_', p, '.h5ad')
  print(paste(" Loading:", path))
  
  # Convert h5ad format to Seurat object
  schard::h5ad2seurat(path)
})

# Assign meaningful names to list elements
names(seu_list) <- paste0("seu.RNA.p", periods)

# Merge all datasets into single Seurat object
print("Merging Seurat objects into a single dataset...")
combined <- merge(
  seu_list[[1]],                    # First dataset as base
  y = seu_list[2:5],               # Remaining datasets to merge
  add.cell.ids = 'prenatal'        # Add prefix to cell barcodes
)

# Clean up memory by removing intermediate objects
rm(seu_list)
gc()  # Force garbage collection

# Cell Type Classification ====================================================
print("Assigning broad cell type categories...")

# Define mapping from specific to broad cell type categories
list.broad.type <- list(
  # Excitatory neurons and precursors
  ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
  
  # Inhibitory neurons and precursors  
  InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
  
  # Glial cells
  Glia = c('Glia'), 
  
  # Non-neural cells
  NNC = c('NNC'), 
  
  # Progenitor
  Progenitor = c('Progenitor')
)

# Initialize broad type classification with existing major class
combined$broad.type <- combined$majorclass_v1

# Map specific cell types to broad categories
for (t in names(list.broad.type)) {
  combined$broad.type[combined$broad.type %in% list.broad.type[[t]]] <- t
}

# Quality Control and Filtering ===============================================
print("Filtering cells based on quality control metrics...")

# Filter cells with poor quality indicators:
# - High mitochondrial content (>15%) indicates dying/stressed cells
# - Low ribosomal content (<0.25%) indicates low metabolic activity
combined <- combined[, combined$percent_mt <= 15 & combined$percent_ribo >= 0.25]

# Cell Downsampling for Computational Efficiency ==============================
print("Performing cell downsampling to balance cell type representation...")

# Set maximum number of cells per cell type to downsample
max_cells <- 3000

# Downsample cells while maintaining cell type diversity
downsampled_cells <- combined@meta.data %>%
  mutate(cell = rownames(.)) %>%              # Add cell barcode as column
  group_by(celltype_v1) %>%                   # Group by detailed cell type
  sample_n(size = min(n(), max_cells)) %>%    # Sample up to max_cells per type
  pull(cell)                                  # Extract cell barcodes

# Create downsampled Seurat object
combined_downsampled <- subset(combined, cells = downsampled_cells)

# Standard Preprocessing Pipeline ==============================================
print("Running standard Seurat preprocessing pipeline...")

# Step 1: Normalize gene expression data
# Uses LogNormalization: log(count/total_count * scale_factor + 1)
combined_downsampled <- NormalizeData(combined_downsampled)

# Step 2: Identify highly variable genes for downstream analysis
# Finds genes with high cell-to-cell variation
combined_downsampled <- FindVariableFeatures(combined_downsampled)

# Step 3: Scale data and regress out unwanted variation
# Centers and scales gene expression values
combined_downsampled <- ScaleData(combined_downsampled)

# Step 4: Perform Principal Component Analysis (PCA)
# Reduces dimensionality while preserving major variation
combined_downsampled <- RunPCA(combined_downsampled)

# Step 5: Generate UMAP embedding for visualization
# Creates 2D representation preserving local neighborhood structure
combined_downsampled <- RunUMAP(combined_downsampled, dims = 1:20)

# Save Processed Data =========================================================
print("Saving processed Seurat object...")

# Save the complete and downsample processed object for future analysis
# Note: Saving the full 'combined' object, not just downsampled version
saveRDS(combined, file = "combined_brainvar_rna_prenatal.rds")
saveRDS(combined_downsampled, file = "combined_downsample_brainvar_rna_prenatal.rds")

print("Analysis pipeline completed successfully!")

# ==============================================================================
# Analysis Complete. Repeat for postnatal as well.
# ==============================================================================