
# ==============================================================================
# DEG Analysis Pipeline for BrainVar Data
# ==============================================================================
# Purpose: Compare gene expression between neurons and non-neuronal cell types
#          in prenatal and postnatal brain development datasets
# Author: Narjes Rohani
# Date: Aug 11, 2025
# R Version: >= 4.0.0
# Analysis: Neuron vs Non-Neuron DEG identification using FindAllMarkers
# ==============================================================================

# Package Management ===========================================================

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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
})

# Global Configuration =========================================================

set.seed(123)
plan("multisession", workers = 2)
options(future.globals.maxSize = 110595407872)

# Main Analysis Function =======================================================

#' Process Seurat object and perform differential expression analysis
#' 
#' This function performs comprehensive cell type classification and differential
#' expression analysis comparing neurons vs non-neuronal cells
#' 
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param group_label Character string identifying the dataset (e.g., "prenatal")
#' @return NULL (saves results as CSV files)
process_seurat <- function(seurat_obj, group_label) {
  
  print(paste("Processing", group_label, "dataset..."))
  print(paste("Initial cell count:", ncol(seurat_obj)))
  print(paste("Initial gene count:", nrow(seurat_obj)))

  # Cell Type Classification =================================================
  print("Step 1: Assigning broad cell type categories...")
  
  # Define hierarchical cell type mapping
  # Maps specific developmental cell types to broader functional categories
  list.broad.type <- list(
    # Excitatory neurons and their developmental precursors
    ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'),
    
    # Inhibitory neurons from different ganglionic eminences
    InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'),
    
    # Glial cell populations (astrocytes, oligodendrocytes, microglia)
    Glia = c('Glia'),
    
    # Non-neural cells (vasculature, immune, meningeal, etc.)
    NNC = c('NNC'),
    
    # Neural progenitor and stem cell populations
    Progenitor = c('Progenitor')
  )
  
  # Initialize broad type classification using existing annotations
  seurat_obj$broad.type <- seurat_obj$majorclass_v1
  
  # Apply hierarchical mapping to create broad categories
  for (cell_type in names(list.broad.type)) {
    matching_cells <- seurat_obj$broad.type %in% list.broad.type[[cell_type]]
    seurat_obj$broad.type[matching_cells] <- cell_type
  }
  
  print(paste("Broad cell types assigned:", 
              paste(unique(seurat_obj$broad.type), collapse = ", ")))
  
  # Quality Control Filtering ================================================
  print("Step 2: Applying quality control filters...")
  
  initial_count <- ncol(seurat_obj)
  
  # Apply stringent quality control criteria
  # Remove cells showing signs of stress, death, or poor capture
  seurat_obj <- seurat_obj[, 
                           seurat_obj$percent_mt <= 15 &      # Mitochondrial gene threshold (stress indicator)
                             seurat_obj$percent_ribo >= 0.25    # Ribosomal gene threshold (activity indicator)
  ]
  
  final_count <- ncol(seurat_obj)
  cells_removed <- initial_count - final_count
  
  print(paste("Quality control applied:"))
  print(paste("  - Cells removed:", cells_removed))
  print(paste("  - Cells retained:", final_count))
  print(paste("  - Retention rate:", round(final_count/initial_count * 100, 1), "%"))
  
  # Neuron vs Non-Neuron Classification ======================================
  print("Step 3: Classifying cells into Neuron vs Non-Neuron categories...")
  
  # Define specific non-neuronal cell types for targeted comparison
  # These are well-characterized, functionally distinct cell populations
  # non_neuron_types <- c(
  #   "VLMC",         # Vascular cells
  #   "SMC",         # Smooth Muscle cells
  #   "Oligo",        # Oligodendrocytes
  #   "Immune",       # Immune cells
  #   "Astrocyte",    # Glial cells providing metabolic and structural support
  #   "Micro",        # Microglia (brain-resident immune cells)
  #   "Endo",         # Endothelial cells (blood-brain barrier)
  #   "OPC"           # Oligodendrocyte precursor cells (myelination)
  # )
  non_neuron_types <- c(

    "Astrocyte",    # Glial cells providing metabolic and structural support
    "Micro",        # Microglia (brain-resident immune cells)
    "Endo",         # Endothelial cells (blood-brain barrier)
    "OPC"           # Oligodendrocyte precursor cells (myelination)
  )
  non_neuron_types <- c(
    "VLMC",         # Vascular cells
    "SMC",         # Smooth Muscle cells
    "Oligo",        # Oligodendrocytes
    "Astrocyte",    # Glial cells providing metabolic and structural support
    "Micro",        # Microglia (brain-resident immune cells)
    "Endo",         # Endothelial cells (blood-brain barrier)
    "OPC"           # Oligodendrocyte precursor cells (myelination)
  )
  # Create binary classification: Neuron vs specific non-neuronal types
  # Excludes "Unknown" cells to ensure clear cell type identity
  exclude_types <- c("Unknown", "Oligo", "Immune")
  exclude_types <- c("Unknown", "VLMC", "SMC", "Oligo", "Immune")
  
  seurat_obj$cell_type_group <- ifelse(
    seurat_obj$celltype_v1 %in% exclude_types,
    NA,
    ifelse(seurat_obj$celltype_v1 %in% non_neuron_types,
           seurat_obj$celltype_v1,  # Keep the specific non-neuronal type name
           "Neuron"                 # All other non-excluded types become "Neuron"
    )
  )
  print(unique(seurat_obj$cell_type_group))
  ## Remove cells with unknown identity for clean comparison
  #valid_cells <- !is.na(seurat_obj$cell_type_group)
  #seurat_obj <- seurat_obj[, valid_cells]
  
  # Display final cell type distribution
  #cell_type_counts <- table(seurat_obj$cell_type_group)
  #print("Final cell type distribution:")
  #print(cell_type_counts)
  
  # Validate sufficient cells for comparison
  if (length(unique(seurat_obj$cell_type_group)) < 2) {
    stop(paste("Insufficient cell type diversity in", group_label, 
               "dataset. Need at least 2 cell types for comparison."))
  }
  
  # Differential Expression Analysis =========================================
  print("Step 4: Performing differential expression analysis...")
  print("  Analysis type: Neuron vs Non-Neuron (all pairwise comparisons)")
  
  # Set cell identities for differential expression analysis
  Idents(seurat_obj) <- seurat_obj$cell_type_group
  
  # Verify identities are properly set
  print(paste("Active identities:", paste(levels(Idents(seurat_obj)), collapse = ", ")))
  
  # Perform comprehensive differential expression analysis
  # FindAllMarkers compares each cell type against all others
  print("  Running FindAllMarkers analysis...")
  
  # This analysis design
  # - Compares each cell type (Neuron, Astrocyte, Micro, OPC)
  # - Uses Wilcoxon rank-sum test (default) for statistical testing
  # - Returns both upregulated and downregulated genes (only.pos = FALSE)

  markers_fine_all <- FindAllMarkers(
    seurat_obj, 
    min.cells.group = 10,
    min.cells.feature = 10,
    only.pos = FALSE,  # Include both up- and down-regulated genes
    group.by = "cell_type_group",    # Group variable for comparison
    verbose = TRUE                   # Show progress during computation
  )
  
  # Add metadata to results for better interpretation
  markers_fine_all$dataset <- group_label
  markers_fine_all$analysis_type <- "neuron_vs_nonneuron"
  markers_fine_all$analysis_date <- Sys.Date()
  
  # Generate output filename with descriptive naming
  output_file <- paste0('/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Enrichment_analysis/RNAmarkers1692025/', group_label, "_markers_neuronvsnn_celltypes_all.csv")
  
  # Save comprehensive results
  write.csv(markers_fine_all, output_file, row.names = FALSE)
  
  print(paste("Analysis completed for", group_label, "dataset!"))
  print(paste("Results saved to:", output_file))
  print(paste("Total DEGs identified:", nrow(markers_fine_all)))
  print(paste("Significant DEGs (p_val_adj < 0.05):", 
              sum(markers_fine_all$p_val_adj < 0.05, na.rm = TRUE)))
  
  # Display summary by cell type
  print("\nDEG summary by cell type:")
  deg_summary <- markers_fine_all %>%
    group_by(cluster) %>%
    summarise(
      total_genes = n(),
      significant_genes = sum(p_val_adj < 0.05, na.rm = TRUE),
      upregulated = sum(avg_log2FC > 0, na.rm = TRUE),
      downregulated = sum(avg_log2FC < 0, na.rm = TRUE),
      .groups = 'drop'
    )
  print(deg_summary)
  
}

# Data Loading and Analysis Execution =========================================
print("Loading preprocessed datasets...")

# Load both developmental stage datasets
# These should be pre-processed Seurat objects with quality control applied
print("Loading postnatal dataset...")
postnatal <- readRDS("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Data/combined_downsample_brainvar_rna_postnatal.rds")

print("Loading prenatal dataset...")
prenatal <- readRDS("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Data/combined_downsample_brainvar_rna_prenatal.rds")

# Execute differential expression analysis for both datasets
print("STARTING DIFFERENTIAL EXPRESSION ANALYSIS")

# Process prenatal dataset
process_seurat(prenatal, "prenatal")

print("\n")

# Process postnatal dataset  
process_seurat(postnatal, "postnatal")


print("ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!")

print("Output files generated:")
print("  - prenatal_markers_neuronvsnn_celltypes_all.csv")
print("  - postnatal_markers_neuronvsnn_celltypes_all.csv")

# ==============================================================================
# Analysis Complete 
# ==============================================================================
