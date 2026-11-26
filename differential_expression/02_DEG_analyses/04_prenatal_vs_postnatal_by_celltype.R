# ==============================================================================
# DEG Analysis Pipeline for BrainVar Data
# ==============================================================================
# Purpose: Find cell-type specific DEGs for prenatal vs postnatal
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

# Install Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}

# Core Package Installation ====================================================

required_packages <- c(
  "ArchR",
  "Seurat",
  "future"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "ArchR") {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Library Loading =============================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(future)
})

# Global Configuration =========================================================

set.seed(123)
plan("multisession", workers = 10)
addArchRThreads(threads = 3)
options(future.globals.maxSize = 110595407872)
addArchRGenome("hg38")

# Set the main directory
setwd("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Enrichment_analysis")

# Read the combined rds data including both prenatal and postnatal samples
seurat_obj <- readRDS("/Users/nrohani/Documents/Sanders_lab/combined_all_brainvar_rna.rds")

# Remove cells with unknown celltypes
seurat_obj <- subset(seurat_obj, subset = celltype_v1 != "Unknown")

# Run the analysis for major cell types
cell_types <- c("ExN", "InN", "NNC", "Glia")
markers_by_celltype <- list()
Idents(seurat_obj) <- seurat_obj$broad.type
for (ct in cell_types) {
  print(paste("Finding markers for broad type:", ct))
  subset_ct <- subset(seurat_obj, idents = ct)
  Idents(subset_ct) <- subset_ct$group
  if (length(unique(subset_ct$group)) == 2) {
    markers <- FindMarkers(subset_ct,min.cells.group = 10,only.pos = FALSE,min.cells.feature = 10, 
                           ident.1 = "prenatal", ident.2 = "postnatal", group.by = "group")
    markers_by_celltype[[ct]] <- markers
    out_file <- paste0("markers_9152025/","markers_", ct, "_prenatal_vs_postnatal.csv")
    write.csv(markers, out_file)
    print(paste("    Saved:", out_file))
  } else {
    print(paste("    Skipping", ct, "— not enough group diversity"))
  }
}
#--------------------------------------------
# Run the analysis for cell subtypes
cell_types <- unique(seurat_obj@meta.data$celltype_v1)
markers_by_celltype <- list()

Idents(seurat_obj) <- seurat_obj$celltype_v1
for (ct in cell_types) {
  print(paste("Finding markers for cell subtype:", ct))
  subset_ct <- subset(seurat_obj, idents = ct)
  Idents(subset_ct) <- subset_ct$group
  if (length(unique(subset_ct$group)) == 2) {
    tryCatch({
      markers <- FindMarkers(
        subset_ct,
        min.cells.group = 10,
        only.pos = FALSE,
        min.cells.feature = 10,
        ident.1 = "prenatal",
        ident.2 = "postnatal",
        group.by = "group"
      )
      
      markers_by_celltype[[ct]] <- markers
      safe_cluster_name <- gsub("/", "_", ct)
      
      out_file <- paste0("markers_9152025/", "markers_", safe_cluster_name, "_prenatal_vs_postnatal.csv")
      write.csv(markers, out_file)
      print(paste("    Saved:", out_file))
    },
    error = function(e) {
      message(paste("Error in cell type:", ct, "->", e$message))
    },
    warning = function(w) {
      message(paste("Warning in cell type:", ct, "->", w$message))
      invokeRestart("muffleWarning")  # suppress repeated warnings
    },
    finally = {
      message(paste("✔ Finished processing cell type:", ct))
    })
    } else {
    print(paste("    Skipping", ct, "— not enough group diversity"))
  }
}
print("\nAnalysis pipeline completed successfully!")

# ==============================================================================
# Analysis Complete
# ==============================================================================