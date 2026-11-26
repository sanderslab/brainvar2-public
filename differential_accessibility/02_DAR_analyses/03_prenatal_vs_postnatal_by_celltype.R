# ==============================================================================
# ATAC Prenatal vs Postnatal Differential Accessibility Analysis.
# ==============================================================================
# Purpose: Differential accessibility analysis comparing prenatal vs postnatal
#          brain development stages in each specific cell type using merged
#          ArchR projects
# Author: Narjes Rohani
# Date: Aug 17, 2025
# R Version: >= 4.0.0
# ==============================================================================

# Library Loading ==============================================================
print("Loading required libraries...")

library(ArchR)   # Single-nuclei ATAC-seq analysis toolkit

# Global Configuration =========================================================
print("Configuring global analysis settings...")

# Set random seed for reproducible results
set.seed(123)

# Set ArchR threading
addArchRThreads(threads = 3)

# Set maximum object size for parallel processing (~103 GB)
options(future.globals.maxSize = 110595407872)

# Configure ArchR genome
addArchRGenome("hg38")

# Data Loading =================================================================
print("Loading preprocessed ArchR projects...")

# Load and add metadata
metadata <- read.delim('/mnt/s3/Sanders_lab/BrainVar2/ArchR_brainvar/ATAC_P3-P15_meta_D240727.tsv',
                       sep = ',', header = TRUE, row.names = 1)

# Function to add metadata to ArchR project
add_metadata_to_archr <- function(archr_proj, metadata, project_name) {
  # Find common cells
  cellNames <- getCellNames(archr_proj)
  common_cells <- intersect(rownames(metadata), cellNames)
  
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
  
  return(archr_proj)
}

# Add metadata to combined project (including both prenatal and postnatal samples)
combined <- loadArchRProject("./ArchR_brainvar_metadata")
combined <- add_metadata_to_archr(combined, metadata, "all")

periods <- getCellColData(combined)$Period

# Create broader developmental stage groupings
group <- ifelse(
  periods >= 3 & periods <= 6, "prenatal",
  ifelse(periods == 7, "transition",
         "postnatal")
)

# Add group to project
combined <- addCellColData(
  ArchRProj = combined,
  data = group,
  cells = getCellNames(combined),
  name = "group",
  force = FALSE
)
valid_cells <- !is.na(group)
combined <- combined[getCellNames(combined)[valid_cells], ]
group <- group[valid_cells]  # update group vector accordingly

# Display group distribution
group_counts <- table(group)
print("Developmental group distribution:")
print(group_counts)

# Define mapping from specific to broad cell type categories
list.broad.type <- list(
  ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'),
  InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'),
  Glia = c('Glia'),
  NNC = c('NNC'),
  Progenitor = c('Progenitor')
)

# Get cell type data and create broad types
rna_classes <- getCellColData(combined, select = "rna_majorclass_v1")

# Create lookup for broad types
lookup <- unlist(
  lapply(names(list.broad.type), function(broad) {
    setNames(rep(broad, length(list.broad.type[[broad]])), list.broad.type[[broad]])
  })
)

# Map to broad types
broad_type <- lookup[as.character(rna_classes$rna_majorclass_v1)]

# Add broad type to project
combined <- addCellColData(
  ArchRProj = combined,
  data = broad_type,
  cells = getCellNames(combined),
  name = "broad.type",
  force = FALSE
)

# Get cells that are either prenatal or postnatal (exclude transition)
cells_to_keep <- getCellNames(combined)[
  getCellColData(combined)$group %in% c("prenatal", "postnatal")
]

# Subset the ArchR project to keep only these cells
combined <- combined[cells_to_keep, ]

# Verify the subsetting worked
table(getCellColData(combined)$group)

# Remove Unknown cells
cell_data <- getCellColData(combined)
unknown_cells <- which(cell_data$rna_celltype_v1 == "Unknown")
if (length(unknown_cells) > 0) {
  valid_cells <- getCellNames(combined)[-unknown_cells]
  combined <- combined[valid_cells, ]
}

# Create a new combined grouping variable
combined_grouping <- paste(
  getCellColData(combined)$broad.type,
  getCellColData(combined)$group,
  sep = "_"
)

# Add this to the project
combined <- addCellColData(
  ArchRProj = combined,
  data = combined_grouping,
  cells = getCellNames(combined),
  name = "celltype_group",
  force = TRUE
)

# Get cells that have broad_type in the specified categories
cells_to_keep <- getCellNames(combined)[
  getCellColData(combined)$broad.type %in% c("ExN", "InN", "Glia", "NNC")
]

# Subset the ArchR project
combined <- combined[cells_to_keep, ]

# Add group coverages for the combined grouping
combined <- addGroupCoverages(
  ArchRProj = combined,
  groupBy = "celltype_group",
  force = TRUE
)

# Use the combined grouping for peak calling
combined <- addReproduciblePeakSet(
  ArchRProj = combined,
  groupBy = "celltype_group",
  pathToMacs2 = findMacs2(),
  force = TRUE
)

# Add peak matrix
combined <- addPeakMatrix(combined, force = TRUE)

saveArchRProject(ArchRProj = combined, outputDirectory = "./ArchR_all_with_peaks")

# First analysis: Broad cell types ============================================
cell_types <- unique(getCellColData(combined)$broad.type)
cell_types <- cell_types[!is.na(cell_types)]
markers_by_celltype <- list()

for (ct in cell_types) {
  print(paste("Finding markers for broad type:", ct))
  
  # Get cells of this type
  cells_of_type <- getCellNames(combined)[getCellColData(combined)$broad.type == ct]
  
  if (length(cells_of_type) > 0) {
    subset_ct <- combined[cells_of_type, ]
    
    # Check group diversity
    groups <- getCellColData(subset_ct)$group
    if (length(unique(groups)) == 2) {
      
      markers <- getMarkerFeatures(
        ArchRProj = subset_ct,
        useMatrix = "PeakMatrix",
        groupBy = "group",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon",
        maxCells = 10000
      )
      
      # Extract markers
      markerList <- getMarkers(markers, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.1")
      
      # Convert to data frame
      markers_df <- data.frame()
      
      for (stage in names(markerList)) {
        if (!is.null(markerList[[stage]]) && nrow(markerList[[stage]]) > 0) {
          temp_df <- as.data.frame(markerList[[stage]])
          temp_df$stage <- stage
          temp_df$peak <- rownames(temp_df)
          markers_df <- rbind(markers_df, temp_df)
        }
      }
      
      # Filter for prenatal markers (equivalent to ident.1 = "prenatal")
      if ("prenatal" %in% markers_df$stage) {
        prenatal_markers <- markers_df[markers_df$stage == "prenatal", ]
        markers_by_celltype[[ct]] <- prenatal_markers
        
        dir.create("ATACmarkers_9152025", showWarnings = FALSE)
        out_file <- paste0("ATACmarkers_9152025/", "markers_", ct, "_prenatal_vs_postnatal.csv")
        write.csv(prenatal_markers, out_file, row.names = FALSE)
        print(paste("    Saved:", out_file))
      }
      
    } else {
      print(paste("    Skipping", ct, "— not enough group diversity"))
    }
  }
}

# Second analysis: Individual cell subtypes ===================================
cell_types <- unique(getCellColData(combined)$rna_celltype_v1)
cell_types <- cell_types[!is.na(cell_types)]
markers_by_celltype <- list()

for (ct in cell_types) {
  message(paste("Finding markers for cell subtype:", ct))
  
  cells_of_type <- getCellNames(combined)[getCellColData(combined)$rna_celltype_v1 == ct]
  
  if (length(cells_of_type) > 0) {
    subset_ct <- combined[cells_of_type, ]
    groups <- getCellColData(subset_ct)$group
    
    if (length(unique(groups)) == 2) {
      
      tryCatch({
        # Handle warnings gracefully
        withCallingHandlers({
          markers <- getMarkerFeatures(
            ArchRProj = subset_ct,
            useMatrix = "PeakMatrix",
            groupBy = "group",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "wilcoxon",
            maxCells = 10000
          )
        }, warning = function(w) {
          message(paste("Warning in cell type:", ct, "->", conditionMessage(w)))
          invokeRestart("muffleWarning")
        })
        
        markerList <- getMarkers(markers, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
        
        safe_cluster_name <- gsub("/", "_", ct)
        out_file <- paste0("ATACmarkers/", "markers_", safe_cluster_name, "_prenatal_vs_postnatal.csv")
        write.csv(data.frame(markerList), out_file, row.names = FALSE)
        message(paste("Saved:", out_file))
        
      }, error = function(e) {
        message(paste("Error in cell type:", ct, "->", e$message))
      }, finally = {
        message(paste("Finished processing cell type:", ct))
      })
      
    } else {
      message(paste("Skipping", ct, "— not enough group diversity"))
    }
  }
}

print('Analysis Finished.')