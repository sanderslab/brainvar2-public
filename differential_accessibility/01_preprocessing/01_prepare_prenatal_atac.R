############################
# 31/01/2025 Narjes Rohani #
# Code for calculating     #
# pseudobulk ATAC &        #
# normalisation.           #
# coverage plot            #
############################
set.seed(123)
# Import necessary packages
library("rtracklayer")  # v1.64.0
library("Seurat")       # v4.4.0
library("tidyverse")    # v2.0.0
library("Signac")       # v1.14.0
library("data.table")   # v1.16.0
library(gtools)         # v3.9.5
library(viridis)        # v0.6.5
library(AnnotationHub)
library("EnsDb.Hsapiens.v86");
library("BSgenome.Hsapiens.UCSC.hg38");
print('Running started...')
##############################################
# Global variables (file locations)
# retrieve gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86);
seqlevelsStyle(annotations) <- "UCSC";
genome(annotations) <- "hg38";
path_atac_data <- "ATAC_data/"
path_fragment_data <- "brainvar_fragments/"

sample_meta_data_path <- "/mnt/s3/Sanders_lab/BrainVar2/Data/brainvar2_single_cell_detailed_metadata_v1.3.tsv"
dataset_meta_data_path <- "/mnt/s3/Sanders_lab/BrainVar2/Data/brainvar2_library_metadata_v1.3.tsv"

print('preparing metadata...')

# Read and process sample metadata
sample_meta <- read_delim(sample_meta_data_path)
sample_meta <- sample_meta[sample_meta$Usable_ATAC == 'Yes',]
dataset_meta <- read_delim(dataset_meta_data_path)
dataset_meta <- dataset_meta[dataset_meta$LibraryModality == 'ATAC' & dataset_meta$Usable == 'Yes',]

# Label datasets
dataset_meta$dataset <- ifelse(dataset_meta$Werling2020 == "Yes", "BrainVar1", "BrainVar2")
dataset_meta$dataset <- ifelse(grepl("^[0-9]", dataset_meta$BioSampleID), "MG", dataset_meta$dataset)
dataset_meta <- dataset_meta[, c("BioSampleID", "dataset")]

# Merge dataset metadata with sample metadata
sample_meta <- merge(sample_meta, dataset_meta, by = "BioSampleID", all.x = TRUE)
sample_meta$dataset <- ifelse(grepl("^[0-9]", sample_meta$BioSampleID), "MG", sample_meta$dataset)

# Initialize list to store Seurat objects
seurat_list <- list()
fragment_list <- list()
print('reading and creating seurat objects...')
# Loop through periods P3 to P15
for (P in paste0("P", 3:6)) {
  print(paste("Processing", P, "..."))
  
  # Define file paths
  peak_matrix_path <- sprintf("%s%s_majorclass_Peakmatrix.rds", path_atac_data, P)
  cell_meta_path <- sprintf("%s%s_majorclass_cellmeta.rds", path_atac_data, P)
  fragment_path <- sprintf("%s%s_atac_fragments.tsv.gz", path_fragment_data, P)
  
  # Read RDS files
  peak_matrix <- readRDS(peak_matrix_path)
  cell_meta <- readRDS(cell_meta_path)
  
  # Filter sample metadata
  sample_meta_subset <- sample_meta[sample_meta$BioSampleID %in% cell_meta$Sample, ]
  sample_meta_subset <- sample_meta_subset[, c("BioSampleID", "dataset", "SexForAnalysis", "Age", "AgeUnits",
                                               "AgeDays", "Period", "Epoch", "AncestryPCA", "PMI_hours",
                                               "pH", "SexReported", "XIST_depth_RNA", "AncestryReported",
                                               "matchingDNALibrary", "CaptureID")]
  
  # Convert cell metadata and recode cell types
  cell_meta <- as.data.frame(cell_meta)
  cell_meta$cell_type <- as.factor(cell_meta$majorclass.atac)
  cell_meta$cell_type <- dplyr::recode(cell_meta$cell_type,
                                       "Astro.OPC.Oligo.gIP" = "Glia",
                                       "Astro/OPC/Oligo/gIP" = "Glia",
                                       "RG" = "Glia",
                                       "Blood" = "NN",
                                       "Micro.Immune" = "NN",
                                       "nIP" = "Progenitor",
                                       "Micro/Immune" = "NN",
                                       "Astro" = "Glia",
                                       "Oligo" = "Glia",
                                       "OPC" = "Glia")
  
  # Merge cell metadata with sample metadata
  rownames_cell_meta <- rownames(cell_meta)
  cell_meta <- merge(cell_meta, sample_meta_subset, by.x = "Sample", by.y = "BioSampleID", all.x = TRUE)
  rownames(cell_meta) <- rownames_cell_meta
  
  # Create Chromatin Assay
  cell_meta_new <- cell_meta[match(colnames(peak_matrix@assays@data@listData[["PeakMatrix"]]), rownames(cell_meta)), ]
  chrom_assay <- CreateChromatinAssay(
    counts = peak_matrix@assays@data$PeakMatrix,
    ranges = peak_matrix@rowRanges,
    sep = c(":", "-"),
    validate.fragments = TRUE,
    min.cells = 0,
    min.features = 0,
    fragments = fragment_path
  )
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = cell_meta)
  peaks.keep <- seqnames(granges(seu)) %in% standardChromosomes(granges(seu))
  seu <- seu[as.vector(peaks.keep), ]
  seu <- AddMetaData(seu, cell_meta_new)
  print('Done!')
  print('TF-IDF normalisation...')
  #temp <- RunTFIDF(seu)
  print('Saving the normalised seurat object for the period...')
  #saveRDS(temp, file = sprintf("./%s_seurat_object.rds", P))
  print('Saved!')
  # Add period identifier
  seu$period_data <- P
  frag <- CreateFragmentObject(
    path = fragment_path,
    cells = rownames(cell_meta)
  )
  # Store in list
  seurat_list[[P]] <- seu
  fragment_list[[P]] <- frag
  
}

print("All periods processed!")

print('Combining periods...')
# Unify peak regions across datasets
combined.peaks <- UnifyPeaks(object.list = seurat_list, mode = "reduce")

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
print('# Compute feature matrices for each dataset')
# Compute feature matrices for each dataset
for (P in names(seurat_list)) {
  seu <- seurat_list[[P]]
  seu.counts <- FeatureMatrix(
    fragments = fragment_list[[P]],
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(seu)
  )
  seu[['re_peaks']] <- CreateChromatinAssay(
        seu.counts,
        fragments = Fragments(seu),
        sep = c(":", "-"),
        genome = "hg38",
        annotation = annotations);
  DefaultAssay(seu) <- "re_peaks";

  
  seurat_list[[P]] <- seu
}
print('Merging all periods...')
# Merge all datasets
combined <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))

print("Merged merged Seurat object created!")
print('Saving the combined object..')
combined <- RunTFIDF(combined)

saveRDS(combined, file = 'combined_prenatal_seurat_object.rds')
print('combined object saved!')
print('Finished!')
