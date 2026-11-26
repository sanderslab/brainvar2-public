process_celltype <- function(cell_type) {
  # ----------- INPUT PATHS -----------
  qnorm_file <- paste0('~/Dropbox/BrainVar2_WD/Tables/caQTL_input/Prenatal/prenatal.', cell_type, '_atac.lcpm.qnorm')
  peakloc_file <- '~/Dropbox/BrainVar2_WD/Tables/prenatal.peak_locs.txt'
  covariate_file <- paste0('~/Dropbox/BrainVar2_WD/Tables/caQTL_input/Prenatal/table.prenatal.', cell_type, '.cov.20250730.txt')
  output_file <- paste0('~/Dropbox/BrainVar2_WD/Tables/caQTL_input/Prenatal/prenatal.', cell_type, '_phenos_tensorqtl.bed')
  
  # ----------- READ DATA -----------
  d0 <- fread(qnorm_file)
  d1 <- fread(peakloc_file)
  covs <- fread(covariate_file)
  
  d0 <- d0 %>% dplyr::rename(ATAC_ID = V1)
  
  merged <- d1 %>%
    right_join(d0, by = "ATAC_ID") %>%
    relocate(`#chr`, start, end, ATAC_ID)
  
  merged$`#chr` <- paste0("chr", merged$`#chr`)
  
  # ----------- SAMPLE MAPPING -----------
  sample_map <- covs %>% dplyr::select(id, BioSample)
  valid_samples <- intersect(sample_map$BioSample, colnames(merged))
  
  rename_map <- sample_map %>%
    filter(BioSample %in% valid_samples) %>%
    distinct()
  
  merged_subset <- merged %>%
    dplyr::select(c("#chr", "start", "end", "ATAC_ID"), all_of(rename_map$BioSample))
  
  bio_cols <- intersect(colnames(merged_subset), rename_map$BioSample)
  id_names <- rename_map$id[match(bio_cols, rename_map$BioSample)]
  
  new_colnames <- colnames(merged_subset)
  new_colnames[match(bio_cols, new_colnames)] <- id_names
  colnames(merged_subset) <- new_colnames
  
  # ----------- FILTER PEAKS -----------
  
  expr_mat <- merged_subset %>%
    dplyr::select(-c(`#chr`, start, end, ATAC_ID))
  
  n_samples <- ncol(expr_mat)
  max_zeros <- n_samples * 0.2
  
  nonzero_filter <- rowSums(expr_mat == 0, na.rm = TRUE) <= max_zeros
  filtered_merged <- merged_subset[nonzero_filter, ]
  
  cat("\n============================\n")
  cat("Cell type:", cell_type, "\n")
  cat("Original peaks:", nrow(merged_subset), "\n")
  cat("Filtered peaks:", nrow(filtered_merged), "\n")
  
  # ----------- Ensure donor order matches -----------
  donor_ids <- covs$id
  final_expr <- filtered_merged %>%
    dplyr::select(`#chr`, start, end, ATAC_ID, all_of(donor_ids))
  
  stopifnot(identical(colnames(final_expr)[5:ncol(final_expr)], donor_ids))
  
  # ----------- SAVE -----------
  fwrite(
    filtered_merged,
    file = output_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
}


cell_types <- c("ExN", "InN", "NNC", "Glia", "Progenitor")

for (ct in cell_types) {
  process_celltype(ct)
}


