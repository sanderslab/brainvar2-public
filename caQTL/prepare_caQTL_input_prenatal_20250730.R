rm(list = ls())
library(tidyverse)
date = '20250730'
cell_types <- c("ExN", "InN", "Glia", "NNC", "Progenitor")

example = read.delim('~/Dropbox/BrainVar2_WD/Data/prenatal/ExN.cov_0425.txt')
duplicated(example$id) %>% table()

example = read.delim('~/Dropbox/BrainVar2_WD/Data/prenatal/InN.cov_0425.txt')
duplicated(example$id) %>% table()

example = read.delim('~/Dropbox/BrainVar2_WD/Data/prenatal/Glia.cov_0425.txt')
duplicated(example$id) %>% table()

example = read.delim('~/Dropbox/BrainVar2_WD/Data/prenatal/NNC.cov_0425.txt')
duplicated(example$id) %>% table()

if(T){
  files <- c(
    '~/Dropbox/BrainVar2_WD/Data/prenatal/ExN.cov_0425.txt',
    '~/Dropbox/BrainVar2_WD/Data/prenatal/InN.cov_0425.txt',
    '~/Dropbox/BrainVar2_WD/Data/prenatal/Glia.cov_0425.txt',
    '~/Dropbox/BrainVar2_WD/Data/prenatal/NNC.cov_0425.txt'
  )
  
  example_all <- files %>%
    lapply(function(file) {
      read.delim(file) %>% select(id, batch, BioSample)
    }) %>%
    bind_rows() %>%
    distinct(id, .keep_all = TRUE)
}

example_all$id %>% duplicated() %>% table()

parse_score_string <- function(x) {
  as.numeric(str_split(str_remove_all(x, "\\[|\\]"), ",")[[1]])
}

for (ct in cell_types) {
  message("Processing: ", ct)
  
  pca_cov_path <- paste0('~/Dropbox/BrainVar2_WD/Tables/prenatal.PCA.', ct, '.20250725.txt')
  dm <- read.delim(pca_cov_path)
  
  genotype_pca_path <- if (ct == "NNC") {
    paste0('~/Dropbox/BrainVar2_WD/Data/prenatal_gPC/brainvar2.022924.joint_called.recalibrated.varQC.AF0.02.r20.1.prenatal.',
           ct, '.v2.PCAscores.txt')
  } else {
    paste0('~/Dropbox/BrainVar2_WD/Data/prenatal_gPC/brainvar2.022924.joint_called.recalibrated.varQC.AF0.02.r20.1.prenatal.',
           ct, '.PCAscores.txt')
  }
  genotype_pca <- read.delim(genotype_pca_path)
  
  genotype_pca_expanded <- genotype_pca %>%
    mutate(score_vec = map(scores, parse_score_string)) %>%
    unnest_wider(score_vec, names_sep = "_") %>%
    rename_with(~ paste0("gt_PC", seq_along(.)), starts_with("score_vec_")) %>%
    mutate(scores = str_remove_all(scores, "\\[|\\]")) %>%
    dplyr::select(-scores) %>%
    mutate(id = s) %>%
    dplyr::select(-s)
  
  dm2 <- dm %>%
    dplyr::select(
      id = matching_WGS_library,
      atac_PC1 = PC1, atac_PC2 = PC2, atac_PC3 = PC3, atac_PC4 = PC4, atac_PC5 = PC5,
      Sex = SexForAnalysis,
      # RIN,
      BioSample = BioSampleID,
      Period, AncestryPCA,
      PMI_hours
    )
  
  print(table(is.na(dm2$AncestryPCA)))
  print(table(dm2$AncestryPCA))
  
  # dm2 = dm2 %>%
  #   mutate(Asian = ifelse(AncestryPCA == 'Asian', 1, 0),
  #          Mix_Am = ifelse(AncestryPCA == 'Admixed American', 1, 0),
  #          AA = ifelse(AncestryPCA == 'African American', 1, 0),
  #          Hispanic = ifelse(AncestryPCA == 'Hispanic', 1, 0),
  #          Caucasian = ifelse(AncestryPCA == 'Caucasian', 1, 0))
  
  dm2 = dm2 %>%
    dplyr::select(-AncestryPCA)
  
  dup_ids <- dm2$id[duplicated(dm2$id)]
  
  if (length(dup_ids) > 0) {
    
    dm2_dup <- dm2 %>% filter(id %in% dup_ids)
    
    dm2_dup_filtered <- dm2_dup %>%
      filter(BioSample %in% example_all$BioSample)
    
    missing_ids <- setdiff(dup_ids, dm2_dup_filtered$id)
    if (length(missing_ids) > 0) {
      message("⚠️ No matching BioSample found for some duplicated IDs: ", paste(unique(missing_ids), collapse = ", "))
    }
    
    dm2_nondup <- dm2 %>% filter(!id %in% dup_ids)
    dm2 <- bind_rows(dm2_nondup, dm2_dup_filtered)
  }
  
  dm2 = dm2 %>%
    left_join(example_all %>% dplyr::select(-BioSample), by = 'id')
  
  dm2 <- dm2 %>%
    mutate(batch = ifelse(id %in% c("321_D2", "793_D2", "175_D1", "515_D2") & is.na(batch), 2, batch))
  
  dm3 <- dm2 %>%
    left_join(genotype_pca_expanded %>% dplyr::select(id, gt_PC1, gt_PC2, gt_PC3),
              by = "id")
  
  dup_ids <- dm3$id[duplicated(dm3$id)]
  if (length(dup_ids) > 0) {
    message("⚠️ Duplicate IDs found in dm3: ", paste(unique(dup_ids), collapse = ", "))
  } else {
    message("✅ No duplicate IDs in dm3.")
  }
  
  output_path <- paste0('~/Dropbox/BrainVar2_WD/Tables/table.prenatal.', ct, '.cov.', date, '.txt')
  write.table(dm3, output_path, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
}
