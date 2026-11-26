library(tidyverse)
library(data.table)

celltypes <- c("NNC")

meta <- read.delim('~/Dropbox/BrainVar2_WD/Data/brainvar2.prenatal.metadata.20250620.txt')

out_dir <- "~/Dropbox/BrainVar2_WD/Tables/caQTL_input/Prenatal/"

exclude_samples <- list(
  NNC = c('H292DFCm2', 'H430DFCm', 'H292DFCm1')
)

for (ct in celltypes) {
  message("Processing: ", ct)
  
  qnorm_path <- paste0("~/Dropbox/BrainVar2_WD/Tables/prenatal.", ct, "_atac.lcpm.qnorm")
  
  if (!file.exists(qnorm_path)) {
    warning("File not found: ", qnorm_path)
    next
  }
  
  df <- fread(qnorm_path)
  
  excluded <- exclude_samples[[ct]]
  samples_filtered <- setdiff(colnames(df)[-1], excluded)
  df_filtered <- df[, c("V1", samples_filtered), with = FALSE]
  
  meta2 <- meta %>%
    filter(BioSampleID %in% samples_filtered) %>%
    dplyr::select(matching_WGS_library) %>%
    distinct() %>%
    mutate(CellType = ct,
           Stage = "Prenatal")
  
  fwrite(
    meta2,
    file = file.path(out_dir, paste0("prenatal.", ct, ".sample_metadata.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  message("✅ Saved ", nrow(meta2), " samples for ", ct)
}








