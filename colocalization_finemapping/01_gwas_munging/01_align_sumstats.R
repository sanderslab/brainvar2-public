# Munge GWAS summary statistics
library(BiocManager)
library(tidyverse)
library(MungeSumstats)

# NOTE: FREQ CANNOT BE FLIPPED FOR MULTIALLELICS, SO IF NOT SUBSETTING TO BIALLELIC, CANNOT TRUST FREQ INFO. 

# To create custom column name mapping file:
# data('sumstatsColHeaders')
# sumstatsColHeaders <- sumstatsColHeaders %>% add_row(Uncorrected = 'IMPINFO', Corrected = 'INFO') %>% as.data.frame

gwas_path <- '~/Desktop/Sanders/coloc_tinkering/data/GWAS'

# gwas_files <- list.files(paste0(gwas_path, 'raw'), full.names = TRUE)

for (f in gwas_files) {
  MungeSumstats::format_sumstats(path = f,
                                 ref_genome = 'GRCh37', # can leave as NULL for automatic detection
                                 # convert_ref_genome = 'GRCh38', # uncomment to convert
                                 bi_allelic_filter = FALSE,
                                 flip_frq_as_biallelic = TRUE,
                                 INFO_filter = 0, # keep all
                                 drop_indels = TRUE,
                                 # mapping_file = sumstatsColHeaders,
                                 save_path = file.path(gwas_path, 'processed', 'ASD', paste0('munged_noINFOFilter_GRCh37_', gsub('.gz', '.txt', basename(f)))),
                                 log_folder = file.path(gwas_path, 'processed', 'logs'),
                                 log_folder_ind = TRUE,
                                 log_mungesumstats_msg = TRUE,
                                 force_new = TRUE)
}
