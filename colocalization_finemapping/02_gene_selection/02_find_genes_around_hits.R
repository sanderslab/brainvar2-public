library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

# Your function
find_genes_around_peaks <- function(sumstats, gtf = 'data/gencode.v47.basic.annotation.gtf.gz', 
                                    range = 2000000, out_intervals, out_genes) {
  sumstats <- read_tsv(sumstats)
  
  sumstats <- sumstats %>%
    dplyr::filter(P < 5e-8) %>%
    distinct(CHR, BP, SNP) %>%
    mutate(
      interval_start = pmax(BP - (range / 2), 0),
      interval_end = BP + (range / 2)
    )
  
  gtf <- rtracklayer::import(gtf) %>%
    as_tibble() %>%
    dplyr::filter(type == 'gene', gene_type %in% c('lncRNA', 'protein_coding')) %>%
    dplyr::select(chr = seqnames, start, end, gene_id, gene_type, gene_name) %>%
    mutate(
      chr = gsub('chr', '', as.character(chr)),
      gene_id = gsub('\\..*', '', gene_id)
    ) %>%
    dplyr::filter(chr %in% as.character(1:22)) %>%
    mutate(chr = as.numeric(chr))
  
  hits_annot <- sumstats %>%
    rowwise() %>%
    mutate(
      genes_in_interval = list(
        gtf %>%
          dplyr::filter(chr == CHR, start >= interval_start, end <= interval_end) %>%
          dplyr::select(gene_name, gene_id) %>%
          distinct()
      )
    ) %>%
    unnest_longer(genes_in_interval) %>%
    group_by(CHR, interval_start, interval_end) %>%
    summarize(
      interval_gene_names = paste(unique(genes_in_interval$gene_name), collapse = ','),
      interval_gene_ids = paste(unique(genes_in_interval$gene_id), collapse = ','),
      .groups = 'drop'
    )
  
  write_tsv(hits_annot, out_intervals)
  
  all_genes <- unique(unlist(strsplit(paste(hits_annot$interval_gene_ids, collapse = ','), ',')))
  fwrite(list(all_genes), out_genes)
}

# --- Automation for all GWAS files in a directory ---
gwas_dir <- "/mnt/s3/GWAS_data/"  # directory containing GWAS files
gtf_file <- "/mnt/s3/reference_genome/gencode.v47lift37.basic.annotation.gtf.gz"
out_dir <- "/mnt/s3/eqtl/eQTL/GWAS_subset/"
# List all GWAS .txt files
gwas_files <- list.files(gwas_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
print('list of gwas')
print(gwas_files)
# Run the function for each file
walk(gwas_files, function(file) {
  file_base <- tools::file_path_sans_ext(basename(file))
  
  out_intervals <- file.path(out_dir, paste0(file_base, ".2Mb_intervals.txt"))
  out_genes <- file.path(out_dir, paste0("all_", file_base, "_genes_in_2Mb_intervals.txt"))
  tryCatch(
  {
  find_genes_around_peaks(
    sumstats = file,
    gtf = gtf_file,
    out_intervals = out_intervals,
    out_genes = out_genes
  )    },
  error = function(e) {
  message("ERROR processing the file")
  })
})
