## TMM normalization and get prepared for eQTL analysis
library(tidyverse)
library(parallel)
library(reshape2)
library(anndata)
library(reticulate)
library(schard)

list.broad.type = list(ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
                       InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
                       NNC = c('NNC'), 
                       Progenitor = c('Progenitor'),
                       Glia = c('Glia'))

list.broad.type = list(ExN = c('ExN-IT', 'ExN-NonIT'), 
                       InN = c('InN-CGE', 'InN-MGE'), 
                       NNC = c('NNC'),
                       Glia = c('Glia'))

sum_expression = NULL
sum_expression_sample = NULL

for(p in 3:15){
  sc.temp = schard::h5ad2seurat(paste0('./single_cell/BrainVar_v1_anno_period_', p,'.h5ad'))
  sc.temp$broad.type = sc.temp$majorclass_v1
  sc.temp$broad.type[sc.temp$majorclass_v1 %in% list.broad.type$ExN] = 'ExN'
  sc.temp$broad.type[sc.temp$majorclass_v1 %in% list.broad.type$InN] = 'InN'
  sc.temp$subjectID = meta.data.library.GEX$SubjectID[match(sc.temp$orig.ident, meta.data.library.GEX$realname.sc)]
  sc.temp = sc.temp[, sc.temp$percent_mt <= 15 & sc.temp$percent_ribo >= 0.25]

  # ExN, Glia, InN, NNC
  cell_type_id = melt(table(sc.temp$orig.ident, sc.temp$broad.type))
  cell_type_id = cell_type_id[cell_type_id$Var2 %in% c('ExN', 'Glia', 'InN', 'NNC'), ]
  cell_type_id$Var2 = factor(cell_type_id$Var2, levels = c('ExN', 'Glia', 'InN', 'NNC'))
  colnames(cell_type_id) = c('subject', 'cell.type', 'num.cells')
  
  # cell_type_id$subjectID = meta.data.library.GEX$SubjectID[match(cell_type_id$subject, meta.data.library.GEX$realname.sc)]
  expression_temp = NULL
  for(i in 1:nrow(cell_type_id)){
    cnt.mtx = sc.temp@assays$RNA@counts[, sc.temp$broad.type == cell_type_id$cell.type[i] & sc.temp$orig.ident == cell_type_id$subject[i]]
    if(is.null(dim(cnt.mtx))){
      tot_counts = cnt.mtx
      n_expressed = as.numeric(cnt.mtx > 0)
      n_cells = 1
    }else{
      tot_counts <- Matrix::rowSums(cnt.mtx)
      n_expressed <- Matrix::rowSums(cnt.mtx > 0)
      n_cells <- ncol(cnt.mtx)
    }
    cnt.df = data.frame(counts = tot_counts, n_expressed = n_expressed, n_cells = n_cells)  %>%
      rownames_to_column('hgnc') %>%
      mutate(perc.exprs = round(n_expressed*100/n_cells, digits = 2), 
             cell.type = cell_type_id$cell.type[i], 
             subject = cell_type_id$subject[i]) %>%
      dplyr::select(cell.type, subject, hgnc, counts, n_expressed, n_cells, perc.exprs)
    cnt.df$period = p 
    cnt.df$subjectID = meta.data.library.GEX$SubjectID[match(cnt.df$subject, meta.data.library.GEX$realname.sc)]
    expression_temp = rbind(expression_temp, cnt.df)
  }
  
  temp <- expression_temp %>% group_by(subjectID, hgnc, cell.type) %>% 
    dplyr::summarize(counts = mean(counts, na.rm=TRUE), 
                     n_expressed = mean(n_expressed, na.rm = TRUE), 
                     n_cells = mean(n_cells, na.rm = TRUE))
  
  sum_expression = rbind(sum_expression, temp)
  
  cat('period ', p, '\n')
}

# 137 sample 4 cell types
sum_expression <- sum_expression %>% 
  group_by(cell.type, subject) %>% 
  mutate(counts_scaled_1M = counts*10^6/sum(counts)) %>% 
  ungroup() %>% 
  as_tibble()

sumstats_gene <- sum_expression %>% group_by(cell.type, hgnc) %>% 
  summarise(counts_per_celltype = sum(counts),
            n_expressed_per_celltype = sum(n_expressed),
            n_cells_per_celltype = sum(n_cells),
            perc_expressed_per_celltype = round(n_expressed_per_celltype*100/n_cells_per_celltype, digits=2),
            n_samples_expressed = sum(counts>0),
            perc_samples_expressed = round(n_samples_expressed*100/length(counts),digits=2),
            mean_cpm = mean(counts_scaled_1M)) %>% 
  group_by(hgnc) %>% 
  mutate(prop_max_cpm = mean_cpm/max(mean_cpm)) %>% 
  ungroup()
# We do not have ambient records

# In total we have 137 subjects
# 1. Expressed in at least 10 individuals
# 2. Expressed with a mean CPM of at least 1
gene_to_keep_per_cell_type <- sumstats_gene %>% 
  filter(n_samples_expressed > 9, mean_cpm > 1) %>% 
  dplyr::select(cell.type, hgnc) %>% 
  unique()

sum_expression <- sum_expression %>% 
  inner_join(.,gene_to_keep_per_cell_type,by=c('cell.type','hgnc'))

(dplyr::count(sum_expression,cell.type,subject) %>% 
    dplyr::select(-subject) %>% 
    unique() %>% 
    arrange(-n))
# A tibble: 4 × 2
# cell.type     n
# <fct>     <int>
# 1 NNC       19766
# 2 Glia      19560
# 3 ExN       19003
# 4 InN       18987

(n_cells_ind_cell_type <- sum_expression %>% 
    dplyr::select(cell.type, subject, n_cells) %>% 
    unique() %>% arrange(n_cells))
# Remove individuals with less than 10 cells for a given cell type.
keep <- filter(n_cells_ind_cell_type,n_cells>=10)
dplyr::count(keep, cell.type) %>% arrange(-n)
# A tibble: 4 × 2
#cell.type     n
#<fct>     <int>
#1 Glia        137
#2 InN         137
#3 ExN         137
#4 NNC         134

sum_expression <- inner_join(sum_expression, keep, by=c('cell.type','subjectID','avg.n_cells'))

meta.data.library = read_tsv('brainvar2_library_metadata_v1.2_XW.tsv')

sum_expression$sample_id = meta.data.library$SubjectID[match(sum_expression$subject, meta.data.library$realname.sc)]
sum_expression$DNA_library = meta.data.library.WGS$LibraryID[match(sum_expression$sample_id, meta.data.library.WGS$SubjectID)]


gtf = readRDS('../Panos_project/new_data_store/gtf.RDS')
# Write pseudobulk expressions
get_fastqtl_pheno <- function(cell){
  #keep all expression for the cell type and patient that have genotype
  bed <- sum_expression %>% dplyr::filter(cell.type == cell) %>% 
    dplyr::select(subject, hgnc, counts) %>% 
    spread(subject, counts) %>% 
    column_to_rownames('hgnc') %>% 
    edgeR::DGEList(counts = .) %>% #TMM normalize the data using edgeR
    edgeR::calcNormFactors(.) %>% 
    edgeR::cpm(.) %>%
    as.data.frame() %>% 
    rownames_to_column('hgnc')
  
  bed <- bed %>% inner_join(.,gtf, by = 'hgnc') %>% 
    dplyr::select(seqnames,TSS_start,TSS_end,gene,everything()) %>% 
    arrange(seqnames,TSS_start)
  subject.names <- colnames(bed)[! colnames(bed)%in%c("seqnames", "TSS_start", "TSS_end", "gene", "hgnc", "ensembl")]
  cat(as.character(cell), ': ', length(subject.names), '\n')
  # subject.names.new = meta.data.library.WGS$LibraryID[match(meta.data.library$SubjectID[ match(subject.names, meta.data.library$realname.sc) ], meta.data.library.WGS$SubjectID )]
  
  colnames(bed)[c(1,2,3,4)] <- c('#Chr','start','end','ID')
  bed <- bed[, c('#Chr','start','end','ID', subject.names)]
  colnames(bed)[-c(1, 2, 3, 4)] = subject.names

  write_tsv(bed, path = paste0('./pseudobulk/',make.names(cell),'.bed'))
  # There is no need to separate the chromosomes for pseudobulk samples
}

cells <- unique(sum_expression$cell.type)
TMM.pseudo.bulk = lapply(cells, get_fastqtl_pheno)
names(TMM.pseudo.bulk) = cells
saveRDS(TMM.pseudo.bulk, file = './pseudobulk/pseudobulk_tmm.RDS')
