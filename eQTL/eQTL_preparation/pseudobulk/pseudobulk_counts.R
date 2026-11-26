# get pseudobulk data with# Preprocess: filter genes and samples
#### single-cell data read-in packages
library(anndata)
library(reticulate)
library(schard) 
library(reshape2)

## Raw Counts =  summation of counts from the same (cell.type * sample)

########## Prenatal #################
list.broad.type = list(ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
                       InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
                       Glia = c('Glia'), 
                       NNC = c('NNC'), 
                       Progenitor = c('Progenitor'))

cell.counts.mat = NULL
meta.pb = NULL
pb.counts = NULL

for(p in 1:6){
  sc.temp = schard::h5ad2seurat(paste0('./single_cell/BrainVar_v1_anno_period_', p,'.h5ad')) # your path to the single-cell data.
  sc.temp$broad.type = sc.temp$majorclass_v1
  for(t in names(list.broad.type)){
    sc.temp$broad.type[sc.temp$broad.type %in% list.broad.type[[t]] ] = t
  }
  sc.temp = sc.temp[, sc.temp$percent_mt <= 15 & sc.temp$percent_ribo >= 0.25]       # filter mt and ribo
  
  cell.tab = table(sc.temp$orig.ident, sc.temp$broad.type)
  cell.counts.mat = rbind(cell.counts.mat, cell.tab)
  
  m.cell.tab = melt(cell.tab)
  colnames(m.cell.tab) = c('subject', 'cell.type', 'cell.number')
  
  pb.names = NULL
  pb.counts.temp = NULL
  for(i in 1:nrow(m.cell.tab)){
    if(m.cell.tab$cell.number[i] > 10){
      pb.counts.temp = cbind(pb.counts.temp, rowSums( sc.temp@assays$RNA$counts[, sc.temp$orig.ident == m.cell.tab$subject[i] & 
                                                                                  sc.temp$broad.type == m.cell.tab$cell.type[i]] ))
      pb.names = c(pb.names, paste0(m.cell.tab$subject[i], '_', m.cell.tab$cell.type[i]))
    }
  }
  colnames(pb.counts.temp) = pb.names
  meta.pb = rbind(meta.pb, m.cell.tab[m.cell.tab$cell.number > 10, ])               # Filter pb.samples with less than 10 cells.
  #save(pb.counts.temp, file = paste0('pseudobulk/pseudobulk_counts_fetal_period_', p,'.rds'))
  #write.table(m.cell.tab[m.cell.tab$cell.number > 10, ], file = paste0('meta_pseudobulk_fetal_period_', p,'.txt'), row.names = FALSE, quote = FALSE)
  
  pb.counts = cbind(pb.counts, pb.counts.temp)
  cat(p, '\n')
}
save(meta.pb, pb.counts, file = 'pseudobulk_fetal.RData')

##### Postnatal ############
list.broad.type = list(ExN = c('ExN-IT', 'ExN-NonIT'), 
                       InN = c('InN-CGE', 'InN-MGE'), 
                       NNC = c('NNC'),
                       Glia = c('Glia'))
cell.counts.mat = NULL
meta.pb = NULL
pb.counts = NULL

for(p in 8:15){
  sc.temp = schard::h5ad2seurat(paste0('./single_cell/BrainVar_v1_anno_period_', p,'.h5ad')) # your path to the single cell data
  sc.temp$broad.type = sc.temp$majorclass_v1
  for(t in names(list.broad.type)){
    sc.temp$broad.type[sc.temp$broad.type %in% list.broad.type[[t]] ] = t
  }
  sc.temp = sc.temp[, sc.temp$percent_mt <= 15 & sc.temp$percent_ribo >= 0.25]        # filter mt and ribo
  
  cell.tab = table(sc.temp$orig.ident, sc.temp$broad.type)
  cell.counts.mat = rbind(cell.counts.mat, cell.tab)
  
  m.cell.tab = melt(cell.tab)
  colnames(m.cell.tab) = c('subject', 'cell.type', 'cell.number')
  
  pb.names = NULL
  pb.counts.temp = NULL
  for(i in 1:nrow(m.cell.tab)){
    if(m.cell.tab$cell.number[i] > 10){
      pb.counts.temp = cbind(pb.counts.temp, rowSums( sc.temp@assays$RNA$counts[, sc.temp$orig.ident == m.cell.tab$subject[i] & 
                                                                                  sc.temp$broad.type == m.cell.tab$cell.type[i]] ))
      pb.names = c(pb.names, paste0(m.cell.tab$subject[i], '_', m.cell.tab$cell.type[i]))
    }
  }
  colnames(pb.counts.temp) = pb.names
  meta.pb = rbind(meta.pb, m.cell.tab[m.cell.tab$cell.number > 10, ])                         # Filter pb.samples with less than 10 cells.
  #save(pb.counts.temp, file = paste0('pseudobulk/pseudobulk_counts_dev_period_', p,'.rds'))
  #write.table(m.cell.tab[m.cell.tab$cell.number > 10, ], file = paste0('meta_pseudobulk_dev_period_', p,'.txt'), row.names = FALSE, quote = FALSE)
  
  pb.counts = cbind(pb.counts, pb.counts.temp)
  cat(p, '\n')
}
save(meta.pb, pb.counts, file = 'pseudobulk_dev.RData')
