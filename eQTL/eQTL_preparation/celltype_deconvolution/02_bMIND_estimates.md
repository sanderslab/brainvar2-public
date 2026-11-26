# Prepare
```
library(readr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(cowplot)

library(schard) # read in files for Seurat
library(Seurat)
```

# Marker Gene Selection 
The differences between prenatal and postnatal single-cell profiles are not ignorable. We are going to select the marker genes, estimate cell type proportions, and estimate cell-type-specific expression from bMIND separately. 
We use `markerpen` to select marker genes. Let's get prepared for marker gene selection: load the package, load and process single-cell data, and load bulk data.

## Load single-cell data
Take fetal data as an example
```
sc.fetal = NULL
counts.mat.fetal = NULL
meta.fetal = NULL
# 10% of the cells
for(p in 3:6){
  sc.temp = schard::h5ad2seurat(paste0('./single_cell/BrainVar_v1_anno_period_', p,'.h5ad'))
  sid = sort(sample(ncol(sc.temp), round(ncol(sc.temp)*0.1)))
  counts.mat.fetal = cbind(counts.mat.fetal, sc.temp@assays$RNA@counts[, sid])
  meta.fetal = rbind(meta.fetal, sc.temp@meta.data[sid, ])
}

sc.fetal = Seurat::CreateSeuratObject(counts = counts.mat.fetal, meta.data = meta.fetal)
selected_c <- WhichCells(sc.fetal, expression = nFeature_RNA > 200)
selected_f <- rownames(sc.fetal)[rowMeans(sc.fetal@assays$RNA$counts > 0) > 0.01]

sc.fetal <- subset(sc.fetal, features=selected_f, cells=selected_c)
sc.fetal <- NormalizeData(sc.fetal, normalization.method = "LogNormalize", scale.factor = 1e6)
VlnPlot(sc.fetal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sc.fetal <- FindVariableFeatures(sc.fetal, selection.method = "vst", nfeatures = 2000)

# filter genes that expressed in at least 1% of the cells
sc.fetal <- sc.fetal[rowMeans(sc.fetal@assays$RNA$counts > 0)> 0.01, ]
# 16174
all.genes <- rownames(sc.fetal)
sc.fetal <- ScaleData(sc.fetal, features = all.genes)
saveRDS(sc.fetal, file = 'data_store/sc_fetal_brain.RDS')
```
merge into 4 major cell types as we did for pseudobulk. 
```
list.broad.type = list(ExN = c('ExN-Newborn', 'ExN-NonIT', 'ExN-IT', 'IT-Immature', 'NonIT-Immature'), 
                       InN = c('InN-CGE', 'InN-CGE-Immature', 'InN-MGE', 'InN-MGE-Immature'), 
                       NNC = c('NNC'), 
                       Progenitor = c('Progenitor'),
                       Glia = c('Glia'))
sc.fetal$broad.type = sc.fetal$majorclass_v1
for(t in names(list.broad.type)){
  sc.fetal$broad.type[sc.fetal$broad.type %in% list.broad.type[[t]] ] = t
}
sc.fetal$broad.type = factor(sc.fetal$broad.type, levels = names(list.broad.type))
Idents(sc.fetal) = sc.fetal$broad.type

# Find a long list of marker genes using Seurat FindAllMarkers()
marker.mat.fetal = FindAllMarkers(sc.fetal, group.by = 'broad.type')
marker.fetal = marker.mat.fetal[marker.mat.fetal$p_val_adj < 0.05 & marker.mat.fetal$avg_log2FC > 1, ]
#        ExN    Progenitor        InN        NNC       Glia 
#        1021       1932          421       2312       1272 
```
## Prepare the package
```
# load markerpen package
install.packages("markerpen")
library(markerpen)
library(biomaRt)

# construct single cell data separately for prenatal and postnatal brains. 
# Discard period 7 (giving birth) 

vis_cor = function(mat_exp, markers){ # correlation visialization
  all_genes = colnames(mat_exp)
  markers = intersect(unlist(markers), all_genes)
  cor_markers = cor(mat_exp[, unlist(markers)])
  p = nrow(cor_markers)
  cols = c("#08306b", "#08519c", "#2171b5", "#6baed6", "#9ecae1", "#c6dbef", "#deebf7",
           "#ffffff",
           "#fcf1f1", "#fae1e1", "#facdcd", "#f49c9c", "#f56566", "#f13a3c", "#d00003")
  
  ncols = length(cols)
  cols = scales::gradient_n_pal(cols, values = (0:ncols) / ncols)((1:100) / 100)
  op = par(mar = c(0, 0, 0, 0))
  image(cor_markers[, p:1], col = cols, breaks = (-50:50) / 50, asp = 1, axes = FALSE)
  par(op)
}
```
## Use markerpen to refine the marker gene list with bulk data
```
load(file = 'bulk/bulk_all_2504.RData')
bulk_dev = bulk_all[, which(colnames(bulk_all) %in% meta.data.library.bulkRNA$LibraryID[meta.data.library.bulkRNA$Period >= 8])]
bulk_fetal = bulk_all[, which(colnames(bulk_all) %in% meta.data.library.bulkRNA$LibraryID[meta.data.library.bulkRNA$Period <= 6])]
gene.names = bulk_all$gene

# Start with fetal brain

bulk.counts = data.matrix(bulk_fetal)
# filter out genes that are expressed in at least 5% of the samples

bulk.counts  = bulk.counts[rowMeans(bulk.counts > 0) > 0.05, ]
# 47559

marker.fetal$gene
# --------------------------------------
# Filter out markers by markerpen
# --------------------------------------
mat_exp = t(bulk.counts)
mat_exp = mat_exp / rowSums(mat_exp) * 10^6
marker.fetal = marker.fetal[marker.fetal$gene %in% colnames(mat_exp), ]

nb.prior = 200
nb.range = 1000
celltype.names = names(list.broad.type)
markers.re = list()
mat_rest = mat_exp[, marker.fetal$gene[marker.fetal$gene %in% colnames(mat_exp)]]
for (ct in celltype.names){
  message(ct)
  markers_ct = marker.fetal[marker.fetal$cluster == ct, 'gene']
  markers_ct = markers_ct[markers_ct %in% colnames(mat_rest)]
  markers_ct = markers_ct[!is.na(apply(mat_rest[,markers_ct], 2, sd)) & apply(mat_rest[,markers_ct], 2, sd)>0]
  
  markers_range = markers_ct[1:min(nb.range, nrow(markers_ct))]  # top 500 as range
  markers_prior = markers_range[1:min(nb.prior, nrow(markers_ct))]  # top 500 as range
  
  re = markerpen::refine_markers(mat_rest, markers_range, markers_prior, 
                                 lambda = 0.4, w = 1.5, alpha=0.01, maxit = 500, verbose = 0)
  mat_rest = mat_rest[, setdiff(colnames(mat_rest), re$markers)]
  markers.re[[ct]] = re$markers
}
str(markers.re)

markers.df = data.frame(gene = unlist(markers.re),
                        celltype = rep(names(markers.re), sapply(markers.re, length)),
                        cluster = rep(names(markers.re), sapply(markers.re, length)))
save(markers.df, file=paste0('bulk/markers_bulk_fetal_0425.rds'))
```
visualization
```
cor_markers = cor(mat_exp[, unlist(markers.re)])
markers_ord = sort_markers(cor_markers, markers.re)
markers_ord = lapply(markers_ord, head, n = 50)
vis_cor(mat_exp = mat_exp, markers = markers_ord)
```
![plot](./figures/fetal_markerpen_cor.png)

# Cell type proportion deconvolution
We choose to use `Bisque` for deconvolution and fetal brain samples as an example. 
```
library(BisqueRNA)
library(Biobase)

# we only care about the expression of marker genes
load('bulk/markers_bulk_fetal_0425.rds')
table(markers.df$celltype)
genes = markers.df$gene
str(genes) # 846

sc.counts.tmp = sc.fetal@assays$RNA$counts[genes,]
bulk.counts.tmp = bulk.counts[genes, ]

# Construct eset object for both bulk and single cell data
bulk.eset = Biobase::ExpressionSet(assayData = bulk.counts.tmp)
sample.ids <- colnames(sc.counts.tmp)
sc.fetal$SubjectID = meta.data.library$SubjectID[match(sc.fetal$library, meta.data.library$BioSampleID)]
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names = F, check.rows = F,
                       stringsAsFactors = F,
                       row.names = sample.ids,
                       SubjectName = sc.fetal$SubjectID,
                       cellType = sc.fetal$broad.type)
sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))
sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData = data.matrix(sc.counts.tmp),
                                  phenoData=sc.pdata)
markers.fetal = unique(c(markers.df))
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers= markers.df, use.overlap = FALSE)
save(res, file = 'bulk_pi_bisque_fetal_0425.RData')

visualization
m.bulk.props = melt(res$bulk.props)
colnames(m.bulk.props) = c('Celltype', 'SampleID', 'Proportion')
ggplot(m.bulk.props, aes(Celltype, Proportion)) + geom_boxplot() + theme_minimal() + ggtitle('Period 3-6')
```
![decov](./figures/fetal_deconv_prop.png)

# bMIND estimation for cell-type-specific expression
## Pre-processing by matching the sample names
```
load('bulk/bulk_all_2504.RData') # bulk_all

meta.data.library = read_tsv('brainvar2_library_metadata_v1.2.tsv') # 842 library
meta.data.library = meta.data.library[meta.data.library$Usable == 'Yes', ] #770 library
meta.data.library.bulkRNA = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'bulkRNA', ] #230
meta.data.library.WGS = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'WGS', ] #237

load('bulk_pi_bisque_fetal_0425.RData')
pi_fetal = t(res$bulk.props)

# The covariates
covariate = df.PC.new[, c('Sex', 'RIN', 'Period', 'batch', 'Asian', 'Mix_Am', 'AA', 'Hispanic', 'Caucasian')]
# change the covariates rownames
rownames(covariate) = meta.data.library.bulkRNA$LibraryID[match(meta.data.library.WGS$SubjectID[match(rownames(covariate), meta.data.library.WGS$LibraryID)], meta.data.library.bulkRNA$SubjectID)]
```

## bMIMD for fetal brain
```
frac = pi_fetal
# remove the outlier samples from bulk eqtl analysis
bulk.id = setdiff(intersect(rownames(frac), rownames(covariate)), c("194_R1", "470_R2", "919_R2", '416_R1'))
frac = frac[bulk.id, ]
bulk = data.matrix(bulk_all[, bulk.id])
rownames(bulk) = bulk_all$gene
cov.sub = data.frame(covariate[bulk.id, ])
# save(bulk, frac, cov.sub, file = 'bmind/bmind_prepare_fetal_0425.RData')

bmind.fetal = bMIND(bulk, frac = frac, covariate = cov.sub)
saveRDS(bmind.fetal, file = 'bmind/bmind_fetal_0425.rds')
```

## bMIND estimates processing for eQTL
Similar to pseudobulk process.
```
sum_expression_fetal <- melt(A.fetal) %>% rename(ensembl = Var1, subject = Var3, cell.type = Var2, counts = value) %>%
  group_by(subject) %>%
  mutate(counts_scaled_1M = counts*10^6/sum(counts)) %>% 
  ungroup() %>% 
  as_tibble()

sumstats_gene <- sum_expression_fetal %>% group_by(cell.type, ensembl) %>% 
  summarise(counts_per_celltype = sum(counts),
            n_samples_expressed = sum(counts>0),
            perc_samples_expressed = round(n_samples_expressed*100/length(counts),digits=2),
            mean_cpm = mean(counts_scaled_1M)) %>% 
  group_by(ensembl) %>% 
  mutate(prop_max_cpm = mean_cpm/max(mean_cpm)) %>% 
  ungroup()

gene_to_keep_per_cell_type <- sumstats_gene %>% 
  filter(n_samples_expressed > 12, mean_cpm > 1) %>% 
  dplyr::select(cell.type, ensembl) %>% 
  unique()

sum_expression_fetal <- sum_expression_fetal %>% 
  inner_join(.,gene_to_keep_per_cell_type,by=c('cell.type','ensembl'))

(dplyr::count(sum_expression_fetal, cell.type, subject) %>% 
    dplyr::select(-subject) %>% 
    unique() %>% 
    arrange(-n))
# cell.type     n
# 1 InN        16308
# 2 Glia       16229
# 3 NNC        16197
# 4 ExN        16021
# 5 Progenitor 15910


sum_expression_fetal$sample_id = meta.data.library$SubjectID[match(sum_expression_fetal$subject, meta.data.library$LibraryID)]
sum_expression_fetal$DNA_library = meta.data.library.WGS$LibraryID[match(sum_expression_fetal$sample_id, meta.data.library.WGS$SubjectID)]
```

