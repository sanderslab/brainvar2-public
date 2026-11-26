# Prepare
```
library(readr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(cowplot)
load('data_store/df_PC.RData')
```

## Genotype PC filters
```
df.PC.new = df.PC[df.PC$remain == 'yes', ]
df.PC.new$Sex = as.numeric(df.PC.new$Sex == 'XY') + 1
df.PC.new$Asian = as.numeric(df.PC.new$Ancestry == 'Asian')
df.PC.new$Mix_Am = as.numeric(df.PC.new$Ancestry == 'Admixed American')
df.PC.new$AA = as.numeric(df.PC.new$Ancestry == 'African American')
df.PC.new$Hispanic = as.numeric(df.PC.new$Ancestry == 'Hispanic')
df.PC.new$Caucasian = as.numeric(df.PC.new$Ancestry == 'Caucasian')
colnames(df.PC.new)[1:3] = paste0('gt_PC', 1:3)

# 222 samples "175_R1" "515_R2" has bulkRNA but no DNA
common.bulk.subject = intersect(rownames(df.PC.new), meta.data.library.bulkRNA$matchingDNALibrary)
df.PC.new = df.PC.new[common.bulk.subject, ]

df.PC.new$pH = meta.data.library.bulkRNA$pH[match(common.bulk.subject, meta.data.library.bulkRNA$matchingDNALibrary)]
df.PC.new$PMI_hours = meta.data.library.bulkRNA$PMI_hours[match(common.bulk.subject, meta.data.library.bulkRNA$matchingDNALibrary)]
df.PC.new$RIN = meta.data.library.bulkRNA$RIN[match(common.bulk.subject, meta.data.library.bulkRNA$matchingDNALibrary)]
df.PC.new$XIST_depth_RNA = meta.data.library.bulkRNA$XIST_depth_RNA[match(common.bulk.subject, meta.data.library.bulkRNA$matchingDNALibrary)]
# Add batch
df.PC.new$batch = as.numeric( RIGHT(meta.data.library.bulkRNA$LibraryID, 1)[match(common.bulk.subject, meta.data.library.bulkRNA$matchingDNALibrary)] )

df.PC.new = df.PC.new[, c(paste0('gt_PC', 1:3), 'Sex', 'RIN', 'Age', 'Period', 'batch', 'Asian', 'Mix_Am', 'AA', 'Hispanic', 
                          'Caucasian', 'PMI_hours', 'pH', 'XIST_depth_RNA')]
```

## Gene expression (TMM normalized)
```
d = read_tsv('bulk/bulk_new_2504.bed')
# We remove poor QC sample: 416_D1
df.PC.new = df.PC.new[rownames(df.PC.new) != '416_D1', ] 
d = d[, colnames(d)!='416_D1']
common.bulk.subject = common.bulk.subject[common.bulk.subject!='416_D1']
pca <- d[-c(1,2,3,4)] %>% t() %>% prcomp(.,scale.=T)

# Elbow Plot from PCs
plot(pca$sdev[1:60], main = 'bulk all') # 10
```
![plot](../../bulk_RNAseq/quality_control/figures/pca_elbow_bulk_all_first_0425.jpg)

```
df.pc.all = as.data.frame(pca$x[, 1:5])
colnames(df.pc.all) = paste0('exprs_PC', 1:5)

df.pc.all = df.pc.all[common.bulk.subject, ]
df.pc.all = cbind(df.pc.all, df.PC.new[rownames(df.pc.all), ])

bulk_all_cor_matrix <- cor(df.pc.all)
bulk_all_cor_matrix['RIN', ] = cor(df.pc.all$RIN[!is.na(df.pc.all$RIN)], df.pc.all[!is.na(df.pc.all$RIN), ])
bulk_all_cor_matrix[, 'RIN'] = cor(df.pc.all$RIN[!is.na(df.pc.all$RIN)], df.pc.all[!is.na(df.pc.all$RIN), ])
bulk_all_cor_matrix['PMI_hours', ] = cor(df.pc.all$PMI_hours[!is.na(df.pc.all$PMI_hours)], df.pc.all[!is.na(df.pc.all$PMI_hours), ])
bulk_all_cor_matrix[, 'PMI_hours'] = cor(df.pc.all$PMI_hours[!is.na(df.pc.all$PMI_hours)], df.pc.all[!is.na(df.pc.all$PMI_hours), ])
bulk_all_cor_matrix['pH', ] = cor(df.pc.all$pH[!is.na(df.pc.all$pH)], df.pc.all[!is.na(df.pc.all$pH), ])
bulk_all_cor_matrix[, 'pH'] = cor(df.pc.all$pH[!is.na(df.pc.all$pH)], df.pc.all[!is.na(df.pc.all$pH), ])


# corrplot(bulk_all_cor_matrix)
# round(bulk_all_cor_matrix, 3)
na.index = which(is.na(bulk_all_cor_matrix), arr.ind = T)
for(i in 1:nrow(na.index)){
  bulk_all_cor_matrix[na.index[i, 1], na.index[i, 2]] <- 0
}

a = heatmap(abs(bulk_all_cor_matrix))
cov_order =rownames(bulk_all_cor_matrix)[a$rowInd]
bulk_all_cor_matrix = bulk_all_cor_matrix[cov_order, cov_order]

#pdf('heatmap_exprs_PCA_bulk_all_corr_0425.pdf', width = 8, height = 9)
corrplot(bulk_all_cor_matrix, title = 'bulk all')
#dev.off()
```
![plot2](../../bulk_RNAseq/quality_control/figures/heatmap_exprs_PCA_bulk_all_corr_0425.pdf)

```
# Significant ones, what we should watching here

pdf('exprs_PCA_bulk_all_0425.pdf', width = 10, height = 8)
gg1 = ggplot(df.pc.all, aes(exprs_PC1, exprs_PC2, col = as.factor(batch))) + geom_point() + theme_minimal() + 
  ggtitle('Batch') + scale_color_brewer(palette = 'Set1')
gg2 = ggplot(df.pc.all, aes(exprs_PC1, exprs_PC2, col = Period)) + geom_point() + theme_minimal() + 
  ggtitle('Period') + scale_color_distiller(palette = 'Spectral')
gg3 = ggplot(df.pc.all, aes(exprs_PC2, exprs_PC5, col = as.factor(batch))) + geom_point() + theme_minimal() + 
  ggtitle('Batch') + scale_color_brewer(palette = 'Set1')
gg4 = ggplot(df.pc.all, aes(exprs_PC2, exprs_PC3, col = RIN)) + geom_point() + theme_minimal() + 
  ggtitle('RIN') + scale_color_distiller(palette = 'Spectral')
plot_grid(gg1, gg2, gg3, gg4)
dev.off()
```

## Brainvar 2 results
```
d = read_tsv('bulk/bulk_new_2504.bed')
temp = d[rownames(df.PC.new)[df.PC.new$batch == 2]] %>% t()
temp = temp[, colSums(temp) > 0]
pca <- temp %>% prcomp(.,scale.=T)

df.pc.all.bv2 = as.data.frame(pca$x[, 1:5])
df.pc.all.bv2$bv2 = (nchar(rownames(temp)) == 7) + 1
df.pc.all.bv2 = cbind(df.pc.all.bv2, df.PC.new[rownames(df.pc.all.bv2), ])

pdf('exprs_PCA_bulk_brainvar2_0425.pdf', width = 8, height = 8)
gg1 = ggplot(df.pc.all.bv2, aes(PC1, PC2, col = as.factor(bv2), label = rownames(df.pc.all.bv2) ) ) + geom_point() + 
  theme_minimal() + scale_color_brewer(palette = 'Set2')  + ggrepel::geom_text_repel() +theme(legend.position="none")
gg2 = ggplot(df.pc.all.bv2, aes(PC2, PC3, col = as.factor(bv2), label = rownames(df.pc.all.bv2) )) + geom_point() + 
  theme_minimal() + scale_color_brewer(palette = 'Set2') + theme(legend.position="none") + ggrepel::geom_text_repel()
gg3 = ggplot(df.pc.all.bv2, aes(PC1, PC2, col = Period, label = rownames(df.pc.all.bv2) )) + geom_point() + 
  theme_minimal() + scale_color_viridis_c()   + ggrepel::geom_text_repel()
gg4 = ggplot(df.pc.all.bv2, aes(PC2, PC3, col = Period, label = rownames(df.pc.all.bv2) )) + geom_point() + 
  theme_minimal() + scale_color_viridis_c()   + ggrepel::geom_text_repel()
plot_grid(gg1, gg2, gg3, gg4)
dev.off()
# 583_D1 and 919_D2 they have some issue. 
temp = d[rownames(df.PC.new)[df.PC.new$batch == 1]] %>% t()
temp = temp[, colSums(temp) > 0]
pca <- temp %>% prcomp(.,scale.=T)

df.pc.all.bv1 = as.data.frame(pca$x[, 1:5])
df.pc.all.bv1 = cbind(df.pc.all.bv1, df.PC.new[rownames(df.pc.all.bv1), ])
pdf('exprs_PCA_bulk_brainvar1_0425.pdf', width = 5.5, height = 10)
gg1 = ggplot(df.pc.all.bv1, aes(PC1, PC2, label = rownames(df.pc.all.bv1), col = Period) ) + geom_point() + 
  theme_minimal() + scale_color_viridis_c()  + ggrepel::geom_text_repel()
gg2 = ggplot(df.pc.all.bv1, aes(PC2, PC3, col = Period, label = rownames(df.pc.all.bv1) )) + geom_point() + 
  theme_minimal() + scale_color_viridis_c()   + ggrepel::geom_text_repel()
plot_grid(gg1, gg2, nrow = 2)
dev.off()
# 194_D1 can be removed

bulk_all_bv1_cor_matrix <- cor(df.pc.all.bv1)
bulk_all_bv1_cor_matrix = bulk_all_bv1_cor_matrix[-13, -13]
bulk_all_bv1_cor_matrix['RIN', ] = cor(df.pc.all.bv1$RIN[!is.na(df.pc.all.bv1$RIN)], df.pc.all.bv1[!is.na(df.pc.all.bv1$RIN), ])
bulk_all_bv1_cor_matrix[, 'RIN'] = cor(df.pc.all.bv1$RIN[!is.na(df.pc.all.bv1$RIN)], df.pc.all.bv1[!is.na(df.pc.all.bv1$RIN), ])
bulk_all_bv1_cor_matrix['PMI_hours', ] = cor(df.pc.all.bv1$PMI_hours[!is.na(df.pc.all.bv1$PMI_hours)], df.pc.all.bv1[!is.na(df.pc.all.bv1$PMI_hours), ])
bulk_all_bv1_cor_matrix[, 'PMI_hours'] = cor(df.pc.all.bv1$PMI_hours[!is.na(df.pc.all.bv1$PMI_hours)], df.pc.all.bv1[!is.na(df.pc.all.bv1$PMI_hours), ])
bulk_all_bv1_cor_matrix['pH', ] = cor(df.pc.all.bv1$pH[!is.na(df.pc.all.bv1$pH)], df.pc.all.bv1[!is.na(df.pc.all.bv1$pH), ])
bulk_all_bv1_cor_matrix[, 'pH'] = cor(df.pc.all.bv1$pH[!is.na(df.pc.all.bv1$pH)], df.pc.all.bv1[!is.na(df.pc.all.bv1$pH), ])


corrplot(bulk_all_bv1_cor_matrix)
round(bulk_all_bv1_cor_matrix, 3)
na.index = which(is.na(bulk_all_bv1_cor_matrix), arr.ind = T)
for(i in 1:nrow(na.index)){
  bulk_all_bv1_cor_matrix[na.index[i, 1], na.index[i, 2]] <- 0
}

a = heatmap(abs(bulk_all_bv1_cor_matrix))
cov_order =rownames(bulk_all_bv1_cor_matrix)[a$rowInd]
bulk_all_bv1_cor_matrix = bulk_all_bv1_cor_matrix[cov_order, cov_order]

pdf('heatmap_exprs_PCA_bulk_bv1_corr.pdf', width = 8, height = 9)
corrplot(bulk_all_bv1_cor_matrix, title = 'bulk all')
dev.off()


sample.remain = setdiff(colnames(d)[-c(1:4)], c("194_D1", "470_D2", "919_D2", '416_D1'))
df.pc.all[c("194_D1", "470_D2", "919_D2"), ]

d = read_tsv('bulk/bulk_new_2504.bed')
pca <- d[sample.remain] %>% t() %>% prcomp(.,scale.=T)
jpeg(filename = 'pca_elbow_bulk_all_remain_0425.jpg', width = 500, height = 300)
plot(pca$sdev[1:60], main = 'bulk all, remaining') # 10
dev.off()

df.pc.all = as.data.frame(pca$x[, 1:5])
colnames(df.pc.all) = paste0('exprs_PC', 1:5)

df.pc.all = cbind(df.pc.all, df.PC.new[rownames(df.pc.all), ])

bulk_all_cor_matrix <- cor(df.pc.all)
bulk_all_cor_matrix['RIN', ] = cor(df.pc.all$RIN[!is.na(df.pc.all$RIN)], df.pc.all[!is.na(df.pc.all$RIN), ])
bulk_all_cor_matrix[, 'RIN'] = cor(df.pc.all$RIN[!is.na(df.pc.all$RIN)], df.pc.all[!is.na(df.pc.all$RIN), ])
bulk_all_cor_matrix['PMI_hours', ] = cor(df.pc.all$PMI_hours[!is.na(df.pc.all$PMI_hours)], df.pc.all[!is.na(df.pc.all$PMI_hours), ])
bulk_all_cor_matrix[, 'PMI_hours'] = cor(df.pc.all$PMI_hours[!is.na(df.pc.all$PMI_hours)], df.pc.all[!is.na(df.pc.all$PMI_hours), ])
bulk_all_cor_matrix['pH', ] = cor(df.pc.all$pH[!is.na(df.pc.all$pH)], df.pc.all[!is.na(df.pc.all$pH), ])
bulk_all_cor_matrix[, 'pH'] = cor(df.pc.all$pH[!is.na(df.pc.all$pH)], df.pc.all[!is.na(df.pc.all$pH), ])


corrplot(bulk_all_cor_matrix)
round(bulk_all_cor_matrix, 3)
na.index = which(is.na(bulk_all_cor_matrix), arr.ind = T)
for(i in 1:nrow(na.index)){
  bulk_all_cor_matrix[na.index[i, 1], na.index[i, 2]] <- 0
}

a = heatmap(abs(bulk_all_cor_matrix))
cov_order =rownames(bulk_all_cor_matrix)[a$rowInd]
bulk_all_cor_matrix = bulk_all_cor_matrix[cov_order, cov_order]

pdf('heatmap_exprs_PCA_bulk_all_corr_remain_0425.pdf', width = 8, height = 9)
corrplot(bulk_all_cor_matrix, title = 'bulk all')
dev.off()
# Significant ones, what we should watching here

pdf('exprs_PCA_bulk_all_remain_0425.pdf', width = 10, height = 8)
gg1 = ggplot(df.pc.all, aes(exprs_PC1, exprs_PC2, col = as.factor(batch))) + geom_point() + theme_minimal() + 
  ggtitle('Batch') + scale_color_brewer(palette = 'Set1')
gg2 = ggplot(df.pc.all, aes(exprs_PC1, exprs_PC2, col = Period)) + geom_point() + theme_minimal() + 
  ggtitle('Period') + scale_color_distiller(palette = 'Spectral')
gg3 = ggplot(df.pc.all, aes(exprs_PC2, exprs_PC5, col = as.factor(batch))) + geom_point() + theme_minimal() + 
  ggtitle('Batch') + scale_color_brewer(palette = 'Set1')
gg4 = ggplot(df.pc.all, aes(exprs_PC2, exprs_PC3, col = RIN)) + geom_point() + theme_minimal() + 
  ggtitle('RIN') + scale_color_distiller(palette = 'Spectral')
plot_grid(gg1, gg2, gg3, gg4)
dev.off()

# 219 sample
# elbow plot for select the best number of expression PCs.
# bulk data Prenatal
d = read_tsv('bulk_prenatal/bulk_pre_new_2504.bed')
pca <- d[intersect(colnames(d), sample.remain)] %>% t() %>% prcomp(.,scale.=T)
jpeg(filename = 'pca_elbow_bulk_pre_remain_0425.jpg', width = 500, height = 300)
plot(pca$sdev[1:60], main = 'bulk prenatal') # 11
dev.off()

# get First few PCs colored by batch
df.pc.pre = as.data.frame(pca$x[, 1:5])
colnames(df.pc.pre) = paste0('exprs_PC', 1:5)
df.pc.pre = cbind(df.pc.pre, df.PC.new[match(rownames(df.pc.pre), rownames(df.PC.new)), ])

df.pc.pre = df.pc.pre[, -20]
bulk_pre_cor_matrix <- cor(df.pc.pre)
bulk_pre_cor_matrix['RIN', ] = cor(df.pc.pre$RIN[!is.na(df.pc.pre$RIN)], df.pc.pre[!is.na(df.pc.pre$RIN), ])
bulk_pre_cor_matrix[, 'RIN'] = cor(df.pc.pre$RIN[!is.na(df.pc.pre$RIN)], df.pc.pre[!is.na(df.pc.pre$RIN), ])
bulk_pre_cor_matrix['PMI_hours', ] = cor(df.pc.pre$PMI_hours[!is.na(df.pc.pre$PMI_hours)], df.pc.pre[!is.na(df.pc.pre$PMI_hours), ])
bulk_pre_cor_matrix[, 'PMI_hours'] = cor(df.pc.pre$PMI_hours[!is.na(df.pc.pre$PMI_hours)], df.pc.pre[!is.na(df.pc.pre$PMI_hours), ])
#bulk_pre_cor_matrix['pH', ] = cor(df.pc.pre$pH[!is.na(df.pc.pre$pH)], df.pc.pre[!is.na(df.pc.pre$pH), ])
#bulk_pre_cor_matrix[, 'pH'] = cor(df.pc.pre$pH[!is.na(df.pc.pre$pH)], df.pc.pre[!is.na(df.pc.pre$pH), ])


corrplot(bulk_pre_cor_matrix)
round(bulk_pre_cor_matrix, 3)
na.index = which(is.na(bulk_pre_cor_matrix), arr.ind = T)
for(i in 1:nrow(na.index)){
  bulk_pre_cor_matrix[na.index[i, 1], na.index[i, 2]] <- 0
}

a = heatmap(abs(bulk_pre_cor_matrix))
cov_order =rownames(bulk_pre_cor_matrix)[a$rowInd]
bulk_pre_cor_matrix = bulk_pre_cor_matrix[cov_order, cov_order]

pdf('heatmap_exprs_PCA_bulk_pre_corr_remain_0425.pdf', width = 8, height = 9)
corrplot(bulk_pre_cor_matrix, title = 'bulk prenatal')
dev.off()

# What happened to pH? we only have 5 available data point
bulk_pre_cor_matrix = bulk_pre_cor_matrix[-1, -1]
a = heatmap(abs(bulk_pre_cor_matrix))
cov_order =rownames(bulk_pre_cor_matrix)[a$rowInd]
bulk_pre_cor_matrix = bulk_pre_cor_matrix[cov_order, cov_order]
pdf('heatmap_exprs_PCA_bulk_pre_corr_no_pH.pdf', width = 8, height = 9)
corrplot(bulk_pre_cor_matrix, title = 'bulk prenatal')
dev.off()


pdf('exprs_PCA_bulk_pre_remain.pdf', width = 10, height = 8)
gg1 = ggplot(df.pc.pre, aes(exprs_PC1, exprs_PC2, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('Batch') + scale_color_brewer(palette = 'Set1') 
gg2 = ggplot(df.pc.pre, aes(exprs_PC2, exprs_PC3, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('Batch') + scale_color_brewer(palette = 'Set1') 
gg3 = ggplot(df.pc.pre, aes(exprs_PC3, Period, col = as.factor(batch))) + geom_point() + 
  scale_color_brewer(palette = 'Set1') + ggtitle('Period') + theme_minimal()
gg4 = ggplot(df.pc.pre, aes(exprs_PC2, RIN, col = as.factor(batch))) + geom_point() + 
  scale_color_brewer(palette = 'Set1') + ggtitle('RIN') + theme_minimal()
plot_grid(gg1, gg2, gg3, gg4)
dev.off()
pdf('exprs_PCA_bulk_pre_Batch_remain_0425.pdf', width = 10, height = 8)
ggplot(df.pc.pre, aes(exprs_PC1, exprs_PC5, col = as.factor(batch), label = rownames(df.pc.pre))) + geom_point() + 
  theme_minimal() + ggtitle('Batch') + scale_color_brewer(palette = 'Set1')  + ggrepel::geom_text_repel()
dev.off()




d = read_tsv('bulk_postnatal/bulk_post_new_2504.bed')
pca <- d[intersect(colnames(d), sample.remain)] %>% t() %>% prcomp(.,scale.=T)
jpeg(filename = 'pca_elbow_bulk_post_remain_0425.jpg', width = 500, height = 300)
plot(pca$sdev[1:60], main = 'bulk postnatal') # 10
dev.off()

# get First few PCs colored by batch
df.pc.post = as.data.frame(pca$x[, 1:5])
colnames(df.pc.post) = paste0('exprs_PC', 1:5)

df.pc.post = cbind(df.pc.post, df.PC.new[match(rownames(df.pc.post), rownames(df.PC.new)), ])

bulk_post_cor_matrix <- cor(df.pc.post)
bulk_post_cor_matrix[, 'RIN'] = cor(df.pc.post$RIN[!is.na(df.pc.post$RIN)], df.pc.post[!is.na(df.pc.post$RIN), ])
bulk_post_cor_matrix['RIN', ] = cor(df.pc.post$RIN[!is.na(df.pc.post$RIN)], df.pc.post[!is.na(df.pc.post$RIN), ])
bulk_post_cor_matrix['PMI_hours', ] = cor(df.pc.post$PMI_hours[!is.na(df.pc.post$PMI_hours)], df.pc.post[!is.na(df.pc.post$PMI_hours), ])
bulk_post_cor_matrix[, 'PMI_hours'] = cor(df.pc.post$PMI_hours[!is.na(df.pc.post$PMI_hours)], df.pc.post[!is.na(df.pc.post$PMI_hours), ])
bulk_post_cor_matrix['pH', ] = cor(df.pc.post$pH[!is.na(df.pc.post$pH)], df.pc.post[!is.na(df.pc.post$pH), ])
bulk_post_cor_matrix[, 'pH'] = cor(df.pc.post$pH[!is.na(df.pc.post$pH)], df.pc.post[!is.na(df.pc.post$pH), ])

na.index = which(is.na(bulk_post_cor_matrix), arr.ind = T)
for(i in 1:nrow(na.index)){
  bulk_post_cor_matrix[na.index[i, 1], na.index[i, 2]] <- 0
}

corrplot(bulk_post_cor_matrix)
round(bulk_post_cor_matrix, 3)

a = heatmap(abs(bulk_post_cor_matrix))
cov_order =rownames(bulk_post_cor_matrix)[a$rowInd]
bulk_post_cor_matrix = bulk_post_cor_matrix[cov_order, cov_order]

pdf('heatmap_exprs_PCA_bulk_post_corr_remain_0425.pdf', width = 8, height = 9)
corrplot(bulk_post_cor_matrix, title = 'bulk postnatal')
dev.off()

pdf('exprs_PCA_bulk_post_remain.pdf', width = 10, height = 8)
gg1 = ggplot(df.pc.post, aes(exprs_PC1, exprs_PC2, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('Batch') + scale_color_brewer(palette = 'Set1')
gg2 = ggplot(df.pc.post, aes(exprs_PC1, Period, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('Period') + scale_color_brewer(palette = 'Set1')
gg3 = ggplot(df.pc.post, aes(exprs_PC3, Period, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('Period') + scale_color_brewer(palette = 'Set1')
gg4 = ggplot(df.pc.post, aes(exprs_PC2, RIN, col = as.factor(batch))) + geom_point() + 
  theme_minimal() + ggtitle('RIN') + scale_color_brewer(palette = 'Set1')
plot_grid(gg1, gg2, gg3, gg4)
dev.off()
pdf('exprs_PCA_bulk_post_Batch_remain_0425.pdf', width = 10, height = 8)
ggplot(df.pc.post, aes(exprs_PC1, exprs_PC2, col = as.factor(batch), label = rownames(df.pc.post))) + geom_point() + 
  theme_minimal() + ggtitle('Batch') + scale_color_brewer(palette = 'Set1')  + ggrepel::geom_text_repel()
dev.off()
```
