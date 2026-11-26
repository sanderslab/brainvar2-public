# Read in files and Data Preparations
```
library(tidyverse)
library(parallel)
library(reshape2)
```

## Read in meta.data
```
# read in meta.data from brainvar2/meta_data/brainvar2_library_metadata_v1.2.tsv
meta.data.library = read_tsv('brainvar2_library_metadata_v1.2.tsv') # 842 library
meta.data.library = meta.data.library[meta.data.library$Usable == 'Yes', ] #770 library
```

```
# bulk
meta.data.library.bulkRNA = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'bulkRNA', ] #230
# genotype
meta.data.library.WGS = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'WGS', ] #237

# single cell
meta.data.library.GEX = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'GEX', ] #153
```
## Genotype and PCA 
```
# Genotype PCA loading
temp = read.table('./your_path/brainvar2.022924.joint_called.recalibrated.varQC.AF.PCAscores.txt', header = T) # 265 
s = temp$s
genotype.PCA = do.call('rbind', lapply(strsplit(gsub('[]]', replacement = '', gsub('[[]', replacement = '', temp$scores)), ','), as.numeric))
rownames(genotype.PCA) = s
colnames(genotype.PCA) = paste0('PC', 1:20)
```
Genotype vcf file has 265 samples
```
s # name of samples in genotype data
# [1] "092_D1"  "100_D1"  "102_D1"  "103_D1"  "105_D1"  "107_D1"  "113_D1"  "114_D1"  "115_D1"  "115_D2"  "116_D1"  "118_D1"  "119_D1"  "121_D1"  "122_D1"  "127_D1" 
```
```
PC.mat = data.frame(genotype.PCA)
PC.mat = PC.mat[intersect(rownames(PC.mat), meta.data.library.WGS$LibraryID), ] #230 samples left
```
**230 usable samples (with name XXX_D1 XXX_D2) for genotype data.**

```
df.PC = data.frame(PC1 = PC.mat$PC1, PC2 = PC.mat$PC2, PC3 = PC.mat$PC3, Sex = NA, Ancestry = NA, RIN = NA, Age = NA, 
                   SubjectID = NA, Period = NA, Epoch = NA)
rownames(df.PC) = rownames(PC.mat)

common.dna.subject = intersect(rownames(df.PC), meta.data.library.WGS$LibraryID)
df.PC[common.dna.subject, 'Sex'] = meta.data.library.WGS$SexForAnalysis[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'Ancestry'] = meta.data.library.WGS$AncestryPCA[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'RIN'] = meta.data.library.WGS$RIN[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'Age'] = meta.data.library.WGS$AgeDays[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'SubjectID'] = meta.data.library.WGS$SubjectID[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'Period'] = meta.data.library.WGS$Period[match(common.dna.subject, meta.data.library.WGS$LibraryID)]
df.PC[common.dna.subject, 'Epoch'] = meta.data.library.WGS$Epoch[match(common.dna.subject, meta.data.library.WGS$LibraryID)]

# remove sample without Period lable
df.PC = df.PC[!is.na(df.PC$Period), ]
# 229 samples

which(table(df.PC$SubjectID) > 1)
# HSB513

save(df.PC, file = './data_store/df_PC.RData')
```

## Bulk data 
read in the htseq_counts data and meta.data for bulk (updated bulk data: 250401_htseq_counts_BrainVar.tsv)
```
# read in the htseq_counts
# bulk_all = read_tsv('bulk/htseq_counts_BrainVar.tsv') # 261 samples
bulk_all = read_tsv('bulk/250401_htseq_counts_BrainVar.tsv') # 261 samples
# bulk meta.data
meta.data.library.bulkRNA = meta.data.library[meta.data.library$Usable == 'Yes' & meta.data.library$LibraryModality == 'bulkRNA', ] #230
```
column names for the bulk data
```
bulk.real.names = data.frame(bulk.data.names = colnames(bulk_all)[-1])
bulk.real.names$libraryID = bulk.real.names$bulk.data.names
```
#### convert names to standard names for RNA-seq data (XXX_R1 or XXX_R2)
```
bulk.real.names$libraryID[216:235]
[1] "5212" "5699" "5800" "5864" "5871" "6041" "6085" "6119" "6154" "6161" "6195" "6197" "8035" "8039" "8050" "8053" "8073" "8079" "8153" "8164"
bulk.real.names$libraryID[216:235] = paste0(bulk.real.names$bulk.data.names[216:235],'_R2')

bulk.real.names$libraryID[236:261]
# [1] "H115" "H175" "H271" "H279" "H321" "H332" "H404" "H420" "H425" "H462" "H470" "H473" "H500" "H505" "H513" "H515" "H568" "H573" "H583" "H675" "H771" "H793" "H817"
# [24] "H818" "H914" "H919"
bulk.real.names$libraryID[236:261] = paste0(gsub('H', replacement = '', bulk.real.names$bulk.data.names[236:261]), '_R2')

bulk.real.names$libraryID[1:215]
# [1] "HSB100" "HSB102" "HSB105" "HSB107" "HSB115" "HSB116" "HSB121" "HSB122" "HSB127" "HSB132" "HSB139" "HSB148" "HSB150" "HSB152" "HSB154" "HSB171" "HSB172" "HSB174"
bulk.real.names$libraryID[1:215] = paste0(gsub('HSB', replacement = '', bulk.real.names$bulk.data.names[1:215]), '_R1')
```
clean up
```
gene.names = bulk_all$gene
colnames(bulk_all) = c('gene', bulk.real.names$libraryID)

# overlap with usable library
bulk_all = bulk_all[, colnames(bulk_all) %in% meta.data.library.bulkRNA$LibraryID] #230 samples
bulk_all$gene = gene.names
# write_tsv(bulk_all, file = 'bulk/counts_htseq_clean.tsv')
write_tsv(bulk_all, file = 'bulk/counts_htseq_clean_2504.tsv')

# How I get v1.3 with $realnamebulk
meta.data.library$realnamebulk = NA
bulk.real.names.sub = bulk.real.names[bulk.real.names$libraryID %in% meta.data.library$LibraryID, ]
meta.data.library$realnamebulk[match(bulk.real.names.sub$libraryID, meta.data.library$LibraryID)] = bulk.real.names.sub$bulk.data.names
```
#### Remove HSB571, due to no period label
```
table(meta.data.subject$Period[match(meta.data.library.bulkRNA$SubjectID, meta.data.subject$SubjectID)])
# 1 sample does not have reported age. HSB571

bulk_all = bulk_all[, -match('571_R1', colnames(bulk_all))]
bulk_dev = bulk_all[, which(colnames(bulk_all) %in% meta.data.library.bulkRNA$LibraryID[meta.data.library.bulkRNA$Period >= 8])]
bulk_fetal = bulk_all[, which(colnames(bulk_all) %in% meta.data.library.bulkRNA$LibraryID[meta.data.library.bulkRNA$Period <= 6])]

# save(bulk_all, file = 'bulk/bulk_all.RData')
save(bulk_all, file = 'bulk/bulk_all_2504.RData')
```
## Match between bulk and genotype samples
```
sp = rownames(PC.mat) # DNA library names
bulk_all = bulk_all[, colnames(bulk_all) %in% c('gene', meta.data.library.bulkRNA$LibraryID[which(meta.data.library.bulkRNA$matchingDNALibrary %in% sp)])]

# 5 samples are removed
# "515_R1" "321_R2" "332_R2" "505_R2" "793_R2"
```

## Pseudobulk 
Please see brainvar2/pseudobulk/pseudobulk_TMM.R and brainvar2/pseudobulk/pseudobulk_counts.R. 
Single-cell samples are one-on-one matched with DNA library. 




