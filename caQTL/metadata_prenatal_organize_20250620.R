library(tidyverse)

m1 = read.delim('~/Dropbox/BrainVar2_WD/Data/brainvar2_single_cell_detailed_metadata_v1.3.tsv')
m1 = m1 %>%
  filter(Usable_ATAC=='Yes')
table(duplicated(m1$SubjectID))

m2 = read.delim('~/Dropbox/BrainVar2_WD/Data/brainvar2_library_metadata_v1.3.tsv')
m2 = m2 %>%
  filter(LibraryModality=='ATAC')
table(duplicated(m2$SubjectID))

m3 = read.delim('~/Dropbox/BrainVar2_WD/Data/brainvar2_subject_metadata_v1.3.tsv')
m3 = m3 %>%
  filter(ConfirmedID_ATAC=='Yes')
table(duplicated(m3$SubjectID))




dt <- readRDS("~/Dropbox/BrainVar2_WD/Data/pseudo_bulk_SE_prenatal.rds")

colData(dt)$CellType <- sapply(strsplit(as.character(colData(dt)$Group), "_"), `[`, 2)

df = colData(dt) %>% as.data.frame()
df = df %>%
  filter(CellType!='None')

df$BioSampleID <- sapply(strsplit(as.character(df$Group), "_"), `[`, 1)




overlap_12 <- intersect(colnames(m1), colnames(m2))
m2_clean <- m2 %>% select(-all_of(overlap_12[overlap_12 != "BioSampleID"]))

# Join m1 and cleaned m2
df_12 <- inner_join(
  m1 %>% filter(BioSampleID %in% df$BioSampleID),
  m2_clean %>% filter(BioSampleID %in% df$BioSampleID),
  by = "BioSampleID"
)

overlap_123 <- intersect(colnames(df_12), colnames(m3))
m3_clean <- m3 %>% select(-all_of(overlap_123[overlap_123 != "SubjectID"]))

df_all <- inner_join(df_12, m3_clean, by = "SubjectID")
df_all = df_all %>%
  filter(Usable_DNA=='Yes')

write.table(df_all,
            '~/Dropbox/BrainVar2_WD/Data/brainvar2.prenatal.metadata.20250620.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')
