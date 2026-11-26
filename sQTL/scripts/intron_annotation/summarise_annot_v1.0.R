### Last updated: 2025.07.22 by SWK
### project: Annotate introns - summary

rm(list=ls())
sessionInfo()
setwd('~/Dropbox/working_BrainVar2_sQTL/Inputs/')

options(stringsAsFactors = F)

library(tidyverse)
library(data.table)
library(scales)
library(extrafont)

#font_import()
loadfonts()

`%notin%` <- Negate(`%in%`)
select <- dplyr::select

date = '250722'



### Summarize GTF-based annotations
## Raw ----------------------------------------------------------------- 
load('leafviz_unfiltered.RData') # 383851 

cat("Merging annotation fields...\n")
key_coords <- paste(all.introns$chr, all.introns$start, all.introns$end)

all.introns$verdict <- unlist(verdict.list)[match(key_coords, unlist(coord.list))]
all.introns$gene <- unlist(gene.list)[match(key_coords, unlist(coord.list))]
all.introns$ensemblID <- unlist(ensemblID.list)[match(key_coords, unlist(coord.list))]
all.introns$transcripts <- unlist(transcripts.list)[match(key_coords, unlist(coord.list))]
all.introns$constitutive.score <- unlist(constitutive.list)[match(key_coords, unlist(coord.list))]

# Replace missing values
all.introns <- all.introns %>%
  mutate(
    gene = ifelse(is.na(gene), ".", gene),
    ensemblID = ifelse(is.na(ensemblID), ".", ensemblID),
    transcripts = ifelse(transcripts == "" | is.na(transcripts), ".", transcripts),
    constitutive.score = signif(constitutive.score, digits = 2)
  )

# Save result
output_path <- paste0("all.introns.", date, ".rds")
saveRDS(all.introns, file = output_path)
cat("✅ Saved summarised intron annotations to", output_path, "\n")


## Filtered ----------------------------------------------------------------- 
rm(list=ls()); date = '250722'

load('leafviz_filtered.RData') # 285600

cat("Merging annotation fields...\n")
key_coords <- paste(all.introns$chr, all.introns$start, all.introns$end)

all.introns$verdict <- unlist(verdict.list)[match(key_coords, unlist(coord.list))]
all.introns$gene <- unlist(gene.list)[match(key_coords, unlist(coord.list))]
all.introns$ensemblID <- unlist(ensemblID.list)[match(key_coords, unlist(coord.list))]
all.introns$transcripts <- unlist(transcripts.list)[match(key_coords, unlist(coord.list))]
all.introns$constitutive.score <- unlist(constitutive.list)[match(key_coords, unlist(coord.list))]

# Replace missing values
all.introns <- all.introns %>%
  mutate(
    gene = ifelse(is.na(gene), ".", gene),
    ensemblID = ifelse(is.na(ensemblID), ".", ensemblID),
    transcripts = ifelse(transcripts == "" | is.na(transcripts), ".", transcripts),
    constitutive.score = signif(constitutive.score, digits = 2)
  ) %>%
  filter(chr != 'chrID') # 285600

# Save result
output_path <- paste0("filtered.introns.", date, ".rds")
saveRDS(all.introns, file = output_path)
cat("✅ Saved summarised intron annotations to", output_path, "\n")



### Plot proportion
rm(list=ls()); date = '250722'


## Load data ----------------------------------------------------------------- 
# Raw
raw_cts = readRDS('all.introns.250722.rds') # 383,851 introns 
length(unique(raw_cts$clusterID)) # 80,199 clusters
length(unique(raw_cts$gene)) # 20,364 genes

# Filtered
filtered_cts = readRDS('filtered.introns.250722.rds') # 285,600 introns
length(unique(filtered_cts$clusterID)) # 55,268 clusters
length(unique(filtered_cts$gene)) # 14,542 genes


## Distribution of verdict ----------------------------------------------------------------- 
get_verdict_summary <- function(df, dataset_name) {
  df %>%
    count(verdict) %>%
    mutate(
      dataset = dataset_name,  # 먼저 넣고
      percent = n / sum(n),
      label = paste0(scales::comma(n), "\n(", percent_format(accuracy = 0.1)(percent), ")")
    )
}

raw_summary <- get_verdict_summary(raw_cts, " Raw introns")
filtered_summary <- get_verdict_summary(filtered_cts, "Filtered introns")

combined_summary <- bind_rows(raw_summary, filtered_summary)
# write.table(combined_summary, file = paste0("../Tables/AnnotationSummary_GTFbased.", date, ".txt"), sep = '\t', row.names = F, col.names = T, quote=F)


## Pie chart ----------------------------------------------------------------- 
ggplot(combined_summary, aes(x = "", y = percent, fill = verdict)) +
  geom_col(width = 1, color = "white") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5,           # 텍스트 크기 조절
    family = "Arial"      # 텍스트 폰트
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_manual(
    values = c(
      "annotated" = "#E74C3C",
      "cryptic_fiveprime" = "#76D7C4",
      "cryptic_threeprime" = "#1ABC9C",
      "cryptic_unanchored" = "#34495E",
      "novel annotated pair" = "#F5B7B1",
      "unknown_strand" = "#AAB7B8"
    ), name = NULL
  ) +
  labs(
    title = "Intron annotation category"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(family = "Arial", size = 22, face = "bold", hjust = 0.5),
    legend.title = element_text(family = "Arial", size = 16),
    legend.text = element_text(family = "Arial", size = 14),
    strip.text = element_text(family = "Arial", size = 18, face = "bold")
  )
# ggsave(paste('../Figures/AnnotationSummary_GTFbased.PiePlot', date, 'pdf', sep = '.'), last_plot(), width = 12, height = 6)









