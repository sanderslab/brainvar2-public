### Last updated: 2025.07.30 by SWK
### project: Phenotype & covariate matrix process

rm(list=ls())
sessionInfo()
setwd('~/Dropbox/working_BrainVar2_sQTL/Inputs/')

options(stringsAsFactors = F)

library(tidyverse)
library(foreach)
library(rstatix)
library(readxl)
library(data.table)
library(sva)
library(MASS)
library(pracma)
library(ggcorrplot)
library(extrafont)

#font_import()
loadfonts()

`%notin%` <- Negate(`%in%`)
select <- dplyr::select

date = '250730'

# Function definition
theme_step1 <- function(base_size = 11, base_family = "",
                        base_line_size = base_size / 22,
                        base_rect_size = base_size / 22) {
  theme(title = element_text(family = 'Arial', size = 18, color = 'black'), text = element_text(family = 'Arial', size = 16, color = 'black'),
        axis.title = element_text(family = 'Arial', size = 18, color = 'black'), axis.text = element_text(family = 'Arial', size = 16, color = 'black'), 
        strip.text = element_text(family = 'Arial', size = 18, color = 'black'),
        panel.grid.major = element_line(color = alpha('lightgrey', .5), size = rel(1)), 
        panel.grid.minor = element_line(color = alpha('lightgrey', .5), size = rel(.5)),
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", colour = NA), axis.line = element_line(colour = "black", size = rel(1)),
        legend.background = element_rect(color = 'black'), legend.title = element_text(family = 'Arial', size = 16),
        legend.text = element_text(family = 'Arial', size = 14),
        legend.direction = "vertical", 
        legend.box = c("horizontal", "vertical"),
        legend.spacing.x = unit(0.1, 'cm'),
        plot.margin = unit(c(0.25, 1, 1, 0.5), 'cm'))
}

#parse and extract genotype PCs
parse_scores <- function(score_str) {
  score_str %>%
    str_remove_all("\\[|\\]") %>%  # remove brackets
    str_split(",") %>%
    unlist() %>%
    as.numeric()
}

extract_gPC_df <- function(filepath) {
  # Read file
  df <- read.delim(filepath)
  
  # Parse and convert scores column
  scores_matrix <- df$scores %>%
    map(parse_scores) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  # Set column names
  colnames(scores_matrix) <- paste0("gt_PC", seq_len(ncol(scores_matrix)))
  
  # Combine with ID
  parsed_df <- bind_cols(df %>% rename(WGS_ID = s) %>% select(WGS_ID), scores_matrix)
  
  return(parsed_df)
}



### Load data
# metadata
meta = read.delim('../Tables/sample_list/SampleMetadata_Annotated.250626.txt') %>% # 230*16
  filter(matchingDNALibrary %notin% c('', '505_D2', '793_D2', '332_D2', '321_D1')) %>% # 5 samples excluded from WGS QC
  filter(matchingDNALibrary %notin% c('175_D1', '515_D2')) %>% # 2 samples excluded (gPC outliers)
  filter(Stage %in% c('Prenatal', 'Postnatal')) %>% # 219 (124 + 95)
  mutate(Stage = factor(Stage, levels = c('Prenatal', 'Postnatal')))

# gPC
gPC_pre <- extract_gPC_df("../Tables/covariates/brainvar2.022924.joint_called.recalibrated.varQC.AF0.02.r20.1.sQTLprenatal.PCAscores.txt") %>%
  select(WGS_ID, gt_PC1, gt_PC2, gt_PC3)
gPC_post <- extract_gPC_df("../Tables/covariates/brainvar2.022924.joint_called.recalibrated.varQC.AF0.02.r20.1.sQTLpostnatal.PCAscores.txt") %>%
  select(WGS_ID, gt_PC1, gt_PC2, gt_PC3)

# normalized counts
cts = fread('leafcutter_perind.counts.nochr.gz.qqnorm_all.gz')
cts <- cts[cts$'#Chr' != '#Chr', ]



### Phenotype process
cts_df <- as.data.frame(cts) #285600*229

valid_samples <- meta$matchingDNALibrary # 219

# Colnames
old_subj <- colnames(cts_df)[5:ncol(cts_df)]
new_subj <- sapply(old_subj, function(x) strsplit(x, split="\\.")[[1]][1])
colnames(cts_df) <- c("#Chr", "start", "end", "ID", new_subj)

# Filter valid samples only
keep_samples <- intersect(valid_samples, new_subj) # 219

keep_idx <- match(keep_samples, colnames(cts_df)) # 219

# Subset
no_combat <- cts_df[, c(1:4, keep_idx)] # 223 (219+4)
no_combat_mat <- as.matrix(no_combat[, 5:ncol(no_combat)])

storage.mode(no_combat_mat) <- "numeric"
rownames(no_combat_mat) <- no_combat$ID

# Match metadata
meta_subset <- meta %>%
  select(matchingDNALibrary, batch, Stage)


## Split by Stage -----------------------------------------------------------------
prenatal_samples <- meta_subset$matchingDNALibrary[meta_subset$Stage == "Prenatal" & !is.na(meta_subset$Stage)] #124
postnatal_samples <- meta_subset$matchingDNALibrary[meta_subset$Stage == "Postnatal" & !is.na(meta_subset$Stage)] #95

# No ComBat + prenatal
no_combat_prenatal <- no_combat_mat[, prenatal_samples, drop=FALSE]
no_combat_prenatal_df <- as.data.frame(no_combat_prenatal) #285600*124
# write.table(no_combat_prenatal_df, "leafcutter_perind.counts.qqnorm_all.NoComBat.prenatal.tsv", quote=F, sep="\t", col.names=T, row.names=T)

no_combat_prenatal_bed <- cbind(cts[,1:4], no_combat_prenatal_df) #285600*128
# write.table(no_combat_prenatal_bed, "leafcutter_perind.counts.qqnorm_all.NoComBat.prenatal.bed", quote=F, sep="\t", col.names=T, row.names=F)

# No ComBat + postnatal
no_combat_postnatal <- no_combat_mat[, postnatal_samples, drop=FALSE]
no_combat_postnatal_df <- as.data.frame(no_combat_postnatal) #285600*95
# write.table(no_combat_postnatal_df, "leafcutter_perind.counts.qqnorm_all.NoComBat.postnatal.tsv", quote=F, sep="\t", col.names=T, row.names=T)

no_combat_postnatal_bed <- cbind(cts[,1:4], no_combat_postnatal_df) #285600*99
# write.table(no_combat_postnatal_bed, "leafcutter_perind.counts.qqnorm_all.NoComBat.postnatal.bed", quote=F, sep="\t", col.names=T, row.names=F)



### PCs from phenotype matrix
no_combat_prenatal = fread("leafcutter_perind.counts.qqnorm_all.NoComBat.prenatal.tsv", data.table=F) # 285600*125
no_combat_postnatal = fread("leafcutter_perind.counts.qqnorm_all.NoComBat.postnatal.tsv", data.table=F) # 285600*95

dataset_list <- list(
  NoComBat_Prenatal = no_combat_prenatal,
  NoComBat_Postnatal = no_combat_postnatal
)

## Run PCA & elbow plot -----------------------------------------------------------------
output_dir <- "../Figures/ElbowPlots/"
dir.create(output_dir, showWarnings=FALSE)

# Each dataset
all_pca_df <- list()

for (name in rev(names(dataset_list))) {
  cat("\n\n=============================\n")
  cat("Dataset:", name, "\n")
  cat("=============================\n")
  
  dat <- dataset_list[[name]]
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  
  # Scale
  dat_mat <- t(scale(t(as.matrix(dat)), center=TRUE, scale=TRUE))
  
  # PCA
  pca <- prcomp(t(dat_mat), scale=FALSE)
  
  # var_exp <- pca$sdev^2
  # var_exp_ratio <- var_exp / sum(var_exp)
  # cum_var_exp <- cumsum(var_exp_ratio)
  # 
  # cat("Cumulative variance explained by PCs:\n")
  # print(round(cum_var_exp[1:10],3))
  # 
  # # SVA
  # mod <- model.matrix(~1, data=as.data.frame(t(dat_mat)))
  # n.sv <- num.sv(dat_mat, mod=mod, method="be")
  # cat("Estimated number of surrogate variables by SVA (method=be):", n.sv, "\n")
  # 
  # # Elbow plot
  # df_elbow <- data.frame(
  #   PC=1:length(var_exp),
  #   Variance=var_exp
  # )
  # 
  # p_elbow <- ggplot(df_elbow, aes(x=PC, y=Variance)) +
  #   geom_point() +
  #   geom_line() +
  #   ggtitle(paste("Elbow Plot:", name)) +
  #   ylab("Eigenvalue") +
  #   xlab("Principal Component") +
  #   theme_minimal(base_size=14) +
  #   geom_vline(xintercept=n.sv, linetype="dashed", color="red") +
  #   annotate("text", x=n.sv+10, y=max(var_exp)*0.9,
  #            label=paste("SVA =", n.sv), color="red") +
  #   scale_x_continuous(breaks = seq(0, max(df_elbow$PC), by = 5))
  # 
  # print(p_elbow)
  # 
  # # save
  # ggsave(
  #   filename=paste0(output_dir, name, "_ElbowPlot.pdf"),
  #   plot=p_elbow,
  #   width=8, height=5
  # )
  
  # PCA1~20 save
  pca_df <- as.data.frame(pca$x[, 1:20])
  pca_df$Sample <- rownames(pca_df)
  pca_df$Dataset <- name
  
  # Metadata merge
  pca_df <- pca_df %>%
    left_join(meta, by=c("Sample"="matchingDNALibrary"))
  
  all_pca_df[[name]] <- pca_df
}


# Combine PCA results
plot_df <- bind_rows(all_pca_df)
# write.table(plot_df, paste('../Tables/PCs_AllCombined', date, 'txt', sep = '.'), sep = '\t', row.names = F, col.names = T, quote=F)


## PCA plot -----------------------------------------------------------------
# Period to discrete factor
plot_df$Period <- factor(plot_df$Period)

# abs_max
abs_max <- plot_df %>%
  select(PC1, PC2) %>%
  summarise(max_val = max(abs(c(PC1, PC2)))) %>%
  pull(max_val)

# Set colors
make_pca_plot <- function(color_var) {
  
  base_plot <- ggplot(plot_df, aes(x = PC1, y = PC2, fill = .data[[color_var]])) +
    geom_point(size = 2, alpha = 0.8, shape = 21) +
    facet_wrap(~ Dataset, scales = "fixed") +
    theme_step1() +
    ggtitle(paste("PCA colored by", color_var)) +
    xlab("PC1") +
    ylab("PC2") +
    coord_cartesian(xlim = c(-abs_max, abs_max),
                    ylim = c(-abs_max, abs_max)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal")
  
  if (color_var == "BrainVar") {
    base_plot + scale_fill_manual(values = c("BrainVar1" = "#95DB90", "BrainVar2" = "#6B508E"))
  } else if (color_var == "Sex") {
    base_plot + scale_fill_manual(values = c("Male" = "#90C0DF", "Female" = "#C593C2"))
  } else if (color_var == "AncestryPCA") {
    base_plot + scale_fill_d3(palette = "category20")
  } else if (color_var == "Period") {
    base_plot + scale_fill_d3(palette = "category20b")
  } else {
    base_plot
  }
}

# Plotting!
p1 <- make_pca_plot("BrainVar")
p2 <- make_pca_plot("Sex")
p3 <- make_pca_plot("AncestryPCA")
p4 <- make_pca_plot("Period")

# Merge plots
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2)
# ggsave(paste("../Figures/PCA_BatchCorrection_Facet", date, "pdf", sep = '.'), last_plot(), width=16, height=13)



### Covariate process
## Distribution of metadata -----------------------------------------------------------------
# batch
table(meta$Stage, meta$batch)

p1 <- ggplot(data = meta) +
  geom_mosaic(aes(weight = 1, x = product(Stage), fill = batch), na.rm = TRUE) +
  labs(title = "Stage × Batch", x = "", y = "") +
  scale_fill_manual(values = c("#95DB90", '#6B508E'), guide = NULL) +
  theme_step1() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5))
p1

# Sex
table(meta$Stage, meta$Sex)

p2 <- ggplot(data = meta) +
  geom_mosaic(aes(weight = 1, x = product(Stage), fill = Sex), na.rm = TRUE) +
  labs(title = "Stage × Sex", x = "", y = "") +
  scale_fill_manual(values = c('#C593C2', '#90C0DF'), guide = NULL) +
  theme_step1() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5))
p2

# Ancestry (PCA)
table(meta$Stage, meta$AncestryPCA)

p3 <- ggplot(data = meta) +
  geom_mosaic(aes(weight = 1, x = product(Stage), fill = AncestryPCA), na.rm = TRUE) +
  labs(title = "Stage × Ancestry (PCA)", x = "", y = "") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00"), guide = NULL) +
  theme_step1() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5))
p3

# PMI hours
p4 <- ggplot(meta, aes(x = PMI_hours)) +
  geom_histogram(binwidth = 5, fill = "#FDB462", color = "black", na.rm = TRUE) +
  facet_wrap(~Stage, scales = "free_y") +
  labs(title = "PMI hours distribution", x = "PMI hours", y = "Count") +
  theme_step1() + 
  coord_cartesian(expand = FALSE)

p4

# RIN
p5 <- ggplot(meta, aes(x = RIN)) +
  geom_histogram(binwidth = 0.5, fill = "#80B1D3", color = "black", na.rm = TRUE) +
  facet_wrap(~Stage, scales = "free_y") +
  labs(title = "RIN distribution", x = "RIN", y = "Count") +
  theme_step1() + 
  coord_cartesian(expand = FALSE)

p5

# Period
p6 <- ggplot(meta, aes(x = Period)) +
  geom_histogram(binwidth = 1, fill = "#B3DE69", color = "black", na.rm = TRUE) +
  facet_wrap(~Stage, scales = "free") +
  labs(title = "Period distribution", x = "Period", y = "Count") +
  scale_x_continuous(breaks = seq(floor(min(meta$Period, na.rm = TRUE)),
                                  ceiling(max(meta$Period, na.rm = TRUE)), by = 1)) +
  theme_step1() + 
  coord_cartesian(expand = FALSE)

p6

# Merge plots
ggarrange(ggarrange(p1, p2, p3, ncol = 3, widths = c(1, .9, 1.2),
                    labels = c('A', 'B', 'C'), font.label = list(size = 28), label.y = 1.01),
          ggarrange(p4, p5, p6, ncol = 3,
                    labels = c('D', 'E', 'F'), font.label = list(size = 28), label.y = 1.01), 
          nrow = 2)
# ggsave(paste('../Figures/BrainVar2_samples_sQTL', date, 'pdf', sep = '.'), last_plot(), width = 18, height = 12)




## Correlations between covariates -----------------------------------------------------------------
all_cov = plot_df %>%
  left_join(gPC_pre %>% rbind(gPC_post), by = c('Sample' = 'WGS_ID')) %>%
  mutate(Sex = ifelse(Sex == "Male", 2, 1),
         batch = ifelse(batch == 'BrainVar1', 1, 2)) %>%
  select(ID = Sample, paste0('PC', 1:8), paste0('gt_PC', 1:3), Sex, RIN, Period, batch, PMI_hours, Stage) 

pre_cov = all_cov %>% filter(Stage == 'Prenatal') %>% select(-Stage) # 124*17
# write.table(pre_cov, paste('../Tables/sQTL_input/Prenatal_cov', date, 'txt', sep = '.'), sep = '\t', row.names = F, col.names = T, quote=F)

post_cov = all_cov %>% filter(Stage == 'Postnatal') %>% select(-Stage) # 95*17
# write.table(post_cov, paste('../Tables/sQTL_input/Postnatal_cov', date, 'txt', sep = '.'), sep = '\t', row.names = F, col.names = T, quote=F)

# Correlation plot
plot_corr_matrix <- function(df, title = "Correlation Matrix") {
  # Numeric only
  df_numeric <- df %>% select(where(is.numeric))
  
  # Correlations
  corr <- cor(df_numeric, use = "pairwise.complete.obs")
  
  # -------- Print high correlations (|r| >= 0.7) --------
  corr_df <- as.data.frame(as.table(corr)) %>%
    filter(Var1 != Var2) %>%  # exclude self-correlation
    filter(abs(Freq) >= 0.7) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(Var1, Var2, Correlation = Freq)
  
  if (nrow(corr_df) > 0) {
    cat("⚠️ High correlation pairs (|r| ≥ 0.7):\n")
    print(corr_df)
  } else {
    cat("✅ No high correlation pairs (|r| ≥ 0.7) found.\n")
  }
  
  # Visualization
  ggcorrplot(
    corr,
    method = "circle",        
    type = "full",           
    lab = FALSE,              
    colors = c("navy", "white", "darkred"),
    outline.color = "white",
    title = title
  ) +
    theme_minimal(base_family = "Arial") +
    theme(
      plot.title = element_text(family = "Arial", size = 22, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = element_text(family = "Arial", angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(family = "Arial", size = 16),
      legend.title = element_text(family = 'Arial', size = 16),
      legend.text = element_text(family = 'Arial', size = 14)
    )
}

plot_corr_matrix(pre_cov, "Prenatal NoComBat (PC8)")
# ✅ No high correlation pairs (|r| ≥ 0.7) found.
# ggsave(paste('../Figures/Cov_PreNoCombatPC8.CorrPlot', date, 'pdf', sep = '.'), last_plot(), width = 9, height = 9)

plot_corr_matrix(post_cov, "Postnatal NoComBat (PC8)")
# ⚠️ High correlation pairs (|r| ≥ 0.7): batch PC1 (0.769)
# ggsave(paste('../Figures/Cov_PostNoCombatPC8.CorrPlot', date, 'pdf', sep = '.'), last_plot(), width = 9, height = 9)









