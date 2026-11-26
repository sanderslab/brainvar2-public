setwd("/Users/nrohani/Library/CloudStorage/Dropbox-Izzy.barraclough/Narjes Rohani/BrainVar_NRohani/Figure_data/Fig2")
rm(list=ls()); gc()
options(stringsAsFactors = F)

############################################
# 1. Load Required Packages
############################################
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(grid)
library(readxl)
require(data.table)

############################################
# 2. Load Data
############################################
file_path <- "BrainVar-LDSC-H2-Selected-Results-2026-04-06.xlsx"

df_prop <- as.data.frame(read_excel(file_path, sheet = "proportion.h2"))
colnames(df_prop)[1] <- "Traits"; rownames(df_prop) <- df_prop[,"Traits"]
df_se   <- as.data.frame(read_excel(file_path, sheet = "proportion.se"))
colnames(df_se)[1]   <- "Traits"; rownames(df_se)   <- df_se[,"Traits"]
df_z    <- as.data.frame(read_excel(file_path, sheet = "Z.regr.coef"))
colnames(df_z)[1]    <- "Traits"; rownames(df_z)    <- df_z[,"Traits"]

############################################
# 3. Cell-type columns — v2 interleaved order
############################################
z_columns <- c(
  "ExN_pre","ExN_post",
  "InN_pre","InN_post",
  "Micro_pre","Micro_post",
  "NNC_pre","NNC_post",
  "Astrocyte_post","Oligo_post","OPC_post"
)

############################################
# 4. Trait order
############################################
desired_order <- c(
  "DD","ASD","ADHD","CEPI","GEPI","TS","AN",
  "PAU","SUD","SCZ","BIP","PANIC","OCD",
  "PTSD","MDD","MS","ALS","FTLD-TDP",
  "PD","ALZ","Height","AOW","IQ","EA"
)

unreliable_traits <- c("DD", "PD", "FTLD-TDP")

df_prop <- df_prop[desired_order, ]
df_se   <- df_se[desired_order, ]
df_z    <- df_z[desired_order, ]

############################################
# 5. Proportion matrix — used for FILL COLOR only
############################################
prop_mat <- as.matrix(df_prop[, c("h2", z_columns)])
rownames(prop_mat) <- df_prop[["Traits"]]
prop_mat[prop_mat < 0] <- 0

############################################
# 6. Z-score matrix — used for SIGNIFICANCE only
#    Exactly as in the original script:
#    - h2 col      : df_prop$Z  (overall heritability Z reported by ldsc)
#    - cell-type cols: taken directly from Z.regr.coef sheet
#    Negative values zeroed out (same as original)
############################################
z_mat <- as.matrix(cbind.data.frame(
  data.frame(h2 = df_prop$Z),
  df_z[, z_columns]
))
z_mat[z_mat < 0] <- 0

############################################
# 7. P-values + Bonferroni
#    N = number of non-negative Z cells (exactly as original)
#    ** : p < 0.01 / N
#    *  : p < 0.05 / N
############################################
pval_mat <- pnorm(-abs(z_mat))   # one-tailed

N       <- sum(z_mat != 0)
bonf_01 <- 0.01 / N
bonf_05 <- 0.05 / N

cat(sprintf("Bonferroni N = %d  |  ** threshold = %.6f  |  * threshold = %.6f\n",
            N, bonf_01, bonf_05))

############################################
# 8. Blue color scale anchored to actual data range
############################################
pos_vals   <- prop_mat[prop_mat > 0]
val_median <- median(pos_vals, na.rm = TRUE)
val_max    <- max(pos_vals,    na.rm = TRUE)

col_fun <- colorRamp2(
  c(0,       0.3,       val_max),
  c("white", "#87CEEB", "#0047AB")
)

############################################
# 9. Row name colours — grey for unreliable traits
############################################
row_name_cols <- ifelse(
  rownames(prop_mat) %in% unreliable_traits,
  "grey60", "black"
)

############################################
# 10. Column / row groups — v2 layout
############################################
column_groups <- factor(
  c("H2",
    "Excitatory","Excitatory",
    "Inhibitory","Inhibitory",
    "Microglia","Microglia",
    "NNC","NNC",
    "Glia","Glia","Glia"),
  levels = c("H2","Excitatory","Inhibitory","Microglia","Glia","NNC")
)

row_groups <- c(rep("grp1",6), rep("grp2",9), rep("grp3",5), rep("grp4",4))

############################################
# 11. Heatmap
############################################
ht <- Heatmap(
  prop_mat,
  name = "Prop. h²",
  col  = col_fun,
  
  rect_gp = gpar(fill = "white", col = "grey90"),
  
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_split    = column_groups,
  row_split       = row_groups,
  
  border    = FALSE,
  border_gp = gpar(col = "grey60", lwd = 1.2),
  
  column_title    = NULL,
  row_title       = NULL,
  
  show_row_names      = TRUE,
  row_names_side      = "left",
  row_names_gp        = gpar(fontsize = 11, col = row_name_cols),
  row_names_max_width = unit(4, "cm"),
  
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 10),
  column_names_rot  = 45,
  
  height = unit(nrow(prop_mat) * 0.65, "cm"),
  width  = unit(ncol(prop_mat) * 0.9,  "cm"),
  
  heatmap_legend_param = list(
    title         = "h² (proportion)",
    legend_height = unit(4, "cm")
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    score <- prop_mat[i, j]
    pval  <- pval_mat[i, j]
    trait <- rownames(prop_mat)[i]
    
    # Greyscale for unreliable rows
    if (trait %in% unreliable_traits) {
      base_fill <- colorRamp2(
        c(0, val_median, val_max),
        c("#F0F0F0", "#BDBDBD", "#969696")
      )(score)
    } else {
      base_fill <- col_fun(score)
    }
    
    # Box size / border driven by nominal p < 0.05
    sq_size    <- ifelse(!is.na(pval) && pval < 0.05, 0.9, 0.7)
    border_col <- ifelse(!is.na(pval) && pval < 0.05, "black", "grey90")
    
    grid.rect(x, y,
              width  = width  * sq_size,
              height = height * sq_size,
              gp     = gpar(fill = base_fill, col = border_col))
    
    # Bonferroni star — single * for all cells with p < 0.05/N
    if (!is.na(pval) && pval < bonf_05) {
      star_col <- ifelse(trait %in% unreliable_traits, "grey40", "black")
      grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold", col = star_col))
    }
  }
)

############################################
# 12. Legends
############################################
legend_pval <- Legend(
  title     = "p-value < 0.05",
  labels    = c("Significant", "Not significant"),
  type      = "points",
  pch       = 15,
  size      = unit(c(6, 3), "mm"),
  legend_gp = gpar(fill = c("grey50","grey50"), col = c("black","grey90"))
)

legend_star <- Legend(
  title     = sprintf("Bonferroni (N = %d)", N),
  labels    = c("* p < 0.05/N"),
  type      = "points",
  pch       = "*",
  size      = unit(4, "mm"),
  legend_gp = gpar(fontsize = 12, col = "black")
)

legend_unreliable <- Legend(
  title     = "Reliability",
  labels    = c("Unreliable (DD, PD, FTLD-TDP)", "Reliable"),
  type      = "points",
  pch       = 15,
  size      = unit(5, "mm"),
  legend_gp = gpar(fill = c("grey70","#6BAED6"), col = "grey60")
)

############################################
# 13. Draw + export
############################################
draw_heatmap <- function() {
  draw(
    ht,
    annotation_legend_list = list(legend_pval, legend_star, legend_unreliable),
    heatmap_legend_side    = "right",
    annotation_legend_side = "bottom",
    padding                = unit(c(12, 15, 15, 10), "mm")
  )
}

draw_heatmap()

pdf("LDSC_heatmap_blue_v7.pdf", width = 10, height = 14)
draw_heatmap()
dev.off()

cat("Saved: LDSC_heatmap_blue_v7.pdf\n")