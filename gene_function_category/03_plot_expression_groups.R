library(ggplot2)
library(dplyr)
library(gghalves)
library(ggbeeswarm)

df <- read.csv("/Users/nrohani/Library/CloudStorage/Dropbox-Izzy.barraclough/Narjes Rohani/BrainVar_NRohani/Figure_data/Fig2/violin_data_for_R_v2.csv")
df$cell_type  <- factor(df$cell_type,  levels = c("ExN", "InN", "Glia", "NNC"))
df$period     <- factor(df$period,     levels = c("Prenatal", "Postnatal"))
df$gene_group <- factor(df$gene_group, levels = c("Transcriptional Regulation", "Neuronal Communication"))

group_colors <- c(
  "Transcriptional Regulation"        = "#C0392B",
  "Neuronal Communication" = "#1A3A6B"
)

y_min <- floor(min(df$log_expr, na.rm = TRUE)) - 1
y_max <- ceiling(max(df$log_expr, na.rm = TRUE)) + 1.5

# --- Compute p-values per facet (cell_type × period) ---
pval_df <- df %>%
  group_by(cell_type, period) %>%
  summarise(
    p_value = wilcox.test(
      log_expr[gene_group == "Transcriptional Regulation"],
      log_expr[gene_group == "Neuronal Communication"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_stars = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    p_numeric = ifelse(
      p_value < 0.001,
      paste0("p=", formatC(p_value, format = "e", digits = 2)),
      paste0("p=", signif(p_value, 3))
    ),
    y_pos   = y_max - 0.6,
    x_start = 1,
    x_end   = 2,
    x_label = 1.5
  )

p <- ggplot(df, aes(x = gene_group, y = log_expr, fill = gene_group, color = gene_group)) +
  
  geom_violin(
    trim      = FALSE,
    alpha     = 0.40,
    linewidth = 0.4,
    scale     = "width"
  ) +
  
  stat_summary(
    fun       = median,
    geom      = "crossbar",
    width     = 0.4,
    linewidth = 0.6,
    color     = "black",
    fatten    = 0
  ) +
  
  geom_beeswarm(
    size        = 0.4,   # <-- smaller dot size
    alpha       = 0.70,
    cex         = 0.8,   # <-- tighter spacing between dots
    priority    = "random",
    shape       = 16,
    color       = "black",
    stroke      = 0,
    show.legend = FALSE
  ) +
  
  # --- Left tick ---
  geom_segment(
    data        = pval_df,
    aes(x = x_start, xend = x_start, y = y_pos - 0.15, yend = y_pos),
    inherit.aes = FALSE, color = "black", linewidth = 0.4
  ) +
  # --- Right tick ---
  geom_segment(
    data        = pval_df,
    aes(x = x_end, xend = x_end, y = y_pos - 0.15, yend = y_pos),
    inherit.aes = FALSE, color = "black", linewidth = 0.4
  ) +
  # --- Horizontal bracket bar ---
  geom_segment(
    data        = pval_df,
    aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
    inherit.aes = FALSE, color = "black", linewidth = 0.4
  ) +
  
  # --- Stars above bracket ---
  geom_text(
    data        = pval_df,
    aes(x = x_label, y = y_pos + 0.10, label = p_stars),
    inherit.aes = FALSE,
    size        = 3.5,
    fontface    = "bold",
    color       = "black"
  ) +
  
  # --- Exact p-value above stars ---
  geom_text(
    data        = pval_df,
    aes(x = x_label, y = y_pos + 0.45, label = p_numeric),
    inherit.aes = FALSE,
    size        = 2.5,
    color       = "grey30"
  ) +
  
  facet_grid(cell_type ~ period, scales = "free_x") +
  
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(ceiling(y_min), floor(y_max), by = 1)
  ) +
  
  labs(
    title = "Gene Expression: Prenatal vs Postnatal",
    x     = NULL,
    y     = "log\u2082(CPM + 1)",
    fill  = "Gene Group",
    color = "Gene Group"
  ) +
  
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "#F0F0F0", color = "grey70"),
    strip.text         = element_text(face = "bold", size = 10),
    axis.text.x        = element_text(angle = 30, hjust = 1, size = 9),
    axis.text.y        = element_text(size = 9),
    axis.title.y       = element_text(size = 10, margin = margin(r = 8)),
    legend.position    = "top",
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 10),
    legend.key.size    = unit(0.9, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    panel.spacing.x    = unit(0.8, "lines"),
    panel.spacing.y    = unit(0.5, "lines"),
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5,
                                      margin = margin(b = 10))
  )

# --- Save PDF ---
pdf_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/output_plots/violin_GR_NC_prenatal_postnatal.pdf"

tryCatch({
  ggsave(
    pdf_path,
    plot   = p,
    width  = 3,
    height = 7  )
  cat("✓ PDF saved to:", pdf_path, "\n")
}, error = function(e) {
  cat("✗ cairo_pdf failed, trying pdf() fallback...\n")
  pdf(pdf_path, width = 8, height = 14)
  print(p)
  dev.off()
  cat("✓ PDF saved via fallback to:", pdf_path, "\n")
})

# --- Save PNG ---
png_path <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/output_plots/violin_GR_NC_prenatal_postnatal.png"

ggsave(
  png_path,
  plot   = p,
  width  = 8,
  height = 14,
  dpi    = 200
)
cat("✓ PNG saved to:", png_path, "\n")