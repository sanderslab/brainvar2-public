
# ============================================================
# SOX5 GRN — SAFE VERSION
#
# Edge colour = Cell type
# Node fill   = P-value
# Node size   = Score (median)
#
# NODE ORDER:
#   1. Same cell types stay together
#   2. Cell-type groups with larger nodes come first
#   3. Within each cell type:
#      largest node -> smallest node
#
# FONT:
#   Standard sans font
#   7 pt throughout
#
# OUTPUT:
#   PDF only
#   No ggsave()
#   No Cairo
#   No showtext
#   No custom fonts
# ============================================================


# ── Libraries ────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)
library(patchwork)
library(scales)
library(grid)


# ============================================================
# 1. LOAD DATA
# ============================================================

df_all <- read_excel(
  "/Users/nrohani/Downloads/phase_grn_summary_all_celltypes-2.xlsx",
  sheet = "Combined"
)


# ============================================================
# 2. FILTER SOX5
# ============================================================

df_all <- df_all %>%
  filter(TF == "SOX5") %>%
  rename(
    tf        = TF,
    pval      = `P-value`,
    phase     = `Development stage`,
    cell_type = `Cell type`
  ) %>%
  mutate(
    median = as.numeric(
      str_extract(
        Score,
        "(?<=median=)[-0-9\\.]+"
      )
    ),
    
    pval = as.numeric(pval),
    
    cell_type = as.character(cell_type),
    
    phase = tolower(
      trimws(as.character(phase))
    )
  ) %>%
  filter(
    !is.na(median),
    !is.na(pval),
    !is.na(cell_type),
    !is.na(`Target gene`)
  )


# ============================================================
# 3. CHECK DATA
# ============================================================

cat("\nNumber of SOX5 rows:", nrow(df_all), "\n\n")

cat("Developmental stages:\n")
print(table(df_all$phase))

cat("\nCell types:\n")
print(sort(unique(df_all$cell_type)))


# Stop if no SOX5 data were found

if (nrow(df_all) == 0) {
  stop("No SOX5 data were found.")
}


# ============================================================
# 4. CELL-TYPE COLOURS
# ============================================================
#
# One colour per cell type.
# Same colour is used in prenatal and postnatal.
#

cell_types <- sort(
  unique(df_all$cell_type)
)


celltype_colours <- setNames(
  hue_pal(
    h = c(15, 375),
    c = 100,
    l = 55
  )(
    length(cell_types)
  ),
  
  cell_types
)


# ============================================================
# 5. P-VALUE COLOURS
# ============================================================

node_colours <- c(
  "< 1e-10" = "#C0392B",
  "< 1e-6"  = "#E8A09A",
  "< 0.001" = "#A0A0A0",
  "< 0.01"  = "#C8C8C8",
  "< 0.05"  = "#E8E8E8"
)


# ============================================================
# 6. GRN FUNCTION
# ============================================================

make_grn <- function(
    df,
    title_label,
    center_colour,
    celltype_colours
) {
  
  
  # ----------------------------------------------------------
  # Safety check
  # ----------------------------------------------------------
  
  if (nrow(df) == 0) {
    stop(
      paste(
        "No data found for:",
        title_label
      )
    )
  }
  
  
  # ==========================================================
  # DETERMINE ORDER OF CELL-TYPE GROUPS
  # ==========================================================
  #
  # Cell type with the largest maximum node is placed first.
  #
  # This means cell types containing the biggest circles
  # appear earlier around the network.
  # ==========================================================
  
  cell_order <- df %>%
    group_by(cell_type) %>%
    summarise(
      max_node = max(
        median,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    arrange(
      desc(max_node)
    ) %>%
    pull(cell_type)
  
  
  # ==========================================================
  # SORT NODES
  # ==========================================================
  #
  # First:
  #   group by cell type
  #
  # Then within each cell type:
  #   largest circle -> smallest circle
  # ==========================================================
  
  df <- df %>%
    mutate(
      cell_type = factor(
        cell_type,
        levels = cell_order
      )
    ) %>%
    arrange(
      cell_type,
      desc(median),
      pval,
      `Target gene`
    )
  
  
  # Rank after sorting
  
  df$rank <- seq_len(
    nrow(df)
  )
  
  
  n <- nrow(df)
  
  
  # ==========================================================
  # CIRCULAR LAYOUT
  # ==========================================================
  #
  # First node starts at 12 o'clock.
  # Nodes proceed clockwise.
  # ==========================================================
  
  df$theta <-
    pi / 2 -
    2 * pi * (df$rank - 1) / n
  
  
  df$x <- cos(df$theta)
  
  df$y <- sin(df$theta)
  
  
  # ==========================================================
  # P-VALUE CATEGORIES
  # ==========================================================
  
  df$col_tier <- cut(
    df$pval,
    
    breaks = c(
      -Inf,
      1e-10,
      1e-6,
      1e-3,
      1e-2,
      Inf
    ),
    
    labels = c(
      "< 1e-10",
      "< 1e-6",
      "< 0.001",
      "< 0.01",
      "< 0.05"
    )
  )
  
  
  # ==========================================================
  # LABEL POSITIONS
  # ==========================================================
  
  df$label_x <- df$x * 1.18
  
  df$label_y <- df$y * 1.18
  
  
  # Right side -> text extends right
  # Left side  -> text extends left
  
  df$label_hjust <- ifelse(
    df$x >= 0,
    0,
    1
  )
  
  
  # ==========================================================
  # EDGE DATA
  # ==========================================================
  
  edges <- data.frame(
    x         = 0,
    y         = 0,
    xend      = df$x,
    yend      = df$y,
    cell_type = df$cell_type
  )
  
  
  # ==========================================================
  # CREATE PLOT
  # ==========================================================
  
  p <- ggplot() +
    
    
    # --------------------------------------------------------
  # EDGES
  #
  # Edge colour = cell type
  # --------------------------------------------------------
  
  geom_segment(
    data = edges,
    
    aes(
      x      = x,
      y      = y,
      xend   = xend,
      yend   = yend,
      colour = cell_type
    ),
    
    linewidth = 0.5,
    alpha     = 0.75,
    lineend   = "round"
  ) +
    
    
    # --------------------------------------------------------
  # TARGET NODES
  #
  # Size = median Score
  # Fill = P-value
  # --------------------------------------------------------
  
  geom_point(
    data = df,
    
    aes(
      x    = x,
      y    = y,
      size = median,
      fill = col_tier
    ),
    
    shape = 21,
    
    colour = "white",
    
    stroke = 0.3
  ) +
    
    
    # --------------------------------------------------------
  # TARGET GENE LABELS
  #
  # 7 pt
  # --------------------------------------------------------
  
  geom_text(
    data = df,
    
    aes(
      x     = label_x,
      y     = label_y,
      label = `Target gene`,
      hjust = label_hjust
    ),
    
    family = "sans",
    
    size = 7 / 2.845,
    
    fontface = "plain",
    
    colour = "#1A1A1A"
  ) +
    
    
    # --------------------------------------------------------
  # CENTER SOX5 NODE
  # --------------------------------------------------------
  
  geom_point(
    aes(
      x = 0,
      y = 0
    ),
    
    size = 10,
    
    shape = 21,
    
    fill = center_colour,
    
    colour = "white",
    
    stroke = 0.5
  ) +
    
    
    # --------------------------------------------------------
  # SOX5 LABEL
  #
  # 7 pt
  # --------------------------------------------------------
  
  annotate(
    "text",
    
    x = 0,
    y = 0,
    
    label = "SOX5",
    
    family = "sans",
    
    size = 7 / 2.845,
    
    fontface = "bold",
    
    colour = "white"
  ) +
    
    
    # ========================================================
  # EDGE COLOUR LEGEND
  #
  # Cell type
  # ========================================================
  
  scale_colour_manual(
    values = celltype_colours,
    
    name = "Cell type",
    
    drop = FALSE,
    
    guide = guide_legend(
      order = 3,
      
      title.position = "top",
      
      title.hjust = 0,
      
      override.aes = list(
        linewidth = 1.5,
        alpha = 1
      )
    )
  ) +
    
    
    # ========================================================
  # NODE FILL LEGEND
  #
  # P-value
  # ========================================================
  
  scale_fill_manual(
    values = node_colours,
    
    name = "P",
    
    drop = FALSE,
    
    guide = guide_legend(
      order = 1,
      
      override.aes = list(
        size = 3,
        shape = 21,
        stroke = 0.3
      ),
      
      title.position = "top",
      
      title.hjust = 0
    )
  ) +
    
    
    # ========================================================
  # NODE SIZE LEGEND
  #
  # Score
  # ========================================================
  
  scale_size_continuous(
    name = "Score",
    
    range = c(
      1,
      11
    ),
    
    breaks = c(
      0.01,
      0.04,
      0.08,
      0.15
    ),
    
    guide = guide_legend(
      order = 2,
      
      override.aes = list(
        fill = "grey60",
        colour = "white",
        shape = 21,
        stroke = 0.3
      ),
      
      title.position = "top",
      
      title.hjust = 0
    )
  ) +
    
    
    # ========================================================
  # COORDINATES
  # ========================================================
  
  coord_fixed(
    xlim = c(
      -1.70,
      1.70
    ),
    
    ylim = c(
      -1.70,
      1.70
    ),
    
    clip = "off"
  ) +
    
    
    # ========================================================
  # TITLE
  # ========================================================
  
  labs(
    title = title_label
  ) +
    
    
    # ========================================================
  # THEME
  #
  # Standard sans font
  # 7 pt throughout
  # ========================================================
  
  theme_void(
    base_family = "sans",
    base_size = 7
  ) +
    
    theme(
      
      text = element_text(
        family = "sans",
        size = 7
      ),
      
      
      # Title
      plot.title = element_text(
        family = "sans",
        size = 7,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 6
        )
      ),
      
      
      # Legend titles
      legend.title = element_text(
        family = "sans",
        size = 7,
        face = "bold"
      ),
      
      
      # Legend text
      legend.text = element_text(
        family = "sans",
        size = 7
      ),
      
      
      legend.position = "right",
      
      legend.direction = "vertical",
      
      legend.box = "vertical",
      
      legend.spacing.y = unit(
        0.25,
        "cm"
      ),
      
      legend.key.size = unit(
        0.4,
        "cm"
      ),
      
      
      # Extra space so labels are not cut off
      plot.margin = margin(
        t = 15,
        r = 25,
        b = 15,
        l = 25
      )
    )
  
  
  return(p)
}


# ============================================================
# 7. SPLIT PRENATAL / POSTNATAL
# ============================================================

df_pre <- df_all %>%
  filter(
    phase == "prenatal"
  )


df_post <- df_all %>%
  filter(
    phase == "postnatal"
  )


cat(
  "\nPrenatal links:",
  nrow(df_pre),
  "\n"
)

cat(
  "Postnatal links:",
  nrow(df_post),
  "\n\n"
)


# Safety checks

if (nrow(df_pre) == 0) {
  stop("No prenatal SOX5 data found.")
}


if (nrow(df_post) == 0) {
  stop("No postnatal SOX5 data found.")
}


# ============================================================
# 8. PRENATAL PLOT
# ============================================================

p_pre <- make_grn(
  df = df_pre,
  
  title_label = "GRN for SOX5 | Prenatal",
  
  center_colour = "#185FA5",
  
  celltype_colours = celltype_colours
)


# ============================================================
# 9. POSTNATAL PLOT
# ============================================================

p_post <- make_grn(
  df = df_post,
  
  title_label = "GRN for SOX5 | Postnatal",
  
  center_colour = "#0F6E56",
  
  celltype_colours = celltype_colours
)


# ============================================================
# 10. COMBINE PLOTS
# ============================================================

combined <- p_pre + p_post +
  
  plot_layout(
    ncol = 2,
    guides = "collect"
  ) &
  
  theme(
    text = element_text(
      family = "sans",
      size = 7
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      family = "sans",
      size = 7,
      face = "bold"
    ),
    
    legend.text = element_text(
      family = "sans",
      size = 7
    )
  )


# ============================================================
# 11. SHOW PLOT IN R
# ============================================================

print(combined)


# ============================================================
# 12. SAVE PDF
#
# Standard R PDF device.
#
# NO:
#   ggsave()
#   cairo_pdf()
#   showtext
#   sysfonts
#   quartz
#   custom fonts
# ============================================================

pdf(
  file = "SOX5_GRN_prenatal_postnatal.pdf",
  
  width = 12,
  
  height = 6,
  
  family = "Helvetica",
  
  useDingbats = FALSE,
  
  onefile = TRUE
)


# IMPORTANT:
# Explicit print is required when writing ggplot/patchwork
# objects to a graphics device.

print(combined)


# Close PDF device

dev.off()


# ============================================================
# 13. CONFIRM OUTPUT
# ============================================================

message(
  "Successfully saved: SOX5_GRN_prenatal_postnatal.pdf"
)
