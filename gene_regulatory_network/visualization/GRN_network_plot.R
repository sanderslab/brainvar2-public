# ============================================================
# AHI1 GRN — ExN (prenatal)
# From: causal_variant_tf_summary.xlsx
# ============================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)

# ── 1. Load Excel ───────────────────────────────────────────
df <- read_excel("/Users/nrohani/Downloads/causal_variant_tf_summary.xlsx")
df <- read_excel("/Users/nrohani/Downloads/AHI1_glia_micro_significant_tf_summary.xlsx")

# ── 2. Filter your data ─────────────────────────────────────
df <- df %>%
  filter(
    `target gene` == "AHI1",
    loci == "6-135497759-C-A",
    `cell type` == "Micro",
    phase == "postnatal"
  )

# ── 3. Extract median score ─────────────────────────────────
df$median <- abs(as.numeric(
  str_extract(df$`mean/median score`, "(?<=median=)[-0-9\\.]+")
))

# Rename for compatibility
df <- df %>%
  rename(
    tf   = TF,
    pval = `p-values`
  )
df$tf <- toupper(df$tf)
df$pval <- as.numeric(df$pval)

# Remove NA just in case
df <- df %>%
  filter(!is.na(median), !is.na(pval))

# Sort by significance
df <- df %>%
  arrange(pval)

df$rank <- seq_len(nrow(df))

# ── 4. Layout (sunflower spiral) ────────────────────────────
golden <- pi * (3 - sqrt(5))
r_min  <- 0.18
r_max  <- 0.92
n      <- nrow(df)

if (n == 1) {
  df$r     <- 0.5
  df$theta <- 0
} else {
  df$r     <- r_min + ((df$rank - 1) / (n - 1)) * (r_max - r_min)
  df$theta <- (df$rank - 1) * golden
}

df$x <- df$r * cos(df$theta)
df$y <- df$r * sin(df$theta)

# ── 5. Aesthetics ───────────────────────────────────────────
df$col_tier <- cut(
  df$pval,
  breaks = c(-Inf, 1e-10, 1e-6, 1e-3, 1e-2, Inf),
  labels = c("p < 1e-10", "p < 1e-6", "p < 0.001", "p < 0.01", "p < 0.05")
)

node_colours <- c(
  "p < 1e-10" = "#C0392B",
  "p < 1e-6"  = "#E8A09A",
  "p < 0.001" = "#A0A0A0",
  "p < 0.01"  = "#C8C8C8",
  "p < 0.05"  = "#E8E8E8"
)

# Bigger median score = bigger node
if (length(unique(df$median)) > 1) {
  df$node_size <- scales::rescale(df$median, to = c(2, 10))
} else {
  df$node_size <- rep(6, nrow(df))
}

df$bold       <- ifelse(df$rank <= 6, "bold", "plain")
df$label_size <- ifelse(df$rank <= 6, 2.2, 1.7)

# ── 6. Build edges ──────────────────────────────────────────
edges <- data.frame(
  x = 0,
  y = 0,
  xend = df$x,
  yend = df$y
)

# ── 7. Plot ─────────────────────────────────────────────────
p <- ggplot() +
  
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "black",
    linewidth = 0.3,
    lineend = "round"
  ) +
  
  geom_point(
    data = df,
    aes(x = x, y = y, size = node_size, fill = col_tier),
    shape = 21,
    colour = "white",
    stroke = 0.3
  ) +
  
  geom_text(
    data = df,
    aes(x = x, y = y, label = tf, fontface = bold),
    size = df$label_size,
    colour = "#1a1a1a"
  ) +
  
  # Center node
  geom_point(
    aes(x = 0, y = 0),
    size = 8,
    shape = 21,
    fill = "#185FA5",
    colour = "white"
  ) +
  
  annotate(
    "text",
    x = 0, y = 0,
    label = "AHI1",
    size = 3.2,
    fontface = "bold",
    colour = "white"
  ) +
  
  scale_fill_manual(values = node_colours) +
  scale_size_identity(guide = "none") +
  
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  
  labs(title = "GRN for AHI1 — Microglia (postnatal)") +
  
  theme_void()

# ── 8. Save ─────────────────────────────────────────────────
ggsave("AHI1_micro_postnatal_GRN.pdf", p, width = 5, height = 5)

message("Saved: AHI1_micro_postnatal_GRN.pdf")

