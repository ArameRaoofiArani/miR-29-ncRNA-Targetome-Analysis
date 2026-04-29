library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)

# Load correlation matrices
cor_hacat <- readRDS("cor_hacat.rds")
cor_tigk  <- readRDS("cor_tigk.rds")

# Define functional categories used in the report
gene_categories <- list(
  ECM_remodeling = c(
    "ADAM12", "LOX", "COL1A1", "FN1", "COL5A1",
    "FBN1", "MMP2", "SPARC", "LAMA2"
  ),
  adhesion_barrier = c(
    "LAMC2", "DSC2", "DSG2", "ITGB1", "ITGA6"
  ),
  growth_factor_signaling = c(
    "WNT4", "PDGFRB", "VEGFA", "SMAD3", "TGFB1", "TGFB2"
  ),
  inflammation_stress = c(
    "TNFAIP3", "REL", "AHR"
  ),
  survival_apoptosis = c(
    "MCL1", "BCL2", "CASP8"
  ),
  epigenetic_regulation = c(
    "HDAC9", "HDAC4", "DNMT3A", "TET2", "SMARCC1"
  )
)

cat_df <- imap_dfr(gene_categories, ~data.frame(
  gene = .x,
  category = .y,
  stringsAsFactors = FALSE
))

category_levels <- c(
  "ECM remodeling",
  "Cell adhesion",
  "Growth factor signaling",
  "Inflammation",
  "Cell survival",
  "Epigenetic regulation"
)

make_bubble_df <- function(cor_mat) {
  as.data.frame(cor_mat) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "lncRNA", values_to = "cor") |>
    mutate(abs_cor = abs(cor)) |>
    left_join(cat_df, by = "gene") |>
    filter(!is.na(category)) |>
    mutate(
      category = recode(
        category,
        ECM_remodeling = "ECM remodeling",
        adhesion_barrier = "Cell adhesion",
        growth_factor_signaling = "Growth factor signaling",
        inflammation_stress = "Inflammation",
        survival_apoptosis = "Cell survival",
        epigenetic_regulation = "Epigenetic regulation"
      )
    )
}

df_hacat <- make_bubble_df(cor_hacat)
df_tigk  <- make_bubble_df(cor_tigk)

# Shared gene order across both plots
combined_df <- bind_rows(
  df_hacat |> mutate(celltype = "HaCaT"),
  df_tigk  |> mutate(celltype = "TIGK")
)

gene_order <- combined_df |>
  group_by(category, gene) |>
  summarise(mean_abs_cor = mean(abs_cor, na.rm = TRUE), .groups = "drop") |>
  mutate(category = factor(category, levels = category_levels)) |>
  arrange(category, desc(mean_abs_cor)) |>
  pull(gene)

format_bubble_df <- function(df) {
  df |>
    mutate(
      lncRNA = factor(lncRNA, levels = c("H19", "LINC00511", "SNHG7")),
      category = factor(category, levels = category_levels),
      gene = factor(gene, levels = rev(gene_order))
    )
}

df_hacat <- format_bubble_df(df_hacat)
df_tigk  <- format_bubble_df(df_tigk)

plot_bubble <- function(df, title_text) {
  ggplot(df, aes(x = lncRNA, y = gene, size = abs_cor, color = cor)) +
    geom_point(alpha = 0.95) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    scale_size(range = c(2.5, 7)) +
    scale_x_discrete(expand = c(0.15, 0)) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0, size = 11, face = "bold"),
      plot.title.position = "plot",
      strip.placement = "outside",
      strip.background = element_rect(fill = "#E8E8E8", color = NA),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 10),
      axis.title = element_blank(),
      panel.spacing.y = grid::unit(0.8, "lines")
    ) +
    labs(
      title = title_text,
      color = "Pearson r",
      size = "|r|"
    )
}

p_hacat <- plot_bubble(
  df_hacat,
  "Candidate lncRNA correlations with miR-29-associated genes in HaCaT"
)

p_tigk <- plot_bubble(
  df_tigk,
  "Candidate lncRNA correlations with miR-29-associated genes in TIGK"
)

ggsave("bubbleplot_HaCaT.pdf", plot = p_hacat, width = 8, height = 9)
ggsave("bubbleplot_TIGK.pdf", plot = p_tigk, width = 8, height = 9)

ggsave("bubbleplot_HaCaT.png", plot = p_hacat, width = 8, height = 9, dpi = 300)
ggsave("bubbleplot_TIGK.png", plot = p_tigk, width = 8, height = 9, dpi = 300)
