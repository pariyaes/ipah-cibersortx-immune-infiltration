# Generate publication-oriented immune-composition and effect-size figures.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap"))

input_path <- file.path(
  PATHS$processed,
  "cibersortx_metadata_integrated_qc_pass.tsv"
)
immune <- read_table_auto(input_path)

immune <- immune |>
  dplyr::mutate(
    Group = factor(Group, levels = GROUP_LEVELS)
  ) |>
  dplyr::arrange(Group, Batch, Sample)

sample_order <- immune$Sample

immune_long <- immune |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(LM22_CELL_TYPES),
    names_to = "Cell",
    values_to = "Fraction"
  ) |>
  dplyr::mutate(
    Sample = factor(Sample, levels = sample_order),
    Cell = factor(Cell, levels = LM22_CELL_TYPES)
  )

# Figure 1: immune composition per sample.
stacked_plot <- ggplot2::ggplot(
  immune_long,
  ggplot2::aes(x = Sample, y = Fraction, fill = Cell)
) +
  ggplot2::geom_col(width = 1) +
  ggplot2::facet_grid(
    ~ Group,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(
    title = "Global immune-cell composition across IPAH and control samples",
    x = NULL,
    y = "Estimated immune-cell fraction",
    fill = "LM22 cell type"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.6, "lines"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(PATHS$figures, "Figure1_Immune_StackedBar_CIBERSORTx.png"),
  stacked_plot,
  width = 16,
  height = 6,
  dpi = 300
)

# Figure 2: IPAH versus Control.
disease_plot <- ggplot2::ggplot(
  immune_long,
  ggplot2::aes(x = Status, y = Fraction, fill = Status)
) +
  ggplot2::geom_violin(trim = FALSE, scale = "width", na.rm = TRUE) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, na.rm = TRUE) +
  ggplot2::facet_wrap(~ Cell, scales = "free_y", ncol = 5) +
  ggplot2::labs(
    title = "Immune-cell infiltration differences between IPAH and controls",
    x = NULL,
    y = "Estimated immune-cell fraction"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(PATHS$figures, "Figure2_IPAH_vs_Control_Violin_Box.png"),
  disease_plot,
  width = 14,
  height = 10,
  dpi = 300
)

# Figure 3: sex comparison within IPAH.
ipah_long <- immune_long |>
  dplyr::filter(Status == "IPAH")

sex_plot <- ggplot2::ggplot(
  ipah_long,
  ggplot2::aes(x = Sex, y = Fraction, fill = Sex)
) +
  ggplot2::geom_violin(trim = FALSE, scale = "width", na.rm = TRUE) +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, na.rm = TRUE) +
  ggplot2::facet_wrap(~ Cell, scales = "free_y", ncol = 5) +
  ggplot2::labs(
    title = "Sex-specific immune-cell infiltration patterns within IPAH",
    x = NULL,
    y = "Estimated immune-cell fraction"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(
  file.path(PATHS$figures, "Figure3_Female_vs_Male_IPAH_Violin_Box.png"),
  sex_plot,
  width = 14,
  height = 10,
  dpi = 300
)

# Figure 4: row-scaled heatmap. Rows = immune-cell types; columns = samples.
immune_matrix <- as.matrix(immune[, LM22_CELL_TYPES])
storage.mode(immune_matrix) <- "numeric"
rownames(immune_matrix) <- immune$Sample
heatmap_matrix <- t(immune_matrix)
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix[!is.finite(heatmap_matrix)] <- 0

annotation_col <- immune |>
  dplyr::select(Status, Sex, Batch)
annotation_col <- as.data.frame(annotation_col)
rownames(annotation_col) <- immune$Sample

pheatmap::pheatmap(
  heatmap_matrix,
  scale = "none",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  main = "Immune-cell infiltration in IPAH and control samples",
  filename = file.path(
    PATHS$figures,
    "Figure4_Immune_Infiltration_Heatmap_CIBERSORTx.png"
  ),
  width = 12,
  height = 8
)

# Figure 5: Cliff's delta ranking.
effect_path <- file.path(PATHS$results, "cliffs_delta_lm22_results.tsv")
if (file.exists(effect_path)) {
  effects <- read_table_auto(effect_path) |>
    dplyr::filter(Analysis_layer == "Primary") |>
    dplyr::mutate(
      Cell = stats::reorder(Cell, Cliffs_delta)
    )

  effect_plot <- ggplot2::ggplot(
    effects,
    ggplot2::aes(x = Cliffs_delta, y = Cell)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = Cliffs_delta, yend = Cell)
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = Magnitude), size = 2.4) +
    ggplot2::facet_wrap(~ Comparison, scales = "free_y") +
    ggplot2::labs(
      title = "Effect-size ranking of LM22 immune-cell fractions",
      x = "Cliff's delta",
      y = NULL,
      shape = "Magnitude"
    ) +
    ggplot2::theme_bw(base_size = 10)

  ggplot2::ggsave(
    file.path(PATHS$figures, "Figure5_Cliffs_Delta_Ranking.png"),
    effect_plot,
    width = 12,
    height = 9,
    dpi = 300
  )
}
