# Compute Cliff's delta for each LM22 cell type and comparison.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tidyr", "tibble"))

input_path <- file.path(
  PATHS$processed,
  "cibersortx_metadata_integrated_qc_pass.tsv"
)
data <- read_table_auto(input_path)

run_effect_size <- function(data, comparison) {
  subset_data <- comparison_subset(data, comparison)
  group_variable <- comparison$group_variable

  dplyr::bind_rows(lapply(LM22_CELL_TYPES, function(cell_type) {
    group_1_values <- subset_data[
      subset_data[[group_variable]] == comparison$group_1,
      cell_type,
      drop = TRUE
    ]
    group_2_values <- subset_data[
      subset_data[[group_variable]] == comparison$group_2,
      cell_type,
      drop = TRUE
    ]

    delta <- cliffs_delta(group_1_values, group_2_values)

    data.frame(
      Cell = cell_type,
      Comparison = comparison$name,
      Analysis_layer = comparison$analysis_layer,
      Group_1 = comparison$group_1,
      Group_2 = comparison$group_2,
      Cliffs_delta = delta,
      Absolute_delta = abs(delta),
      Direction = ifelse(
        is.na(delta), NA_character_,
        ifelse(delta > 0, "Higher_in_Group_1",
               ifelse(delta < 0, "Lower_in_Group_1", "No_difference"))
      ),
      Magnitude = cliffs_magnitude(delta),
      stringsAsFactors = FALSE
    )
  }))
}

effect_results <- dplyr::bind_rows(lapply(
  PRIMARY_COMPARISONS,
  function(comparison) run_effect_size(data, comparison)
)) |>
  dplyr::arrange(Analysis_layer, Comparison, dplyr::desc(Absolute_delta))

write_tsv_safe(
  effect_results,
  file.path(PATHS$results, "cliffs_delta_lm22_results.tsv")
)

if (file.exists(file.path(PATHS$results, "wilcoxon_lm22_results.tsv"))) {
  statistical_results <- read_table_auto(
    file.path(PATHS$results, "wilcoxon_lm22_results.tsv")
  )
  integrated_results <- statistical_results |>
    dplyr::left_join(
      effect_results |>
        dplyr::select(
          Cell, Comparison, Cliffs_delta, Absolute_delta,
          Direction, Magnitude
        ),
      by = c("Cell", "Comparison")
    )
  write_tsv_safe(
    integrated_results,
    file.path(PATHS$results, "integrated_statistics_and_effect_sizes.tsv")
  )
}
