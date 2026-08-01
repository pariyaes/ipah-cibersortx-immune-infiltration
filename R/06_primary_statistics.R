# Primary two-sided Wilcoxon comparisons across LM22 cell fractions.
#
# Benjamini-Hochberg correction is applied separately within each
# biological comparison across the 22 LM22 cell types.

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
assert_required_columns(
  data,
  c("Sample", "Status", "Sex", "Batch", LM22_CELL_TYPES),
  "integrated CIBERSORTx data"
)

run_comparison <- function(data, comparison) {
  subset_data <- comparison_subset(data, comparison)
  group_variable <- comparison$group_variable

  rows <- lapply(LM22_CELL_TYPES, function(cell_type) {
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

    data.frame(
      Cell = cell_type,
      Comparison = comparison$name,
      Analysis_layer = comparison$analysis_layer,
      Group_1 = comparison$group_1,
      Group_2 = comparison$group_2,
      Group_1_n = sum(is.finite(group_1_values)),
      Group_2_n = sum(is.finite(group_2_values)),
      Group_1_median = median(group_1_values, na.rm = TRUE),
      Group_2_median = median(group_2_values, na.rm = TRUE),
      Difference_in_medians =
        median(group_1_values, na.rm = TRUE) -
        median(group_2_values, na.rm = TRUE),
      P_value = safe_wilcox(group_1_values, group_2_values),
      stringsAsFactors = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  result$FDR <- stats::p.adjust(result$P_value, method = "BH")
  result
}

all_results <- dplyr::bind_rows(lapply(
  PRIMARY_COMPARISONS,
  function(comparison) run_comparison(data, comparison)
)) |>
  dplyr::arrange(Analysis_layer, Comparison, FDR, P_value)

write_tsv_safe(
  all_results,
  file.path(PATHS$results, "wilcoxon_lm22_results.tsv")
)

significant <- all_results |>
  dplyr::filter(!is.na(FDR), FDR < 0.05)
write_tsv_safe(
  significant,
  file.path(PATHS$results, "wilcoxon_lm22_fdr_lt_0.05.tsv")
)

message("Completed ", length(PRIMARY_COMPARISONS),
        " comparisons across ", length(LM22_CELL_TYPES), " cell types.")
