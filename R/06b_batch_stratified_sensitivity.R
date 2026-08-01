# Batch-stratified sensitivity analysis.
#
# Tests are run independently within each observed batch. The script reports
# within-batch direction and FDR but does not pool batch-specific p-values and
# does not apply ComBat or regress batch out of CIBERSORTx fractions.

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

run_batch_comparison <- function(data, comparison, batch_value) {
  batch_data <- data[data$Batch == batch_value, , drop = FALSE]
  subset_data <- comparison_subset(batch_data, comparison)
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

    if (length(group_1_values) == 0 || length(group_2_values) == 0) {
      return(NULL)
    }

    median_difference <-
      median(group_1_values, na.rm = TRUE) -
      median(group_2_values, na.rm = TRUE)

    data.frame(
      Cell = cell_type,
      Comparison = comparison$name,
      Batch = as.character(batch_value),
      Group_1_n = sum(is.finite(group_1_values)),
      Group_2_n = sum(is.finite(group_2_values)),
      Difference_in_medians = median_difference,
      Direction = ifelse(
        median_difference > 0, "Higher_in_Group_1",
        ifelse(median_difference < 0, "Lower_in_Group_1", "No_difference")
      ),
      P_value = safe_wilcox(group_1_values, group_2_values),
      stringsAsFactors = FALSE
    )
  })

  result <- dplyr::bind_rows(rows)
  if (nrow(result) > 0) {
    result$FDR_within_batch <- stats::p.adjust(result$P_value, method = "BH")
  }
  result
}

batch_results <- list()
for (comparison in PRIMARY_COMPARISONS[1:2]) {
  for (batch_value in unique(data$Batch)) {
    batch_results[[paste(comparison$name, batch_value, sep = "__")]] <-
      run_batch_comparison(data, comparison, batch_value)
  }
}

batch_results <- dplyr::bind_rows(batch_results)

direction_summary <- batch_results |>
  dplyr::group_by(Comparison, Cell) |>
  dplyr::summarise(
    Batches_evaluated = dplyr::n(),
    Positive_directions = sum(Difference_in_medians > 0, na.rm = TRUE),
    Negative_directions = sum(Difference_in_medians < 0, na.rm = TRUE),
    Direction_consistent =
      Positive_directions == Batches_evaluated |
      Negative_directions == Batches_evaluated,
    .groups = "drop"
  )

write_tsv_safe(
  batch_results,
  file.path(PATHS$results, "batch_stratified_wilcoxon_results.tsv")
)
write_tsv_safe(
  direction_summary,
  file.path(PATHS$results, "batch_direction_consistency.tsv")
)
