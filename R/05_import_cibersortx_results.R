# Import CIBERSORTx LM22 output, apply sample-level QC, and merge metadata.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tidyr", "tibble"))

cibersortx_path <- file.path(PATHS$public, "CIBERSORTx_Job8_Adjusted.tsv")
metadata_path <- file.path(PATHS$processed, "metadata_cibersortx_clean.tsv")

if (!file.exists(metadata_path)) {
  stop("Run R/02_prepare_metadata.R before this script.")
}

cibersortx <- read_table_auto(cibersortx_path)
metadata <- read_table_auto(metadata_path)

assert_required_columns(
  cibersortx,
  c("Mixture", LM22_CELL_TYPES, CIBERSORTX_QC_COLUMNS),
  "CIBERSORTx output"
)
assert_required_columns(
  metadata,
  c("Sample", "Status", "Sex", "Batch", "Group"),
  "metadata"
)

colnames(cibersortx)[colnames(cibersortx) == "Mixture"] <- "Sample"
assert_unique_values(cibersortx$Sample, "CIBERSORTx Sample")
assert_unique_values(metadata$Sample, "metadata Sample")

fraction_matrix <- as.matrix(cibersortx[, LM22_CELL_TYPES])
storage.mode(fraction_matrix) <- "numeric"

fraction_sums <- rowSums(fraction_matrix, na.rm = TRUE)
if (any(abs(fraction_sums - 1) > 1e-6)) {
  warning("Some LM22 fraction rows do not sum to 1 within tolerance.")
}

cibersortx <- cibersortx |>
  dplyr::mutate(
    Fraction_sum = fraction_sums,
    CIBERSORTx_QC_pass = .data[["P-value"]] < CIBERSORTX_PVALUE_CUTOFF
  )

integrated <- cibersortx |>
  dplyr::inner_join(metadata, by = "Sample")

missing_metadata <- setdiff(cibersortx$Sample, metadata$Sample)
missing_cibersortx <- setdiff(metadata$Sample, cibersortx$Sample)

if (length(missing_metadata) > 0) {
  stop("CIBERSORTx samples without metadata: ",
       paste(missing_metadata, collapse = ", "))
}
if (length(missing_cibersortx) > 0) {
  stop("Metadata samples without CIBERSORTx output: ",
       paste(missing_cibersortx, collapse = ", "))
}

integrated_qc <- integrated |>
  dplyr::filter(CIBERSORTx_QC_pass)

long <- integrated_qc |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(LM22_CELL_TYPES),
    names_to = "Cell",
    values_to = "Fraction"
  )

write_tsv_safe(
  integrated,
  file.path(PATHS$processed, "cibersortx_metadata_integrated_all.tsv")
)
write_tsv_safe(
  integrated_qc,
  file.path(PATHS$processed, "cibersortx_metadata_integrated_qc_pass.tsv")
)
write_tsv_safe(
  long,
  file.path(PATHS$processed, "cibersortx_metadata_long_qc_pass.tsv")
)

qc_summary <- data.frame(
  Metric = c(
    "Samples_in_CIBERSORTx_output",
    "Samples_passing_CIBERSORTx_P_value",
    "Minimum_fraction_sum",
    "Maximum_fraction_sum",
    "Median_correlation",
    "Median_RMSE"
  ),
  Value = c(
    nrow(integrated),
    nrow(integrated_qc),
    min(integrated$Fraction_sum),
    max(integrated$Fraction_sum),
    median(integrated$Correlation),
    median(integrated$RMSE)
  )
)
write_tsv_safe(qc_summary, file.path(PATHS$results, "cibersortx_qc_summary.tsv"))

message("Integrated ", nrow(integrated_qc),
        " QC-passing samples with metadata.")
