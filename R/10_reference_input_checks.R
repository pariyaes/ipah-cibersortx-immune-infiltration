# Reproduce basic checks on the included public-size inputs.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tibble"))

metadata <- read_table_auto(
  file.path(PATHS$public, "sample_metadata_cibersortx.tsv")
)
cibersortx <- read_table_auto(
  file.path(PATHS$public, "CIBERSORTx_Job8_Adjusted.tsv")
)

assert_required_columns(
  cibersortx,
  c("Mixture", LM22_CELL_TYPES, CIBERSORTX_QC_COLUMNS),
  "CIBERSORTx output"
)

fraction_matrix <- as.matrix(cibersortx[, LM22_CELL_TYPES])
storage.mode(fraction_matrix) <- "numeric"
fraction_sums <- rowSums(fraction_matrix)

checks <- data.frame(
  Check = c(
    "Metadata samples",
    "CIBERSORTx samples",
    "LM22 cell types",
    "Samples with P-value < 0.05",
    "Minimum fraction sum",
    "Maximum fraction sum",
    "Metadata samples absent from CIBERSORTx",
    "CIBERSORTx samples absent from metadata"
  ),
  Value = c(
    nrow(metadata),
    nrow(cibersortx),
    length(LM22_CELL_TYPES),
    sum(cibersortx[["P-value"]] < CIBERSORTX_PVALUE_CUTOFF),
    min(fraction_sums),
    max(fraction_sums),
    length(setdiff(metadata$Sample, cibersortx$Mixture)),
    length(setdiff(cibersortx$Mixture, metadata$Sample))
  )
)

print(checks)
write_tsv_safe(
  checks,
  file.path(PATHS$results, "reference_input_checks.tsv")
)
