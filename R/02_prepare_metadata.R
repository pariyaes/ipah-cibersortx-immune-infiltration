# Clean and validate sample metadata for the CIBERSORTx branch.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tibble"))

input_path <- file.path(PATHS$public, "sample_metadata_cibersortx.tsv")
output_path <- file.path(PATHS$processed, "metadata_cibersortx_clean.tsv")

metadata <- read_table_auto(input_path)

rename_map <- c(
  Gender = "Sex",
  gender = "Sex",
  batch = "Batch",
  sample = "Sample",
  status = "Status"
)

for (old_name in names(rename_map)) {
  if (old_name %in% colnames(metadata) &&
      !rename_map[[old_name]] %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == old_name] <- rename_map[[old_name]]
  }
}

assert_required_columns(
  metadata,
  c("Sample", "Status", "Sex", "Batch"),
  "metadata"
)
assert_unique_values(metadata$Sample, "metadata$Sample")

metadata <- metadata |>
  dplyr::mutate(
    Sample = as.character(Sample),
    Status = factor(Status, levels = c("Control", "IPAH")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    Batch = factor(Batch),
    Group = factor(
      paste(as.character(Status), as.character(Sex), sep = "_"),
      levels = GROUP_LEVELS
    ),
    Group_Batch = paste(Group, Batch, sep = "_Batch")
  ) |>
  dplyr::arrange(Status, Sex, Batch, Sample)

if (anyNA(metadata$Status)) stop("Unexpected Status value in metadata.")
if (anyNA(metadata$Sex)) stop("Unexpected Sex value in metadata.")

write_tsv_safe(metadata, output_path)

summary_table <- metadata |>
  dplyr::count(Status, Sex, Batch, name = "n")
write_tsv_safe(
  summary_table,
  file.path(PATHS$results, "metadata_group_counts.tsv")
)

message("Saved cleaned metadata: ", output_path)
print(summary_table)
