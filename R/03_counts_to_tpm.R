# Convert raw gene counts to TPM using featureCounts gene lengths.
#
# This script calculates TPM at Ensembl-gene level, maps Ensembl IDs to
# HGNC symbols, removes unmapped genes, and sums TPM values for duplicated
# symbols. The result is suitable for CIBERSORTx formatting.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tibble", "AnnotationDbi", "org.Hs.eg.db"))

counts_path <- file.path(PATHS$processed, "final_counts_matrix.csv")
lengths_path <- file.path(PATHS$processed, "gene_lengths_from_featurecounts.tsv")
metadata_path <- file.path(PATHS$processed, "metadata_cibersortx_clean.tsv")
output_path <- file.path(PATHS$processed, "TPM_for_CIBERSORTx.tsv")

if (!file.exists(counts_path)) stop("Missing count matrix: ", counts_path)
if (!file.exists(lengths_path)) stop("Missing gene-length table: ", lengths_path)
if (!file.exists(metadata_path)) stop("Missing cleaned metadata: ", metadata_path)

counts_data <- utils::read.csv(
  counts_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_column <- intersect(c("Geneid", "gene_id", "GeneID"), colnames(counts_data))[1]
if (is.na(gene_column)) stop("Count matrix has no recognized gene-ID column.")

gene_ids <- as.character(counts_data[[gene_column]])
count_matrix <- as.matrix(counts_data[, setdiff(colnames(counts_data), gene_column)])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- gene_ids

metadata <- read_table_auto(metadata_path)
assert_required_columns(metadata, c("Sample"), "metadata")

missing_samples <- setdiff(metadata$Sample, colnames(count_matrix))
if (length(missing_samples) > 0) {
  stop("Count matrix is missing metadata samples: ",
       paste(missing_samples, collapse = ", "))
}
count_matrix <- count_matrix[, metadata$Sample, drop = FALSE]

lengths <- read_table_auto(lengths_path)
assert_required_columns(lengths, c("gene_id", "gene_length_bp"), "gene lengths")
lengths$gene_id <- as.character(lengths$gene_id)
assert_unique_values(lengths$gene_id, "gene-length table")

gene_length_bp <- lengths$gene_length_bp[
  match(rownames(count_matrix), lengths$gene_id)
]

valid <- is.finite(gene_length_bp) & gene_length_bp > 0
if (!all(valid)) {
  message("Removing ", sum(!valid), " genes without valid positive lengths.")
  count_matrix <- count_matrix[valid, , drop = FALSE]
  gene_length_bp <- gene_length_bp[valid]
}

length_kb <- gene_length_bp / 1000
rpk <- sweep(count_matrix, 1, length_kb, "/")
per_million_scaling <- colSums(rpk, na.rm = TRUE) / 1e6

if (any(!is.finite(per_million_scaling) | per_million_scaling <= 0)) {
  stop("One or more samples have an invalid TPM scaling factor.")
}

tpm_ensembl <- sweep(rpk, 2, per_million_scaling, "/")
ensembl_no_version <- strip_ensembl_version(rownames(tpm_ensembl))

gene_symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = ensembl_no_version,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
gene_symbols <- unname(gene_symbols)

keep_symbol <- !is.na(gene_symbols) & gene_symbols != ""
tpm_symbol <- rowsum(
  tpm_ensembl[keep_symbol, , drop = FALSE],
  group = gene_symbols[keep_symbol],
  reorder = TRUE
)

tpm_output <- data.frame(
  GeneSymbol = rownames(tpm_symbol),
  tpm_symbol,
  check.names = FALSE
)
write_tsv_safe(tpm_output, output_path)

column_sums <- colSums(tpm_symbol)
qc <- data.frame(
  Sample = names(column_sums),
  TPM_sum = unname(column_sums),
  Deviation_from_1e6 = unname(column_sums - 1e6)
)
write_tsv_safe(qc, file.path(PATHS$results, "tpm_column_sum_qc.tsv"))

message("Saved TPM matrix: ", output_path)
message("Dimensions: ", nrow(tpm_symbol), " gene symbols x ",
        ncol(tpm_symbol), " samples")
