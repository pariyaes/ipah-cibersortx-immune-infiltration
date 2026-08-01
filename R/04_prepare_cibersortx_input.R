# Validate and export the CIBERSORTx mixture file.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tibble"))

input_path <- file.path(PATHS$processed, "TPM_for_CIBERSORTx.tsv")
output_path <- file.path(PATHS$processed, "CIBERSORTx_mixture_LM22.tsv")

tpm <- read_table_auto(input_path)
gene_column <- intersect(c("GeneSymbol", "gene_symbol", "Gene", "Symbol"), colnames(tpm))[1]
if (is.na(gene_column)) stop("No gene-symbol column found in TPM matrix.")

gene_symbols <- as.character(tpm[[gene_column]])
expression <- as.matrix(tpm[, setdiff(colnames(tpm), gene_column)])
storage.mode(expression) <- "numeric"

valid_gene <- !is.na(gene_symbols) & gene_symbols != ""
expression <- expression[valid_gene, , drop = FALSE]
gene_symbols <- gene_symbols[valid_gene]

if (any(expression < 0, na.rm = TRUE)) {
  stop("TPM matrix contains negative values.")
}
if (any(!is.finite(expression))) {
  stop("TPM matrix contains missing or non-finite values.")
}

if (anyDuplicated(gene_symbols)) {
  message("Summing TPM values for duplicated gene symbols.")
  expression <- rowsum(expression, group = gene_symbols, reorder = TRUE)
  gene_symbols <- rownames(expression)
}

nonzero <- rowSums(expression) > 0
expression <- expression[nonzero, , drop = FALSE]
gene_symbols <- gene_symbols[nonzero]

output <- data.frame(
  GeneSymbol = gene_symbols,
  expression,
  check.names = FALSE
)

readr::write_tsv(output, output_path, na = "")

qc <- data.frame(
  Metric = c("Genes", "Samples", "Duplicated_symbols", "All_zero_genes_removed"),
  Value = c(
    nrow(expression),
    ncol(expression),
    0,
    sum(!nonzero)
  )
)
write_tsv_safe(qc, file.path(PATHS$results, "cibersortx_input_qc.tsv"))

message("Saved CIBERSORTx mixture file: ", output_path)
message("Upload this file to CIBERSORTx using LM22, 100 permutations,")
message("quantile normalization disabled for RNA-seq, and the project batch-correction setting.")
