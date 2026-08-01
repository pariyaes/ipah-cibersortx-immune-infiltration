# Build a genes-by-samples count matrix from featureCounts-style files.
#
# Expected manifest columns:
# Sample, File, Dataset, Status, Sex, Batch, Include_in_CIBERSORTx
#
# Each input file should contain Geneid, Length, and one sample-count column.

if (!file.exists("R/00_config.R")) {
  stop("Run this script from the repository root.")
}
source("R/00_config.R")
source("R/00_utils.R")
assert_packages(c("readr", "dplyr", "tibble"))

manifest_path <- file.path(PATHS$raw, "sample_manifest.tsv")
output_counts <- file.path(PATHS$processed, "final_counts_matrix.csv")
output_lengths <- file.path(PATHS$processed, "gene_lengths_from_featurecounts.tsv")
output_metadata <- file.path(PATHS$processed, "sample_metadata_all.tsv")

manifest <- read_table_auto(manifest_path)
assert_required_columns(
  manifest,
  c("Sample", "File", "Dataset", "Status", "Sex", "Batch"),
  "sample manifest"
)
assert_unique_values(manifest$Sample, "Sample manifest")

annotation_columns <- c(
  "Geneid", "GeneID", "gene_id", "Chr", "Start", "End",
  "Strand", "Length"
)

read_featurecounts_sample <- function(path, sample_id) {
  if (!grepl("^(/|[A-Za-z]:)", path)) {
    path <- file.path(PROJECT_ROOT, path)
  }
  if (!file.exists(path)) stop("Count file not found for ", sample_id, ": ", path)

  data <- utils::read.delim(
    path,
    header = TRUE,
    comment.char = "#",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  gene_column <- intersect(c("Geneid", "GeneID", "gene_id"), colnames(data))[1]
  length_column <- intersect(c("Length", "gene_length", "GeneLength"), colnames(data))[1]

  if (is.na(gene_column)) stop("No gene-ID column found in: ", path)
  if (is.na(length_column)) stop("No gene-length column found in: ", path)

  sample_matches <- colnames(data)[grepl(sample_id, colnames(data), fixed = TRUE)]
  sample_matches <- setdiff(sample_matches, annotation_columns)

  if (length(sample_matches) == 1) {
    count_column <- sample_matches
  } else {
    candidate_columns <- setdiff(colnames(data), annotation_columns)
    numeric_candidates <- candidate_columns[
      vapply(data[candidate_columns], is.numeric, logical(1))
    ]
    if (length(numeric_candidates) == 0) {
      stop("No numeric count column found in: ", path)
    }
    count_column <- tail(numeric_candidates, 1)
  }

  result <- data.frame(
    gene_id = as.character(data[[gene_column]]),
    gene_length_bp = suppressWarnings(as.numeric(data[[length_column]])),
    count = suppressWarnings(as.numeric(data[[count_column]])),
    stringsAsFactors = FALSE
  )

  result <- result[
    !is.na(result$gene_id) & result$gene_id != "" &
      !is.na(result$count) & result$count >= 0,
    ,
    drop = FALSE
  ]

  if (anyDuplicated(result$gene_id)) {
    stop("Duplicated gene IDs found in count file: ", path)
  }

  result
}

count_vectors <- list()
length_vectors <- list()

for (index in seq_len(nrow(manifest))) {
  sample_id <- manifest$Sample[index]
  message("Reading ", sample_id, " (", index, "/", nrow(manifest), ")")
  sample_data <- read_featurecounts_sample(manifest$File[index], sample_id)

  count_vectors[[sample_id]] <- stats::setNames(
    sample_data$count,
    sample_data$gene_id
  )
  length_vectors[[sample_id]] <- stats::setNames(
    sample_data$gene_length_bp,
    sample_data$gene_id
  )
}

all_genes <- Reduce(union, lapply(count_vectors, names))
count_matrix <- matrix(
  0L,
  nrow = length(all_genes),
  ncol = length(count_vectors),
  dimnames = list(all_genes, names(count_vectors))
)

for (sample_id in names(count_vectors)) {
  values <- count_vectors[[sample_id]]
  count_matrix[names(values), sample_id] <- as.integer(round(values))
}

gene_lengths <- vapply(all_genes, function(gene_id) {
  observed <- unique(na.omit(vapply(
    length_vectors,
    function(lengths) unname(lengths[gene_id]),
    numeric(1)
  )))
  if (length(observed) > 1) {
    warning("Inconsistent featureCounts lengths for ", gene_id,
            "; using the first observed value.")
  }
  if (length(observed) == 0) NA_real_ else observed[1]
}, numeric(1))

if (anyNA(gene_lengths)) {
  warning(sum(is.na(gene_lengths)), " genes have no valid gene length.")
}

counts_output <- data.frame(
  Geneid = rownames(count_matrix),
  count_matrix,
  check.names = FALSE
)
utils::write.csv(counts_output, output_counts, row.names = FALSE, quote = TRUE)

length_output <- data.frame(
  gene_id = names(gene_lengths),
  gene_length_bp = unname(gene_lengths)
)
write_tsv_safe(length_output, output_lengths)
write_tsv_safe(manifest, output_metadata)

message("Saved count matrix: ", output_counts)
message("Dimensions: ", nrow(count_matrix), " genes x ",
        ncol(count_matrix), " samples")
