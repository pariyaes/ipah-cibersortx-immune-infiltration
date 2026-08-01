# Reusable validation and I/O helpers.

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap"
)

assert_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R package(s): ", paste(missing, collapse = ", "),
      "\nInstall them before running the workflow."
    )
  }
}

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  first_line <- readLines(path, n = 1, warn = FALSE)
  delimiter <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
  readr::read_delim(
    path,
    delim = delimiter,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
}

assert_required_columns <- function(data, columns, object_name = "data") {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) {
    stop(object_name, " is missing required column(s): ",
         paste(missing, collapse = ", "))
  }
}

assert_unique_values <- function(values, label) {
  duplicated_values <- unique(values[duplicated(values)])
  if (length(duplicated_values) > 0) {
    stop(label, " contains duplicated value(s): ",
         paste(head(duplicated_values, 10), collapse = ", "))
  }
}

strip_ensembl_version <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

write_tsv_safe <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(data, path, na = "")
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)

  greater <- sum(vapply(x, function(value) sum(value > y), numeric(1)))
  less <- sum(vapply(x, function(value) sum(value < y), numeric(1)))
  (greater - less) / (length(x) * length(y))
}

cliffs_magnitude <- function(delta) {
  if (is.na(delta)) return(NA_character_)
  absolute <- abs(delta)
  if (absolute < 0.147) return("negligible")
  if (absolute < 0.330) return("small")
  if (absolute < 0.474) return("medium")
  "large"
}

comparison_subset <- function(data, comparison) {
  keep <- eval(comparison$subset_expression)
  subset_data <- data[keep, , drop = FALSE]
  group_values <- subset_data[[comparison$group_variable]]
  subset_data[group_values %in% c(comparison$group_1, comparison$group_2), ,
              drop = FALSE]
}

safe_wilcox <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  suppressWarnings(
    stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)$p.value
  )
}
