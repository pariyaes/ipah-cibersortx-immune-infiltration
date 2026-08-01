# Install the packages used by this repository.
# Run manually; project scripts never install packages automatically.

cran_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioconductor_packages <- c("AnnotationDbi", "org.Hs.eg.db")
missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
