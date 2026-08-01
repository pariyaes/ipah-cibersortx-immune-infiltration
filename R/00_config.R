# Project-wide configuration
# Run scripts from the repository root.

PROJECT_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

PATHS <- list(
  raw = file.path(PROJECT_ROOT, "data", "raw"),
  public = file.path(PROJECT_ROOT, "data", "public"),
  processed = file.path(PROJECT_ROOT, "data", "processed"),
  results = file.path(PROJECT_ROOT, "results"),
  figures = file.path(PROJECT_ROOT, "figures")
)

invisible(lapply(PATHS[c("raw", "processed", "results", "figures")], dir.create,
                 recursive = TRUE, showWarnings = FALSE))

LM22_CELL_TYPES <- c(
  "B cells naive",
  "B cells memory",
  "Plasma cells",
  "T cells CD8",
  "T cells CD4 naive",
  "T cells CD4 memory resting",
  "T cells CD4 memory activated",
  "T cells follicular helper",
  "T cells regulatory (Tregs)",
  "T cells gamma delta",
  "NK cells resting",
  "NK cells activated",
  "Monocytes",
  "Macrophages M0",
  "Macrophages M1",
  "Macrophages M2",
  "Dendritic cells resting",
  "Dendritic cells activated",
  "Mast cells resting",
  "Mast cells activated",
  "Eosinophils",
  "Neutrophils"
)

CIBERSORTX_QC_COLUMNS <- c("P-value", "Correlation", "RMSE")
CIBERSORTX_PVALUE_CUTOFF <- 0.05

GROUP_LEVELS <- c(
  "Control_Female",
  "Control_Male",
  "IPAH_Female",
  "IPAH_Male"
)

PRIMARY_COMPARISONS <- list(
  list(
    name = "IPAH_vs_Control",
    subset_expression = quote(rep(TRUE, nrow(data))),
    group_variable = "Status",
    group_1 = "IPAH",
    group_2 = "Control",
    analysis_layer = "Primary"
  ),
  list(
    name = "Female_vs_Male_within_IPAH",
    subset_expression = quote(data$Status == "IPAH"),
    group_variable = "Sex",
    group_1 = "Female",
    group_2 = "Male",
    analysis_layer = "Primary"
  ),
  list(
    name = "Female_IPAH_vs_Female_Control",
    subset_expression = quote(data$Sex == "Female"),
    group_variable = "Status",
    group_1 = "IPAH",
    group_2 = "Control",
    analysis_layer = "Secondary"
  )
)
