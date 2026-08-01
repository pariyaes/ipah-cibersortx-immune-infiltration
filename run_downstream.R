# Run the downstream analysis from the retained CIBERSORTx output.
# Execute from the repository root:
# source("run_downstream.R")

source("R/02_prepare_metadata.R")
source("R/05_import_cibersortx_results.R")
source("R/06_primary_statistics.R")
source("R/06b_batch_stratified_sensitivity.R")
source("R/07_cliffs_delta.R")
source("R/08_generate_figures.R")
source("R/09_session_info.R")
source("R/10_reference_input_checks.R")

message("Downstream CIBERSORTx workflow completed.")
