# Run count-matrix and TPM preparation.
#
# Before running:
# 1. Copy data/templates/sample_manifest_template.tsv to
#    data/raw/sample_manifest.tsv and replace the example rows.
# 2. Place the featureCounts files at the paths listed in the manifest.
#
# The CIBERSORTx web step is not automated by this repository.

source("R/01_merge_featurecounts_files.R")
source("R/02_prepare_metadata.R")
source("R/03_counts_to_tpm.R")
source("R/04_prepare_cibersortx_input.R")

message("Preprocessing completed. Upload the mixture file to CIBERSORTx.")
