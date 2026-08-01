# Public-release checklist

Before publishing this repository:

- [x] Confirm supervisor/coauthor permission to release the CIBERSORTx-derived output.
- [x] Confirm the associated manuscript is publicly available.
- [x] Run `source("run_downstream.R")` in R and inspect all warnings.
- [ ] Run preprocessing against the original featureCounts files and compare matrix dimensions.
- [ ] Verify the gene-length source and TPM column sums.
- [ ] Confirm the CIBERSORTx settings recorded in `docs/methods.md`.
- [ ] Remove any local absolute paths, tokens, or credentials.
- [ ] Confirm that no private patient data, reports, or unpublished candidate tables are present.
- [ ] Add the final GitHub URL to `CITATION.cff`.
- [ ] Pin the repository on the GitHub profile only after the workflow has been run successfully.
