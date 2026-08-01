# Code provenance and scope

This repository is a cleaned reconstruction of the analysis workflow from retained project records, matrices, metadata, CIBERSORTx output, and final methodological decisions.

## Directly supported by retained files

- 64-sample final count matrix dimensions and sample accessions
- 48-sample GSE303084 metadata with Status, Sex, and Batch
- 48-sample TPM matrix
- 48-sample CIBERSORTx LM22 output
- CIBERSORTx quality metrics
- stacked-bar, violin/box, heatmap, Wilcoxon/FDR, and effect-size workflow decisions

## Reconstructed reusable components

The original individual featureCounts files and the exact original count-merging script were not included in the retained upload. Therefore:

- `01_merge_featurecounts_files.R` is a transparent reusable reconstruction for the known featureCounts file structure;
- `03_counts_to_tpm.R` reconstructs TPM normalization using the featureCounts `Length` field and Ensembl-to-symbol mapping;
- these scripts should be tested against the user's original local files before public release.

The downstream scripts are aligned with the retained metadata and CIBERSORTx output included in this repository.

## Excluded material

The repository intentionally excludes:

- AI conversation transcripts
- preliminary or duplicated code blocks
- screenshots embedded in Word documents
- unpublished narrative interpretations
- LM22 signature files
- full count and TPM matrices
- patient-identifiable clinical information
