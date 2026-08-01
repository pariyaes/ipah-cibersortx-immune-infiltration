# IPAH RNA-seq Count-Matrix Preparation and CIBERSORTx Immune-Infiltration Analysis

This repository contains a cleaned, reproducible R workflow for:

1. assembling featureCounts outputs into a genes-by-samples count matrix;
2. validating sample metadata;
3. converting raw counts to TPM using featureCounts gene lengths;
4. preparing an RNA-seq mixture file for CIBERSORTx;
5. processing CIBERSORTx LM22 output;
6. comparing immune-cell fractions by disease status and sex;
7. quantifying Cliff's delta effect sizes; and
8. generating stacked composition plots, violin–box plots, a row-scaled heatmap, and an effect-size ranking plot.

The workflow was developed for idiopathic pulmonary arterial hypertension (IPAH) transcriptomic data and emphasizes sex-stratified immune phenotyping.

## Repository scope

The broader count matrix contained **63,187 genes across 64 samples** from **GSE126262** and **GSE303084**. The immune-deconvolution branch retained **48 GSE303084 samples** with complete disease-status, sex, and batch metadata. The retained TPM matrix contained **34,976 gene symbols across 48 samples**.

This repository includes the cleaned 48-sample metadata and the small CIBERSORTx LM22 output so that the downstream analysis can be reproduced. The full count and TPM matrices are intentionally excluded because they are large derived files and should be regenerated from public source data.

## Scientific workflow

```text
featureCounts files
        |
        v
combined raw-count matrix
        |
        v
featureCounts gene lengths + Ensembl-to-symbol mapping
        |
        v
TPM matrix for selected samples
        |
        v
CIBERSORTx web analysis
LM22 | 100 permutations | QN disabled for RNA-seq
project batch-correction setting enabled
        |
        v
22 immune-cell fractions + P-value + correlation + RMSE
        |
        v
metadata integration: Status + Sex + Batch
        |
        +--> two-sided Wilcoxon tests
        |    - IPAH vs Control
        |    - Female vs Male within IPAH
        |    - Female IPAH vs Female Control (secondary)
        |    - BH correction across 22 cell types
        |
        +--> Cliff's delta effect sizes
        |
        +--> batch-stratified direction sensitivity analysis
        |
        v
figures and tabular results
```

## Key methodological decisions

- **Primary deconvolution:** CIBERSORTx using the LM22 leukocyte signature.
- **RNA-seq settings:** 100 permutations and quantile normalization disabled.
- **Sample quality:** CIBERSORTx P-value `< 0.05`; the retained project output contains 48 passing samples.
- **Primary tests:** two-sided Wilcoxon rank-sum tests.
- **Multiple testing:** Benjamini–Hochberg correction separately within each biological comparison across the 22 LM22 cell types.
- **Effect size:** Cliff's delta with direction and conventional magnitude categories.
- **Batch handling:** batch-corrected CIBERSORTx output plus within-batch sensitivity analysis; the scripts do not apply ComBat or regress batch out of the estimated fractions.
- **Heatmap orientation:** immune-cell types are rows and samples are columns; values are scaled within each immune-cell row.

## Files

```text
R/
├── 00_config.R
├── 00_dependencies.R
├── 00_utils.R
├── 01_merge_featurecounts_files.R
├── 02_prepare_metadata.R
├── 03_counts_to_tpm.R
├── 04_prepare_cibersortx_input.R
├── 05_import_cibersortx_results.R
├── 06_primary_statistics.R
├── 06b_batch_stratified_sensitivity.R
├── 07_cliffs_delta.R
├── 08_generate_figures.R
├── 09_session_info.R
└── 10_reference_input_checks.R

data/
├── public/
│   ├── sample_metadata_cibersortx.tsv
│   ├── CIBERSORTx_Job8_Adjusted.tsv
│   └── data_dimensions_summary.tsv
├── example/
│   └── tpm_format_example.tsv
├── templates/
│   ├── sample_manifest_template.tsv
│   └── gene_annotation_template.tsv
├── raw/
└── processed/

docs/
├── methods.md
├── data_dictionary.md
├── code_provenance.md
├── input_validation_report.md
└── public_release_checklist.md
```

## Reproduce the downstream analysis

Install dependencies:

```r
source("R/00_dependencies.R")
```

Run the complete downstream workflow:

```r
source("run_downstream.R")
```

This creates:

```text
results/
├── metadata_group_counts.tsv
├── cibersortx_qc_summary.tsv
├── wilcoxon_lm22_results.tsv
├── wilcoxon_lm22_fdr_lt_0.05.tsv
├── batch_stratified_wilcoxon_results.tsv
├── batch_direction_consistency.tsv
├── cliffs_delta_lm22_results.tsv
├── integrated_statistics_and_effect_sizes.tsv
├── reference_input_checks.tsv
└── sessionInfo.txt

figures/
├── Figure1_Immune_StackedBar_CIBERSORTx.png
├── Figure2_IPAH_vs_Control_Violin_Box.png
├── Figure3_Female_vs_Male_IPAH_Violin_Box.png
├── Figure4_Immune_Infiltration_Heatmap_CIBERSORTx.png
└── Figure5_Cliffs_Delta_Ranking.png
```

## Reconstruct the count and TPM matrices

Copy the manifest template:

```bash
cp data/templates/sample_manifest_template.tsv data/raw/sample_manifest.tsv
```

Replace the example rows with real sample IDs and featureCounts file paths, then run:

```r
source("run_preprocessing.R")
```

The CIBERSORTx web step must be performed separately. Do not upload or redistribute the LM22 signature matrix in this repository.

## Data and privacy

The sample identifiers are public GEO accessions. No directly identifying clinical information is included. The full count and TPM matrices are not distributed here.

Before making this repository public, confirm that the corresponding manuscript has been published and that the project supervisor or collaborators permit release of derived output files.

## Associated publication

Pourdadashi A, **Eskandari P**, Ghayour Najafabadi Z, Tajari A, Heydarzadeh S, Panahi M.  
*Integrated transcriptomic and machine learning analysis identifies female-associated candidate biomarkers in idiopathic pulmonary arterial hypertension.*  
**Computational Biology and Chemistry.** 2026; Article 109276.  
DOI: `10.1016/j.compbiolchem.2026.109276`

## Author

**Pariya Eskandari, MSc**  
Molecular & Human Genetics Researcher  
Scientific portfolio: `https://pariyaes.github.io/`

## License

The source code is released under the MIT License. Public-dataset accessions and CIBERSORTx-derived outputs remain subject to their original data-source and platform terms.
