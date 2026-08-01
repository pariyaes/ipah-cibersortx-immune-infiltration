# Data dictionary

## `data/public/sample_metadata_cibersortx.tsv`

| Column | Meaning |
|---|---|
| `Sample` | Public GEO sample accession |
| `Status` | `Control` or `IPAH` |
| `Sex` | `Female` or `Male` |
| `Batch` | Internal sequencing batch recorded for GSE303084 |
| `Dataset` | GEO series accession |
| `Group` | Combined `Status_Sex` label |

## `data/public/CIBERSORTx_Job8_Adjusted.tsv`

| Column group | Meaning |
|---|---|
| `Mixture` | Sample accession |
| 22 LM22 columns | Estimated relative immune-cell fractions |
| `P-value` | CIBERSORTx permutation P-value |
| `Correlation` | CIBERSORTx mixture-fit correlation |
| `RMSE` | CIBERSORTx root-mean-square error |

The 22 fraction columns sum to approximately 1 for each sample.

## Preprocessing inputs not distributed

### Sample manifest

Required columns:

- `Sample`
- `File`
- `Dataset`
- `Status`
- `Sex`
- `Batch`
- `Include_in_CIBERSORTx`

### featureCounts files

Expected annotation columns include:

- `Geneid`
- `Chr`
- `Start`
- `End`
- `Strand`
- `Length`

The final numeric count column should correspond to one sample.

## Generated files

- `final_counts_matrix.csv`: raw integer counts, genes by samples
- `gene_lengths_from_featurecounts.tsv`: gene ID and featureCounts length
- `TPM_for_CIBERSORTx.tsv`: gene-symbol TPM matrix
- `CIBERSORTx_mixture_LM22.tsv`: validated web-upload file
- `cibersortx_metadata_integrated_qc_pass.tsv`: wide immune fractions plus metadata
- `cibersortx_metadata_long_qc_pass.tsv`: long-format immune fractions plus metadata
