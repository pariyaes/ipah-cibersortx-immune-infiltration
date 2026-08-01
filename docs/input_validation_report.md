# Input validation report

This report was generated from the retained project files before packaging the repository.

## Dimensions

- Combined count matrix: **63,187 genes × 64 samples**
- GSE126262 samples in count matrix: **16**
- GSE303084 samples in count matrix: **48**
- TPM matrix used for CIBERSORTx: **34,976 gene symbols × 48 samples**
- CIBERSORTx output: **48 samples × 22 LM22 immune-cell fractions**

## Identifier checks

- Duplicated count-matrix gene IDs: **0**
- Duplicated count-matrix sample IDs: **0**
- Duplicated TPM gene symbols: **0**
- Duplicated TPM sample IDs: **0**
- Metadata samples missing from TPM: **0**
- Metadata samples missing from CIBERSORTx output: **0**

## CIBERSORTx quality checks

- Samples with CIBERSORTx P-value < 0.05: **48/48**
- Minimum LM22 fraction sum: **1.000000000000**
- Maximum LM22 fraction sum: **1.000000000000**
- Median correlation: **0.5511**
- Median RMSE: **0.8340**

## Metadata balance

| Status | Sex | Batch | n |
|---|---|---:|---:|
| Control | Female | 1 | 10 |
| Control | Female | 2 | 3 |
| Control | Male | 1 | 5 |
| Control | Male | 2 | 6 |
| IPAH | Female | 1 | 10 |
| IPAH | Female | 2 | 3 |
| IPAH | Male | 1 | 5 |
| IPAH | Male | 2 | 6 |

The disease-status and sex distributions are matched across the two recorded batches in the retained 48-sample metadata.
