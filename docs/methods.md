# Methods

## Count-matrix preparation

Individual featureCounts-style sample files were imported using the `Geneid` field as the row identifier and the sample-specific count column as the expression value. The featureCounts `Length` field was retained for TPM normalization. Sample files were merged into a genes-by-samples matrix, and duplicated gene identifiers were treated as an error. The completed project matrix contained 63,187 gene rows and 64 samples across GSE126262 and GSE303084.

The repository script validates file existence, sample uniqueness, gene-ID uniqueness, nonnegative counts, and consistency of gene lengths across files.

## Sample selection and metadata

The CIBERSORTx branch used 48 samples from GSE303084. Metadata variables were standardized to:

- `Sample`
- `Status`: Control or IPAH
- `Sex`: Female or Male
- `Batch`: internal sequencing batch
- `Group`: combined disease-status and sex label

The 48 samples are balanced by disease status within each sex and batch stratum.

## TPM normalization

For each Ensembl gene and sample:

1. featureCounts gene length was converted from base pairs to kilobases;
2. reads per kilobase (RPK) were calculated;
3. each sample's RPK values were divided by the sum of RPK values in millions;
4. Ensembl version suffixes were removed for annotation;
5. Ensembl gene IDs were mapped to HGNC symbols using `org.Hs.eg.db`;
6. unmapped symbols were excluded; and
7. TPM values for duplicated gene symbols were summed.

The resulting project TPM matrix contained 34,976 gene symbols and 48 samples.

## CIBERSORTx

The TPM matrix was formatted as a tab-delimited gene-symbol-by-sample mixture matrix. Immune-cell deconvolution was performed through the CIBERSORTx web platform using:

- LM22 signature matrix;
- 100 permutations;
- quantile normalization disabled for RNA-seq; and
- the project batch-correction setting enabled.

The output contained relative fractions for 22 LM22 immune-cell populations together with permutation P-value, correlation, and RMSE. Samples with CIBERSORTx P-value below 0.05 were retained. All 48 retained project samples passed this criterion.

## Statistical analysis

Two-sided Wilcoxon rank-sum tests were performed for:

1. IPAH versus Control;
2. Female versus Male within IPAH; and
3. Female IPAH versus Female Control as a secondary comparison.

For each comparison, Benjamini–Hochberg correction was applied across the 22 LM22 cell types. The output reports group sizes, group medians, difference in medians, raw P-value, and FDR.

Cliff's delta was calculated for each cell type and comparison. Positive values indicate higher estimated fractions in Group 1; negative values indicate lower estimated fractions in Group 1.

A batch-stratified sensitivity analysis repeated the primary comparisons separately within each batch. It reports within-batch direction and significance without pooling batch-specific P-values and without applying ComBat or post hoc regression to the CIBERSORTx fractions.

## Visualization

- A stacked bar plot displays the relative LM22 composition of each sample.
- Violin plots with embedded boxplots display distributions for IPAH versus Control and Female versus Male within IPAH.
- The heatmap contains immune-cell types as rows and samples as columns; each immune-cell row is standardized to a Z-score before Euclidean/complete-linkage clustering.
- A Cliff's delta plot ranks immune-cell populations by effect size.

The stacked composition plot and heatmap are descriptive visualizations and are not used as independent statistical evidence.
