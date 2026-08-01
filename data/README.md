# Data directories

- `public/`: small metadata and CIBERSORTx output retained for downstream reproducibility.
- `example/`: format-only subset of the TPM matrix; not suitable for biological inference or CIBERSORTx deconvolution.
- `templates/`: input templates for rebuilding the count and TPM matrices.
- `raw/`: local source files; ignored by Git.
- `processed/`: generated matrices and integrated tables; ignored by Git except `.gitkeep`.

The full count and TPM matrices are not included in the repository archive.
