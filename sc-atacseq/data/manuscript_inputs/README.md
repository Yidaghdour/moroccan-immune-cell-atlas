# Manuscript analysis inputs

This directory contains the compact, non-donor-level inputs required to reproduce the Figure 2D TFBS analysis.

## Files

- `Gene_Overlaps_RNA-ATAC_TFAnalysis.csv` is the submitted Figure 2D gene-set table. It contains 110 cell-type/sex/gene rows and 83 distinct genes. The file is semicolon-delimited; its empty fifth column is retained from the submitted table.
- `figure2d_reference_expression.tsv.gz` contains one median `mean_log_cpm_overall` value per gene. It is used to retain reference promoters with an available expression estimate when the standardized promoter-background table is built. It contains no donor-level values or differential-test statistics.

The expression lookup can be regenerated from the full stratified RNA summary with:

```bash
make figure2d_reference_expression
```

The implementation is `code/scripts/11_build_figure2d_reference_expression.py`. The full RNA summary is not versioned because it is a large generated analysis table and is not an input to the Figure 2D enrichment model beyond this compact lookup.

## SHA-256 checksums

```text
ed2c0e7b74671221e478093a3b83f95e485133e627a113ceeb72d2091665d532  Gene_Overlaps_RNA-ATAC_TFAnalysis.csv
25aec9ab0cb2b9f4334977e082bf5ce835b99bba9446840787a7c85afec62bee  figure2d_reference_expression.tsv.gz
```
