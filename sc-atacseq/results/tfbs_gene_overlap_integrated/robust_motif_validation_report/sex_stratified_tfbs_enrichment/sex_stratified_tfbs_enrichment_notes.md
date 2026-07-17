# Cell-type- and sex-stratified TFBS enrichment

## Support-sensitivity design

- Every full-rank TF-by-cell-type-by-sex model was fitted before support filtering.
- The multiplicity universe was fixed at 3768 hypotheses in every analysis.
- Sixty-four configurations crossed minimum target-region counts of 5, 10, 15, or 20; minimum motif-present and motif-absent target counts of 1, 2, 3, or 5; and corresponding background counts of 1, 5, 10, or 25.
- Within each configuration, failed or ineligible hypotheses were assigned P=1 before global Benjamini-Hochberg correction.
- Threshold-robust enrichment required adjusted odds ratio >1 and global FDR <=0.05 in all 64 configurations.

## Result counts

- Successfully fitted models before support filtering: 3768.
- Threshold-robust motif-stratum combinations: 24.
- Combinations significant in some but not all configurations: 3.
- Significant calls per configuration ranged from 24 to 27.

## Threshold-sensitive combinations

- B cells | M | TFDP1: significant in 8/64 configurations.
- Naive CD8 T cells | F | IRF6: significant in 8/64 configurations.
- Naive CD8 T cells | F | SMAD5: significant in 8/64 configurations.
