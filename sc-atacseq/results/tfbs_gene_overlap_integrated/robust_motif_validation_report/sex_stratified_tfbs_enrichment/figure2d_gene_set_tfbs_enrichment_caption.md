# Figure 2D gene-set TFBS enrichment

## Suggested caption

TFBS enrichment across cell-type- and sex-stratified concordant RNA-ATAC gene sets. For each Unstimulated gene set shown in Figure 2D, reference promoters (TSS +/-2 kb) and linked ATAC peaks were compared with non-target background regions using JASPAR 2024 and HOCOMOCO v14 consensus TF-symbol annotations. Binomial models adjusted for regulatory-region class, GC fraction, CpG density, and sequence length. A 64-configuration support analysis varied target-set and motif-state count requirements while retaining a fixed multiple-testing universe. Dot color represents the adjusted log2 odds ratio, dot area represents the largest global BH FDR across the configurations, and black outlines identify combinations with FDR <=0.05 and odds ratio >1 in every configuration. Closely related motifs, particularly IRF motifs, cannot distinguish the exact TF family member or establish TF occupancy.

## Main descriptive result

Female lymphoid gene sets show a recurrent IRF-like TFBS signature. Female B-cell targets additionally support MAFB and STAT3 motifs, whereas the Male Memory T-cell set supports IRF7 together with ZBED2, HIC1, and HAND1 motifs. No motif passed global FDR in the Female Monocyte set. These are gene-set-specific enrichment patterns, not formal Female-versus-Male effects.

## Displayed motifs

IRF7, IRF2, IRF9, IRF5, IRF4, IRF8, IRF6, STAT3, MAFB, ZBED2, HIC1, HAND1.

## Significant combinations

- CD4 Naive | Female: IRF7; adjusted OR=4.43, worst-case global FDR=0.00145.
- CD4 Naive | Female: IRF2; adjusted OR=3.91, worst-case global FDR=0.00317.
- CD4 Naive | Female: IRF9; adjusted OR=4.05, worst-case global FDR=0.0033.
- CD4 Naive | Female: IRF5; adjusted OR=4.05, worst-case global FDR=0.0081.
- CD4 Naive | Female: IRF6; adjusted OR=4.61, worst-case global FDR=0.0304.
- CD4 Naive | Female: IRF4; adjusted OR=3.04, worst-case global FDR=0.0387.
- CD4 Naive | Female: IRF8; adjusted OR=2.94, worst-case global FDR=0.0456.
- Cytotoxic T/NK | Female: IRF9; adjusted OR=6.79, worst-case global FDR=0.00651.
- Cytotoxic T/NK | Female: IRF2; adjusted OR=5.98, worst-case global FDR=0.0112.
- Cytotoxic T/NK | Female: IRF8; adjusted OR=5.35, worst-case global FDR=0.0149.
- Cytotoxic T/NK | Female: IRF7; adjusted OR=5.67, worst-case global FDR=0.0221.
- B cells | Female: IRF5; adjusted OR=4.97, worst-case global FDR=0.00645.
- B cells | Female: MAFB; adjusted OR=3.73, worst-case global FDR=0.0134.
- B cells | Female: IRF4; adjusted OR=3.62, worst-case global FDR=0.0304.
- B cells | Female: STAT3; adjusted OR=3.35, worst-case global FDR=0.041.
- Memory T | Female: IRF7; adjusted OR=3.37, worst-case global FDR=0.0081.
- Memory T | Female: IRF2; adjusted OR=3.16, worst-case global FDR=0.0121.
- Memory T | Female: IRF5; adjusted OR=3.52, worst-case global FDR=0.0149.
- Memory T | Female: IRF9; adjusted OR=2.84, worst-case global FDR=0.0287.
- Memory T | Female: IRF4; adjusted OR=2.84, worst-case global FDR=0.0334.
- Memory T | Male: ZBED2; adjusted OR=9.77, worst-case global FDR=2.78e-05.
- Memory T | Male: IRF7; adjusted OR=3.88, worst-case global FDR=0.00645.
- Memory T | Male: HIC1; adjusted OR=3.84, worst-case global FDR=0.0149.
- Memory T | Male: HAND1; adjusted OR=3.20, worst-case global FDR=0.041.
