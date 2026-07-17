# Manuscript Methods: Figure 2D TFBS enrichment summary

## Manuscript-ready Methods text

### Target gene sets and regulatory regions

Transcription factor binding-site enrichment was evaluated for the pre-specified unstimulated gene sets represented in Figure 2D. The input table (`Gene_Overlaps_RNA-ATAC_TFAnalysis.csv`) contained 110 cell-type/sex/gene rows representing 83 distinct genes. Gene symbols, major cell-type labels, stimulation status, and sex were standardized before interval construction. All submitted genes were retained during region construction, without an additional RNA or ATAC effect-size or FDR filter.

For each submitted gene, a reference promoter was defined as the hg38 TSS +/-2 kb using the TSS annotation exported from ArchR [4]. Filtered peak-to-gene links were joined to the submitted genes by major cell type and gene symbol. Linked ATAC peaks were divided into three mutually exclusive classes according to the absolute distance between the peak center and the linked-gene TSS: TSS-proximal (`<=2 kb`), proximal (`>2-50 kb`), and distal (`>50 kb`, with candidate links limited to 250 kb during peak-to-gene linking). Repeated promoter coordinates and repeated peak identifiers were deduplicated within each analysis stratum.

### Background regions

Promoter backgrounds were sampled from the ArchR-derived hg38 reference-gene set after excluding all submitted target genes. One promoter per reference gene was retained, defined as TSS +/-2 kb on chromosomes 1-22, X, or Y. Up to 5,000 candidate background promoters were sampled with random seed 1729, and promoters with available RNA abundance estimates were retained in the analysis table.

Linked-peak backgrounds were drawn from the filtered peak-to-gene-link universe. Background rows were required to occur on chromosomes 1-22, X, or Y, to be linked to a gene outside the submitted target list, and not to reuse a target peak. Background peaks were classified with the same `<=2 kb`, `>2-50 kb`, and `>50 kb` distance rules and deduplicated by major cell type, region class, and peak identifier. For each Figure 2D stratum, promoter targets were compared with the common non-target promoter background, whereas linked target peaks were compared with non-target linked peaks from the same major cell type and distance class.

### Motif scanning and cross-database consensus

Target and background sequences were extracted from `BSgenome.Hsapiens.UCSC.hg38`. Motifs were scanned with `motifmatchr::matchMotifs` v1.28.0 [3] at a per-PWM match threshold of `P < 5 x 10^-5`. The screen used 879 JASPAR 2024 CORE vertebrate non-redundant position-frequency matrices [1] and 1,595 HOCOMOCO v14 H14CORE matrices [2]. Motif names were normalized to TF symbols. The primary binary motif indicator for a region and TF required at least one JASPAR motif and at least one HOCOMOCO motif assigned to the same TF symbol in that region. The 471 normalized TF symbols occurring among target-region hits from both databases were carried forward for testing.

### Cell-type- and sex-stratified enrichment model

Each submitted major-cell-type/sex combination was analyzed as a separate stratum (for example, the female-monocyte gene set). For a given stratum, the foreground comprised the deduplicated regulatory intervals assigned to its submitted genes: fixed reference-promoter intervals spanning TSS +/-2 kb and filtered ATAC peaks linked to those genes. Each interval remained a separate model observation. Linked peaks were subdivided by their absolute distance from the linked-gene TSS into TSS-proximal (`<=2 kb`), proximal (`>2-50 kb`), and distal (`>50 kb`) classes. The model therefore contained four region classes: one fixed promoter class and three linked-ATAC-peak classes. Although the promoter and TSS-proximal-peak classes both cover TSS-adjacent sequence, they represent different intervals: the promoter is a predefined 4-kb reference window, whereas the TSS-proximal interval is an observed ATAC peak that passed the peak-to-gene-link filters.

For every TF symbol and eligible cell-type/sex stratum, all foreground intervals and their corresponding non-target background intervals were analyzed together in one binomial generalized linear model. Target-region status (foreground = 1; background = 0) was the outcome, and consensus TF motif presence was the predictor of interest. Region class was included as a categorical covariate to account for differences among promoter and linked-peak compartments; GC fraction, CpG density, and log-transformed sequence length were included as additional covariates. Continuous sequence covariates were standardized, and HC3 heteroscedasticity-consistent standard errors were used. Thus, each TF model produced one region-class-adjusted enrichment estimate for the cell-type/sex stratum, rather than a separate adjusted estimate for each region class.

To avoid selecting a single motif-count cutoff, every full-rank TF-by-cell-type-by-sex model was first fitted without a foreground-size or motif-count filter. Model-support sensitivity was then evaluated across 64 configurations formed by crossing minimum foreground sizes of 5, 10, 15, or 20 regions; minimum motif-present and motif-absent foreground counts of 1, 2, 3, or 5; and corresponding background counts of 1, 5, 10, or 25. The same fitted coefficient and nominal P value were used whenever a model met a configuration's support requirements. Within every configuration, failed or ineligible hypotheses were assigned `P=1`, and Benjamini-Hochberg correction [5] was applied over the same fixed universe of 3,768 TF-by-cell-type-by-sex hypotheses. This prevented the multiple-testing family from changing with the support rule.

A motif-stratum combination was retained for the manuscript summary only when its covariate-adjusted odds ratio was greater than 1 and its global FDR was `<=0.05` in all 64 configurations. The largest FDR across the configurations was retained as the reported FDR. Thus, support requirements were assessed as a multiverse sensitivity analysis rather than used to select one primary `10/2/5` cutoff. The odds-ratio requirement restricted reporting to motif overrepresentation and did not impose a minimum effect-size threshold.

### Figure construction

The Figure 2D summary displays the six cell-type/sex strata eligible under all foreground-size configurations: female Monocytes, female CD4 Naive T cells, female Cytotoxic T/NK cells, female B cells, female Memory T cells, and male Memory T cells. Motifs were included in the displayed matrix when they met the all-configuration enrichment criterion in at least one displayed stratum; estimates for those motifs were then shown across all six strata. Dot color represents the covariate-adjusted log2 odds ratio. Dot area increases with `-log10` of the largest global BH FDR across the 64 configurations, with FDR values below `10^-6` assigned the same maximum size. Black outlines identify combinations with global FDR `<=0.05` and a covariate-adjusted odds ratio greater than 1 in every configuration. Twenty-four motif-stratum combinations met this rule; three additional combinations were significant in only 8 of 64 configurations and were not outlined or used to select displayed motifs. The target gene and regulatory-region counts shown above each column were calculated from the deduplicated foreground used by the model.

An extended family-context panel retained the same threshold-robust motifs and added manually selected representatives of AP-1/NF-kB (`FOS`, `JUN`, `NFKB1`, and `REL`), BACH (`BACH2`), RUNX (`RUNX1`), CTCF (`CTCF`), FOX (`FOXA2`), and HIF (`HIF1A`) programs. These additional motifs were selected by TF family rather than statistical significance and were displayed across the same six strata regardless of significance to provide a common comparison set. The same color, size, and black-outline encodings were used; an additional motif without a black outline did not meet the all-configuration enrichment rule.

## Reproducibility record

### Script sequence

The complete sequence is executed with:

```bash
make figure2d_tfbs
```

| Stage | Script | Main output used by the Figure 2D workflow |
|---|---|---|
| Build target promoters and linked peaks | `code/scripts/12_build_tfbs_gene_overlap_regions.py` | `results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz` |
| Scan target sequences | `code/scripts/13_scan_tfbs_gene_overlap_motifs.R` | `results/tfbs_gene_overlap/tfbs_motif_hits.tsv.gz` |
| Build promoter background | `code/scripts/14_prepare_reference_promoter_background.py` | `results/tfbs_gene_overlap_integrated/background_reference_promoters.tsv.gz` |
| Build linked-peak background | `code/scripts/15_prepare_linked_peak_background.py` | `results/tfbs_gene_overlap_integrated/background_linked_peaks.tsv.gz` |
| Scan promoter-background sequences | `code/scripts/13_scan_tfbs_gene_overlap_motifs.R`, with the promoter-background region table supplied through `--regions` | `results/tfbs_gene_overlap_integrated/background_reference_promoter_motif_hits.tsv.gz` |
| Scan linked-peak-background sequences | `code/scripts/13_scan_tfbs_gene_overlap_motifs.R`, with the linked-peak-background region table supplied through `--regions` | `results/tfbs_gene_overlap_integrated/background_linked_peak_motif_hits.tsv.gz` |
| Extract foreground/background sequences | `code/scripts/16_prepare_standard_motif_enrichment_fastas.R` | FASTA files under `standard_motif_enrichment/fastas/` |
| Build standardized region tables and JASPAR/HOCOMOCO presence matrices | `code/scripts/17_run_robust_motif_enrichment.py` | `peak_analysis_regions.tsv.gz`, `promoter_analysis_regions.tsv.gz`, motif-presence matrices, and `tested_tf_symbols.tsv` |
| Fit all full-rank Figure 2D cell-type/sex models | `code/scripts/18_run_sex_stratified_tfbs_enrichment.py` | `sex_stratified_tfbs_enrichment_all_estimable.tsv.gz` and target coverage table |
| Evaluate the 64 support configurations | `code/scripts/19_assess_tfbs_support_sensitivity.py` | Final enrichment table, configuration table, and threshold-robust calls |
| Construct the manuscript summary | `code/scripts/20_build_figure2d_tfbs_summary.py` | `figure2d_gene_set_tfbs_enrichment.pdf`, PNG, values table, and caption |
| Construct the extended family-context panel | `code/scripts/21_build_figure2d_tfbs_extended_families.py` | `figure2d_gene_set_tfbs_enrichment_extended_families.pdf`, PNG, values table, and notes |

The same `motifmatchr` scan settings were applied to target promoters, target linked peaks, non-target promoters, and non-target linked peaks. The publication target invokes region construction with `--skip-differential-metadata`; differential-expression and differential-accessibility tables are therefore neither joined nor used by this workflow. Target inclusion and the Figure 2D enrichment model use the submitted gene set and the region/motif variables described above.

### Principal input files

| Role | File |
|---|---|
| Submitted Figure 2D gene sets | `data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv` |
| Gene-level reference-expression lookup used for promoter-background eligibility | `data/manuscript_inputs/figure2d_reference_expression.tsv.gz` |
| Reference TSS coordinates | `results/annotations/reference_tss.bed.gz` |
| Filtered peak-to-gene links | `results/annotations/peak_gene_links.filtered.tsv.gz` |
| JASPAR motif collection | `references/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt` |
| HOCOMOCO motif collection | `references/motifs/HOCOMOCO_v14_H14CORE_jaspar_explicit.txt` |
| hg38 genome sequence | `BSgenome.Hsapiens.UCSC.hg38` |

### Principal parameter values

| Stage | Parameter | Value |
|---|---|---|
| Target promoter | Window | TSS +/-2 kb |
| Linked-peak classes | TSS-proximal; proximal; distal | `<=2 kb`; `>2-50 kb`; `>50 kb and <=250 kb` |
| Input peak-to-gene links | Retention rule inherited from linking stage | `abs(Pearson r) >= 0.30`; within-cell-type BH FDR `<=0.05` |
| Promoter background | Maximum candidates and random seed | 5,000; 1729 |
| Linked-peak background | Sampling | All eligible non-target links; no downsampling |
| Motif scanning | Per-PWM match threshold | `P < 5 x 10^-5` |
| Motif indicator | Cross-database requirement | At least one JASPAR and one HOCOMOCO hit for the same TF symbol in the same region |
| Enrichment model | Model and covariance | Binomial GLM; HC3 robust covariance |
| Enrichment model | Covariates | Region class, standardized GC fraction, CpG density, and log sequence length |
| Support sensitivity | Minimum foreground regions | 5, 10, 15, or 20 |
| Support sensitivity | Minimum motif-present and motif-absent foreground counts | 1, 2, 3, or 5 for each state |
| Support sensitivity | Minimum motif-present and motif-absent background counts | 1, 5, 10, or 25 for each state |
| Multiple testing | Family | Fixed set of 3,768 hypotheses in each configuration; ineligible/failed hypotheses assigned `P=1` before BH correction |
| Display rule | Robust enrichment | Global FDR `<=0.05` and odds ratio `>1` in all 64 configurations |

### Analysis dimensions retained in the generated files

| Quantity | Value |
|---|---:|
| Submitted gene-stratum rows | 110 |
| Distinct submitted genes | 83 |
| Submitted cell-type/sex strata | 8 |
| Strata eligible under all foreground-size configurations | 6 |
| Unique target promoters in the standardized analysis table | 83 |
| Target linked peak-by-cell-type rows | 185 |
| Non-target reference promoters | 2,987 |
| Non-target linked peak-by-cell-type rows | 10,594 |
| TF symbols occurring among target-region hits from both databases | 471 |
| Fixed TF-by-cell-type-by-sex hypothesis universe | 3,768 |
| Successfully fitted TF-by-cell-type-by-sex models | 3,768 |
| Support configurations | 64 |
| Motif-stratum combinations significant in all configurations | 24 |
| Motif-stratum combinations significant in only a subset of configurations | 3 |

### Foreground represented in the final summary

| Figure 2D stratum | Target genes | Promoters | Linked peaks `<=2 kb` | Linked peaks `>2-50 kb` | Linked peaks `>50 kb` | Total regulatory regions |
|---|---:|---:|---:|---:|---:|---:|
| Monocytes, Female | 34 | 34 | 5 | 17 | 3 | 59 |
| CD4 Naive T cells, Female | 17 | 17 | 9 | 18 | 5 | 49 |
| Cytotoxic T/NK cells, Female | 8 | 8 | 3 | 13 | 0 | 24 |
| B cells, Female | 11 | 11 | 6 | 22 | 5 | 44 |
| Memory T cells, Female | 19 | 19 | 7 | 28 | 6 | 60 |
| Memory T cells, Male | 15 | 15 | 10 | 18 | 4 | 47 |

### Software and database versions

| Component | Version |
|---|---|
| R | 4.4.3 |
| ArchR | 1.0.3 (source of the reference TSS annotation) |
| motifmatchr | 1.28.0 |
| TFBSTools | 1.44.0 |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 |
| GenomicRanges | 1.58.0 |
| Biostrings | 2.74.1 |
| data.table | 1.18.2.1 |
| Python | 3.12.13 |
| NumPy | 1.26.4 |
| pandas | 2.3.3 |
| SciPy | 1.17.1 |
| statsmodels | 0.14.6 |
| Matplotlib | 3.10.8 |
| JASPAR | 2024 CORE vertebrate non-redundant; 879 PFMs in the scanned bundle |
| HOCOMOCO | v14 H14CORE; 1,595 PFMs in the scanned bundle |

## References

1. Rauluseviciute I, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. *Nucleic Acids Research*. 2024;52:D174-D182. [doi:10.1093/nar/gkad1059](https://doi.org/10.1093/nar/gkad1059).
2. Vorontsov IE, et al. HOCOMOCO in 2024: a rebuild of the curated collection of binding models for human and mouse transcription factors. *Nucleic Acids Research*. 2024;52:D154-D163. [doi:10.1093/nar/gkad1077](https://doi.org/10.1093/nar/gkad1077).
3. Schep AN. motifmatchr: fast motif matching in R. Bioconductor. [doi:10.18129/B9.bioc.motifmatchr](https://doi.org/10.18129/B9.bioc.motifmatchr).
4. Granja JM, et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nature Genetics*. 2021;53:403-411. [doi:10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6).
5. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*. 1995;57:289-300. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).
