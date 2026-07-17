# ATAC and Figure 2D manuscript-companion reproduction

## Scope

This workflow reproduces the analysis described in:

- `docs/manuscript_methods_atac_processing_label_transfer.md`
- `docs/manuscript_methods_figure2d_tfbs.md`

It covers consensus-peak requantification, ATAC object construction, PEAKVI, mosaic construction, MultiVI, reference reconciliation, RNA-to-ATAC label transfer, ArchR peak annotation, donor-matched peak-to-gene links, Figure 2D TFBS region construction, motif scanning, enrichment modelling, and the final figure. It does not execute differential-accessibility testing. The Figure 2D gene list is treated as a pre-specified input and no RNA or ATAC differential-test threshold is added by this workflow.

Run all commands from the repository root.

## Software environments

Create the four workflow environments:

```bash
conda env create -f code/envs/scvi.yml
conda env create -f code/envs/signac.yml
conda env create -f code/envs/rbio.yml
conda env create -f code/envs/tfbs.yml
```

Install ArchR and its hg38 annotation dependencies once in `rbio`:

```bash
make install_archr_refs
```

The default Make variables are `CONDA_ENV=scvi`, `SIGNAC_ENV=signac`, `RBIO_ENV=rbio`, and `TFBS_ENV=tfbs`. They can be overridden on the command line.

## Tracked publication inputs

The repository includes:

- the submitted Figure 2D gene-set table under `data/manuscript_inputs/`
- the compact gene-level expression lookup used to define the eligible promoter background
- donor-assignment tables for libraries `1_001`-`1_016` under `data/metadata/V/`
- JASPAR 2024 and HOCOMOCO v14 motif bundles under `references/motifs/`
- the hg38 reference gene/TSS tables and consensus-peak annotations under `results/annotations/`
- compact motif-presence matrices and the Figure 2D reference outputs

Verify immutable inputs with:

```bash
shasum -a 256 -c docs/manuscript_input_checksums.sha256
make validate_manuscript_companion
```

## External data layout

The full ATAC reconstruction requires controlled or large inputs that are not stored in Git. Stage them at these paths:

```text
data/raw/atac_cellranger/1_001/outs/{peaks.bed,fragments.tsv.gz,fragments.tsv.gz.tbi,singlecell.csv}
...
data/raw/atac_cellranger/1_017/outs/{peaks.bed,fragments.tsv.gz,fragments.tsv.gz.tbi,singlecell.csv}
data/raw/reference/merged_QC_split_SCT_merged_harmony-batch_Dims22_slim.h5ad
data/raw/reference/Morocco_scRNA-seq.h5ad
data/raw/reference/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5
```

The exact required filenames are also checked by the individual CLI programs before computation begins. Large generated `h5ad`, `h5mu`, model, count-matrix, and complete peak-link files remain outside Git; see `docs/data_availability.md`.

## ATAC processing and label transfer

To start from the 17 Cell Ranger ATAC output directories:

```bash
make manuscript_atac_from_cellranger
```

This first runs `code/scripts/01_build_consensus_peaks.R`, then executes these stages sequentially:

```text
make atac
make peakvi
make mosaic
make multivi
make filter_reference_expression
make transfer
make annotate_peaks
make peak_gene_links
```

To start from already requantified matrices under `data/processed/atac_requantified/`, run:

```bash
make manuscript_atac
```

The principal terminal products are:

```text
results/final/tala_integrated_multivi_labeled.h5mu
results/final/label_transfer_summary.tsv
results/annotations/consensus_peak_annotations.tsv.gz
results/annotations/peak_gene_links.filtered.tsv.gz
```

The consensus-requantification CLI defaults reproduce the original R execution, including the original GRanges interpretation of BED starts. `--bed-starts-zero-based true` is available for a coordinate-corrected sensitivity reconstruction, but it does not recreate the count matrix used for the reported analysis.

## Figure 2D TFBS workflow

After `results/annotations/peak_gene_links.filtered.tsv.gz` is available, run the full sequence:

```bash
make figure2d_tfbs
```

The target performs these steps in order:

1. Build submitted-gene promoters and filtered linked-peak foregrounds, plus non-target promoter and linked-peak backgrounds.
2. Scan every foreground and background region with the same JASPAR/HOCOMOCO motif settings.
3. Extract hg38 sequences and calculate sequence-composition covariates.
4. Build the cross-database TF-symbol consensus motif-presence matrices.
5. Fit every full-rank cell-type- and sex-stratified TFBS enrichment model.
6. Evaluate all models across the 64 count-support configurations.
7. Render the primary Figure 2D-aligned panel and the extended family-context panel from threshold-robust results.

To reproduce only the final models and figure from the tracked standardized region tables and motif matrices:

```bash
make figure2d_tfbs_from_matrices
```

The primary and extended outputs are:

```text
results/tfbs_gene_overlap_integrated/robust_motif_validation_report/sex_stratified_tfbs_enrichment/figure2d_gene_set_tfbs_enrichment.pdf
results/tfbs_gene_overlap_integrated/robust_motif_validation_report/sex_stratified_tfbs_enrichment/figure2d_gene_set_tfbs_enrichment_extended_families.pdf
```

The extended panel contains the same threshold-robust motifs as the primary panel and adds manually selected representatives of AP-1/NF-kB, BACH, RUNX, CTCF, FOX, and HIF programs. These context motifs are displayed whether or not they pass the all-configuration enrichment rule and therefore must not be interpreted as positive evidence without a black ring.

The tracked reference result profile contains 110 submitted gene-stratum rows, 83 distinct genes, 471 cross-database TF symbols, 3,768 TF-by-stratum hypotheses, 3,768 successfully fitted models, and 24 motif-stratum combinations meeting the global FDR and odds-ratio rule in all 64 support configurations. `make validate_manuscript_companion` checks these values and the foreground/background dimensions.

Cache locations are not hardcoded by the scripts. On systems without writable default cache directories, supply an optional command prefix, for example `make RUN_ENV='XDG_CACHE_HOME=/tmp/tala-cache MPLCONFIGDIR=/tmp/tala-cache NUMBA_CACHE_DIR=/tmp/tala-cache' compile`.

## Syntax and command checks

```bash
make compile
make -n manuscript_atac
make -n figure2d_tfbs
```

`make -n` prints the exact command sequence without executing the computationally intensive model-training and motif-scanning stages.
