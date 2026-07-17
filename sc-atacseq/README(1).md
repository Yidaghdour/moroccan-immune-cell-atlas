# ATAC and Figure 2D manuscript companion

This branch contains the minimal code, inputs, compact checkpoints, and instructions required to reproduce two manuscript Methods workflows:

1. scATAC-seq consensus-peak construction, PEAKVI embedding, RNA-ATAC mosaic integration with MultiVI, RNA-to-ATAC label transfer, peak annotation, and donor-matched peak-to-gene linking.
2. Figure 2D cell-type- and sex-stratified TFBS enrichment over submitted-gene promoters and linked ATAC peaks, including a fixed-universe 64-configuration support-sensitivity analysis and primary plus extended family-context figures.

Differential-accessibility testing, exploratory analyses, notebooks, historical reports, and unrelated manuscript figures are intentionally excluded.

## Start here

Read [the reproduction guide](docs/manuscript_companion_atac_tfbs_reproduction.md), then inspect the manuscript-ready Methods:

- [ATAC processing and label transfer](docs/manuscript_methods_atac_processing_label_transfer.md)
- [Figure 2D TFBS enrichment](docs/manuscript_methods_figure2d_tfbs.md)

Create the environments listed in the reproduction guide and run:

```bash
make help
make compile
make validate_manuscript_companion
```

## Workflow entry points

From staged requantified ATAC matrices:

```bash
make manuscript_atac
```

From Cell Ranger ATAC outputs:

```bash
make manuscript_atac_from_cellranger
```

From filtered peak-gene links through the complete Figure 2D TFBS analysis:

```bash
make figure2d_tfbs
```

From the tracked standardized region tables and consensus motif matrices:

```bash
make figure2d_tfbs_from_matrices
```

## Repository contents

- `code/`: only scripts invoked by the two workflows
- `code/envs/`: four workflow-specific conda environments
- `data/metadata/`: donor-assignment tables consumed during ATAC object construction
- `data/manuscript_inputs/`: submitted Figure 2D genes and compact promoter-background expression lookup
- `references/motifs/`: JASPAR 2024 and HOCOMOCO v14 motif bundles
- `results/annotations/`: compact hg38 annotation checkpoints
- `results/tfbs_gene_overlap*/`: target-region, motif-matrix, model-result, and final-figure checkpoints

Raw data, count matrices, trained models, full `h5ad`/`h5mu` objects, complete peak-link tables, and generated motif-scan intermediates are not stored in Git. Their required paths are documented in [data availability](docs/data_availability.md).
