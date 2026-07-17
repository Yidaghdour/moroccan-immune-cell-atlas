# Manuscript Methods: scATAC-seq processing and RNA-to-ATAC label transfer

## Manuscript-ready Methods text

### Consensus peak set and scATAC-seq count matrix construction

Peak calls from the 17 scATAC-seq libraries (`1_001`-`1_017`) were concatenated, coordinate-sorted, and merged with a maximum inter-peak gap of 50 bp using `bedtools merge -d 50` [9] to generate a common peak universe. Fragment counts for each library were then requantified over this common set of intervals with the `FeatureMatrix` function in Signac v1.16.0 [1]. Cell barcodes identified as cells by Cell Ranger (`is__cell_barcode == 1`) were retained during requantification. The resulting sparse peak-by-cell matrices, barcode tables, and common peak coordinates were used as the inputs to the downstream Python workflow.

The 17 requantified matrices were read into AnnData and concatenated across libraries. Donor assignments were imported from the library-specific demultiplexing tables. Cells assigned as doublets or unassigned were removed from the study libraries; library `1_017`, which supplied the public paired multiome bridge, was assigned the identifier `Public_Donor_1`. Raw peak counts were retained in the `counts` layer. Cells with at least 1,000 detected peaks were retained, and peaks detected in at least 10 cells were kept for downstream analysis.

### PEAKVI embedding and clustering

The filtered scATAC-seq count matrix was modelled with PEAKVI [2] using scvi-tools v1.4.0.post1 [4], with `library_batch` supplied as the batch covariate. The model was trained for a maximum of 500 epochs with early stopping and random seed 0. The PEAKVI latent representation was stored as `X_peakvi`. A 30-nearest-neighbor graph was constructed from this latent representation with Scanpy [6], followed by UMAP [7] and Leiden clustering [8] at resolution 1.0.

### Mosaic RNA-ATAC integration with MultiVI

A mosaic integration object was assembled from three components: our PEAKVI-processed scATAC-seq data, the corresponding scRNA-seq reference, and a public paired 10x Genomics PBMC multiome dataset. For the paired bridge, gene-expression features were read from the public 10x feature-barcode matrix and accessibility profiles were taken from library `1_017`, which had been requantified over the same consensus peak set as our ATAC libraries. RNA and ATAC bridge cells were matched by cell barcode and combined into a paired multiome object. Library `1_017` was excluded from our unpaired ATAC component after construction of the bridge.

The paired bridge, our RNA-only cells, and our ATAC-only cells were organized with `scvi.data.organize_multiome_anndatas`. Raw RNA and ATAC count matrices were stored as separate modalities in a MuData object [10]. MultiVI [3] was configured with modality as the batch key and trained for a maximum of 500 epochs with early stopping and random seed 0. The shared latent representation was stored as `X_multivi`. A 30-nearest-neighbor graph was calculated in the MultiVI latent space using cosine distance, followed by UMAP [7] and Leiden clustering [8] at resolution 1.0.

### Reference reconciliation and label transfer

Before label transfer, integrated expression cells were reconciled against the final Morocco scRNA-seq reference (`Morocco_scRNA-seq.h5ad`). Modality suffixes were removed from integrated cell identifiers, and expression cells absent from the Morocco reference were excluded; accessibility and paired bridge cells were retained. Curated `Cluster_celltypes` and `Sub_Cluster_celltypes` labels were then copied from barcode-matched Morocco reference cells into the integrated RNA modality.

Major cell-type labels were transferred from seeded RNA cells to ATAC cells by two complementary procedures. First, a distance-weighted 50-nearest-neighbor classifier was fitted in `X_multivi` using the scikit-learn default Euclidean distance. Second, labels were propagated iteratively over the cosine-based MultiVI connectivity graph. Self-edges were removed before graph propagation, edge weights and the purity of previously assigned labels contributed to each vote, and only cells with purity at least 0.8 were permitted to vote in subsequent iterations. Propagation was run for at most 50 iterations.

Confidence thresholds were calibrated separately for each method using a 20% holdout of the seeded RNA cells. Twenty-six equally spaced candidate thresholds from 0.50 to 0.99 were evaluated, and the lowest threshold attaining at least 90% holdout precision was selected. The final major-cell-type column was populated from thresholded KNN assignments and then completed with thresholded graph assignments where the KNN result was missing. Agreement between the two thresholded methods was recorded separately.

Subcluster labels were transferred after major-cell-type assignment. B-cell subclusters were trained and transferred within the `B cells` lineage. T/NK subclusters were trained and transferred within `Naive CD4 T cells`, `Naive CD8 T cells`, `Memory T cells`, and `Cytotoxic T/NK cells`. The same KNN, graph-propagation, holdout-calibration, and KNN-first combination procedures were applied within each lineage.

### Consensus peak annotation

Consensus peaks were annotated against hg38 using ArchR v1.0.3 [5]. Genome and gene annotations were generated with `createGenomeAnnotation` and `createGeneAnnotation`, respectively. A one-base TSS was derived from each ArchR gene interval. Each peak was represented by its center coordinate and assigned to the closest annotated TSS with `GenomicRanges::distanceToNearest`, ignoring strand for the distance calculation. The exported annotation table contained peak coordinates, peak width, nearest gene identifier and symbol, gene strand, TSS coordinate, and absolute peak-center-to-TSS distance. Reference gene and TSS BED files derived from the same annotation were exported for downstream analyses.

### Donor-matched peak-to-gene links

Peak-to-gene links were calculated within each transferred major cell type using matched donor-level RNA and ATAC pseudobulks. Candidate pairs consisted of each peak and its nearest annotated gene when the peak center was within 250 kb of that gene's TSS and both features were present in the integrated object. For each donor and major cell type, raw RNA counts and raw ATAC peak counts were summed across cells and transformed as `log1p(count / library total x 10,000)`. A donor-cell-type profile was retained when it contained at least 20 RNA cells and 20 ATAC cells, and a major cell type required at least 10 eligible shared donors.

For each candidate pair, Pearson correlation was calculated between donor-level peak accessibility and expression of the nearest gene. Two-sided correlation P values were calculated from the Pearson t statistic with `n - 2` degrees of freedom and adjusted within each major cell type by the Benjamini-Hochberg method. Filtered links required an absolute Pearson correlation of at least 0.30 and FDR of at most 0.05. Both positive and negative correlations were retained. Peak-to-gene correlations were calculated at the major-cell-type level using all eligible donor profiles rather than in separate sex or stimulation strata.

## Reproducibility record

### Executed Make targets

```bash
make manuscript_atac
```

`make manuscript_atac` executes `make atac`, `make peakvi`, `make mosaic`, `make multivi`, `make filter_reference_expression`, `make transfer`, `make annotate_peaks`, and `make peak_gene_links` sequentially. `make manuscript_atac_from_cellranger` additionally constructs and requantifies the consensus peak set. The Python targets use `CONDA_ENV=scvi`, Signac requantification uses `SIGNAC_ENV=signac`, and the R annotation target uses `RBIO_ENV=rbio` by default.

### Target-to-file map

| Stage or Make target | Implementation | Principal input | Principal output |
|---|---|---|---|
| `make consensus_peaks` | `code/scripts/01_build_consensus_peaks.R` | Cell Ranger `peaks.bed`, `fragments.tsv.gz`, `fragments.tsv.gz.tbi`, and `singlecell.csv` files for libraries `1_001`-`1_017` under `data/raw/atac_cellranger/` | `data/processed/atac_requantified/consensus_peak_universe.bed` and library-specific requantified matrices |
| `make atac` | `code/scripts/02_build_atac_anndata.py` | `data/processed/atac_requantified/*_q_matrix.mtx`, `*_q_barcodes.tsv`, `*_q_peaks.bed`; `data/metadata/V/Pool*/donor_ids.tsv` | `results/intermediate/filtered_atac.h5ad` |
| `make peakvi` | `code/scripts/03_train_peakvi.py` | `results/intermediate/filtered_atac.h5ad` | `results/intermediate/peakvi_17.h5ad`; `results/models/all_batches_17/` |
| `make mosaic` | `code/scripts/04_build_mosaic.py` | `results/intermediate/peakvi_17.h5ad`; `data/raw/reference/merged_QC_split_SCT_merged_harmony-batch_Dims22_slim.h5ad`; `data/raw/reference/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5` | `results/intermediate/mosaic.h5ad` |
| `make multivi` | `code/scripts/05_train_multivi.py` | `results/intermediate/mosaic.h5ad` | `results/intermediate/tala_integrated_multivi.h5mu`; `results/models/muon_mvi_mosaic/` |
| `make filter_reference_expression` | `code/scripts/06_filter_integrated_expression_reference.py` | `results/intermediate/tala_integrated_multivi.h5mu`; `data/raw/reference/Morocco_scRNA-seq.h5ad` | `results/intermediate/tala_integrated_multivi_morocco_expr_only.h5mu`; `results/intermediate/tala_integrated_multivi_morocco_expr_only_filter_summary.tsv`; `results/intermediate/tala_integrated_multivi_non_morocco_expression_cells.tsv.gz` |
| `make transfer` | `code/scripts/07_transfer_labels.py` | `results/intermediate/tala_integrated_multivi_morocco_expr_only.h5mu`; `data/raw/reference/Morocco_scRNA-seq.h5ad` | `results/final/tala_integrated_multivi_labeled.h5mu`; `results/final/label_transfer_summary.tsv` |
| `make annotate_peaks` | `code/scripts/09_annotate_peaks_archr.R` | `data/processed/atac_requantified/consensus_peak_universe.bed` | `results/annotations/consensus_peak_annotations.tsv.gz`; `results/annotations/reference_genes.bed.gz`; `results/annotations/reference_tss.bed.gz` |
| `make peak_gene_links` | `code/scripts/10_peak_gene_links.py` | `results/final/tala_integrated_multivi_labeled.h5mu`; `results/annotations/consensus_peak_annotations.tsv.gz` | `results/annotations/peak_gene_links.tsv.gz`; `results/annotations/peak_gene_links.filtered.tsv.gz`; `results/annotations/peak_gene_link_summary.tsv` |

### Principal parameter values

| Stage | Parameter | Value |
|---|---|---:|
| Consensus peaks | Maximum merge gap | 50 bp |
| ATAC filtering | Minimum detected peaks per cell | 1,000 |
| ATAC filtering | Minimum cells per peak | 10 |
| PEAKVI | Batch key | `library_batch` |
| PEAKVI | Maximum epochs | 500, with early stopping |
| PEAKVI graph | Neighbors | 30 |
| PEAKVI clustering | Leiden resolution | 1.0 |
| MultiVI | Batch key | `modality` |
| MultiVI | Maximum epochs | 500, with early stopping |
| MultiVI graph | Neighbors and metric | 30; cosine |
| MultiVI clustering | Leiden resolution | 1.0 |
| Model training | Random seed | 0 |
| KNN transfer | Neighbors and weighting | 50; distance-weighted |
| KNN transfer | Distance | Euclidean, the scikit-learn default |
| Transfer calibration | Holdout and target precision | 20%; 0.90 |
| Graph transfer | Maximum iterations | 50 |
| Graph transfer | Minimum purity permitted to vote | 0.80 |
| Peak-gene candidates | Maximum TSS distance | 250 kb |
| Peak-gene pseudobulk | Minimum RNA and ATAC cells | 20 per donor-cell-type modality |
| Peak-gene correlation | Minimum shared donors | 10 |
| Filtered peak-gene link | Correlation and FDR | `abs(Pearson r) >= 0.30`; BH FDR `<= 0.05` |

### Data-dependent transfer thresholds

The stored transfer summary records 317,336 major-label RNA seeds and 167,801 ATAC query cells. Selected KNN and graph thresholds were 0.50 and 0.50 for major labels, 0.6372 and 0.6764 for B-cell subclusters, and 0.5784 and 0.7156 for T/NK subclusters, respectively. The corresponding values are retained in `results/final/label_transfer_summary.tsv`.

### Software version provenance

| Component | Version used or retained provenance |
|---|---|
| Signac | 1.16.0; package built under R 4.3.1 in the dn177 upstream peak-requantification library |
| bedtools | Used for `sort` and `merge -d 50`; the executable version was not retained in the repository |
| scvi-tools | 1.4.0.post1, recorded independently in both saved PEAKVI and MultiVI checkpoints |
| Retained dn177 Python stack | Python 3.13.3; AnnData 0.11.4; Scanpy 1.11.1; MuData 0.3.1; NumPy 2.2.5; pandas 2.2.3; SciPy 1.15.2; scikit-learn 1.5.2; PyTorch 2.5.1.post306; umap-learn 0.5.7; leidenalg 0.10.2; igraph 0.11.8 |
| MuData and muon | MuData 0.3.1 in the retained project stack; muon 0.1.7 in the local reproduction environment; the original muon version was not stored in the MultiVI checkpoint |
| ArchR annotation environment | R 4.4.3; ArchR 1.0.3; GenomicRanges 1.58.0; IRanges 2.40.1; S4Vectors 0.44.0; GenomicFeatures 1.58.0; AnnotationDbi 1.68.0 |
| hg38 annotation resources | BSgenome.Hsapiens.UCSC.hg38 1.4.5; TxDb.Hsapiens.UCSC.hg38.knownGene 3.20.0; org.Hs.eg.db 3.20.0 |

The scvi-tools version is taken from the saved model metadata. Supporting Python versions describe the retained project environment used for data handling and transfer; model checkpoints do not store every supporting package version.

## References

1. Stuart T, et al. Single-cell chromatin state analysis with Signac. *Nature Methods*. 2021;18:1333-1341. [doi:10.1038/s41592-021-01282-5](https://doi.org/10.1038/s41592-021-01282-5).
2. Ashuach T, Reidenbach DA, Gayoso A, Yosef N. PeakVI: A deep generative model for single-cell chromatin accessibility analysis. *Cell Reports Methods*. 2022;2:100182. [doi:10.1016/j.crmeth.2022.100182](https://doi.org/10.1016/j.crmeth.2022.100182).
3. Ashuach T, Gabitto MI, Koodli RV, Saldi GA, Jordan MI, Yosef N. MultiVI: deep generative model for the integration of multimodal data. *Nature Methods*. 2023;20:1222-1231. [doi:10.1038/s41592-023-01909-9](https://doi.org/10.1038/s41592-023-01909-9).
4. Gayoso A, et al. A Python library for probabilistic analysis of single-cell omics data. *Nature Biotechnology*. 2022;40:163-166. [doi:10.1038/s41587-021-01206-w](https://doi.org/10.1038/s41587-021-01206-w).
5. Granja JM, et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nature Genetics*. 2021;53:403-411. [doi:10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6).
6. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*. 2018;19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0).
7. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. 2018. [arXiv:1802.03426](https://doi.org/10.48550/arXiv.1802.03426).
8. Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*. 2019;9:5233. [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z).
9. Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26:841-842. [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033).
10. Bredikhin D, Kats I, Stegle O. MUON: multimodal omics analysis framework. *Genome Biology*. 2022;23:42. [doi:10.1186/s13059-021-02577-8](https://doi.org/10.1186/s13059-021-02577-8).
