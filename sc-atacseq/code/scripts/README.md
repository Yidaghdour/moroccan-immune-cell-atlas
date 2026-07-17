# Workflow scripts

Only scripts used by the two manuscript workflows are retained. Numbering follows execution order without gaps; `_common.py` contains shared ATAC-loading helpers and is not an executable stage.

1. `01_build_consensus_peaks.R`: merge Cell Ranger peaks and requantify fragments with Signac.
2. `02_build_atac_anndata.py`: combine requantified libraries, attach donors, remove doublets/unassigned cells, and filter cells and peaks.
3. `03_train_peakvi.py`: train PEAKVI and calculate the latent graph, UMAP, and Leiden clusters.
4. `04_build_mosaic.py`: combine our RNA and ATAC data with the paired public multiome bridge.
5. `05_train_multivi.py`: train MultiVI and calculate the shared latent graph.
6. `06_filter_integrated_expression_reference.py`: reconcile integrated expression cells against the final Morocco reference.
7. `07_transfer_labels.py`: transfer major and lineage-restricted labels by KNN and graph propagation.
8. `08_install_archr_hg38_refs.R`: install the ArchR and hg38 annotation prerequisites.
9. `09_annotate_peaks_archr.R`: annotate consensus peaks with closest-gene and TSS information.
10. `10_peak_gene_links.py`: calculate donor-matched, major-cell-type peak-to-gene correlations.
11. `11_build_figure2d_reference_expression.py`: regenerate the compact promoter-background expression lookup.
12. `12_build_tfbs_gene_overlap_regions.py`: build target promoters and linked-peak intervals.
13. `13_scan_tfbs_gene_overlap_motifs.R`: scan target or background regions with JASPAR and HOCOMOCO motifs.
14. `14_prepare_reference_promoter_background.py`: build the non-target promoter background.
15. `15_prepare_linked_peak_background.py`: build cell-type- and distance-matched linked-peak backgrounds.
16. `16_prepare_standard_motif_enrichment_fastas.R`: extract hg38 sequences for sequence covariates.
17. `17_run_robust_motif_enrichment.py`: standardize regions and build cross-database motif-presence matrices.
18. `18_run_sex_stratified_tfbs_enrichment.py`: fit all full-rank cell-type/sex TFBS models.
19. `19_assess_tfbs_support_sensitivity.py`: evaluate the fixed hypothesis universe across 64 support configurations.
20. `20_build_figure2d_tfbs_summary.py`: build the threshold-robust manuscript figure and caption.
21. `21_build_figure2d_tfbs_extended_families.py`: add selected AP-1/NF-kB, BACH, RUNX, CTCF, FOX, and HIF context to the threshold-robust figure.
22. `22_validate_manuscript_companion.py`: validate inputs, checkpoints, model dimensions, sensitivity results, and figure files.
