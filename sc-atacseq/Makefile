CONDA_ENV ?= scvi
RBIO_ENV ?= rbio
SIGNAC_ENV ?= signac
TFBS_ENV ?= tfbs

RUN_ENV ?=
PYTHON_RUN = $(RUN_ENV) conda run -n $(CONDA_ENV) python
R_RUN = $(RUN_ENV) conda run -n $(RBIO_ENV) Rscript
SIGNAC_R_RUN = $(RUN_ENV) conda run -n $(SIGNAC_ENV) Rscript
TFBS_PYTHON_RUN = $(RUN_ENV) conda run -n $(TFBS_ENV) python

FILTERED_MULTIVI = results/intermediate/tala_integrated_multivi_morocco_expr_only.h5mu
FIGURE2D_GENE_LIST ?= data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv
FIGURE2D_REFERENCE_EXPRESSION ?= data/manuscript_inputs/figure2d_reference_expression.tsv.gz
FIGURE2D_ROBUST_DIR ?= results/tfbs_gene_overlap_integrated/robust_motif_validation
FIGURE2D_OUTPUT_DIR ?= results/tfbs_gene_overlap_integrated/robust_motif_validation_report/sex_stratified_tfbs_enrichment

.PHONY: help consensus_peaks atac peakvi mosaic multivi filter_reference_expression transfer install_archr_refs annotate_peaks peak_gene_links manuscript_atac manuscript_atac_from_cellranger figure2d_reference_expression figure2d_tfbs_regions figure2d_tfbs_scans figure2d_tfbs_fastas figure2d_tfbs_standardize figure2d_tfbs_enrichment figure2d_tfbs_sensitivity figure2d_tfbs_figure figure2d_tfbs_extended_figure figure2d_tfbs figure2d_tfbs_from_matrices manuscript_companion validate_manuscript_companion compile docs

help:
	@echo "ATAC processing and label transfer:"
	@echo "  make manuscript_atac - run from staged requantified ATAC matrices"
	@echo "  make manuscript_atac_from_cellranger - include consensus-peak requantification"
	@echo ""
	@echo "Figure 2D TFBS enrichment:"
	@echo "  make figure2d_tfbs - run from filtered peak-gene links through the primary and extended figures"
	@echo "  make figure2d_tfbs_from_matrices - rebuild the model and both figures from tracked matrices"
	@echo ""
	@echo "Validation:"
	@echo "  make compile"
	@echo "  make validate_manuscript_companion"

consensus_peaks:
	$(SIGNAC_R_RUN) code/scripts/01_build_consensus_peaks.R

atac:
	$(PYTHON_RUN) code/scripts/02_build_atac_anndata.py

peakvi:
	$(PYTHON_RUN) code/scripts/03_train_peakvi.py

mosaic:
	$(PYTHON_RUN) code/scripts/04_build_mosaic.py

multivi:
	$(PYTHON_RUN) code/scripts/05_train_multivi.py

filter_reference_expression:
	$(PYTHON_RUN) code/scripts/06_filter_integrated_expression_reference.py --output $(FILTERED_MULTIVI)

transfer:
	$(PYTHON_RUN) code/scripts/07_transfer_labels.py --input $(FILTERED_MULTIVI)

install_archr_refs:
	$(R_RUN) code/scripts/08_install_archr_hg38_refs.R

annotate_peaks:
	$(R_RUN) code/scripts/09_annotate_peaks_archr.R

peak_gene_links:
	$(PYTHON_RUN) code/scripts/10_peak_gene_links.py

manuscript_atac:
	$(MAKE) atac
	$(MAKE) peakvi
	$(MAKE) mosaic
	$(MAKE) multivi
	$(MAKE) filter_reference_expression
	$(MAKE) transfer
	$(MAKE) annotate_peaks
	$(MAKE) peak_gene_links

manuscript_atac_from_cellranger:
	$(MAKE) consensus_peaks
	$(MAKE) manuscript_atac

figure2d_reference_expression:
	$(TFBS_PYTHON_RUN) code/scripts/11_build_figure2d_reference_expression.py --output $(FIGURE2D_REFERENCE_EXPRESSION)

figure2d_tfbs_regions:
	$(TFBS_PYTHON_RUN) code/scripts/12_build_tfbs_gene_overlap_regions.py --gene-list $(FIGURE2D_GENE_LIST) --skip-differential-metadata
	$(TFBS_PYTHON_RUN) code/scripts/14_prepare_reference_promoter_background.py --target-genes $(FIGURE2D_GENE_LIST)
	$(TFBS_PYTHON_RUN) code/scripts/15_prepare_linked_peak_background.py --target-genes $(FIGURE2D_GENE_LIST)

figure2d_tfbs_scans:
	$(R_RUN) code/scripts/13_scan_tfbs_gene_overlap_motifs.R --regions results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz --output results/tfbs_gene_overlap/tfbs_motif_hits.tsv.gz
	$(R_RUN) code/scripts/13_scan_tfbs_gene_overlap_motifs.R --regions results/tfbs_gene_overlap_integrated/background_reference_promoters.tsv.gz --output results/tfbs_gene_overlap_integrated/background_reference_promoter_motif_hits.tsv.gz
	$(R_RUN) code/scripts/13_scan_tfbs_gene_overlap_motifs.R --regions results/tfbs_gene_overlap_integrated/background_linked_peaks.tsv.gz --output results/tfbs_gene_overlap_integrated/background_linked_peak_motif_hits.tsv.gz

figure2d_tfbs_fastas:
	$(R_RUN) code/scripts/16_prepare_standard_motif_enrichment_fastas.R

figure2d_tfbs_standardize:
	$(TFBS_PYTHON_RUN) code/scripts/17_run_robust_motif_enrichment.py --reference-expression $(FIGURE2D_REFERENCE_EXPRESSION) --output-dir $(FIGURE2D_ROBUST_DIR)

figure2d_tfbs_enrichment:
	$(TFBS_PYTHON_RUN) code/scripts/18_run_sex_stratified_tfbs_enrichment.py --robust-dir $(FIGURE2D_ROBUST_DIR) --output-dir $(FIGURE2D_OUTPUT_DIR)

figure2d_tfbs_sensitivity:
	$(TFBS_PYTHON_RUN) code/scripts/19_assess_tfbs_support_sensitivity.py --input-dir $(FIGURE2D_OUTPUT_DIR)

figure2d_tfbs_figure:
	$(TFBS_PYTHON_RUN) code/scripts/20_build_figure2d_tfbs_summary.py --input-dir $(FIGURE2D_OUTPUT_DIR)

figure2d_tfbs_extended_figure:
	$(TFBS_PYTHON_RUN) code/scripts/21_build_figure2d_tfbs_extended_families.py --input-dir $(FIGURE2D_OUTPUT_DIR)

figure2d_tfbs:
	$(MAKE) figure2d_tfbs_regions
	$(MAKE) figure2d_tfbs_scans
	$(MAKE) figure2d_tfbs_fastas
	$(MAKE) figure2d_tfbs_standardize
	$(MAKE) figure2d_tfbs_enrichment
	$(MAKE) figure2d_tfbs_sensitivity
	$(MAKE) figure2d_tfbs_figure
	$(MAKE) figure2d_tfbs_extended_figure

figure2d_tfbs_from_matrices:
	$(MAKE) figure2d_tfbs_enrichment
	$(MAKE) figure2d_tfbs_sensitivity
	$(MAKE) figure2d_tfbs_figure
	$(MAKE) figure2d_tfbs_extended_figure

manuscript_companion:
	$(MAKE) manuscript_atac
	$(MAKE) figure2d_tfbs

validate_manuscript_companion:
	$(TFBS_PYTHON_RUN) code/scripts/22_validate_manuscript_companion.py

compile:
	python -m py_compile \
		code/scripts/_common.py \
		code/scripts/02_build_atac_anndata.py \
		code/scripts/03_train_peakvi.py \
		code/scripts/04_build_mosaic.py \
		code/scripts/05_train_multivi.py \
		code/scripts/06_filter_integrated_expression_reference.py \
		code/scripts/07_transfer_labels.py \
		code/scripts/10_peak_gene_links.py \
		code/scripts/11_build_figure2d_reference_expression.py \
		code/scripts/12_build_tfbs_gene_overlap_regions.py \
		code/scripts/14_prepare_reference_promoter_background.py \
		code/scripts/15_prepare_linked_peak_background.py \
		code/scripts/17_run_robust_motif_enrichment.py \
		code/scripts/18_run_sex_stratified_tfbs_enrichment.py \
		code/scripts/19_assess_tfbs_support_sensitivity.py \
		code/scripts/20_build_figure2d_tfbs_summary.py \
		code/scripts/21_build_figure2d_tfbs_extended_families.py \
		code/scripts/22_validate_manuscript_companion.py
	$(R_RUN) -e "invisible(parse(file='code/scripts/01_build_consensus_peaks.R'))"
	$(R_RUN) -e "invisible(parse(file='code/scripts/08_install_archr_hg38_refs.R'))"
	$(R_RUN) -e "invisible(parse(file='code/scripts/09_annotate_peaks_archr.R'))"
	$(R_RUN) -e "invisible(parse(file='code/scripts/13_scan_tfbs_gene_overlap_motifs.R'))"
	$(R_RUN) -e "invisible(parse(file='code/scripts/16_prepare_standard_motif_enrichment_fastas.R'))"

docs:
	@echo "See docs/manuscript_companion_atac_tfbs_reproduction.md"
