#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


GENE_LIST = Path("data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv")
EXPRESSION_LOOKUP = Path("data/manuscript_inputs/figure2d_reference_expression.tsv.gz")
JASPAR = Path("references/motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
HOCOMOCO = Path("references/motifs/HOCOMOCO_v14_H14CORE_jaspar_explicit.txt")
ROBUST_DIR = Path("results/tfbs_gene_overlap_integrated/robust_motif_validation")
FIGURE_DIR = Path(
    "results/tfbs_gene_overlap_integrated/robust_motif_validation_report/"
    "sex_stratified_tfbs_enrichment"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the tracked manuscript inputs and published Figure 2D result profile."
    )
    parser.add_argument(
        "--mode",
        choices=["inputs", "published"],
        default="published",
        help="Validate only immutable inputs or also validate the tracked publication outputs.",
    )
    return parser.parse_args()


class Validator:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.passes: list[str] = []

    def require_file(self, path: Path) -> bool:
        if not path.is_file() or path.stat().st_size == 0:
            self.failures.append(f"Missing or empty file: {path}")
            return False
        return True

    def equal(self, observed: object, expected: object, label: str) -> None:
        if observed != expected:
            self.failures.append(f"{label}: expected {expected!r}, observed {observed!r}")
        else:
            self.passes.append(f"{label}: {observed}")

    def truth(self, condition: bool, label: str) -> None:
        if condition:
            self.passes.append(label)
        else:
            self.failures.append(label)


def count_motifs(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def validate_inputs(validator: Validator) -> None:
    required = [
        GENE_LIST,
        EXPRESSION_LOOKUP,
        JASPAR,
        HOCOMOCO,
        Path("results/annotations/reference_genes.bed.gz"),
        Path("results/annotations/reference_tss.bed.gz"),
        Path("results/annotations/consensus_peak_annotations.tsv.gz"),
        Path("results/final/label_transfer_summary.tsv"),
    ]
    if not all(validator.require_file(path) for path in required):
        return

    genes = pd.read_csv(GENE_LIST, sep=";")
    genes = genes.dropna(subset=["Celltype", "Condition", "Sex", "Gene"])
    validator.equal(genes.shape[0], 110, "Submitted gene-stratum rows")
    validator.equal(genes["Gene"].astype(str).str.upper().nunique(), 83, "Distinct submitted genes")
    validator.equal(
        genes[["Celltype", "Sex"]].drop_duplicates().shape[0],
        8,
        "Submitted cell-type/sex strata",
    )
    validator.truth(
        genes["Condition"].astype(str).eq("Unstimulated").all(),
        "All submitted Figure 2D rows are Unstimulated",
    )

    expression = pd.read_csv(EXPRESSION_LOOKUP, sep="\t")
    validator.equal(expression.shape[0], 36_601, "Reference-expression genes")
    validator.truth(
        expression["gene_symbol"].astype(str).is_unique,
        "Reference-expression gene symbols are unique",
    )
    target_genes = set(genes["Gene"].astype(str).str.upper())
    expression_genes = set(expression["gene_symbol"].astype(str).str.upper())
    validator.truth(
        target_genes.issubset(expression_genes),
        "Every submitted gene has a reference-expression value",
    )
    validator.equal(count_motifs(JASPAR), 879, "JASPAR PFMs")
    validator.equal(count_motifs(HOCOMOCO), 1_595, "HOCOMOCO PFMs")


def validate_published_outputs(validator: Validator) -> None:
    files = {
        "peak_table": ROBUST_DIR / "peak_analysis_regions.tsv.gz",
        "promoter_table": ROBUST_DIR / "promoter_analysis_regions.tsv.gz",
        "tested_tfs": ROBUST_DIR / "tested_tf_symbols.tsv",
        "all_estimable": FIGURE_DIR / "sex_stratified_tfbs_enrichment_all_estimable.tsv.gz",
        "results": FIGURE_DIR / "sex_stratified_tfbs_enrichment_results.tsv.gz",
        "coverage": FIGURE_DIR / "sex_stratified_tfbs_target_coverage.tsv",
        "significant": FIGURE_DIR / "sex_stratified_tfbs_enrichment_significant.tsv",
        "configurations": FIGURE_DIR / "tfbs_support_sensitivity_configurations.tsv",
        "sensitive": FIGURE_DIR / "tfbs_support_sensitive_calls.tsv",
        "values": FIGURE_DIR / "figure2d_gene_set_tfbs_enrichment_values.tsv",
        "caption": FIGURE_DIR / "figure2d_gene_set_tfbs_enrichment_caption.md",
        "pdf": FIGURE_DIR / "figure2d_gene_set_tfbs_enrichment.pdf",
        "png": FIGURE_DIR / "figure2d_gene_set_tfbs_enrichment.png",
        "extended_values": FIGURE_DIR
        / "figure2d_gene_set_tfbs_enrichment_extended_families_values.tsv",
        "extended_notes": FIGURE_DIR
        / "figure2d_gene_set_tfbs_enrichment_extended_families_notes.md",
        "extended_pdf": FIGURE_DIR
        / "figure2d_gene_set_tfbs_enrichment_extended_families.pdf",
        "extended_png": FIGURE_DIR
        / "figure2d_gene_set_tfbs_enrichment_extended_families.png",
    }
    if not all(validator.require_file(path) for path in files.values()):
        return

    peaks = pd.read_csv(files["peak_table"], sep="\t")
    promoters = pd.read_csv(files["promoter_table"], sep="\t")
    tested_tfs = pd.read_csv(files["tested_tfs"], sep="\t")
    all_estimable = pd.read_csv(files["all_estimable"], sep="\t")
    results = pd.read_csv(files["results"], sep="\t")
    coverage = pd.read_csv(files["coverage"], sep="\t")
    significant = pd.read_csv(files["significant"], sep="\t")
    configurations = pd.read_csv(files["configurations"], sep="\t")
    sensitive = pd.read_csv(files["sensitive"], sep="\t")
    extended_values = pd.read_csv(files["extended_values"], sep="\t")

    validator.equal(int(peaks["target"].eq(1).sum()), 185, "Target linked peak-by-cell-type rows")
    validator.equal(int(peaks["target"].eq(0).sum()), 10_594, "Background linked peak-by-cell-type rows")
    validator.equal(int(promoters["target"].eq(1).sum()), 83, "Target promoters")
    validator.equal(int(promoters["target"].eq(0).sum()), 2_987, "Background promoters")
    validator.equal(tested_tfs["tf_symbol"].nunique(), 471, "Cross-database TF symbols")
    validator.equal(all_estimable.shape[0], 3_768, "Fixed TF-by-stratum hypothesis universe")
    validator.equal(
        int(all_estimable["model_status"].eq("ok").sum()),
        3_768,
        "Successfully fitted models before support filtering",
    )
    validator.equal(configurations.shape[0], 64, "Support-sensitivity configurations")
    validator.equal(significant.shape[0], 24, "Threshold-robust motif-stratum combinations")
    validator.equal(sensitive.shape[0], 3, "Threshold-sensitive motif-stratum combinations")
    validator.truth(
        significant["support_configurations_significant"].eq(64).all(),
        "Every reported enrichment is significant in all support configurations",
    )
    validator.truth(
        significant["maximum_support_grid_fdr"].le(0.05).all(),
        "Every reported enrichment has worst-case global FDR <= 0.05",
    )
    validator.equal(
        int(coverage["support_grid_status"].eq("eligible in all target-size configurations").sum()),
        6,
        "Figure 2D strata eligible in all target-size configurations",
    )
    validator.equal(extended_values.shape[0], 126, "Extended-panel motif-stratum estimates")
    validator.equal(
        extended_values["tf_symbol"].nunique(),
        21,
        "Extended-panel individual TF motifs",
    )
    validator.equal(
        int(extended_values["significant_enrichment"].sum()),
        24,
        "Threshold-robust combinations in the extended panel",
    )
    validator.equal(
        int(
            extended_values.loc[
                extended_values["display_group"].eq("Selected context"),
                "significant_enrichment",
            ].sum()
        ),
        0,
        "Threshold-robust combinations in the selected context block",
    )
    validator.truth(
        extended_values["display_fdr"].equals(
            extended_values["maximum_support_grid_fdr"]
        ),
        "Extended-panel dot sizes use worst-case support-grid FDR",
    )

    with files["pdf"].open("rb") as handle:
        validator.truth(handle.read(4) == b"%PDF", "Figure PDF has a valid header")
    with files["png"].open("rb") as handle:
        validator.truth(
            handle.read(8) == b"\x89PNG\r\n\x1a\n",
            "Figure PNG has a valid header",
        )
    with files["extended_pdf"].open("rb") as handle:
        validator.truth(
            handle.read(4) == b"%PDF",
            "Extended-family PDF has a valid header",
        )
    with files["extended_png"].open("rb") as handle:
        validator.truth(
            handle.read(8) == b"\x89PNG\r\n\x1a\n",
            "Extended-family PNG has a valid header",
        )


def main() -> None:
    args = parse_args()
    validator = Validator()
    validate_inputs(validator)
    if args.mode == "published":
        validate_published_outputs(validator)

    for message in validator.passes:
        print(f"PASS: {message}")
    for message in validator.failures:
        print(f"FAIL: {message}", file=sys.stderr)
    if validator.failures:
        raise SystemExit(1)
    print(f"Validation completed: {len(validator.passes)} checks passed")


if __name__ == "__main__":
    main()
