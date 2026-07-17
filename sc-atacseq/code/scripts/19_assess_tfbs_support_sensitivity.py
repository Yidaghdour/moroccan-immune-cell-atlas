#!/usr/bin/env python3
from __future__ import annotations

import argparse
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_REGION_MINIMA = (5, 10, 15, 20)
TARGET_STATE_MINIMA = (1, 2, 3, 5)
BACKGROUND_STATE_MINIMA = (1, 5, 10, 25)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate TFBS enrichment across a grid of motif-count support rules and "
            "retain results that are significant under every configuration."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=(
            "results/tfbs_gene_overlap_integrated/robust_motif_validation_report/"
            "sex_stratified_tfbs_enrichment"
        ),
    )
    return parser.parse_args()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p_values = np.clip(np.asarray(values, dtype=float), 0.0, 1.0)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    return restored


def require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {sorted(missing)}")


def evaluate_support_grid(
    results: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = {
        "major_label",
        "sex",
        "tf_symbol",
        "target_regions",
        "background_regions",
        "target_consensus_hits",
        "background_consensus_hits",
        "adjusted_odds_ratio",
        "adjusted_p",
        "global_fdr",
        "significant_enrichment",
        "model_status",
    }
    require_columns(results, required, "All-estimable result table")
    if results.duplicated(["major_label", "sex", "tf_symbol"]).any():
        raise ValueError("TF-by-cell-type-by-sex hypotheses are not unique")

    table = results.copy()
    table["target_consensus_absent"] = (
        table["target_regions"] - table["target_consensus_hits"]
    )
    table["background_consensus_absent"] = (
        table["background_regions"] - table["background_consensus_hits"]
    )
    model_ok = table["model_status"].eq("ok") & table["adjusted_p"].notna()
    odds_positive = table["adjusted_odds_ratio"].gt(1).fillna(False)

    n_hypotheses = table.shape[0]
    eligible_count = np.zeros(n_hypotheses, dtype=np.int16)
    significant_count = np.zeros(n_hypotheses, dtype=np.int16)
    minimum_fdr = np.ones(n_hypotheses, dtype=float)
    maximum_fdr = np.zeros(n_hypotheses, dtype=float)
    configuration_rows: list[dict[str, object]] = []

    configurations = list(
        product(
            TARGET_REGION_MINIMA,
            TARGET_STATE_MINIMA,
            BACKGROUND_STATE_MINIMA,
        )
    )
    for index, (minimum_regions, minimum_target_state, minimum_background_state) in enumerate(
        configurations,
        start=1,
    ):
        eligible = (
            model_ok
            & table["target_regions"].ge(minimum_regions)
            & table["target_consensus_hits"].ge(minimum_target_state)
            & table["target_consensus_absent"].ge(minimum_target_state)
            & table["background_consensus_hits"].ge(minimum_background_state)
            & table["background_consensus_absent"].ge(minimum_background_state)
        )
        multiplicity_p = np.ones(n_hypotheses, dtype=float)
        multiplicity_p[eligible.to_numpy()] = table.loc[eligible, "adjusted_p"].to_numpy(
            dtype=float
        )
        adjusted = bh_adjust(multiplicity_p)
        significant = eligible.to_numpy() & (adjusted <= 0.05) & odds_positive.to_numpy()

        eligible_count += eligible.to_numpy(dtype=np.int16)
        significant_count += significant.astype(np.int16)
        minimum_fdr = np.minimum(minimum_fdr, adjusted)
        maximum_fdr = np.maximum(maximum_fdr, adjusted)
        configuration_rows.append(
            {
                "configuration": index,
                "minimum_target_regions": minimum_regions,
                "minimum_target_motif_present": minimum_target_state,
                "minimum_target_motif_absent": minimum_target_state,
                "minimum_background_motif_present": minimum_background_state,
                "minimum_background_motif_absent": minimum_background_state,
                "fixed_hypothesis_universe": n_hypotheses,
                "eligible_models": int(eligible.sum()),
                "significant_enrichments": int(significant.sum()),
            }
        )

    n_configurations = len(configurations)
    table = table.rename(
        columns={
            "global_fdr": "inclusive_global_fdr",
            "significant_enrichment": "inclusive_significant_enrichment",
        }
    )
    table["support_configurations_tested"] = n_configurations
    table["support_configurations_eligible"] = eligible_count
    table["support_configurations_significant"] = significant_count
    table["support_significance_fraction"] = significant_count / n_configurations
    table["minimum_support_grid_fdr"] = minimum_fdr
    table["maximum_support_grid_fdr"] = maximum_fdr
    table["global_fdr"] = maximum_fdr
    table["significant_enrichment"] = (
        model_ok
        & odds_positive
        & pd.Series(significant_count == n_configurations, index=table.index)
    )
    table["support_classification"] = np.select(
        [
            table["significant_enrichment"],
            table["support_configurations_significant"].gt(0),
        ],
        ["threshold-robust", "threshold-sensitive"],
        default="not significant in grid",
    )

    table = table.sort_values(
        [
            "significant_enrichment",
            "support_configurations_significant",
            "global_fdr",
            "major_label",
            "sex",
            "tf_symbol",
        ],
        ascending=[False, False, True, True, True, True],
    ).reset_index(drop=True)
    configuration_table = pd.DataFrame(configuration_rows)
    sensitive = table.loc[table["support_classification"].eq("threshold-sensitive")].copy()
    return table, configuration_table, sensitive


def write_notes(
    results: pd.DataFrame,
    configurations: pd.DataFrame,
    sensitive: pd.DataFrame,
    path: Path,
) -> None:
    robust = results.loc[results["significant_enrichment"]]
    lines = [
        "# Cell-type- and sex-stratified TFBS enrichment",
        "",
        "## Support-sensitivity design",
        "",
        "- Every full-rank TF-by-cell-type-by-sex model was fitted before support filtering.",
        f"- The multiplicity universe was fixed at {results.shape[0]} hypotheses in every analysis.",
        "- Sixty-four configurations crossed minimum target-region counts of 5, 10, 15, or 20; minimum motif-present and motif-absent target counts of 1, 2, 3, or 5; and corresponding background counts of 1, 5, 10, or 25.",
        "- Within each configuration, failed or ineligible hypotheses were assigned P=1 before global Benjamini-Hochberg correction.",
        "- Threshold-robust enrichment required adjusted odds ratio >1 and global FDR <=0.05 in all 64 configurations.",
        "",
        "## Result counts",
        "",
        f"- Successfully fitted models before support filtering: {int(results['model_status'].eq('ok').sum())}.",
        f"- Threshold-robust motif-stratum combinations: {robust.shape[0]}.",
        f"- Combinations significant in some but not all configurations: {sensitive.shape[0]}.",
        f"- Significant calls per configuration ranged from {int(configurations['significant_enrichments'].min())} to {int(configurations['significant_enrichments'].max())}.",
    ]
    if not sensitive.empty:
        lines.extend(["", "## Threshold-sensitive combinations", ""])
        for row in sensitive.itertuples(index=False):
            lines.append(
                f"- {row.major_label} | {row.sex} | {row.tf_symbol}: significant in "
                f"{row.support_configurations_significant}/{row.support_configurations_tested} configurations."
            )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.input_dir)
    all_models_path = output_dir / "sex_stratified_tfbs_enrichment_all_estimable.tsv.gz"
    coverage_path = output_dir / "sex_stratified_tfbs_target_coverage.tsv"
    results = pd.read_csv(all_models_path, sep="\t")
    coverage = pd.read_csv(coverage_path, sep="\t")

    results, configurations, sensitive = evaluate_support_grid(results)
    robust = results.loc[results["significant_enrichment"]].copy()
    maximum_target_regions = max(TARGET_REGION_MINIMA)
    coverage["support_grid_status"] = np.where(
        coverage["target_regions"].ge(maximum_target_regions),
        "eligible in all target-size configurations",
        "not eligible in all target-size configurations",
    )

    results.to_csv(
        output_dir / "sex_stratified_tfbs_enrichment_results.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    robust.to_csv(
        output_dir / "sex_stratified_tfbs_enrichment_significant.tsv",
        sep="\t",
        index=False,
    )
    coverage.to_csv(coverage_path, sep="\t", index=False)
    configurations.to_csv(
        output_dir / "tfbs_support_sensitivity_configurations.tsv",
        sep="\t",
        index=False,
    )
    sensitive.to_csv(
        output_dir / "tfbs_support_sensitive_calls.tsv",
        sep="\t",
        index=False,
    )
    write_notes(
        results,
        configurations,
        sensitive,
        output_dir / "sex_stratified_tfbs_enrichment_notes.md",
    )
    print(
        f"Retained {robust.shape[0]} threshold-robust enrichments; "
        f"{sensitive.shape[0]} additional calls were configuration-sensitive."
    )


if __name__ == "__main__":
    main()
