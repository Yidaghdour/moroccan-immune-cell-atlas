#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact
import statsmodels.api as sm


CELL_ORDER = [
    "B cells",
    "Monocytes",
    "Memory T cells",
    "Naive CD4 T cells",
    "Naive CD8 T cells",
    "Cytotoxic T/NK cells",
]
REGION_ORDER = [
    "promoter",
    "tss_proximal_peak",
    "proximal_peak",
    "distal_peak",
]
SEX_ORDER = ["F", "M"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit all full-rank cell-type- and sex-stratified TFBS enrichment models "
            "over submitted target promoters and linked peaks."
        )
    )
    parser.add_argument(
        "--target-regions",
        default="results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz",
    )
    parser.add_argument(
        "--robust-dir",
        default="results/tfbs_gene_overlap_integrated/robust_motif_validation",
    )
    parser.add_argument(
        "--output-dir",
        default=(
            "results/tfbs_gene_overlap_integrated/robust_motif_validation_report/"
            "sex_stratified_tfbs_enrichment"
        ),
    )
    return parser.parse_args()


def bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    valid = numeric.notna()
    output = pd.Series(np.nan, index=values.index, dtype=float)
    if not valid.any():
        return output
    p_values = numeric.loc[valid].clip(0, 1).to_numpy()
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    output.loc[valid] = restored
    return output


def odds_ratio_haldane(a: int, b: int, c: int, d: int) -> float:
    return float(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))


def standardized_design(table: pd.DataFrame) -> np.ndarray:
    columns: list[np.ndarray] = [np.ones(table.shape[0], dtype=float)]
    region_dummies = pd.get_dummies(
        table["region_class"].astype(str),
        prefix="region",
        drop_first=True,
        dtype=float,
    )
    columns.extend(region_dummies[column].to_numpy(dtype=float) for column in region_dummies.columns)
    for column in ["gc_fraction", "cpg_density", "log_sequence_length"]:
        values = pd.to_numeric(table[column], errors="coerce")
        values = values.fillna(values.median()).fillna(0).to_numpy(dtype=float)
        scale = float(np.std(values, ddof=0))
        if scale < 1e-8:
            continue
        candidate = (values - float(np.mean(values))) / scale
        current = np.column_stack(columns)
        if np.linalg.matrix_rank(np.column_stack([current, candidate])) > np.linalg.matrix_rank(current):
            columns.append(candidate)
    return np.column_stack(columns)


def validate_matrix_alignment(table: pd.DataFrame, matrix: sparse.csr_matrix, label: str) -> None:
    expected = np.arange(table.shape[0], dtype=int)
    observed = table["analysis_index"].to_numpy(dtype=int)
    if not np.array_equal(expected, observed):
        raise ValueError(f"{label} analysis_index is not aligned to table row order")
    if matrix.shape[0] != table.shape[0]:
        raise ValueError(f"{label} motif matrix row count does not match its region table")


def load_inputs(args: argparse.Namespace) -> dict[str, object]:
    robust_dir = Path(args.robust_dir)
    targets = pd.read_csv(args.target_regions, sep="\t", compression="infer")
    peaks = pd.read_csv(robust_dir / "peak_analysis_regions.tsv.gz", sep="\t").sort_values(
        "analysis_index"
    ).reset_index(drop=True)
    promoters = pd.read_csv(
        robust_dir / "promoter_analysis_regions.tsv.gz", sep="\t"
    ).sort_values("analysis_index").reset_index(drop=True)
    peak_matrix = sparse.load_npz(robust_dir / "peak_motif_presence_consensus.npz").tocsr()
    promoter_matrix = sparse.load_npz(
        robust_dir / "promoter_motif_presence_consensus.npz"
    ).tocsr()
    tf_symbols = pd.read_csv(robust_dir / "tested_tf_symbols.tsv", sep="\t")[
        "tf_symbol"
    ].astype(str).tolist()

    validate_matrix_alignment(peaks, peak_matrix, "Peak")
    validate_matrix_alignment(promoters, promoter_matrix, "Promoter")
    if peak_matrix.shape[1] != len(tf_symbols) or promoter_matrix.shape[1] != len(tf_symbols):
        raise ValueError("Motif matrix columns do not match the tested TF-symbol table")

    required = {
        "major_label",
        "stimulation",
        "sex",
        "target_gene",
        "region_class",
        "chrom",
        "start",
        "end",
        "peak_id",
    }
    missing = required.difference(targets.columns)
    if missing:
        raise ValueError(f"Target-region table is missing columns: {sorted(missing)}")
    if set(targets["stimulation"].dropna().astype(str)) != {"Unstim"}:
        raise ValueError("This analysis expects the submitted Unstimulated target-region set")

    return {
        "targets": targets,
        "peaks": peaks,
        "promoters": promoters,
        "peak_matrix": peak_matrix,
        "promoter_matrix": promoter_matrix,
        "tf_symbols": tf_symbols,
    }


def build_stratum_data(
    target_group: pd.DataFrame,
    major_label: str,
    peaks: pd.DataFrame,
    promoters: pd.DataFrame,
    peak_matrix: sparse.csr_matrix,
    promoter_matrix: sparse.csr_matrix,
) -> tuple[pd.DataFrame, sparse.csr_matrix, dict[str, int]]:
    peak_lookup = (
        peaks.drop_duplicates(["major_label", "peak_id"])
        .set_index(["major_label", "peak_id"])["analysis_index"]
    )
    promoter_lookup = (
        promoters.drop_duplicates(["chrom", "start", "end"])
        .set_index(["chrom", "start", "end"])["analysis_index"]
    )
    parts: list[pd.DataFrame] = []
    matrices: list[sparse.csr_matrix] = []
    target_counts = {region_class: 0 for region_class in REGION_ORDER}

    for region_class in REGION_ORDER:
        region_targets = target_group.loc[target_group["region_class"].eq(region_class)].copy()
        if region_targets.empty:
            continue
        if region_class == "promoter":
            region_targets = region_targets.drop_duplicates(["chrom", "start", "end"])
            indices = [
                promoter_lookup.get((row.chrom, int(row.start), int(row.end)), np.nan)
                for row in region_targets.itertuples(index=False)
            ]
            indices = np.asarray(indices, dtype=float)
            if np.isnan(indices).any():
                raise ValueError(
                    f"Some {major_label} target promoters were absent from the robust region table"
                )
            target_indices = indices.astype(int)
            background = promoters.loc[promoters["target"].eq(0)].copy()
            background_indices = background["analysis_index"].to_numpy(dtype=int)
            target_table = promoters.iloc[target_indices].copy()
            matrix = promoter_matrix
        else:
            region_targets = region_targets.drop_duplicates("peak_id")
            indices = [
                peak_lookup.get((major_label, str(peak_id)), np.nan)
                for peak_id in region_targets["peak_id"]
            ]
            indices = np.asarray(indices, dtype=float)
            if np.isnan(indices).any():
                raise ValueError(
                    f"Some {major_label} {region_class} targets were absent from the robust region table"
                )
            target_indices = indices.astype(int)
            background = peaks.loc[
                peaks["target"].eq(0)
                & peaks["major_label"].eq(major_label)
                & peaks["region_class"].eq(region_class)
            ].copy()
            background_indices = background["analysis_index"].to_numpy(dtype=int)
            target_table = peaks.iloc[target_indices].copy()
            matrix = peak_matrix

        target_table["target"] = 1
        background["target"] = 0
        target_counts[region_class] = int(target_table.shape[0])
        parts.extend([target_table, background])
        matrices.extend([matrix[target_indices, :], matrix[background_indices, :]])

    if not parts:
        raise ValueError(f"No target regions were available for {major_label}")
    table = pd.concat(parts, ignore_index=True, sort=False)
    motif_matrix = sparse.vstack(matrices, format="csr")
    table["log_sequence_length"] = np.log1p(
        pd.to_numeric(table["sequence_length"], errors="coerce")
    )
    if table[["gc_fraction", "cpg_density", "log_sequence_length"]].isna().any().any():
        raise ValueError("Sequence covariates were missing after target/background assembly")
    return table, motif_matrix, target_counts


def fit_stratum(
    table: pd.DataFrame,
    motif_matrix: sparse.csr_matrix,
    tf_symbols: list[str],
    major_label: str,
    sex: str,
) -> pd.DataFrame:
    target = table["target"].to_numpy(dtype=int)
    target_rows = np.flatnonzero(target == 1)
    background_rows = np.flatnonzero(target == 0)
    target_count = np.asarray(motif_matrix[target_rows, :].sum(axis=0)).ravel().astype(int)
    background_count = np.asarray(motif_matrix[background_rows, :].sum(axis=0)).ravel().astype(int)
    n_target = target_rows.size
    n_background = background_rows.size
    base = standardized_design(table)
    rows: list[dict[str, object]] = []

    for motif_index, tf_symbol in enumerate(tf_symbols):
        a = int(target_count[motif_index])
        b = n_target - a
        c = int(background_count[motif_index])
        d = n_background - c
        raw_odds_ratio = odds_ratio_haldane(a, b, c, d)
        raw_p = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        record: dict[str, object] = {
            "major_label": major_label,
            "sex": sex,
            "tf_symbol": tf_symbol,
            "target_regions": n_target,
            "background_regions": n_background,
            "target_consensus_hits": a,
            "background_consensus_hits": c,
            "raw_consensus_odds_ratio": raw_odds_ratio,
            "raw_consensus_fisher_p": raw_p,
            "adjusted_log_odds": np.nan,
            "adjusted_log2_odds": np.nan,
            "adjusted_odds_ratio": np.nan,
            "adjusted_ci_low": np.nan,
            "adjusted_ci_high": np.nan,
            "adjusted_p": np.nan,
            "covariance_estimator": pd.NA,
            "model_status": "not fitted",
        }
        motif = motif_matrix[:, motif_index].toarray().ravel().astype(float)
        if np.unique(motif).size < 2:
            record["model_status"] = "constant motif indicator"
            rows.append(record)
            continue
        design = np.column_stack([base, motif])
        if np.linalg.matrix_rank(design) <= np.linalg.matrix_rank(base):
            record["model_status"] = "motif indicator collinear with covariates"
            rows.append(record)
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = sm.GLM(target.astype(float), design, family=sm.families.Binomial())
                fit = model.fit(maxiter=100, disp=0, cov_type="HC3")
            if hasattr(fit, "converged") and not bool(fit.converged):
                raise RuntimeError("GLM did not converge")
            coefficient = float(fit.params[-1])
            standard_error = float(fit.bse[-1])
            p_value = float(fit.pvalues[-1])
            if not all(np.isfinite([coefficient, standard_error, p_value])):
                raise FloatingPointError("non-finite model estimate")
            record.update(
                {
                    "adjusted_log_odds": coefficient,
                    "adjusted_log2_odds": coefficient / np.log(2),
                    "adjusted_odds_ratio": float(np.exp(np.clip(coefficient, -30, 30))),
                    "adjusted_ci_low": float(
                        np.exp(np.clip(coefficient - 1.96 * standard_error, -30, 30))
                    ),
                    "adjusted_ci_high": float(
                        np.exp(np.clip(coefficient + 1.96 * standard_error, -30, 30))
                    ),
                    "adjusted_p": p_value,
                    "covariance_estimator": "HC3",
                    "model_status": "ok",
                }
            )
        except Exception as exc:  # pragma: no cover - data-dependent numerical failure
            record["model_status"] = f"failed: {type(exc).__name__}"
        rows.append(record)
    return pd.DataFrame(rows)


def run_analysis(
    inputs: dict[str, object],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    targets = inputs["targets"]
    peaks = inputs["peaks"]
    promoters = inputs["promoters"]
    peak_matrix = inputs["peak_matrix"]
    promoter_matrix = inputs["promoter_matrix"]
    tf_symbols = inputs["tf_symbols"]
    results: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []

    for major_label in CELL_ORDER:
        for sex in SEX_ORDER:
            target_group = targets.loc[
                targets["major_label"].eq(major_label) & targets["sex"].eq(sex)
            ].copy()
            if target_group.empty:
                coverage_rows.append(
                    {
                        "major_label": major_label,
                        "sex": sex,
                        "submitted_genes": 0,
                        "target_regions": 0,
                        **{f"target_{region}": 0 for region in REGION_ORDER},
                        "analysis_status": "no submitted target set",
                    }
                )
                continue
            table, motif_matrix, target_counts = build_stratum_data(
                target_group,
                major_label,
                peaks,
                promoters,
                peak_matrix,
                promoter_matrix,
            )
            n_target = int(table["target"].sum())
            coverage_rows.append(
                {
                    "major_label": major_label,
                    "sex": sex,
                    "submitted_genes": int(target_group["target_gene"].nunique()),
                    "target_regions": n_target,
                    **{
                        f"target_{region}": int(target_counts[region])
                        for region in REGION_ORDER
                    },
                    "analysis_status": "modeled",
                }
            )
            results.append(
                fit_stratum(
                    table,
                    motif_matrix,
                    tf_symbols,
                    major_label,
                    sex,
                )
            )

    result = pd.concat(results, ignore_index=True)
    result["within_stratum_fdr"] = result.groupby(["major_label", "sex"])[
        "adjusted_p"
    ].transform(bh_adjust)
    result["global_fdr"] = bh_adjust(result["adjusted_p"].fillna(1.0))
    result["significant_enrichment"] = (
        result["model_status"].eq("ok")
        & result["global_fdr"].le(0.05)
        & result["adjusted_odds_ratio"].gt(1)
    )
    coverage = pd.DataFrame(coverage_rows)
    return result, coverage


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs(args)
    results, coverage = run_analysis(inputs, args)
    results = results.sort_values(
        ["significant_enrichment", "global_fdr", "major_label", "sex", "tf_symbol"],
        ascending=[False, True, True, True, True],
    )
    results.to_csv(
        output_dir / "sex_stratified_tfbs_enrichment_all_estimable.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    coverage.to_csv(
        output_dir / "sex_stratified_tfbs_target_coverage.tsv",
        sep="\t",
        index=False,
    )
    print(
        f"Wrote {int(results['model_status'].eq('ok').sum())} estimable TF-by-stratum models "
        f"to {output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
