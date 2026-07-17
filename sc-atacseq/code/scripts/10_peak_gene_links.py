#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as t_dist


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build donor-matched pseudobulk peak-to-gene links from the final labeled MuData object.",
    )
    parser.add_argument(
        "--input",
        default="results/final/tala_integrated_multivi_labeled.h5mu",
        help="Final labeled MuData object.",
    )
    parser.add_argument(
        "--peak-annotations",
        default="results/annotations/consensus_peak_annotations.tsv.gz",
        help="Peak annotation table from the ArchR annotation stage.",
    )
    parser.add_argument(
        "--output-tsv",
        default="results/annotations/peak_gene_links.tsv.gz",
        help="Full peak-to-gene correlation table.",
    )
    parser.add_argument(
        "--filtered-output-tsv",
        default="results/annotations/peak_gene_links.filtered.tsv.gz",
        help="Filtered high-confidence peak-to-gene links, including both positive and negative associations.",
    )
    parser.add_argument(
        "--summary-tsv",
        default="results/annotations/peak_gene_link_summary.tsv",
        help="Summary table per major label.",
    )
    parser.add_argument(
        "--rna-label-col",
        default="morocco_Cluster_celltypes",
        help="RNA label column used to define major-cell pseudobulks.",
    )
    parser.add_argument(
        "--atac-label-col",
        default="morocco_Cluster_celltypes_consensus",
        help="ATAC label column used to define major-cell pseudobulks.",
    )
    parser.add_argument(
        "--donor-col",
        default="donor_id",
        help="Shared donor identifier column.",
    )
    parser.add_argument(
        "--rna-modality",
        default="expression",
        help="Modality value used for RNA cells.",
    )
    parser.add_argument(
        "--atac-modality",
        default="accessibility",
        help="Modality value used for ATAC cells.",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=250_000,
        help="Maximum distance-to-TSS retained from the peak annotation table.",
    )
    parser.add_argument(
        "--min-rna-cells",
        type=int,
        default=20,
        help="Minimum RNA cells per donor and label.",
    )
    parser.add_argument(
        "--min-atac-cells",
        type=int,
        default=20,
        help="Minimum ATAC cells per donor and label.",
    )
    parser.add_argument(
        "--min-shared-donors",
        type=int,
        default=10,
        help="Minimum donor count required to compute links for a major label.",
    )
    parser.add_argument(
        "--min-corr",
        type=float,
        default=0.3,
        help="Minimum absolute Pearson correlation retained in the filtered output.",
    )
    parser.add_argument(
        "--max-fdr",
        type=float,
        default=0.05,
        help="Maximum Benjamini-Hochberg FDR retained in the filtered output.",
    )
    return parser.parse_args()


def ensure_parent_dir(path: str | Path) -> None:
    Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=np.float64)
    n = p_values.size
    if n == 0:
        return p_values
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return out


def log_cpm(matrix: np.ndarray, scale_to: float = 1e4) -> np.ndarray:
    row_sums = matrix.sum(axis=1, keepdims=True)
    safe = np.where(row_sums > 0, row_sums, 1.0)
    return np.log1p((matrix / safe) * scale_to)


def pseudobulk_sum(adata, cell_groups: list[np.ndarray]) -> np.ndarray:
    aggregates: list[np.ndarray] = []
    for cell_idx in cell_groups:
        subset = adata[cell_idx, :].X
        if sparse.issparse(subset):
            vector = np.asarray(subset.sum(axis=0)).ravel()
        else:
            vector = np.asarray(subset).sum(axis=0).ravel()
        aggregates.append(vector.astype(np.float32, copy=False))
    return np.vstack(aggregates) if aggregates else np.zeros((0, adata.n_vars), dtype=np.float32)


def compute_rowwise_pearson(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = x.shape[1]
    x_mean = x.mean(axis=1, keepdims=True)
    y_mean = y.mean(axis=1, keepdims=True)
    x_centered = x - x_mean
    y_centered = y - y_mean
    numerator = np.sum(x_centered * y_centered, axis=1)
    denominator = np.sqrt(np.sum(x_centered ** 2, axis=1) * np.sum(y_centered ** 2, axis=1))
    corr = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=np.float64),
        where=denominator > 0,
    )

    if n <= 2:
        return corr, np.full(corr.shape, np.nan, dtype=np.float64)

    clipped = np.clip(corr, -0.999999, 0.999999)
    t_stat = clipped * np.sqrt((n - 2) / (1 - clipped**2))
    p_values = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)
    p_values = np.where(np.isnan(corr), np.nan, p_values)
    return corr, p_values


def main() -> None:
    args = parse_args()

    mdata = md.read_h5mu(args.input, backed="r")
    obs = mdata.obs.copy()

    required_obs = [args.donor_col, args.rna_label_col, args.atac_label_col, "modality"]
    missing_obs = [column for column in required_obs if column not in obs.columns]
    if missing_obs:
        raise KeyError(f"Missing required columns in mdata.obs: {missing_obs}")

    peak_annotations = pd.read_csv(args.peak_annotations, sep="\t", compression="infer")
    peak_annotations = peak_annotations.dropna(subset=["closest_gene_symbol", "distance_to_tss"]).copy()
    peak_annotations["distance_to_tss"] = peak_annotations["distance_to_tss"].astype(int)
    peak_annotations = peak_annotations.loc[peak_annotations["distance_to_tss"] <= args.max_distance].copy()

    rna_var_names = pd.Index(mdata.mod["rna"].var_names.astype(str))
    atac_var_names = pd.Index(mdata.mod["atac"].var_names.astype(str))
    peak_annotations = peak_annotations.loc[
        peak_annotations["peak_id"].isin(atac_var_names) & peak_annotations["closest_gene_symbol"].isin(rna_var_names)
    ].copy()
    peak_annotations = peak_annotations.drop_duplicates(subset=["peak_id", "closest_gene_symbol"]).reset_index(drop=True)
    if peak_annotations.empty:
        raise ValueError("No candidate peak-gene pairs remained after filtering against the h5mu feature spaces.")

    rna_var_index = pd.Series(np.arange(len(rna_var_names)), index=rna_var_names)
    atac_var_index = pd.Series(np.arange(len(atac_var_names)), index=atac_var_names)
    peak_annotations["peak_idx"] = peak_annotations["peak_id"].map(atac_var_index).astype(int)
    peak_annotations["gene_idx"] = peak_annotations["closest_gene_symbol"].map(rna_var_index).astype(int)

    donor_series = obs[args.donor_col].astype(str)
    modality_series = obs["modality"].astype(str)
    rna_mask = modality_series.eq(args.rna_modality)
    atac_mask = modality_series.eq(args.atac_modality)

    rna_labels = obs[args.rna_label_col].astype("object")
    atac_labels = obs[args.atac_label_col].astype("object")
    major_labels = sorted(set(rna_labels.loc[rna_mask].dropna().astype(str)) & set(atac_labels.loc[atac_mask].dropna().astype(str)))

    all_link_tables: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []

    for label in major_labels:
        rna_group_mask = rna_mask & rna_labels.astype(str).eq(label)
        atac_group_mask = atac_mask & atac_labels.astype(str).eq(label)

        rna_counts = donor_series.loc[rna_group_mask].value_counts()
        atac_counts = donor_series.loc[atac_group_mask].value_counts()
        shared_donors = sorted(
            donor
            for donor in set(rna_counts.index) & set(atac_counts.index)
            if rna_counts[donor] >= args.min_rna_cells and atac_counts[donor] >= args.min_atac_cells
        )

        summary_row = {
            "major_label": label,
            "shared_donors": len(shared_donors),
            "rna_cells": int(rna_group_mask.sum()),
            "atac_cells": int(atac_group_mask.sum()),
            "candidate_pairs": 0,
            "tested_pairs": 0,
            "filtered_links": 0,
            "positive_filtered_links": 0,
            "negative_filtered_links": 0,
        }

        if len(shared_donors) < args.min_shared_donors:
            summary_rows.append(summary_row)
            continue

        shared_set = set(shared_donors)
        rna_obs_idx = np.flatnonzero(rna_group_mask.to_numpy() & donor_series.isin(shared_set).to_numpy())
        atac_obs_idx = np.flatnonzero(atac_group_mask.to_numpy() & donor_series.isin(shared_set).to_numpy())

        rna_group_indices: list[np.ndarray] = []
        atac_group_indices: list[np.ndarray] = []
        for donor in shared_donors:
            donor_rna_idx = np.flatnonzero(rna_group_mask.to_numpy() & donor_series.eq(donor).to_numpy())
            donor_atac_idx = np.flatnonzero(atac_group_mask.to_numpy() & donor_series.eq(donor).to_numpy())
            rna_group_indices.append(donor_rna_idx)
            atac_group_indices.append(donor_atac_idx)

        rna_bulk = pseudobulk_sum(mdata.mod["rna"], rna_group_indices)
        atac_bulk = pseudobulk_sum(mdata.mod["atac"], atac_group_indices)
        rna_bulk = log_cpm(rna_bulk)
        atac_bulk = log_cpm(atac_bulk)

        label_pairs = peak_annotations.copy()
        label_pairs["major_label"] = label
        label_pairs["shared_donors"] = len(shared_donors)
        summary_row["candidate_pairs"] = int(label_pairs.shape[0])

        peak_matrix = atac_bulk[:, label_pairs["peak_idx"].to_numpy()].T.astype(np.float64, copy=False)
        gene_matrix = rna_bulk[:, label_pairs["gene_idx"].to_numpy()].T.astype(np.float64, copy=False)

        corr, p_values = compute_rowwise_pearson(peak_matrix, gene_matrix)
        label_pairs["pearson_r"] = corr
        label_pairs["p_value"] = p_values
        valid_mask = np.isfinite(corr) & np.isfinite(p_values)
        label_pairs = label_pairs.loc[valid_mask].copy()
        if label_pairs.empty:
            summary_rows.append(summary_row)
            continue

        label_pairs["fdr"] = bh_fdr(label_pairs["p_value"].to_numpy())
        label_pairs["atac_median_logcpm"] = np.median(peak_matrix[valid_mask], axis=1)
        label_pairs["rna_median_logcpm"] = np.median(gene_matrix[valid_mask], axis=1)
        label_pairs["abs_pearson_r"] = label_pairs["pearson_r"].abs()
        label_pairs["link_direction"] = np.where(
            label_pairs["pearson_r"] >= 0,
            "positive",
            "negative",
        )
        label_pairs = label_pairs.sort_values(
            ["major_label", "fdr", "abs_pearson_r", "distance_to_tss"],
            ascending=[True, True, False, True],
        ).reset_index(drop=True)

        summary_row["tested_pairs"] = int(label_pairs.shape[0])
        filtered = label_pairs.loc[
            (label_pairs["abs_pearson_r"] >= args.min_corr) & (label_pairs["fdr"] <= args.max_fdr)
        ].copy()
        summary_row["filtered_links"] = int(filtered.shape[0])
        summary_row["positive_filtered_links"] = int(filtered["pearson_r"].ge(0).sum())
        summary_row["negative_filtered_links"] = int(filtered["pearson_r"].lt(0).sum())
        summary_rows.append(summary_row)
        all_link_tables.append(label_pairs)

    if not all_link_tables:
        raise ValueError("No peak-gene links were computed for any major label under the current donor/cell thresholds.")

    full_links = pd.concat(all_link_tables, ignore_index=True)
    filtered_links = full_links.loc[
        (full_links["abs_pearson_r"] >= args.min_corr) & (full_links["fdr"] <= args.max_fdr)
    ].copy()
    summary_df = pd.DataFrame(summary_rows).sort_values("major_label").reset_index(drop=True)

    output_columns = [
        "major_label",
        "peak_id",
        "chrom",
        "start",
        "end",
        "peak_center",
        "peak_width",
        "closest_gene_id",
        "closest_gene_symbol",
        "closest_gene_strand",
        "closest_tss",
        "distance_to_tss",
        "shared_donors",
        "pearson_r",
        "abs_pearson_r",
        "link_direction",
        "p_value",
        "fdr",
        "atac_median_logcpm",
        "rna_median_logcpm",
    ]

    ensure_parent_dir(args.output_tsv)
    full_links.loc[:, output_columns].to_csv(args.output_tsv, sep="\t", index=False, compression="infer")
    filtered_links.loc[:, output_columns].to_csv(args.filtered_output_tsv, sep="\t", index=False, compression="infer")
    summary_df.to_csv(args.summary_tsv, sep="\t", index=False)

    print(f"Wrote full peak-gene links to {Path(args.output_tsv).resolve()}")
    print(f"Wrote filtered peak-gene links to {Path(args.filtered_output_tsv).resolve()}")
    print(f"Wrote peak-gene summary to {Path(args.summary_tsv).resolve()}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
