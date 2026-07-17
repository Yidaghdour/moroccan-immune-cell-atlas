#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier

from _common import materialize_if_gz

DEFAULT_MISSING_LABELS = {"n/a", "Public_n/a", "NA", ""}
DEFAULT_MAJOR_LABEL_KEY = "Cluster_celltypes"
DEFAULT_SUBCLUSTER_LABEL_KEY = "Sub_Cluster_celltypes"
DEFAULT_B_MAJOR_LABELS = ["B cells"]
DEFAULT_T_MAJOR_LABELS = [
    "Naive CD4 T cells",
    "Naive CD8 T cells",
    "Memory T cells",
    "Cytotoxic T/NK cells",
]


def ensure_parent_dir(path: str | Path) -> None:
    Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def atomic_h5mu_write(mdata: md.MuData, output_path: str | Path) -> None:
    output_path = Path(output_path).expanduser().resolve()
    ensure_parent_dir(output_path)
    with tempfile.NamedTemporaryFile(
        prefix=f"{output_path.stem}_",
        suffix=output_path.suffix,
        dir=output_path.parent,
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)

    try:
        mdata.write_h5mu(temp_path)
        os.replace(temp_path, output_path)
    finally:
        temp_path.unlink(missing_ok=True)


def valid_label_mask(values: pd.Series, missing_labels: set[str] | None = None) -> pd.Series:
    missing_labels = set(missing_labels or DEFAULT_MISSING_LABELS)
    as_string = values.astype("string")
    return values.notna() & ~as_string.fillna("").isin(missing_labels)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Transfer Morocco-derived RNA labels to ATAC cells using two methods: KNN in X_multivi and graph propagation on connectivities.",
    )
    parser.add_argument("--input", default="results/intermediate/tala_integrated_multivi.h5mu", help="Input h5mu with X_multivi and connectivities.")
    parser.add_argument("--morocco-h5ad", default="data/raw/reference/Morocco_scRNA-seq.h5ad", help="Morocco RNA h5ad used to seed labels on expression cells.")
    parser.add_argument("--output", default="results/final/tala_integrated_multivi_labeled.h5mu", help="Output h5mu with seed and transfer columns.")
    parser.add_argument("--summary-tsv", default="results/final/label_transfer_summary.tsv", help="Summary TSV for method metrics and method agreement.")
    parser.add_argument("--latent-key", default="X_multivi", help="Embedding key in mdata.obsm for KNN transfer.")
    parser.add_argument("--graph-key", default="connectivities", help="Connectivity graph key in mdata.obsp for graph propagation.")
    parser.add_argument("--label-key", action="append", dest="label_keys", help="Morocco label column to import and transfer. Repeat to transfer multiple labels.")
    parser.add_argument("--major-label-key", default=DEFAULT_MAJOR_LABEL_KEY, help="Morocco major label column to transfer globally.")
    parser.add_argument("--subcluster-label-key", default=DEFAULT_SUBCLUSTER_LABEL_KEY, help="Morocco subcluster label column to transfer within major lineages.")
    parser.add_argument("--reference-modality", nargs="+", default=["expression"], help="Modalities used as seeded RNA reference cells.")
    parser.add_argument("--query-modality", nargs="+", default=["accessibility"], help="Modalities that should receive transferred labels.")
    parser.add_argument("--missing-label", action="append", default=None, help="String value treated as missing. Repeat to add more.")
    parser.add_argument("--n-neighbors", type=int, default=50, help="Neighbor count for KNN transfer.")
    parser.add_argument("--holdout-frac", type=float, default=0.2, help="Fraction of seeded RNA cells used for holdout evaluation.")
    parser.add_argument("--target-precision", type=float, default=0.9, help="Target holdout precision used to choose a confidence threshold.")
    parser.add_argument("--graph-max-iter", type=int, default=50, help="Maximum iterations for graph propagation.")
    parser.add_argument("--graph-min-purity-to-vote", type=float, default=0.8, help="Minimum purity required for a graph-labeled cell to vote in later iterations.")
    parser.add_argument("--b-major-label", action="append", dest="b_major_labels", help="Major labels considered part of the B-cell lineage. Repeat to add more.")
    parser.add_argument("--t-major-label", action="append", dest="t_major_labels", help="Major labels considered part of the T/NK lineage. Repeat to add more.")
    parser.add_argument("--seed-prefix", default="morocco", help="Prefix for seed columns imported from the Morocco object.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed.")
    return parser.parse_args()


def infer_label_keys(obs: pd.DataFrame) -> list[str]:
    preferred = ["Cluster_celltypes", "Sub_Cluster_celltypes"]
    return [column for column in preferred if column in obs.columns]


def base_label_name(seed_prefix: str, label_key: str) -> str:
    return f"{seed_prefix}_{label_key}"


def make_morocco_join_key(index: pd.Index) -> pd.Index:
    return (
        index.astype(str)
        .str.replace("_expression", "", regex=False)
        .str.replace("_accessibility", "", regex=False)
        .str.replace("_paired", "", regex=False)
    )


def attach_morocco_labels(
    obs: pd.DataFrame,
    morocco_obs: pd.DataFrame,
    label_keys: list[str],
    seed_prefix: str,
    reference_modalities: list[str],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    obs = obs.copy()
    join_key = make_morocco_join_key(obs.index)
    ref_mask = obs["modality"].astype(str).isin(reference_modalities)
    join_series = pd.Series(join_key, index=obs.index)

    summary_rows: list[dict[str, object]] = []
    matched_mask = ref_mask & join_series.isin(morocco_obs.index)

    for label_key in label_keys:
        seed_col = base_label_name(seed_prefix, label_key)
        obs[seed_col] = pd.Series(pd.NA, index=obs.index, dtype="object")
        if label_key in morocco_obs.columns:
            matched_index = join_series.loc[matched_mask]
            obs.loc[matched_mask, seed_col] = morocco_obs.loc[matched_index, label_key].to_numpy()

        summary_rows.append(
            {
                "label_key": label_key,
                "lineage": "all",
                "method": "seed_import",
                "reference_cells": int(ref_mask.sum()),
                "matched_seed_cells": int(obs[seed_col].notna().sum()),
                "query_cells": np.nan,
                "holdout_accuracy": np.nan,
                "selected_threshold": np.nan,
                "selected_precision": np.nan,
                "selected_coverage": np.nan,
                "transferred_cells_before_threshold": np.nan,
                "transferred_cells_after_threshold": np.nan,
                "pairwise_agreement": np.nan,
            }
        )

    return obs, summary_rows


def add_summary_row(
    summary_rows: list[dict[str, object]],
    *,
    label_key: str,
    method: str,
    reference_cells: int,
    matched_seed_cells: int,
    query_cells: int | float | None,
    lineage: str = "all",
    holdout_accuracy: float | None = np.nan,
    selected_threshold: float | None = np.nan,
    selected_precision: float | None = np.nan,
    selected_coverage: float | None = np.nan,
    transferred_cells_before_threshold: int | float | None = np.nan,
    transferred_cells_after_threshold: int | float | None = np.nan,
    pairwise_agreement: float | None = np.nan,
) -> None:
    summary_rows.append(
        {
            "label_key": label_key,
            "lineage": lineage,
            "method": method,
            "reference_cells": reference_cells,
            "matched_seed_cells": matched_seed_cells,
            "query_cells": query_cells,
            "holdout_accuracy": holdout_accuracy,
            "selected_threshold": selected_threshold,
            "selected_precision": selected_precision,
            "selected_coverage": selected_coverage,
            "transferred_cells_before_threshold": transferred_cells_before_threshold,
            "transferred_cells_after_threshold": transferred_cells_after_threshold,
            "pairwise_agreement": pairwise_agreement,
        }
    )


def combine_prefer_left(obs: pd.DataFrame, left_col: str, right_col: str, out_col: str) -> pd.DataFrame:
    obs[out_col] = obs[left_col].astype("object")
    missing = obs[out_col].isna()
    obs.loc[missing, out_col] = obs.loc[missing, right_col].astype("object")
    return obs


def combine_from_columns(obs: pd.DataFrame, source_cols: list[str], out_col: str) -> pd.DataFrame:
    obs[out_col] = pd.Series(pd.NA, index=obs.index, dtype="object")
    for column in source_cols:
        fill_mask = obs[out_col].isna()
        obs.loc[fill_mask, out_col] = obs.loc[fill_mask, column].astype("object")
    return obs


def add_method_comparison_from_cols(
    obs: pd.DataFrame,
    left_col: str,
    right_col: str,
    comparison_col: str,
    query_mask: pd.Series,
) -> tuple[pd.DataFrame, float]:
    left = obs.loc[query_mask, left_col].astype("object")
    right = obs.loc[query_mask, right_col].astype("object")

    comparison = pd.Series("both_missing", index=left.index, dtype="object")
    comparison[left.notna() & right.notna() & (left == right)] = "agree"
    comparison[left.notna() & right.notna() & (left != right)] = "disagree"
    comparison[left.notna() & right.isna()] = "knn_only"
    comparison[left.isna() & right.notna()] = "graph_only"
    obs.loc[query_mask, comparison_col] = comparison

    comparable = left.notna() & right.notna()
    pairwise_agreement = float((left[comparable] == right[comparable]).mean()) if comparable.any() else np.nan
    return obs, pairwise_agreement


def make_lineage_seed_column(seed_col: str, lineage_name: str) -> str:
    return f"{seed_col}_{lineage_name}_seed"


def build_lineage_seed_labels(
    obs: pd.DataFrame,
    *,
    seed_col: str,
    major_seed_col: str,
    lineage_major_labels: list[str],
    lineage_name: str,
) -> tuple[pd.DataFrame, str]:
    lineage_seed_col = make_lineage_seed_column(seed_col, lineage_name)
    obs[lineage_seed_col] = pd.Series(pd.NA, index=obs.index, dtype="object")
    lineage_seed_mask = obs[major_seed_col].astype(str).isin(lineage_major_labels) & obs[seed_col].notna()
    obs.loc[lineage_seed_mask, lineage_seed_col] = obs.loc[lineage_seed_mask, seed_col].astype("object")
    return obs, lineage_seed_col


def choose_knn_confidence_threshold(
    embedding: np.ndarray,
    labels: np.ndarray,
    n_neighbors: int,
    holdout_frac: float,
    target_precision: float,
    seed: int,
) -> tuple[float | None, dict[str, float]]:
    label_counts = pd.Series(labels).value_counts()
    if len(labels) < 10 or label_counts.shape[0] < 2:
        return None, {"holdout_accuracy": np.nan, "threshold": np.nan, "coverage": np.nan, "precision": np.nan}

    stratify = labels if (label_counts >= 2).all() else None
    train_x, test_x, train_y, test_y = train_test_split(
        embedding,
        labels,
        test_size=holdout_frac,
        random_state=seed,
        stratify=stratify,
    )
    fitted_neighbors = max(1, min(n_neighbors, len(train_y)))
    model = KNeighborsClassifier(n_neighbors=fitted_neighbors, weights="distance")
    model.fit(train_x, train_y)
    proba = model.predict_proba(test_x)
    pred = model.classes_[proba.argmax(axis=1)]
    conf = proba.max(axis=1)

    threshold_grid = np.linspace(0.5, 0.99, 26)
    best_threshold = None
    best_precision = np.nan
    best_coverage = np.nan

    for threshold in threshold_grid:
        keep_mask = conf >= threshold
        if not keep_mask.any():
            continue
        precision = float((pred[keep_mask] == test_y[keep_mask]).mean())
        coverage = float(keep_mask.mean())
        if precision >= target_precision:
            best_threshold = float(threshold)
            best_precision = precision
            best_coverage = coverage
            break

    metrics = {
        "holdout_accuracy": float((pred == test_y).mean()),
        "threshold": np.nan if best_threshold is None else best_threshold,
        "coverage": best_coverage,
        "precision": best_precision,
    }
    return best_threshold, metrics


def weighted_graph_propagation(
    conn,
    seed_labels: pd.Series,
    max_iter: int,
    min_purity_to_vote: float,
) -> tuple[pd.Categorical, np.ndarray, np.ndarray]:
    labels = seed_labels.astype("category")
    categories = labels.cat.categories
    codes = labels.cat.codes.to_numpy()
    n_cells = conn.shape[0]
    n_labels = len(categories)

    conn_no_self = conn.copy().tocsr()
    conn_no_self.setdiag(0)
    conn_no_self.eliminate_zeros()

    final_codes = codes.copy()
    purity = np.full(n_cells, np.nan, dtype=np.float32)
    iter_assigned = np.full(n_cells, -1, dtype=np.int32)

    labeled_mask = final_codes != -1
    purity[labeled_mask] = 1.0
    iter_assigned[labeled_mask] = 0

    for iteration in range(1, max_iter + 1):
        voters_mask = (final_codes != -1) & (purity >= min_purity_to_vote)
        voters_idx = np.where(voters_mask)[0]
        if voters_idx.size == 0:
            break

        indicator = coo_matrix(
            (
                purity[voters_mask].astype(np.float32),
                (voters_idx, final_codes[voters_mask]),
            ),
            shape=(n_cells, n_labels),
            dtype=np.float32,
        ).tocsr()

        label_scores = conn_no_self @ indicator
        indptr = label_scores.indptr
        indices = label_scores.indices
        data = label_scores.data

        current_unlabeled_idx = np.where(final_codes == -1)[0]
        new_labels = 0

        for cell_idx in current_unlabeled_idx:
            start, end = indptr[cell_idx], indptr[cell_idx + 1]
            if start == end:
                continue
            row_scores = data[start:end]
            total_score = row_scores.sum()
            if total_score <= 0:
                continue
            max_position = row_scores.argmax()
            majority_label = indices[start + max_position]
            majority_score = row_scores[max_position]

            final_codes[cell_idx] = majority_label
            purity[cell_idx] = majority_score / total_score
            iter_assigned[cell_idx] = iteration
            new_labels += 1

        if new_labels == 0:
            break

    final_labels = pd.Categorical.from_codes(final_codes, categories=categories)
    return final_labels, purity, iter_assigned


def choose_graph_threshold(
    conn,
    seed_labels: pd.Series,
    holdout_frac: float,
    target_precision: float,
    max_iter: int,
    min_purity_to_vote: float,
    seed: int,
) -> tuple[float | None, dict[str, float]]:
    valid_mask = seed_labels.notna()
    valid_index = np.where(valid_mask.to_numpy())[0]
    labels = seed_labels.loc[valid_mask].astype(str).to_numpy()
    label_counts = pd.Series(labels).value_counts()
    if valid_index.size < 10 or label_counts.shape[0] < 2:
        return None, {"holdout_accuracy": np.nan, "threshold": np.nan, "coverage": np.nan, "precision": np.nan}

    stratify = labels if (label_counts >= 2).all() else None
    train_idx, test_idx = train_test_split(
        valid_index,
        test_size=holdout_frac,
        random_state=seed,
        stratify=stratify,
    )

    hidden_seed_labels = seed_labels.copy()
    hidden_seed_labels.iloc[test_idx] = np.nan

    propagated, purity, _ = weighted_graph_propagation(
        conn=conn,
        seed_labels=hidden_seed_labels,
        max_iter=max_iter,
        min_purity_to_vote=min_purity_to_vote,
    )

    eval_mask = pd.Series(False, index=seed_labels.index)
    eval_mask.iloc[test_idx] = True
    predicted_mask = eval_mask.to_numpy() & ~pd.isna(propagated) & ~np.isnan(purity)
    if predicted_mask.sum() == 0:
        return None, {"holdout_accuracy": np.nan, "threshold": np.nan, "coverage": 0.0, "precision": np.nan}

    pred = pd.Series(propagated.astype(object), index=seed_labels.index).iloc[test_idx].to_numpy()
    true = seed_labels.iloc[test_idx].astype(str).to_numpy()
    conf = purity[test_idx]
    valid_pred_mask = ~pd.isna(pred) & ~np.isnan(conf)
    pred = pred[valid_pred_mask]
    true = true[valid_pred_mask]
    conf = conf[valid_pred_mask]

    threshold_grid = np.linspace(0.5, 0.99, 26)
    best_threshold = None
    best_precision = np.nan
    best_coverage = np.nan

    for threshold in threshold_grid:
        keep_mask = conf >= threshold
        if not keep_mask.any():
            continue
        precision = float((pred[keep_mask] == true[keep_mask]).mean())
        coverage = float(keep_mask.mean())
        if precision >= target_precision:
            best_threshold = float(threshold)
            best_precision = precision
            best_coverage = coverage
            break

    metrics = {
        "holdout_accuracy": float((pred == true).mean()) if len(pred) else np.nan,
        "threshold": np.nan if best_threshold is None else best_threshold,
        "coverage": best_coverage,
        "precision": best_precision,
    }
    return best_threshold, metrics


def write_knn_results(
    obs: pd.DataFrame,
    seed_col: str,
    query_mask: pd.Series,
    embedding: np.ndarray,
    n_neighbors: int,
    holdout_frac: float,
    target_precision: float,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    valid_mask = valid_label_mask(obs[seed_col], DEFAULT_MISSING_LABELS)
    train_x = embedding[valid_mask.to_numpy()]
    train_y = obs.loc[valid_mask, seed_col].astype(str).to_numpy()

    threshold, metrics = choose_knn_confidence_threshold(
        embedding=train_x,
        labels=train_y,
        n_neighbors=n_neighbors,
        holdout_frac=holdout_frac,
        target_precision=target_precision,
        seed=seed,
    )

    fitted_neighbors = max(1, min(n_neighbors, len(train_y)))
    model = KNeighborsClassifier(n_neighbors=fitted_neighbors, weights="distance")
    model.fit(train_x, train_y)

    raw_col = f"{seed_col}_knn_transfer"
    conf_col = f"{seed_col}_knn_confidence"
    thr_col = f"{seed_col}_knn_transfer_thresholded"
    source_col = f"{seed_col}_knn_source"

    obs[raw_col] = obs[seed_col].astype("object")
    obs[conf_col] = np.nan
    obs[thr_col] = obs[seed_col].astype("object")
    obs[source_col] = np.where(valid_mask, "seed", "missing")
    if query_mask.sum() == 0:
        metrics["transferred_after_threshold"] = 0
        metrics["transferred_before_threshold"] = 0
        return obs, metrics

    query_x = embedding[query_mask.to_numpy()]
    query_proba = model.predict_proba(query_x)
    query_pred = model.classes_[query_proba.argmax(axis=1)]
    query_conf = query_proba.max(axis=1)

    obs.loc[query_mask, raw_col] = query_pred
    obs.loc[query_mask, conf_col] = query_conf

    if threshold is None:
        keep_mask = np.ones_like(query_conf, dtype=bool)
    else:
        keep_mask = query_conf >= threshold
    obs.loc[query_mask, thr_col] = pd.Series(
        np.where(keep_mask, query_pred, np.nan),
        index=obs.index[query_mask],
        dtype="object",
    )
    obs.loc[query_mask, source_col] = np.where(keep_mask, "predicted", "predicted_low_confidence")
    metrics["transferred_after_threshold"] = int(keep_mask.sum())
    metrics["transferred_before_threshold"] = int(query_mask.sum())
    return obs, metrics


def write_graph_results(
    obs: pd.DataFrame,
    seed_col: str,
    query_mask: pd.Series,
    conn,
    max_iter: int,
    min_purity_to_vote: float,
    holdout_frac: float,
    target_precision: float,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    seed_labels = obs[seed_col].astype("object")
    threshold, metrics = choose_graph_threshold(
        conn=conn,
        seed_labels=seed_labels,
        holdout_frac=holdout_frac,
        target_precision=target_precision,
        max_iter=max_iter,
        min_purity_to_vote=min_purity_to_vote,
        seed=seed,
    )
    propagated, purity, iter_assigned = weighted_graph_propagation(
        conn=conn,
        seed_labels=seed_labels,
        max_iter=max_iter,
        min_purity_to_vote=min_purity_to_vote,
    )

    raw_col = f"{seed_col}_graph_transfer"
    conf_col = f"{seed_col}_graph_purity"
    thr_col = f"{seed_col}_graph_transfer_thresholded"
    iter_col = f"{seed_col}_graph_iter"
    source_col = f"{seed_col}_graph_source"

    propagated_series = pd.Series(propagated.astype(object), index=obs.index)
    purity_series = pd.Series(purity, index=obs.index)
    iter_series = pd.Series(iter_assigned, index=obs.index)

    obs[raw_col] = seed_labels
    obs[conf_col] = purity_series.astype("float32")
    obs[thr_col] = seed_labels
    obs[iter_col] = iter_series.astype("int32")
    obs[source_col] = np.where(seed_labels.notna(), "seed", "missing")
    if query_mask.sum() == 0:
        metrics["transferred_before_threshold"] = 0
        metrics["transferred_after_threshold"] = 0
        return obs, metrics

    obs.loc[query_mask, raw_col] = propagated_series.loc[query_mask]
    if threshold is None:
        keep_mask = pd.Series(True, index=obs.index[query_mask])
    else:
        keep_mask = purity_series.loc[query_mask] >= threshold
    obs.loc[query_mask, thr_col] = propagated_series.loc[query_mask].where(keep_mask, np.nan)
    obs.loc[query_mask, source_col] = np.where(keep_mask.to_numpy(), "predicted", "predicted_low_confidence")

    metrics["transferred_before_threshold"] = int(propagated_series.loc[query_mask].notna().sum())
    metrics["transferred_after_threshold"] = int(obs.loc[query_mask, thr_col].notna().sum())
    return obs, metrics


def add_method_comparison(obs: pd.DataFrame, seed_col: str, query_mask: pd.Series) -> tuple[pd.DataFrame, float]:
    knn_col = f"{seed_col}_knn_transfer_thresholded"
    graph_col = f"{seed_col}_graph_transfer_thresholded"
    comparison_col = f"{seed_col}_knn_vs_graph_thresholded"

    left = obs.loc[query_mask, knn_col].astype("object")
    right = obs.loc[query_mask, graph_col].astype("object")

    comparison = pd.Series("both_missing", index=left.index, dtype="object")
    comparison[left.notna() & right.notna() & (left == right)] = "agree"
    comparison[left.notna() & right.notna() & (left != right)] = "disagree"
    comparison[left.notna() & right.isna()] = "knn_only"
    comparison[left.isna() & right.notna()] = "graph_only"
    obs.loc[query_mask, comparison_col] = comparison

    comparable = left.notna() & right.notna()
    pairwise_agreement = float((left[comparable] == right[comparable]).mean()) if comparable.any() else np.nan
    return obs, pairwise_agreement


def main() -> None:
    args = parse_args()

    mdata = md.read_h5mu(args.input)
    with materialize_if_gz(args.morocco_h5ad) as morocco_path:
        morocco = ad.read_h5ad(morocco_path, backed="r")

    if args.latent_key not in mdata.obsm:
        raise KeyError(f"Embedding {args.latent_key!r} was not found in mdata.obsm.")
    if args.graph_key not in mdata.obsp:
        raise KeyError(f"Graph {args.graph_key!r} was not found in mdata.obsp.")
    if "modality" not in mdata.obs.columns:
        raise KeyError("mdata.obs['modality'] is required to define seeded RNA and ATAC query cells.")

    morocco_obs = morocco.obs.copy()
    requested_label_keys = args.label_keys or infer_label_keys(morocco_obs)
    import_label_keys = list(
        dict.fromkeys(
            requested_label_keys + [args.major_label_key, args.subcluster_label_key]
        )
    )
    if not import_label_keys:
        raise ValueError("No Morocco label keys were provided and no default Morocco label columns were found.")
    for required_key in [args.major_label_key, args.subcluster_label_key]:
        if required_key not in morocco_obs.columns:
            raise KeyError(f"Required Morocco label column {required_key!r} was not found.")

    missing_labels = set(DEFAULT_MISSING_LABELS)
    if args.missing_label:
        missing_labels.update(args.missing_label)

    b_major_labels = args.b_major_labels or DEFAULT_B_MAJOR_LABELS
    t_major_labels = args.t_major_labels or DEFAULT_T_MAJOR_LABELS

    obs = mdata.obs.copy()
    obs, summary_rows = attach_morocco_labels(
        obs=obs,
        morocco_obs=morocco_obs,
        label_keys=import_label_keys,
        seed_prefix=args.seed_prefix,
        reference_modalities=args.reference_modality,
    )

    embedding = np.asarray(mdata.obsm[args.latent_key])
    conn = mdata.obsp[args.graph_key].tocsr()
    query_mask = obs["modality"].astype(str).isin(args.query_modality)
    if query_mask.sum() == 0:
        raise ValueError("No ATAC query cells matched the requested query modalities.")

    major_seed_col = base_label_name(args.seed_prefix, args.major_label_key)
    subcluster_seed_col = base_label_name(args.seed_prefix, args.subcluster_label_key)
    for seed_col in [major_seed_col, subcluster_seed_col]:
        obs.loc[~valid_label_mask(obs[seed_col], missing_labels), seed_col] = np.nan

    major_seeded_cells = int(obs[major_seed_col].notna().sum())
    if major_seeded_cells == 0:
        raise ValueError(f"No Morocco seed labels were imported for {args.major_label_key}.")

    obs, knn_metrics = write_knn_results(
        obs=obs,
        seed_col=major_seed_col,
        query_mask=query_mask,
        embedding=embedding,
        n_neighbors=args.n_neighbors,
        holdout_frac=args.holdout_frac,
        target_precision=args.target_precision,
        seed=args.seed,
    )
    add_summary_row(
        summary_rows,
        label_key=args.major_label_key,
        lineage="all",
        method="knn",
        reference_cells=major_seeded_cells,
        matched_seed_cells=major_seeded_cells,
        query_cells=int(query_mask.sum()),
        holdout_accuracy=knn_metrics["holdout_accuracy"],
        selected_threshold=knn_metrics["threshold"],
        selected_precision=knn_metrics["precision"],
        selected_coverage=knn_metrics["coverage"],
        transferred_cells_before_threshold=knn_metrics["transferred_before_threshold"],
        transferred_cells_after_threshold=knn_metrics["transferred_after_threshold"],
    )

    obs, graph_metrics = write_graph_results(
        obs=obs,
        seed_col=major_seed_col,
        query_mask=query_mask,
        conn=conn,
        max_iter=args.graph_max_iter,
        min_purity_to_vote=args.graph_min_purity_to_vote,
        holdout_frac=args.holdout_frac,
        target_precision=args.target_precision,
        seed=args.seed,
    )
    add_summary_row(
        summary_rows,
        label_key=args.major_label_key,
        lineage="all",
        method="graph",
        reference_cells=major_seeded_cells,
        matched_seed_cells=major_seeded_cells,
        query_cells=int(query_mask.sum()),
        holdout_accuracy=graph_metrics["holdout_accuracy"],
        selected_threshold=graph_metrics["threshold"],
        selected_precision=graph_metrics["precision"],
        selected_coverage=graph_metrics["coverage"],
        transferred_cells_before_threshold=graph_metrics["transferred_before_threshold"],
        transferred_cells_after_threshold=graph_metrics["transferred_after_threshold"],
    )

    major_knn_col = f"{major_seed_col}_knn_transfer_thresholded"
    major_graph_col = f"{major_seed_col}_graph_transfer_thresholded"
    major_compare_col = f"{major_seed_col}_knn_vs_graph_thresholded"
    major_consensus_col = f"{major_seed_col}_consensus"
    obs, pairwise_agreement = add_method_comparison_from_cols(
        obs,
        left_col=major_knn_col,
        right_col=major_graph_col,
        comparison_col=major_compare_col,
        query_mask=query_mask,
    )
    add_summary_row(
        summary_rows,
        label_key=args.major_label_key,
        lineage="all",
        method="knn_vs_graph_thresholded",
        reference_cells=major_seeded_cells,
        matched_seed_cells=major_seeded_cells,
        query_cells=int(query_mask.sum()),
        pairwise_agreement=pairwise_agreement,
    )
    obs = combine_prefer_left(obs, major_knn_col, major_graph_col, major_consensus_col)

    subcluster_seeded_cells = int(obs[subcluster_seed_col].notna().sum())
    lineage_specs = [
        ("b_lineage", b_major_labels),
        ("t_lineage", t_major_labels),
    ]
    lineage_knn_cols: list[str] = []
    lineage_graph_cols: list[str] = []
    lineage_query_union = pd.Series(False, index=obs.index)

    for lineage_name, lineage_major_labels in lineage_specs:
        obs, lineage_seed_col = build_lineage_seed_labels(
            obs=obs,
            seed_col=subcluster_seed_col,
            major_seed_col=major_seed_col,
            lineage_major_labels=lineage_major_labels,
            lineage_name=lineage_name,
        )
        lineage_seeded_cells = int(obs[lineage_seed_col].notna().sum())
        lineage_query_mask = query_mask & obs[major_consensus_col].astype(str).isin(lineage_major_labels)
        lineage_query_union = lineage_query_union | lineage_query_mask

        if lineage_seeded_cells == 0:
            add_summary_row(
                summary_rows,
                label_key=args.subcluster_label_key,
                lineage=lineage_name,
                method="seed_import",
                reference_cells=subcluster_seeded_cells,
                matched_seed_cells=subcluster_seeded_cells,
                query_cells=int(lineage_query_mask.sum()),
            )
            continue

        obs, knn_metrics = write_knn_results(
            obs=obs,
            seed_col=lineage_seed_col,
            query_mask=lineage_query_mask,
            embedding=embedding,
            n_neighbors=args.n_neighbors,
            holdout_frac=args.holdout_frac,
            target_precision=args.target_precision,
            seed=args.seed,
        )
        add_summary_row(
            summary_rows,
            label_key=args.subcluster_label_key,
            lineage=lineage_name,
            method="knn",
            reference_cells=subcluster_seeded_cells,
            matched_seed_cells=lineage_seeded_cells,
            query_cells=int(lineage_query_mask.sum()),
            holdout_accuracy=knn_metrics["holdout_accuracy"],
            selected_threshold=knn_metrics["threshold"],
            selected_precision=knn_metrics["precision"],
            selected_coverage=knn_metrics["coverage"],
            transferred_cells_before_threshold=knn_metrics["transferred_before_threshold"],
            transferred_cells_after_threshold=knn_metrics["transferred_after_threshold"],
        )

        obs, graph_metrics = write_graph_results(
            obs=obs,
            seed_col=lineage_seed_col,
            query_mask=lineage_query_mask,
            conn=conn,
            max_iter=args.graph_max_iter,
            min_purity_to_vote=args.graph_min_purity_to_vote,
            holdout_frac=args.holdout_frac,
            target_precision=args.target_precision,
            seed=args.seed,
        )
        add_summary_row(
            summary_rows,
            label_key=args.subcluster_label_key,
            lineage=lineage_name,
            method="graph",
            reference_cells=subcluster_seeded_cells,
            matched_seed_cells=lineage_seeded_cells,
            query_cells=int(lineage_query_mask.sum()),
            holdout_accuracy=graph_metrics["holdout_accuracy"],
            selected_threshold=graph_metrics["threshold"],
            selected_precision=graph_metrics["precision"],
            selected_coverage=graph_metrics["coverage"],
            transferred_cells_before_threshold=graph_metrics["transferred_before_threshold"],
            transferred_cells_after_threshold=graph_metrics["transferred_after_threshold"],
        )

        lineage_knn_col = f"{lineage_seed_col}_knn_transfer_thresholded"
        lineage_graph_col = f"{lineage_seed_col}_graph_transfer_thresholded"
        lineage_compare_col = f"{lineage_seed_col}_knn_vs_graph_thresholded"
        obs, pairwise_agreement = add_method_comparison_from_cols(
            obs,
            left_col=lineage_knn_col,
            right_col=lineage_graph_col,
            comparison_col=lineage_compare_col,
            query_mask=lineage_query_mask,
        )
        add_summary_row(
            summary_rows,
            label_key=args.subcluster_label_key,
            lineage=lineage_name,
            method="knn_vs_graph_thresholded",
            reference_cells=subcluster_seeded_cells,
            matched_seed_cells=lineage_seeded_cells,
            query_cells=int(lineage_query_mask.sum()),
            pairwise_agreement=pairwise_agreement,
        )

        lineage_knn_cols.append(lineage_knn_col)
        lineage_graph_cols.append(lineage_graph_col)

    if lineage_knn_cols:
        combined_subcluster_knn_col = f"{subcluster_seed_col}_lineage_split_knn_transfer_thresholded"
        combined_subcluster_graph_col = f"{subcluster_seed_col}_lineage_split_graph_transfer_thresholded"
        combined_subcluster_compare_col = f"{subcluster_seed_col}_lineage_split_knn_vs_graph_thresholded"
        combined_subcluster_consensus_col = f"{subcluster_seed_col}_consensus"

        obs = combine_from_columns(obs, lineage_knn_cols, combined_subcluster_knn_col)
        obs = combine_from_columns(obs, lineage_graph_cols, combined_subcluster_graph_col)
        obs, pairwise_agreement = add_method_comparison_from_cols(
            obs,
            left_col=combined_subcluster_knn_col,
            right_col=combined_subcluster_graph_col,
            comparison_col=combined_subcluster_compare_col,
            query_mask=lineage_query_union,
        )
        add_summary_row(
            summary_rows,
            label_key=args.subcluster_label_key,
            lineage="lineage_split_combined",
            method="knn_vs_graph_thresholded",
            reference_cells=subcluster_seeded_cells,
            matched_seed_cells=subcluster_seeded_cells,
            query_cells=int(lineage_query_union.sum()),
            pairwise_agreement=pairwise_agreement,
        )
        obs = combine_prefer_left(
            obs,
            combined_subcluster_knn_col,
            combined_subcluster_graph_col,
            combined_subcluster_consensus_col,
        )

    mdata.obs = obs
    atomic_h5mu_write(mdata, args.output)

    summary_df = pd.DataFrame(summary_rows)
    ensure_parent_dir(args.summary_tsv)
    summary_df.to_csv(args.summary_tsv, sep="\t", index=False)
    print(f"Wrote labeled MuData to {Path(args.output).resolve()}")
    print(f"Wrote transfer summary to {Path(args.summary_tsv).resolve()}")


if __name__ == "__main__":
    main()
