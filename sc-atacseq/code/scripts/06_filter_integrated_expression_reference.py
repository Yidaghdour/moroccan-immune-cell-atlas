#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import mudata as md
import pandas as pd

from _common import ensure_parent_dir, materialize_if_gz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Restrict the integrated MultiVI MuData object to expression cells whose "
            "barcodes exist in the final Morocco scRNA-seq reference, while preserving "
            "all non-expression cells."
        )
    )
    parser.add_argument(
        "--input",
        default="results/intermediate/tala_integrated_multivi.h5mu",
        help="Integrated MultiVI MuData object before Morocco label transfer.",
    )
    parser.add_argument(
        "--morocco-h5ad",
        default="data/raw/reference/Morocco_scRNA-seq.h5ad",
        help="Final Morocco scRNA-seq reference used to define the keep set.",
    )
    parser.add_argument(
        "--output",
        default="results/intermediate/tala_integrated_multivi_morocco_expr_only.h5mu",
        help="Filtered MuData output path.",
    )
    parser.add_argument(
        "--summary-tsv",
        default="results/intermediate/tala_integrated_multivi_morocco_expr_only_filter_summary.tsv",
        help="Summary TSV describing kept and dropped cells by modality.",
    )
    parser.add_argument(
        "--dropped-expression-output",
        default="results/intermediate/tala_integrated_multivi_non_morocco_expression_cells.tsv.gz",
        help="TSV.gz list of expression cells removed by the reference filter.",
    )
    parser.add_argument(
        "--reference-modality",
        nargs="+",
        default=["expression"],
        help="Modalities that should be checked against the Morocco reference.",
    )
    return parser.parse_args()


def make_morocco_join_key(index: pd.Index) -> pd.Index:
    return (
        index.astype(str)
        .str.replace("_expression", "", regex=False)
        .str.replace("_accessibility", "", regex=False)
        .str.replace("_paired", "", regex=False)
    )


def main() -> None:
    args = parse_args()

    mdata = md.read_h5mu(args.input, backed="r")
    with materialize_if_gz(args.morocco_h5ad) as morocco_path:
        morocco = ad.read_h5ad(morocco_path, backed="r")

    if "modality" not in mdata.obs.columns:
        raise KeyError("mdata.obs['modality'] is required to filter expression cells against the Morocco reference.")

    obs = mdata.obs.copy()
    obs["modality"] = obs["modality"].astype(str)
    join_key = pd.Series(make_morocco_join_key(obs.index), index=obs.index)
    morocco_index = pd.Index(morocco.obs_names.astype(str))

    ref_mask = obs["modality"].isin(args.reference_modality)
    matched_ref_mask = ref_mask & join_key.isin(morocco_index)
    keep_mask = (~ref_mask) | matched_ref_mask
    dropped_expr = obs.loc[ref_mask & ~matched_ref_mask, ["modality"]].copy()
    dropped_expr.insert(0, "join_key", join_key.loc[dropped_expr.index].to_numpy())
    dropped_expr = dropped_expr.reset_index(names="cell_id")

    summary_rows = []
    modalities = obs["modality"].value_counts().sort_index().index.tolist()
    for modality in modalities:
        original_count = int(obs["modality"].eq(modality).sum())
        kept_count = int((obs["modality"].eq(modality) & keep_mask).sum())
        dropped_count = original_count - kept_count
        summary_rows.append(
            {
                "modality": modality,
                "original_cells": original_count,
                "kept_cells": kept_count,
                "dropped_cells": dropped_count,
            }
        )

    summary_rows.extend(
        [
            {
                "modality": "Morocco_reference",
                "original_cells": int(morocco.n_obs),
                "kept_cells": int(matched_ref_mask.sum()),
                "dropped_cells": int(morocco.n_obs - matched_ref_mask.sum()),
            },
            {
                "modality": "Total",
                "original_cells": int(obs.shape[0]),
                "kept_cells": int(keep_mask.sum()),
                "dropped_cells": int((~keep_mask).sum()),
            },
        ]
    )
    summary = pd.DataFrame(summary_rows)

    ensure_parent_dir(args.output)
    ensure_parent_dir(args.summary_tsv)
    ensure_parent_dir(args.dropped_expression_output)
    filtered = mdata[keep_mask.to_numpy()]
    filtered.write_h5mu(args.output)
    summary.to_csv(args.summary_tsv, sep="\t", index=False)
    dropped_expr.to_csv(args.dropped_expression_output, sep="\t", index=False)
    mdata.file.close()

    print(f"Wrote filtered MuData to {Path(args.output).resolve()}")
    print(f"Wrote filter summary to {Path(args.summary_tsv).resolve()}")
    print(f"Wrote dropped expression cell list to {Path(args.dropped_expression_output).resolve()}")
    print()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
