#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collapse the full major-cell-type RNA summary to the gene-level expression "
            "lookup used to define the Figure 2D promoter background."
        )
    )
    parser.add_argument(
        "--input",
        default="results/dge/stratified/rna_pseudobulk_dge_major_stratified.tsv.gz",
    )
    parser.add_argument(
        "--output",
        default="data/manuscript_inputs/figure2d_reference_expression.tsv.gz",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    table = pd.read_csv(
        args.input,
        sep="\t",
        compression="infer",
        usecols=["gene_symbol", "grouping", "mean_log_cpm_overall"],
    )
    table = table.loc[table["grouping"].astype(str).eq("major")].copy()
    table["gene_symbol"] = table["gene_symbol"].astype(str).str.strip().str.upper()
    table["mean_log_cpm_overall"] = pd.to_numeric(
        table["mean_log_cpm_overall"], errors="coerce"
    )
    table = table.dropna(subset=["gene_symbol", "mean_log_cpm_overall"])
    lookup = (
        table.groupby("gene_symbol", as_index=False, sort=True)["mean_log_cpm_overall"]
        .median()
        .sort_values("gene_symbol")
        .reset_index(drop=True)
    )
    if lookup.empty:
        raise ValueError("No major-cell-type expression values were available")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    lookup.to_csv(
        output,
        sep="\t",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )
    print(
        f"Collapsed {table.shape[0]} major-cell-type rows to {lookup.shape[0]} genes; "
        f"wrote {output.resolve()}"
    )


if __name__ == "__main__":
    main()
