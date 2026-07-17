#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ALLOWED_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a reproducible reference promoter background for TFBS enrichment.")
    parser.add_argument("--reference-genes", default="results/annotations/reference_genes.bed.gz")
    parser.add_argument(
        "--target-genes",
        default="data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv",
    )
    parser.add_argument("--output", default="results/tfbs_gene_overlap_integrated/background_reference_promoters.tsv.gz")
    parser.add_argument("--flank-bp", type=int, default=2000)
    parser.add_argument("--max-background-promoters", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ref = pd.read_csv(
        args.reference_genes,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
        compression="infer",
    )
    ref = ref.loc[ref["chrom"].isin(ALLOWED_CHROMS)].copy()
    ref[["gene_symbol", "gene_id"]] = ref["name"].astype(str).str.split("|", n=1, expand=True)
    ref["gene_symbol_upper"] = ref["gene_symbol"].astype(str).str.upper()
    ref = ref.drop_duplicates("gene_symbol_upper").copy()

    targets = pd.read_csv(args.target_genes, sep=";")
    target_genes = set(targets["Gene"].dropna().astype(str).str.upper())
    ref = ref.loc[~ref["gene_symbol_upper"].isin(target_genes)].copy()

    if args.max_background_promoters > 0 and ref.shape[0] > args.max_background_promoters:
        ref = ref.sample(n=args.max_background_promoters, random_state=args.seed).sort_values(
            ["chrom", "start", "end", "gene_symbol_upper"]
        )

    tss = ref["end"].where(ref["strand"].eq("-"), ref["start"]).astype(int)
    out = pd.DataFrame(
        {
            "region_id": (
                "BG_PROMOTER__"
                + ref["gene_symbol_upper"].astype(str)
                + "__"
                + ref["gene_id"].fillna("").astype(str)
                + "__"
                + ref["chrom"].astype(str)
                + "_"
                + (tss - args.flank_bp).clip(lower=0).astype(str)
                + "_"
                + (tss + args.flank_bp + 1).astype(str)
            ),
            "chrom": ref["chrom"].astype(str),
            "start": (tss - args.flank_bp).clip(lower=0).astype(int),
            "end": (tss + args.flank_bp + 1).astype(int),
            "gene_symbol": ref["gene_symbol"].astype(str),
            "gene_id": ref["gene_id"].fillna("").astype(str),
            "strand": ref["strand"].astype(str),
            "region_class": "reference_promoter_background",
            "promoter_flank_bp": args.flank_bp,
        }
    )
    out = out.loc[out["end"].gt(out["start"])].drop_duplicates("region_id")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False, compression="gzip")
    print(f"Wrote {out.shape[0]} reference promoter background regions to {output}")


if __name__ == "__main__":
    main()
