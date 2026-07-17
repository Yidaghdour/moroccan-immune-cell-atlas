#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


ALLOWED_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare a non-target linked-peak background for TFBS enrichment."
    )
    parser.add_argument("--peak-links", default="results/annotations/peak_gene_links.filtered.tsv.gz")
    parser.add_argument("--target-regions", default="results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz")
    parser.add_argument(
        "--target-genes",
        default="data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv",
    )
    parser.add_argument("--output", default="results/tfbs_gene_overlap_integrated/background_linked_peaks.tsv.gz")
    parser.add_argument("--promoter-flank-bp", type=int, default=2000)
    parser.add_argument("--proximal-bp", type=int, default=50000)
    parser.add_argument("--max-per-celltype-class", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def normalize_gene_symbol(value: object) -> str:
    return str(value).strip().upper()


def safe_id_part(value: object) -> str:
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return text or "NA"


def classify_peak(distance_to_tss: object, promoter_flank_bp: int, proximal_bp: int) -> str:
    if pd.isna(distance_to_tss):
        return "peak_unknown_distance"
    distance = abs(float(distance_to_tss))
    if distance <= promoter_flank_bp:
        return "tss_proximal_peak"
    if distance <= proximal_bp:
        return "proximal_peak"
    return "distal_peak"


def main() -> None:
    args = parse_args()
    target_genes = pd.read_csv(args.target_genes, sep=";")["Gene"].dropna().map(normalize_gene_symbol)
    target_genes = set(target_genes)

    target_regions = pd.read_csv(
        args.target_regions,
        sep="\t",
        compression="infer",
        usecols=["peak_id", "region_class"],
    )
    target_peak_ids = set(target_regions["peak_id"].dropna().astype(str))

    links = pd.read_csv(args.peak_links, sep="\t", compression="infer")
    links = links.loc[links["chrom"].isin(ALLOWED_CHROMS)].copy()
    links["target_gene"] = links["closest_gene_symbol"].map(normalize_gene_symbol)
    links["region_class"] = links["distance_to_tss"].map(
        lambda value: classify_peak(value, args.promoter_flank_bp, args.proximal_bp)
    )
    links = links.loc[
        links["region_class"].isin(["tss_proximal_peak", "proximal_peak", "distal_peak"])
        & ~links["target_gene"].isin(target_genes)
        & ~links["peak_id"].astype(str).isin(target_peak_ids)
    ].copy()
    links = links.drop_duplicates(["major_label", "region_class", "peak_id"])

    if args.max_per_celltype_class > 0:
        links = (
            links.groupby(["major_label", "region_class"], group_keys=False, dropna=False)
            .apply(
                lambda group: group.sample(
                    n=min(args.max_per_celltype_class, group.shape[0]),
                    random_state=args.seed,
                )
            )
            .reset_index(drop=True)
        )

    links = links.sort_values(["major_label", "region_class", "chrom", "start", "end", "peak_id"]).reset_index(drop=True)
    region_id = (
        "BG_LINK__"
        + links["major_label"].map(safe_id_part)
        + "__"
        + links["region_class"].map(safe_id_part)
        + "__"
        + links["peak_id"].map(safe_id_part)
    )
    out = pd.DataFrame(
        {
            "region_id": region_id,
            "chrom": links["chrom"].astype(str),
            "start": links["start"].astype(int),
            "end": links["end"].astype(int),
            "major_label": links["major_label"].astype(str),
            "region_class": links["region_class"].astype(str),
            "peak_id": links["peak_id"].astype(str),
            "closest_gene_symbol": links["closest_gene_symbol"].astype(str),
            "distance_to_tss": links["distance_to_tss"],
            "link_fdr": links["fdr"],
            "background_source": "filtered_non_target_peak_gene_links",
        }
    )
    out = out.loc[out["end"].gt(out["start"])].drop_duplicates("region_id")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False, compression="gzip")
    print(f"Wrote {out.shape[0]} linked-peak background regions to {output}")
    print(
        out.groupby(["major_label", "region_class"], dropna=False)["region_id"]
        .nunique()
        .rename("background_regions")
        .reset_index()
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
