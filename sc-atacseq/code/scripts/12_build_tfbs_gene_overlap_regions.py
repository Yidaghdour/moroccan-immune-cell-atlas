#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import unicodedata
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


CELLTYPE_MAP = {
    "b cells": "B cells",
    "monocytes": "Monocytes",
    "memory t cells": "Memory T cells",
    "cytotoxic t/nk cells": "Cytotoxic T/NK cells",
    "cd4 naive": "Naive CD4 T cells",
    "cd4 naive t cells": "Naive CD4 T cells",
    "naive cd4 t cells": "Naive CD4 T cells",
    "cd8 naive": "Naive CD8 T cells",
    "cd8 naive t cells": "Naive CD8 T cells",
    "naive cd8 t cells": "Naive CD8 T cells",
}

CONDITION_MAP = {
    "unstimulated": "Unstim",
    "unstim": "Unstim",
    "lps": "LPS",
}

SEX_MAP = {
    "female": "F",
    "f": "F",
    "male": "M",
    "m": "M",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build promoter and linked-peak intervals for TFBS scans from the "
            "joint RNA-ATAC target gene list."
        )
    )
    parser.add_argument(
        "--gene-list",
        default="data/manuscript_inputs/Gene_Overlaps_RNA-ATAC_TFAnalysis.csv",
    )
    parser.add_argument("--tss-bed", default="results/annotations/reference_tss.bed.gz")
    parser.add_argument("--peak-links", default="results/annotations/peak_gene_links.filtered.tsv.gz")
    parser.add_argument("--atac-da", default="results/da/stratified/atac_pseudobulk_da_major_stratified.tsv.gz")
    parser.add_argument("--rna-dge", default="results/dge/stratified/rna_pseudobulk_dge_major_stratified.tsv.gz")
    parser.add_argument("--output-regions", default="results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz")
    parser.add_argument("--output-bed", default="results/tfbs_gene_overlap/tfbs_target_regions.bed")
    parser.add_argument("--output-summary", default="results/tfbs_gene_overlap/tfbs_target_region_summary.tsv")
    parser.add_argument("--promoter-flank-bp", type=int, default=2000)
    parser.add_argument("--proximal-bp", type=int, default=50000)
    parser.add_argument(
        "--skip-differential-metadata",
        action="store_true",
        help=(
            "Build regions without joining the external RNA differential-expression and "
            "ATAC differential-accessibility tables. This is the publication Figure 2D mode."
        ),
    )
    parser.add_argument(
        "--override-stimulation",
        choices=["Unstim", "LPS"],
        default=None,
        help="Evaluate all input target rows against this stimulation stratum instead of the CSV Condition column.",
    )
    return parser.parse_args()


def ensure_parent(path: str | Path) -> None:
    Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def strip_accents(value: str) -> str:
    normalized = unicodedata.normalize("NFKD", value)
    return "".join(ch for ch in normalized if not unicodedata.combining(ch))


def norm_key(value: object) -> str:
    text = strip_accents(str(value).strip())
    text = re.sub(r"\s+", " ", text)
    return text.lower()


def normalize_gene_symbol(value: object) -> str:
    return str(value).strip().upper()


def normalize_targets(path: str | Path, override_stimulation: str | None = None) -> pd.DataFrame:
    raw = pd.read_csv(path, sep=";", engine="python")
    raw = raw.loc[:, [column for column in raw.columns if column and not column.startswith("Unnamed")]].copy()
    required = {"Celltype", "Condition", "Sex", "Gene"}
    missing = required.difference(raw.columns)
    if missing:
        raise KeyError(f"Gene list is missing required columns: {sorted(missing)}")

    targets = raw.dropna(subset=["Celltype", "Condition", "Sex", "Gene"]).copy()
    targets["input_celltype"] = targets["Celltype"].astype(str).str.strip()
    targets["input_condition"] = targets["Condition"].astype(str).str.strip()
    targets["input_sex"] = targets["Sex"].astype(str).str.strip()
    targets["target_gene"] = targets["Gene"].map(normalize_gene_symbol)
    targets["major_label"] = targets["input_celltype"].map(lambda x: CELLTYPE_MAP.get(norm_key(x), x))
    targets["stimulation"] = targets["input_condition"].map(lambda x: CONDITION_MAP.get(norm_key(x), x))
    if override_stimulation is not None:
        targets["stimulation"] = override_stimulation
    targets["sex"] = targets["input_sex"].map(lambda x: SEX_MAP.get(norm_key(x), x))
    targets = targets[
        [
            "input_celltype",
            "input_condition",
            "input_sex",
            "major_label",
            "stimulation",
            "sex",
            "target_gene",
        ]
    ].drop_duplicates()
    targets["target_key"] = (
        targets["major_label"].astype(str)
        + "|"
        + targets["stimulation"].astype(str)
        + "|"
        + targets["sex"].astype(str)
        + "|"
        + targets["target_gene"].astype(str)
    )
    return targets.sort_values(["major_label", "stimulation", "sex", "target_gene"]).reset_index(drop=True)


def read_tss(path: str | Path) -> pd.DataFrame:
    columns = ["chrom", "start", "end", "name", "score", "strand"]
    tss = pd.read_csv(path, sep="\t", header=None, names=columns, compression="infer")
    tss["target_gene"] = tss["name"].astype(str).str.split("|", regex=False).str[0].map(normalize_gene_symbol)
    tss["gene_id"] = tss["name"].astype(str).str.split("|", regex=False).str[1]
    tss = tss.dropna(subset=["chrom", "start", "end", "target_gene"]).copy()
    tss["start"] = tss["start"].astype(int)
    tss["end"] = tss["end"].astype(int)
    return tss


def add_target_rna(regions: pd.DataFrame, rna_dge_path: str | Path) -> pd.DataFrame:
    usecols = [
        "gene_symbol",
        "grouping",
        "group_label",
        "major_label",
        "stimulation",
        "sex",
        "stratum_label",
        "cells",
        "donors",
        "rural_donors",
        "urban_donors",
        "logfc_rural_vs_urban",
        "p_value",
        "fdr",
        "mean_log_cpm_rural",
        "mean_log_cpm_urban",
        "mean_log_cpm_overall",
        "detect_frac_rural",
        "detect_frac_urban",
        "detect_frac_max",
        "significant_fdr",
        "significant_fdr_logfc",
    ]
    rna = pd.read_csv(rna_dge_path, sep="\t", usecols=usecols, compression="infer")
    rna = rna.loc[rna["grouping"].astype(str).eq("major")].copy()
    rna["target_gene"] = rna["gene_symbol"].map(normalize_gene_symbol)
    rna = rna.drop(columns=["gene_symbol", "grouping", "group_label"])
    rename = {
        "cells": "target_rna_cells",
        "donors": "target_rna_donors",
        "rural_donors": "target_rna_rural_donors",
        "urban_donors": "target_rna_urban_donors",
        "stratum_label": "target_rna_stratum_label",
        "logfc_rural_vs_urban": "target_rna_logfc_rural_vs_urban",
        "p_value": "target_rna_p_value",
        "fdr": "target_rna_fdr",
        "mean_log_cpm_rural": "target_rna_mean_log_cpm_rural",
        "mean_log_cpm_urban": "target_rna_mean_log_cpm_urban",
        "mean_log_cpm_overall": "target_rna_mean_log_cpm_overall",
        "detect_frac_rural": "target_rna_detect_frac_rural",
        "detect_frac_urban": "target_rna_detect_frac_urban",
        "detect_frac_max": "target_rna_detect_frac_max",
        "significant_fdr": "target_rna_significant_fdr",
        "significant_fdr_logfc": "target_rna_significant_fdr_logfc",
    }
    rna = rna.rename(columns=rename)
    join_cols = ["major_label", "stimulation", "sex", "target_gene"]
    return regions.merge(rna, on=join_cols, how="left")


def add_peak_da(regions: pd.DataFrame, atac_da_path: str | Path) -> pd.DataFrame:
    usecols = [
        "peak_id",
        "grouping",
        "group_label",
        "major_label",
        "stimulation",
        "sex",
        "stratum_label",
        "cells",
        "donors",
        "rural_donors",
        "urban_donors",
        "logfc_rural_vs_urban",
        "p_value",
        "fdr",
        "mean_log_cpm_rural",
        "mean_log_cpm_urban",
        "mean_log_cpm_overall",
        "detect_frac_rural",
        "detect_frac_urban",
        "detect_frac_max",
        "significant_fdr",
        "significant_fdr_logfc",
        "nearest_gene",
        "nearest_gene_distance_to_tss",
    ]
    atac = pd.read_csv(atac_da_path, sep="\t", usecols=usecols, compression="infer")
    atac = atac.loc[atac["grouping"].astype(str).eq("major")].copy()
    atac = atac.drop(columns=["grouping", "group_label"])
    rename = {
        "cells": "atac_cells",
        "donors": "atac_donors",
        "rural_donors": "atac_rural_donors",
        "urban_donors": "atac_urban_donors",
        "stratum_label": "atac_stratum_label",
        "logfc_rural_vs_urban": "atac_logfc_rural_vs_urban",
        "p_value": "atac_p_value",
        "fdr": "atac_fdr",
        "mean_log_cpm_rural": "atac_mean_log_cpm_rural",
        "mean_log_cpm_urban": "atac_mean_log_cpm_urban",
        "mean_log_cpm_overall": "atac_mean_log_cpm_overall",
        "detect_frac_rural": "atac_detect_frac_rural",
        "detect_frac_urban": "atac_detect_frac_urban",
        "detect_frac_max": "atac_detect_frac_max",
        "significant_fdr": "atac_significant_fdr",
        "significant_fdr_logfc": "atac_significant_fdr_logfc",
        "nearest_gene": "atac_nearest_gene",
        "nearest_gene_distance_to_tss": "atac_nearest_gene_distance_to_tss",
    }
    atac = atac.rename(columns=rename)
    join_cols = ["major_label", "stimulation", "sex", "peak_id"]
    return regions.merge(atac, on=join_cols, how="left")


def add_empty_differential_metadata(regions: pd.DataFrame) -> pd.DataFrame:
    out = regions.copy()
    numeric_columns = [
        "atac_logfc_rural_vs_urban",
        "target_rna_logfc_rural_vs_urban",
    ]
    boolean_columns = [
        "atac_significant_fdr_logfc",
        "target_rna_significant_fdr_logfc",
    ]
    for column in numeric_columns:
        out[column] = np.nan
    for column in boolean_columns:
        out[column] = False
    return out


def classify_peak(distance_to_tss: object, promoter_flank_bp: int, proximal_bp: int) -> str:
    if pd.isna(distance_to_tss):
        return "peak_unknown_distance"
    distance = abs(float(distance_to_tss))
    if distance <= promoter_flank_bp:
        return "tss_proximal_peak"
    if distance <= proximal_bp:
        return "proximal_peak"
    return "distal_peak"


def safe_id_part(value: object) -> str:
    text = str(value)
    text = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return text or "NA"


def build_promoter_regions(targets: pd.DataFrame, tss: pd.DataFrame, promoter_flank_bp: int) -> pd.DataFrame:
    promoters = targets.merge(tss, on="target_gene", how="inner")
    if promoters.empty:
        return pd.DataFrame()
    promoters["region_class"] = "promoter"
    promoters["region_source"] = "reference_tss"
    promoters["promoter_flank_bp"] = promoter_flank_bp
    promoters["chrom"] = promoters["chrom"].astype(str)
    promoters["start"] = (promoters["start"].astype(int) - promoter_flank_bp).clip(lower=0)
    promoters["end"] = promoters["end"].astype(int) + promoter_flank_bp
    promoters["peak_id"] = pd.NA
    promoters["link_shared_donors"] = pd.NA
    promoters["link_pearson_r"] = pd.NA
    promoters["link_abs_pearson_r"] = pd.NA
    promoters["link_direction"] = pd.NA
    promoters["link_p_value"] = pd.NA
    promoters["link_fdr"] = pd.NA
    promoters["link_atac_median_logcpm"] = pd.NA
    promoters["link_rna_median_logcpm"] = pd.NA
    promoters["link_distance_to_tss"] = pd.NA
    return promoters


def build_peak_regions(targets: pd.DataFrame, peak_links_path: str | Path, promoter_flank_bp: int, proximal_bp: int) -> pd.DataFrame:
    links = pd.read_csv(peak_links_path, sep="\t", compression="infer")
    links["target_gene"] = links["closest_gene_symbol"].map(normalize_gene_symbol)
    links = links.rename(
        columns={
            "shared_donors": "link_shared_donors",
            "pearson_r": "link_pearson_r",
            "abs_pearson_r": "link_abs_pearson_r",
            "link_direction": "link_direction",
            "p_value": "link_p_value",
            "fdr": "link_fdr",
            "atac_median_logcpm": "link_atac_median_logcpm",
            "rna_median_logcpm": "link_rna_median_logcpm",
            "distance_to_tss": "link_distance_to_tss",
        }
    )
    join_cols = ["major_label", "target_gene"]
    peaks = targets.merge(links, on=join_cols, how="inner", suffixes=("", "_link"))
    if peaks.empty:
        return pd.DataFrame()
    peaks["region_class"] = peaks["link_distance_to_tss"].map(
        lambda value: classify_peak(value, promoter_flank_bp=promoter_flank_bp, proximal_bp=proximal_bp)
    )
    peaks["region_source"] = "filtered_peak_gene_link"
    peaks["promoter_flank_bp"] = promoter_flank_bp
    peaks["gene_id"] = peaks["closest_gene_id"]
    peaks["name"] = peaks["closest_gene_symbol"].astype(str) + "|" + peaks["closest_gene_id"].astype(str)
    peaks["score"] = 0
    peaks["strand"] = peaks["closest_gene_strand"]
    return peaks


def make_region_ids(regions: pd.DataFrame) -> pd.DataFrame:
    out = regions.copy()
    parts = []
    for idx, row in out.reset_index(drop=True).iterrows():
        base = [
            f"R{idx + 1:06d}",
            safe_id_part(row["region_class"]),
            safe_id_part(row["major_label"]),
            safe_id_part(row["stimulation"]),
            safe_id_part(row["sex"]),
            safe_id_part(row["target_gene"]),
        ]
        if pd.notna(row.get("peak_id")):
            base.append(safe_id_part(row["peak_id"]))
        else:
            base.append(f"{row['chrom']}_{int(row['start'])}_{int(row['end'])}")
        parts.append("__".join(base))
    out.insert(0, "region_id", parts)
    return out


def median_or_nan(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    if not numeric.notna().any():
        return float("nan")
    return float(numeric.median())


def main() -> None:
    args = parse_args()
    targets = normalize_targets(args.gene_list, override_stimulation=args.override_stimulation)
    tss = read_tss(args.tss_bed)

    promoter_regions = build_promoter_regions(targets, tss, args.promoter_flank_bp)
    peak_regions = build_peak_regions(targets, args.peak_links, args.promoter_flank_bp, args.proximal_bp)
    if promoter_regions.empty and peak_regions.empty:
        raise ValueError("No promoter or peak regions were built from the supplied target list.")

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="The behavior of DataFrame concatenation with empty or all-NA entries is deprecated",
            category=FutureWarning,
        )
        regions = pd.concat([promoter_regions, peak_regions], ignore_index=True, sort=False)
    regions["chrom"] = regions["chrom"].astype(str)
    regions["start"] = regions["start"].astype(int)
    regions["end"] = regions["end"].astype(int)
    regions = regions.loc[regions["end"] > regions["start"]].copy()

    if args.skip_differential_metadata:
        regions = add_empty_differential_metadata(regions)
    else:
        regions = add_peak_da(regions, args.atac_da)
        regions = add_target_rna(regions, args.rna_dge)
    if args.skip_differential_metadata:
        regions["link_atac_target_logfc_direction_consistent"] = pd.NA
    else:
        atac_logfc = pd.to_numeric(regions["atac_logfc_rural_vs_urban"], errors="coerce")
        target_logfc = pd.to_numeric(regions["target_rna_logfc_rural_vs_urban"], errors="coerce")
        link_corr = pd.to_numeric(regions["link_pearson_r"], errors="coerce")
        direction_consistent = pd.Series(
            np.sign(atac_logfc * target_logfc * link_corr) > 0,
            index=regions.index,
            dtype="object",
        )
        direction_consistent.loc[regions["region_class"].eq("promoter")] = pd.NA
        regions["link_atac_target_logfc_direction_consistent"] = direction_consistent

    key_cols = [
        "region_class",
        "chrom",
        "start",
        "end",
        "major_label",
        "stimulation",
        "sex",
        "target_gene",
        "peak_id",
        "link_pearson_r",
        "link_fdr",
    ]
    regions = regions.drop_duplicates(subset=[column for column in key_cols if column in regions.columns]).copy()
    regions = regions.sort_values(
        ["major_label", "stimulation", "sex", "target_gene", "region_class", "chrom", "start", "end"]
    ).reset_index(drop=True)
    regions = make_region_ids(regions)

    bed = regions[["chrom", "start", "end", "region_id", "target_gene", "strand"]].copy()
    bed["score"] = 0
    bed = bed[["chrom", "start", "end", "region_id", "score", "strand"]]

    summary = (
        regions.groupby(["major_label", "stimulation", "sex", "region_class"], dropna=False)
        .agg(
            regions=("region_id", "nunique"),
            target_genes=("target_gene", "nunique"),
            peaks=("peak_id", lambda x: x.dropna().nunique()),
            linked_regions=("link_fdr", lambda x: x.notna().sum()),
            median_link_abs_pearson=("link_abs_pearson_r", median_or_nan),
            median_link_shared_donors=("link_shared_donors", median_or_nan),
            target_rna_significant=("target_rna_significant_fdr_logfc", lambda x: int(pd.Series(x).fillna(False).sum())),
            atac_da_significant=("atac_significant_fdr_logfc", lambda x: int(pd.Series(x).fillna(False).sum())),
        )
        .reset_index()
        .sort_values(["major_label", "stimulation", "sex", "region_class"])
    )

    ensure_parent(args.output_regions)
    ensure_parent(args.output_bed)
    ensure_parent(args.output_summary)
    regions.to_csv(args.output_regions, sep="\t", index=False, compression="infer")
    bed.to_csv(args.output_bed, sep="\t", index=False, header=False)
    summary.to_csv(args.output_summary, sep="\t", index=False)

    missing_promoter_genes = sorted(set(targets["target_gene"]) - set(tss["target_gene"]))
    print(f"Target strata rows: {targets.shape[0]}")
    print(f"Unique target genes: {targets['target_gene'].nunique()}")
    print(f"Regions written: {regions.shape[0]}")
    print(summary.to_string(index=False))
    if missing_promoter_genes:
        print("Genes without a reference TSS:", ", ".join(missing_promoter_genes))
    print(f"Wrote {Path(args.output_regions).resolve()}")
    print(f"Wrote {Path(args.output_bed).resolve()}")
    print(f"Wrote {Path(args.output_summary).resolve()}")


if __name__ == "__main__":
    main()
