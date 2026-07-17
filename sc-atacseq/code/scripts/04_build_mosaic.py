#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import scanpy as sc
import scvi

from _common import (
    ATAC_METADATA_COLUMNS,
    DEFAULT_BRIDGE_LIBRARY,
    RNA_METADATA_COLUMNS,
    ensure_parent_dir,
    materialize_if_gz,
    pool_to_library_batch,
    sanitize_obs_names,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the mosaic AnnData used for MultiVI.")
    parser.add_argument("--atac-h5ad", default="results/intermediate/peakvi_17.h5ad", help="ATAC h5ad containing library_batch and donor_id.")
    parser.add_argument("--rna-h5ad", default="data/raw/reference/merged_QC_split_SCT_merged_harmony-batch_Dims22_slim.h5ad", help="RNA h5ad with cell type labels.")
    parser.add_argument("--bridge-multiome-h5", default="data/raw/reference/10k_PBMC_Multiome_nextgem_Chromium_X_filtered_feature_bc_matrix.h5", help="Public 10x multiome bridge file.")
    parser.add_argument("--output", default="results/intermediate/mosaic.h5ad", help="Output mosaic h5ad.")
    parser.add_argument("--bridge-library", default=DEFAULT_BRIDGE_LIBRARY, help="ATAC library used as the paired bridge.")
    return parser.parse_args()


def prepare_bridge(pbmc_multi: ad.AnnData, atac_adata: ad.AnnData, bridge_library: str) -> ad.AnnData:
    bridge_atac = atac_adata[atac_adata.obs["library_batch"].astype(str) == bridge_library].copy()
    if bridge_atac.n_obs == 0:
        raise ValueError(f"No ATAC cells found for bridge library {bridge_library}.")
    bridge_atac.var["modality"] = "Peaks"

    feature_column = "feature_types" if "feature_types" in pbmc_multi.var.columns else "modality"
    bridge_rna = pbmc_multi[:, pbmc_multi.var[feature_column].astype(str) == "Gene Expression"].copy()
    bridge_rna.var["modality"] = "Gene Expression"

    common_cells = bridge_rna.obs_names.intersection(bridge_atac.obs_names)
    if len(common_cells) == 0:
        raise ValueError("No shared cell barcodes found between bridge RNA and bridge ATAC.")

    bridge_rna = bridge_rna[common_cells].copy()
    bridge_atac = bridge_atac[common_cells].copy()
    bridge_final = ad.concat([bridge_rna, bridge_atac], axis=1, uns_merge="unique")
    bridge_final.obs["batch_id"] = "Bridge_Multiome"
    bridge_final.obs["donor_id"] = "Public_Donor_1"
    bridge_final.obs["library_batch"] = "Public_10x_Multiome"
    return bridge_final


def harmonize_metadata(tala_rna: ad.AnnData, user_atac: ad.AnnData, bridge_final: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData, ad.AnnData]:
    tala_rna = tala_rna.copy()
    user_atac = user_atac.copy()
    bridge_final = bridge_final.copy()

    tala_rna.var["modality"] = "Gene Expression"
    tala_rna.obs["batch_id"] = "User_RNA"
    user_atac.var["modality"] = "Peaks"
    user_atac.obs["batch_id"] = "User_ATAC"

    for column in RNA_METADATA_COLUMNS:
        if column not in bridge_final.obs.columns:
            bridge_final.obs[column] = "Public_n/a"

    if "Pool" in tala_rna.obs.columns:
        tala_rna.obs["library_batch"] = pool_to_library_batch(tala_rna.obs["Pool"])
    elif "library_batch" not in tala_rna.obs.columns:
        raise KeyError("RNA AnnData must contain either 'Pool' or 'library_batch' in .obs.")

    keep_rna = [column for column in RNA_METADATA_COLUMNS if column in tala_rna.obs.columns]
    tala_rna.obs = tala_rna.obs[keep_rna + ["batch_id", "library_batch"]].copy()

    keep_atac = [column for column in ATAC_METADATA_COLUMNS if column in user_atac.obs.columns]
    user_atac.obs = user_atac.obs[keep_atac + ["batch_id"]].copy()

    donor_covariates = ["Batch", "Lifestyle", "Sex", "Age", "Ethnicity", "Stimulation"]
    available_covariates = [column for column in donor_covariates if column in tala_rna.obs.columns]
    summary_df = (
        tala_rna.obs[["donor_id", "library_batch"] + available_covariates]
        .dropna(subset=["donor_id", "library_batch"])
        .groupby(["donor_id", "library_batch"], dropna=False)
        .first()
        .reset_index()
    )
    summary_df["gsi"] = summary_df["donor_id"].astype(str) + "_" + summary_df["library_batch"].astype(str)
    summary_df = summary_df.set_index("gsi")

    user_atac.obs["gsi"] = user_atac.obs["donor_id"].astype(str) + "_" + user_atac.obs["library_batch"].astype(str)
    if available_covariates:
        user_atac.obs = user_atac.obs.join(summary_df[available_covariates], on="gsi")

    for adata in (bridge_final, tala_rna, user_atac):
        sanitize_obs_names(adata)

    return tala_rna, user_atac, bridge_final


def main() -> None:
    args = parse_args()

    with materialize_if_gz(args.atac_h5ad) as atac_path, materialize_if_gz(args.rna_h5ad) as rna_path, materialize_if_gz(
        args.bridge_multiome_h5
    ) as bridge_path:
        atac_adata = sc.read_h5ad(atac_path)
        tala_rna = sc.read_h5ad(rna_path)
        pbmc_multi = sc.read_10x_h5(bridge_path, gex_only=False)

    bridge_final = prepare_bridge(pbmc_multi, atac_adata, args.bridge_library)
    user_atac = atac_adata[atac_adata.obs["library_batch"].astype(str) != args.bridge_library].copy()
    tala_rna, user_atac, bridge_final = harmonize_metadata(tala_rna, user_atac, bridge_final)

    adata_mvi = scvi.data.organize_multiome_anndatas(
        multi_anndata=bridge_final,
        rna_anndata=tala_rna,
        atac_anndata=user_atac,
    )
    adata_mvi.obs = adata_mvi.obs.copy()
    adata_mvi.obs.fillna("n/a", inplace=True)
    sanitize_obs_names(adata_mvi)

    ensure_parent_dir(args.output)
    adata_mvi.write_h5ad(args.output)
    print(f"Wrote mosaic object to {Path(args.output).resolve()}")
    print(adata_mvi.obs["modality"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
