#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import muon as mu
import scanpy as sc
import scvi
import torch
from mudata import MuData
from scipy.sparse import csr_matrix

from _common import ensure_parent_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train MultiVI from the mosaic object and save a MuData artifact.")
    parser.add_argument("--input", default="results/intermediate/mosaic.h5ad", help="Input mosaic h5ad.")
    parser.add_argument("--output", default="results/intermediate/tala_integrated_multivi.h5mu", help="Output h5mu path.")
    parser.add_argument("--model-dir", default="results/models/muon_mvi_mosaic", help="Directory to save the trained MultiVI model.")
    parser.add_argument("--batch-key", default="modality", help="Batch key passed to MultiVI.")
    parser.add_argument("--max-epochs", type=int, default=500, help="Maximum training epochs.")
    parser.add_argument("--n-neighbors", type=int, default=30, help="Neighbor count for latent graph construction.")
    parser.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden clustering resolution.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed.")
    return parser.parse_args()


def _counts_matrix(adata):
    for layer_name in ("counts", "UMI"):
        if layer_name in adata.layers:
            return csr_matrix(adata.layers[layer_name])
    return csr_matrix(adata.X)


def main() -> None:
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")

    adata_mvi = sc.read_h5ad(args.input)
    rna_mask = adata_mvi.var["modality"].isin(["Gene Expression", "expression"])
    atac_mask = adata_mvi.var["modality"].isin(["Peaks", "accessibility"])
    if rna_mask.sum() == 0 or atac_mask.sum() == 0:
        raise ValueError("Could not identify both RNA and ATAC features from adata.var['modality'].")

    adata_rna = adata_mvi[:, rna_mask].copy()
    adata_atac = adata_mvi[:, atac_mask].copy()
    adata_rna = adata_rna[adata_mvi.obs_names, :].copy()
    adata_atac = adata_atac[adata_mvi.obs_names, :].copy()
    adata_rna.X = _counts_matrix(adata_rna)
    adata_atac.X = _counts_matrix(adata_atac)

    mdata = MuData({"rna": adata_rna, "atac": adata_atac})
    mdata.obs = adata_mvi.obs.copy()

    scvi.model.MULTIVI.setup_mudata(
        mdata,
        batch_key=args.batch_key,
        rna_layer=None,
        atac_layer=None,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "atac",
        },
    )
    model = scvi.model.MULTIVI(
        mdata,
        n_genes=mdata.mod["rna"].n_vars,
        n_regions=mdata.mod["atac"].n_vars,
    )
    model.train(max_epochs=args.max_epochs, early_stopping=True)

    mdata.obsm["X_multivi"] = model.get_latent_representation()
    mu.pp.neighbors(mdata, use_rep="X_multivi", n_neighbors=args.n_neighbors, metric="cosine")
    mu.tl.umap(mdata)
    mu.tl.leiden(mdata, resolution=args.leiden_resolution, key_added="leiden_multivi")

    Path(args.model_dir).mkdir(parents=True, exist_ok=True)
    model.save(args.model_dir, save_anndata=True, overwrite=True)
    ensure_parent_dir(args.output)
    mdata.write_h5mu(args.output)
    print(f"Wrote MultiVI MuData to {Path(args.output).resolve()}")
    print(f"Saved MultiVI model to {Path(args.model_dir).resolve()}")


if __name__ == "__main__":
    main()
