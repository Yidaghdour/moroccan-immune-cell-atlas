#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc
import scvi
import torch

from _common import ensure_parent_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train PEAKVI on the filtered ATAC object.")
    parser.add_argument("--input", default="results/intermediate/filtered_atac.h5ad", help="Input ATAC h5ad.")
    parser.add_argument("--output", default="results/intermediate/peakvi_17.h5ad", help="Output h5ad with latent space and clustering.")
    parser.add_argument("--model-dir", default="results/models/all_batches_17", help="Directory to store the trained PEAKVI model.")
    parser.add_argument("--batch-key", default="library_batch", help="Batch covariate for PEAKVI.")
    parser.add_argument("--counts-layer", default="counts", help="Layer containing raw counts.")
    parser.add_argument("--max-epochs", type=int, default=500, help="Maximum training epochs.")
    parser.add_argument("--n-neighbors", type=int, default=30, help="Neighbor count for graph construction.")
    parser.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden clustering resolution.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")

    adata = sc.read_h5ad(args.input)
    layer = args.counts_layer if args.counts_layer in adata.layers else None

    scvi.model.PEAKVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=args.batch_key,
    )
    model = scvi.model.PEAKVI(adata)
    model.train(max_epochs=args.max_epochs, early_stopping=True)

    adata.obsm["X_peakvi"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_peakvi", n_neighbors=args.n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.leiden_resolution)

    ensure_parent_dir(args.output)
    adata.write_h5ad(args.output)
    Path(args.model_dir).mkdir(parents=True, exist_ok=True)
    model.save(args.model_dir, overwrite=True)
    print(f"Wrote integrated ATAC object to {Path(args.output).resolve()}")
    print(f"Saved PEAKVI model to {Path(args.model_dir).resolve()}")


if __name__ == "__main__":
    main()
