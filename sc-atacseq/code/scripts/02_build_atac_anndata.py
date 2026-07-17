#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import scanpy as sc

from _common import DEFAULT_LIBRARIES, ensure_parent_dir, prepare_atac_library, sanitize_obs_names, set_counts_layer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a single ATAC AnnData from requantified peak matrices and donor assignments.",
    )
    parser.add_argument("--matrix-dir", default="data/processed/atac_requantified", help="Directory containing *_q_matrix.mtx / *_q_barcodes.tsv / *_q_peaks.bed files.")
    parser.add_argument("--donor-root", default="data/metadata/V", help="Directory containing Pool*/donor_ids.tsv files.")
    parser.add_argument("--output", default="results/intermediate/filtered_atac.h5ad", help="Output h5ad path.")
    parser.add_argument("--libraries", nargs="+", default=DEFAULT_LIBRARIES, help="Library IDs to include.")
    parser.add_argument("--min-peaks-per-cell", type=int, default=1000, help="Minimum detected peaks per cell.")
    parser.add_argument("--min-cells-per-peak", type=int, default=10, help="Minimum cells per peak.")
    parser.add_argument("--keep-unassigned", action="store_true", help="Keep donor assignments marked as unassigned.")
    parser.add_argument("--exclude-public-bridge", action="store_true", help="Exclude library 1_017 from the concatenated object.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    libraries = list(args.libraries)
    if args.exclude_public_bridge:
        libraries = [library for library in libraries if library != "1_017"]

    adatas: dict[str, ad.AnnData] = {}
    for library_id in libraries:
        adata = prepare_atac_library(
            matrix_dir=args.matrix_dir,
            donor_root=args.donor_root,
            library_id=library_id,
            include_public_bridge=not args.exclude_public_bridge,
            drop_doublets=True,
            drop_unassigned=not args.keep_unassigned,
        )
        sanitize_obs_names(adata)
        adatas[library_id] = adata

    combined = ad.concat(adatas, label="library_key", keys=list(adatas.keys()), merge="unique")
    sanitize_obs_names(combined)
    set_counts_layer(combined)
    sc.pp.calculate_qc_metrics(combined, inplace=True)
    sc.pp.filter_cells(combined, min_genes=args.min_peaks_per_cell)
    sc.pp.filter_genes(combined, min_cells=args.min_cells_per_peak)

    ensure_parent_dir(args.output)
    combined.write_h5ad(args.output)
    print(f"Wrote {combined.n_obs} cells x {combined.n_vars} peaks to {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
