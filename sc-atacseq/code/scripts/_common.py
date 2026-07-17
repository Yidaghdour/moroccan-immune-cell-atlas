from __future__ import annotations

import gzip
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


DEFAULT_LIBRARIES = [f"1_{i:03d}" for i in range(1, 18)]
DEFAULT_BRIDGE_LIBRARY = "1_017"
POOL_TO_LIBRARY = {str(i): f"1_{i:03d}" for i in range(1, 17)}
RNA_METADATA_COLUMNS = [
    "donor_id",
    "predicted.celltype.l2",
    "predicted.celltype.l1",
    "Cluster_celltypes",
    "Sub_Cluster_celltypes",
    "Sample.ID",
    "Pool",
    "Batch",
    "Lifestyle",
    "Sex",
    "Age",
    "Ethnicity",
    "Stimulation",
]
ATAC_METADATA_COLUMNS = ["donor_id", "library_batch"]


def ensure_parent_dir(path: str | Path) -> None:
    Path(path).expanduser().resolve().parent.mkdir(parents=True, exist_ok=True)


def resolve_input_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.exists():
        return candidate
    gz_candidate = Path(f"{candidate}.gz")
    if gz_candidate.exists():
        return gz_candidate
    raise FileNotFoundError(f"Could not find input file {candidate} or {gz_candidate}.")


@contextmanager
def materialize_if_gz(path: str | Path):
    resolved = resolve_input_path(path)
    if resolved.suffix != ".gz":
        yield resolved
        return

    original_path = Path(path)
    suffix = "".join(original_path.suffixes) or ".tmp"
    with tempfile.NamedTemporaryFile(prefix="tala_atac_", suffix=suffix, delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        with gzip.open(resolved, "rb") as src, open(temp_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        yield temp_path
    finally:
        temp_path.unlink(missing_ok=True)


def normalize_pool_value(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text or None


def pool_to_library_batch(pool_series: pd.Series) -> pd.Series:
    return pool_series.map(lambda value: POOL_TO_LIBRARY.get(normalize_pool_value(value)))


def library_to_pool_name(library_id: str) -> str:
    numeric = int(library_id.split("_")[1])
    return f"Pool{numeric}"


def load_requantified_atac_matrix(
    matrix_path: str | Path,
    barcodes_path: str | Path,
    peaks_path: str | Path,
) -> ad.AnnData:
    matrix = sc.read_mtx(str(resolve_input_path(matrix_path)))
    barcodes_df = pd.read_csv(resolve_input_path(barcodes_path), header=None, names=["barcode"], sep="\t")
    peaks_df = pd.read_csv(
        resolve_input_path(peaks_path),
        sep="\t",
        header=None,
        names=["chr", "start", "end"],
    )

    peak_names = (
        peaks_df["chr"]
        + "-"
        + peaks_df["start"].astype(str)
        + "-"
        + peaks_df["end"].astype(str)
    )

    adata = ad.AnnData(matrix.T)
    adata.obs_names = barcodes_df["barcode"].astype(str).values
    adata.var_names = peak_names.astype(str).values
    adata.var["chr"] = peaks_df["chr"].astype(str).values
    adata.var["start"] = peaks_df["start"].astype(int).values
    adata.var["end"] = peaks_df["end"].astype(int).values
    adata.var["feature_types"] = "Peaks"
    adata.var["modality"] = "Peaks"
    return adata


def read_donor_assignments(donor_root: str | Path, library_id: str) -> pd.DataFrame:
    donor_root = Path(donor_root)
    donor_path = donor_root / library_to_pool_name(library_id) / "donor_ids.tsv"
    donor_df = pd.read_csv(resolve_input_path(donor_path), sep="\t", index_col=0)
    donor_df.index = donor_df.index.astype(str)
    return donor_df


def attach_donor_assignments(
    adata: ad.AnnData,
    donor_df: pd.DataFrame,
    drop_doublets: bool = True,
    drop_unassigned: bool = False,
) -> ad.AnnData:
    common_barcodes = adata.obs_names.intersection(donor_df.index)
    if len(common_barcodes) == 0:
        raise ValueError("No overlapping barcodes found between ATAC matrix and donor assignments.")

    adata = adata[common_barcodes].copy()
    donor_subset = donor_df.loc[adata.obs_names].copy()
    for column in donor_subset.columns:
        adata.obs[column] = donor_subset[column].values

    if drop_doublets and "donor_id" in adata.obs.columns:
        adata = adata[adata.obs["donor_id"].astype(str) != "doublet"].copy()
    if drop_unassigned and "donor_id" in adata.obs.columns:
        adata = adata[adata.obs["donor_id"].astype(str) != "unassigned"].copy()

    return adata


def prepare_atac_library(
    matrix_dir: str | Path,
    donor_root: str | Path,
    library_id: str,
    include_public_bridge: bool = True,
    drop_doublets: bool = True,
    drop_unassigned: bool = False,
) -> ad.AnnData:
    matrix_dir = Path(matrix_dir)
    adata = load_requantified_atac_matrix(
        matrix_dir / f"{library_id}_q_matrix.mtx",
        matrix_dir / f"{library_id}_q_barcodes.tsv",
        matrix_dir / f"{library_id}_q_peaks.bed",
    )
    adata.obs["library_batch"] = library_id

    if library_id == DEFAULT_BRIDGE_LIBRARY:
        if include_public_bridge:
            adata.obs["donor_id"] = "Public_Donor_1"
        return adata

    donor_df = read_donor_assignments(donor_root, library_id)
    return attach_donor_assignments(
        adata,
        donor_df=donor_df,
        drop_doublets=drop_doublets,
        drop_unassigned=drop_unassigned,
    )


def set_counts_layer(adata: ad.AnnData, layer_name: str = "counts") -> None:
    adata.layers[layer_name] = adata.X.copy()


def sanitize_obs_names(adata: ad.AnnData) -> None:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
