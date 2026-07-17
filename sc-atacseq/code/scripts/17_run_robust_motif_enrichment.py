#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from scipy.stats import fisher_exact


PEAK_CLASSES = ["tss_proximal_peak", "proximal_peak", "distal_peak"]
CLASS_LABELS = {
    "promoter": "Reference promoter (TSS +/-2 kb)",
    "tss_proximal_peak": "Linked ATAC peak <=2 kb",
    "proximal_peak": "Linked ATAC peak 2-50 kb",
    "distal_peak": "Linked ATAC peak >50 kb",
}
ALLOWED_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
STAT_IRF = {
    "IRF1",
    "IRF2",
    "IRF3",
    "IRF4",
    "IRF5",
    "IRF6",
    "IRF7",
    "IRF8",
    "IRF9",
    "STAT1",
    "STAT2",
    "STAT3",
    "STAT4",
    "STAT6",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run covariate-adjusted and repeatedly matched motif enrichment tests using "
            "mutually exclusive linked-peak distance bins and a separate promoter analysis."
        )
    )
    parser.add_argument("--target-regions", default="results/tfbs_gene_overlap/tfbs_target_regions.tsv.gz")
    parser.add_argument("--peak-links", default="results/annotations/peak_gene_links.filtered.tsv.gz")
    parser.add_argument("--target-hits", default="results/tfbs_gene_overlap/tfbs_motif_hits.tsv.gz")
    parser.add_argument(
        "--background-peak-regions",
        default="results/tfbs_gene_overlap_integrated/background_linked_peaks.tsv.gz",
    )
    parser.add_argument(
        "--background-peak-hits",
        default="results/tfbs_gene_overlap_integrated/background_linked_peak_motif_hits.tsv.gz",
    )
    parser.add_argument(
        "--background-promoters",
        default="results/tfbs_gene_overlap_integrated/background_reference_promoters.tsv.gz",
    )
    parser.add_argument(
        "--background-promoter-hits",
        default="results/tfbs_gene_overlap_integrated/background_reference_promoter_motif_hits.tsv.gz",
    )
    parser.add_argument(
        "--reference-expression",
        "--rna-dge",
        dest="reference_expression",
        default="data/manuscript_inputs/figure2d_reference_expression.tsv.gz",
        help=(
            "Gene-level expression lookup used to define the eligible promoter background. "
            "A full stratified DGE table is also accepted through the legacy --rna-dge alias."
        ),
    )
    parser.add_argument(
        "--fasta-dir",
        default="results/tfbs_gene_overlap_integrated/standard_motif_enrichment/fastas",
    )
    parser.add_argument(
        "--output-dir",
        default="results/tfbs_gene_overlap_integrated/robust_motif_validation",
    )
    parser.add_argument("--match-repeats", type=int, default=100)
    parser.add_argument("--controls-per-target", type=int, default=1)
    parser.add_argument("--min-target-hits", type=int, default=2)
    parser.add_argument("--min-background-hits", type=int, default=5)
    parser.add_argument("--chunksize", type=int, default=500_000)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def bh_fdr(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    valid = numeric.notna()
    out = pd.Series(np.nan, index=values.index, dtype=float)
    if not valid.any():
        return out
    p = numeric.loc[valid].clip(lower=0.0, upper=1.0).to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[valid] = restored
    return out


def parse_fasta_features(path: Path) -> dict[tuple[str, int, int], tuple[float, float, int]]:
    features: dict[tuple[str, int, int], tuple[float, float, int]] = {}
    name: str | None = None
    chunks: list[str] = []

    def commit() -> None:
        if name is None:
            return
        parts = name.split("__")
        if len(parts) < 3:
            raise ValueError(f"Cannot parse FASTA coordinates from {name!r} in {path}")
        chrom, start_text, end_text = parts[-3:]
        sequence = "".join(chunks).upper()
        valid = sum(sequence.count(base) for base in "ACGT")
        gc = (sequence.count("G") + sequence.count("C")) / valid if valid else np.nan
        cpg = sequence.count("CG") / max(valid - 1, 1) if valid else np.nan
        features[(chrom, int(start_text), int(end_text))] = (gc, cpg, len(sequence))

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    commit()
    return features


def load_sequence_features(fasta_dir: Path) -> dict[tuple[str, int, int], tuple[float, float, int]]:
    paths = [
        fasta_dir / "linked_pooled_all.target.fa",
        fasta_dir / "linked_pooled_all.background.fa",
        fasta_dir / "promoter_all.target.fa",
        fasta_dir / "promoter_all.background.fa",
    ]
    features: dict[tuple[str, int, int], tuple[float, float, int]] = {}
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
        features.update(parse_fasta_features(path))
    return features


def add_sequence_features(
    table: pd.DataFrame,
    features: dict[tuple[str, int, int], tuple[float, float, int]],
) -> pd.DataFrame:
    out = table.copy()
    keys = list(zip(out["chrom"].astype(str), out["start"].astype(int), out["end"].astype(int)))
    values = [features.get(key, (np.nan, np.nan, np.nan)) for key in keys]
    out[["gc_fraction", "cpg_density", "sequence_length"]] = pd.DataFrame(values, index=out.index)
    missing = int(out["gc_fraction"].isna().sum())
    if missing:
        raise ValueError(f"Sequence features were unavailable for {missing} regions")
    return out


def normalize_tf_symbols(motif_name: object, motif_db: object) -> list[str]:
    name = str(motif_name).strip()
    if str(motif_db).startswith("HOCOMOCO"):
        name = name.split(".", 1)[0]
    symbols: list[str] = []
    for token in name.replace("::", "|").split("|"):
        token = re.sub(r"\(.*$", "", token).strip().upper()
        if re.fullmatch(r"[A-Z][A-Z0-9_-]{1,30}", token):
            symbols.append(token)
    return sorted(set(symbols))


def build_motif_catalog(target_hits_path: Path) -> tuple[pd.DataFrame, list[str]]:
    identities = pd.read_csv(
        target_hits_path,
        sep="\t",
        compression="infer",
        usecols=["motif_db", "motif_id", "motif_name"],
    ).drop_duplicates()
    rows: list[dict[str, str]] = []
    for row in identities.itertuples(index=False):
        db_short = "JASPAR" if str(row.motif_db).startswith("JASPAR") else "HOCOMOCO"
        for symbol in normalize_tf_symbols(row.motif_name, row.motif_db):
            rows.append(
                {
                    "motif_db": str(row.motif_db),
                    "motif_id": str(row.motif_id),
                    "motif_name": str(row.motif_name),
                    "database": db_short,
                    "tf_symbol": symbol,
                }
            )
    catalog = pd.DataFrame(rows).drop_duplicates()
    database_sets = catalog.groupby("database")["tf_symbol"].apply(set).to_dict()
    shared = sorted(database_sets.get("JASPAR", set()) & database_sets.get("HOCOMOCO", set()))
    catalog = catalog.loc[catalog["tf_symbol"].isin(shared)].copy()
    return catalog, shared


def load_presence_matrices(
    hit_paths: list[Path],
    region_ids: pd.Series,
    catalog: pd.DataFrame,
    tf_symbols: list[str],
    chunksize: int,
) -> dict[str, sparse.csr_matrix]:
    region_map = {value: idx for idx, value in enumerate(region_ids.astype(str))}
    tf_map = {value: idx for idx, value in enumerate(tf_symbols)}
    map_table = catalog[["motif_db", "motif_id", "database", "tf_symbol"]].drop_duplicates()
    rows: dict[str, list[np.ndarray]] = {"JASPAR": [], "HOCOMOCO": []}
    cols: dict[str, list[np.ndarray]] = {"JASPAR": [], "HOCOMOCO": []}
    for path in hit_paths:
        for chunk in pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            usecols=["region_id", "motif_db", "motif_id"],
            chunksize=chunksize,
        ):
            chunk["region_index"] = chunk["region_id"].astype(str).map(region_map)
            chunk = chunk.loc[chunk["region_index"].notna()].copy()
            if chunk.empty:
                continue
            expanded = chunk.merge(map_table, on=["motif_db", "motif_id"], how="inner")
            if expanded.empty:
                continue
            expanded["tf_index"] = expanded["tf_symbol"].map(tf_map)
            for database in ("JASPAR", "HOCOMOCO"):
                keep = expanded["database"].eq(database)
                if keep.any():
                    rows[database].append(expanded.loc[keep, "region_index"].to_numpy(dtype=np.int32))
                    cols[database].append(expanded.loc[keep, "tf_index"].to_numpy(dtype=np.int32))

    matrices: dict[str, sparse.csr_matrix] = {}
    shape = (len(region_ids), len(tf_symbols))
    for database in ("JASPAR", "HOCOMOCO"):
        if rows[database]:
            row = np.concatenate(rows[database])
            col = np.concatenate(cols[database])
            matrix = sparse.coo_matrix((np.ones(row.size, dtype=np.uint8), (row, col)), shape=shape).tocsr()
            matrix.sum_duplicates()
            matrix.data[:] = 1
        else:
            matrix = sparse.csr_matrix(shape, dtype=np.uint8)
        matrices[database] = matrix
    matrices["consensus"] = matrices["JASPAR"].multiply(matrices["HOCOMOCO"]).astype(np.uint8)
    return matrices


def classify_peak(distance: object) -> str:
    value = abs(float(distance))
    if value <= 2_000:
        return "tss_proximal_peak"
    if value <= 50_000:
        return "proximal_peak"
    return "distal_peak"


def build_peak_table(
    target_regions_path: Path,
    background_regions_path: Path,
    peak_links_path: Path,
    sequence_features: dict[tuple[str, int, int], tuple[float, float, int]],
) -> pd.DataFrame:
    target = pd.read_csv(target_regions_path, sep="\t", compression="infer")
    target = target.loc[
        target["region_class"].isin(PEAK_CLASSES) & target["chrom"].isin(ALLOWED_CHROMS)
    ].copy()
    target = target.sort_values(["major_label", "peak_id", "region_id"]).drop_duplicates(
        ["major_label", "peak_id"], keep="first"
    )
    target = target.rename(
        columns={
            "link_distance_to_tss": "distance_to_tss",
            "link_shared_donors": "shared_donors",
            "link_abs_pearson_r": "abs_pearson_r",
            "link_atac_median_logcpm": "atac_median_logcpm",
        }
    )
    target["target"] = 1
    target["gene_symbol"] = target["target_gene"].astype(str)

    background = pd.read_csv(background_regions_path, sep="\t", compression="infer")
    background = background.drop(
        columns=[
            "peak_width",
            "distance_to_tss",
            "shared_donors",
            "abs_pearson_r",
            "atac_median_logcpm",
            "closest_gene_symbol",
        ],
        errors="ignore",
    )
    links = pd.read_csv(
        peak_links_path,
        sep="\t",
        compression="infer",
        usecols=[
            "major_label",
            "peak_id",
            "peak_width",
            "distance_to_tss",
            "shared_donors",
            "abs_pearson_r",
            "atac_median_logcpm",
            "closest_gene_symbol",
        ],
    ).drop_duplicates(["major_label", "peak_id"])
    background = background.merge(links, on=["major_label", "peak_id"], how="left", validate="one_to_one")
    background["target"] = 0
    background["gene_symbol"] = background["closest_gene_symbol"]
    target_major = set(target["major_label"].astype(str))
    background = background.loc[
        background["region_class"].isin(PEAK_CLASSES)
        & background["major_label"].astype(str).isin(target_major)
        & background["chrom"].isin(ALLOWED_CHROMS)
    ].copy()

    columns = [
        "region_id",
        "chrom",
        "start",
        "end",
        "major_label",
        "region_class",
        "peak_id",
        "gene_symbol",
        "target",
        "peak_width",
        "distance_to_tss",
        "shared_donors",
        "abs_pearson_r",
        "atac_median_logcpm",
    ]
    combined = pd.concat([target[columns], background[columns]], ignore_index=True, sort=False)
    eligible_strata = set(
        combined.loc[combined["target"].eq(1), ["major_label", "region_class"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    combined = combined.loc[
        [
            (str(major_label), str(region_class)) in eligible_strata
            for major_label, region_class in zip(combined["major_label"], combined["region_class"])
        ]
    ].reset_index(drop=True)
    combined["expected_region_class"] = combined["distance_to_tss"].map(classify_peak)
    mismatch = ~combined["region_class"].eq(combined["expected_region_class"])
    if mismatch.any():
        raise ValueError(f"Found {int(mismatch.sum())} peak rows with inconsistent distance bins")
    combined = add_sequence_features(combined, sequence_features)
    combined["log_peak_width"] = np.log1p(pd.to_numeric(combined["peak_width"], errors="coerce"))
    combined["log_abs_distance"] = np.log1p(
        pd.to_numeric(combined["distance_to_tss"], errors="coerce").abs()
    )
    combined["analysis_index"] = np.arange(combined.shape[0], dtype=int)
    return combined


def build_expression_lookup(path: Path) -> pd.Series:
    columns = pd.read_csv(path, sep="\t", compression="infer", nrows=0).columns
    required = {"gene_symbol", "mean_log_cpm_overall"}
    missing = required.difference(columns)
    if missing:
        raise ValueError(f"Reference-expression table is missing columns: {sorted(missing)}")
    usecols = ["gene_symbol", "mean_log_cpm_overall"]
    if "grouping" in columns:
        usecols.append("grouping")
    dge = pd.read_csv(path, sep="\t", compression="infer", usecols=usecols)
    if "grouping" in dge.columns:
        dge = dge.loc[dge["grouping"].astype(str).eq("major")].copy()
    dge["gene_key"] = dge["gene_symbol"].astype(str).str.upper()
    return dge.groupby("gene_key")["mean_log_cpm_overall"].median()


def build_promoter_table(
    target_regions_path: Path,
    background_regions_path: Path,
    expression: pd.Series,
    sequence_features: dict[tuple[str, int, int], tuple[float, float, int]],
) -> pd.DataFrame:
    target = pd.read_csv(target_regions_path, sep="\t", compression="infer")
    target = target.loc[target["region_class"].eq("promoter") & target["chrom"].isin(ALLOWED_CHROMS)].copy()
    target = target.sort_values(["chrom", "start", "end", "region_id"]).drop_duplicates(
        ["chrom", "start", "end"], keep="first"
    )
    target["gene_symbol"] = target["target_gene"].astype(str)
    target["target"] = 1
    target["region_class"] = "promoter"

    background = pd.read_csv(background_regions_path, sep="\t", compression="infer")
    background = background.loc[background["chrom"].isin(ALLOWED_CHROMS)].copy()
    background = background.rename(columns={"region_class": "background_region_class"})
    background["region_class"] = "promoter"
    background["target"] = 0

    columns = ["region_id", "chrom", "start", "end", "gene_symbol", "region_class", "target"]
    combined = pd.concat([target[columns], background[columns]], ignore_index=True, sort=False)
    combined = add_sequence_features(combined, sequence_features)
    combined["gene_key"] = combined["gene_symbol"].astype(str).str.upper()
    combined["median_rna_logcpm"] = combined["gene_key"].map(expression)
    missing_target_expression = combined["target"].eq(1) & combined["median_rna_logcpm"].isna()
    if missing_target_expression.any():
        missing = ", ".join(sorted(combined.loc[missing_target_expression, "gene_symbol"].astype(str)))
        raise ValueError(f"Target promoters lacked RNA abundance estimates: {missing}")
    combined = combined.loc[combined["target"].eq(1) | combined["median_rna_logcpm"].notna()].reset_index(drop=True)
    combined["analysis_index"] = np.arange(combined.shape[0], dtype=int)
    return combined


def standardized_design(table: pd.DataFrame, covariates: list[str], categorical: list[str]) -> np.ndarray:
    columns: list[np.ndarray] = [np.ones(table.shape[0], dtype=float)]
    for name in categorical:
        dummies = pd.get_dummies(table[name].astype(str), prefix=name, drop_first=True, dtype=float)
        for column in dummies.columns:
            columns.append(dummies[column].to_numpy(dtype=float))
    for name in covariates:
        values = pd.to_numeric(table[name], errors="coerce")
        median = float(values.median()) if values.notna().any() else 0.0
        values = values.fillna(median).to_numpy(dtype=float)
        scale = float(np.std(values, ddof=0))
        if scale < 1e-8:
            continue
        candidate = (values - float(np.mean(values))) / scale
        current = np.column_stack(columns)
        if np.linalg.matrix_rank(np.column_stack([current, candidate])) > np.linalg.matrix_rank(current):
            columns.append(candidate)
    return np.column_stack(columns)


def odds_ratio_haldane(a: int, b: int, c: int, d: int) -> float:
    return float(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))


def presence_counts(presence: sparse.csr_matrix, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    target_rows = np.flatnonzero(target == 1)
    background_rows = np.flatnonzero(target == 0)
    target_count = np.asarray(presence[target_rows, :].sum(axis=0)).ravel().astype(int)
    background_count = np.asarray(presence[background_rows, :].sum(axis=0)).ravel().astype(int)
    return target_count, background_count


def fit_adjusted_models(
    table: pd.DataFrame,
    matrices: dict[str, sparse.csr_matrix],
    tf_symbols: list[str],
    region_class: str,
    covariates: list[str],
    categorical: list[str],
    cluster_column: str | None,
    min_target_hits: int,
    min_background_hits: int,
) -> pd.DataFrame:
    index = table.index.to_numpy(dtype=int)
    y = table["target"].to_numpy(dtype=float)
    base = standardized_design(table, covariates, categorical)
    clusters = table[cluster_column].astype(str).to_numpy() if cluster_column else None
    consensus = matrices["consensus"][index, :].tocsr()
    jaspar = matrices["JASPAR"][index, :].tocsr()
    hocomoco = matrices["HOCOMOCO"][index, :].tocsr()
    consensus_target, consensus_bg = presence_counts(consensus, y)
    jaspar_target, jaspar_bg = presence_counts(jaspar, y)
    hocomoco_target, hocomoco_bg = presence_counts(hocomoco, y)
    n_target = int(y.sum())
    n_background = int((1 - y).sum())
    rows: list[dict[str, object]] = []

    for motif_index, tf_symbol in enumerate(tf_symbols):
        a = int(consensus_target[motif_index])
        b = n_target - a
        c = int(consensus_bg[motif_index])
        d = n_background - c
        raw_or = odds_ratio_haldane(a, b, c, d)
        raw_p = float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue)
        result: dict[str, object] = {
            "region_class": region_class,
            "region_label": CLASS_LABELS[region_class],
            "tf_symbol": tf_symbol,
            "target_regions": n_target,
            "background_regions": n_background,
            "target_consensus_hits": a,
            "background_consensus_hits": c,
            "target_jaspar_hits": int(jaspar_target[motif_index]),
            "background_jaspar_hits": int(jaspar_bg[motif_index]),
            "target_hocomoco_hits": int(hocomoco_target[motif_index]),
            "background_hocomoco_hits": int(hocomoco_bg[motif_index]),
            "raw_consensus_odds_ratio": raw_or,
            "raw_consensus_fisher_p": raw_p,
            "adjusted_log_odds": np.nan,
            "adjusted_odds_ratio": np.nan,
            "adjusted_ci_low": np.nan,
            "adjusted_ci_high": np.nan,
            "adjusted_p": np.nan,
            "covariance_estimator": pd.NA,
            "model_status": "insufficient motif-positive regions",
        }
        if a < min_target_hits or c < min_background_hits or b < 1 or d < 5:
            result["model_status"] = "insufficient motif-presence variation"
            rows.append(result)
            continue
        motif = consensus[:, motif_index].toarray().ravel().astype(float)
        design = np.column_stack([base, motif])
        if np.linalg.matrix_rank(design) <= np.linalg.matrix_rank(base):
            result["model_status"] = "motif indicator collinear with adjustment covariates"
            rows.append(result)
            continue
        try:
            if clusters is None:
                model = sm.GLM(y, design, family=sm.families.Binomial())
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = model.fit(maxiter=100, disp=0, cov_type="HC3")
                covariance_estimator = "HC3"
            else:
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        model = sm.GEE(
                            y,
                            design,
                            groups=clusters,
                            family=sm.families.Binomial(),
                            cov_struct=sm.cov_struct.Independence(),
                        )
                        fit = model.fit(maxiter=100)
                    probe = [float(fit.params[-1]), float(fit.bse[-1]), float(fit.pvalues[-1])]
                    if not all(np.isfinite(probe)):
                        raise FloatingPointError("non-finite GEE estimate")
                    covariance_estimator = "binomial GEE, genomic-peak clusters"
                except Exception:
                    model = sm.GLM(y, design, family=sm.families.Binomial())
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fit = model.fit(
                            maxiter=100,
                            disp=0,
                            cov_type="HC3",
                        )
                    covariance_estimator = "HC3 fallback after failed peak-cluster GEE"
            coefficient = float(fit.params[-1])
            standard_error = float(fit.bse[-1])
            p_value = float(fit.pvalues[-1])
            if not all(np.isfinite([coefficient, standard_error, p_value])):
                raise FloatingPointError("non-finite final model estimate")
            result.update(
                {
                    "adjusted_log_odds": coefficient,
                    "adjusted_odds_ratio": float(np.exp(np.clip(coefficient, -30, 30))),
                    "adjusted_ci_low": float(np.exp(np.clip(coefficient - 1.96 * standard_error, -30, 30))),
                    "adjusted_ci_high": float(np.exp(np.clip(coefficient + 1.96 * standard_error, -30, 30))),
                    "adjusted_p": p_value,
                    "covariance_estimator": covariance_estimator,
                    "model_status": "ok",
                }
            )
        except Exception as exc:  # pragma: no cover - data-dependent numerical fallback
            result["model_status"] = f"failed: {type(exc).__name__}"
        rows.append(result)
    return pd.DataFrame(rows)


def match_control_indices(
    table: pd.DataFrame,
    covariates: list[str],
    categorical: list[str],
    controls_per_target: int,
    repeats: int,
    seed: int,
) -> tuple[list[np.ndarray], pd.DataFrame]:
    rng = np.random.default_rng(seed)
    selected_by_repeat: list[list[int]] = [[] for _ in range(repeats)]
    diagnostics: list[dict[str, object]] = []
    group_values = table[categorical].astype(str).agg("|".join, axis=1) if categorical else pd.Series("all", index=table.index)

    for group_name, group_index in group_values.groupby(group_values).groups.items():
        group = table.loc[group_index]
        target = group.loc[group["target"].eq(1)]
        background = group.loc[group["target"].eq(0)]
        if target.empty or background.empty:
            continue
        values = group[covariates].apply(pd.to_numeric, errors="coerce")
        values = values.fillna(values.median()).fillna(0.0)
        mean = values.mean(axis=0)
        scale = values.std(axis=0, ddof=0).replace(0, 1.0)
        standardized = (values - mean) / scale
        target_x = standardized.loc[target.index].to_numpy(dtype=float)
        background_x = standardized.loc[background.index].to_numpy(dtype=float)
        background_indices = background.index.to_numpy(dtype=int)
        effective_controls = min(controls_per_target, background.shape[0] // target.shape[0])
        if effective_controls < 1:
            continue
        target_copies = np.repeat(np.arange(target.shape[0], dtype=int), effective_controls)
        base_cost = np.square(target_x[target_copies, None, :] - background_x[None, :, :]).sum(axis=2)

        for repeat in range(repeats):
            jitter = rng.gumbel(loc=0.0, scale=0.01, size=base_cost.shape)
            _, selected_columns = linear_sum_assignment(base_cost + jitter)
            selected_by_repeat[repeat].extend(background_indices[selected_columns].astype(int).tolist())

        first_selected = np.asarray(selected_by_repeat[0], dtype=int)
        first_selected = first_selected[np.isin(first_selected, background.index.to_numpy(dtype=int))]
        for covariate in covariates:
            target_values = pd.to_numeric(target[covariate], errors="coerce")
            background_values = pd.to_numeric(background[covariate], errors="coerce")
            matched_values = pd.to_numeric(table.loc[first_selected, covariate], errors="coerce")
            pooled_sd = float(np.sqrt((target_values.var(ddof=1) + background_values.var(ddof=1)) / 2))
            matched_sd = float(np.sqrt((target_values.var(ddof=1) + matched_values.var(ddof=1)) / 2))
            diagnostics.append(
                {
                    "group": group_name,
                    "covariate": covariate,
                    "target_regions": target.shape[0],
                    "background_regions": background.shape[0],
                    "matched_background_regions": matched_values.shape[0],
                    "smd_before": (float(target_values.mean()) - float(background_values.mean())) / max(pooled_sd, 1e-8),
                    "smd_after": (float(target_values.mean()) - float(matched_values.mean())) / max(matched_sd, 1e-8),
                }
            )
    return [np.asarray(values, dtype=int) for values in selected_by_repeat], pd.DataFrame(diagnostics)


def run_matching_sensitivity(
    table: pd.DataFrame,
    consensus: sparse.csr_matrix,
    tf_symbols: list[str],
    region_class: str,
    covariates: list[str],
    categorical: list[str],
    controls_per_target: int,
    repeats: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected, diagnostics = match_control_indices(
        table,
        covariates,
        categorical,
        controls_per_target,
        repeats,
        seed,
    )
    target_index = table.index[table["target"].eq(1)].to_numpy(dtype=int)
    target_count = np.asarray(consensus[target_index, :].sum(axis=0)).ravel().astype(int)
    n_target = len(target_index)
    rows: list[pd.DataFrame] = []
    for repeat, background_index in enumerate(selected):
        background_count = np.asarray(consensus[background_index, :].sum(axis=0)).ravel().astype(int)
        n_background = len(background_index)
        records: list[dict[str, object]] = []
        for motif_index, tf_symbol in enumerate(tf_symbols):
            a = int(target_count[motif_index])
            b = n_target - a
            c = int(background_count[motif_index])
            d = n_background - c
            records.append(
                {
                    "region_class": region_class,
                    "repeat": repeat + 1,
                    "tf_symbol": tf_symbol,
                    "target_hits": a,
                    "background_hits": c,
                    "matched_background_regions": n_background,
                    "matched_log2_odds_ratio": float(np.log2(odds_ratio_haldane(a, b, c, d))),
                    "matched_fisher_p": float(fisher_exact([[a, b], [c, d]], alternative="greater").pvalue),
                }
            )
        rows.append(pd.DataFrame(records))
    replicates = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if replicates.empty:
        return replicates, replicates, diagnostics
    summary = (
        replicates.groupby(["region_class", "tf_symbol"], as_index=False)
        .agg(
            match_repeats=("repeat", "nunique"),
            median_matched_log2_or=("matched_log2_odds_ratio", "median"),
            q05_matched_log2_or=("matched_log2_odds_ratio", lambda x: x.quantile(0.05)),
            q95_matched_log2_or=("matched_log2_odds_ratio", lambda x: x.quantile(0.95)),
            fraction_repeats_enriched=("matched_log2_odds_ratio", lambda x: float((x > 0).mean())),
            fraction_repeats_nominal_p_lt_0_05=("matched_fisher_p", lambda x: float((x < 0.05).mean())),
            median_matched_p=("matched_fisher_p", "median"),
        )
    )
    return replicates, summary, diagnostics


def plot_results(results: pd.DataFrame, matching: pd.DataFrame, output_dir: Path) -> None:
    merged = results.copy()
    eligible = merged.loc[
        merged["model_status"].eq("ok")
        & merged["adjusted_odds_ratio"].gt(1)
        & merged["fraction_repeats_enriched"].ge(0.8)
    ].copy()
    eligible["plot_rank"] = -np.log10(eligible["analysis_fdr"].clip(lower=1e-300))
    selected = (
        eligible.sort_values(["region_class", "plot_rank"], ascending=[True, False])
        .groupby("region_class", group_keys=False)
        .head(12)
    )
    selected = pd.concat(
        [selected, merged.loc[merged["tf_symbol"].isin(STAT_IRF)]], ignore_index=True
    ).drop_duplicates(["region_class", "tf_symbol"])
    if selected.empty:
        return
    order = (
        selected.groupby("tf_symbol")["plot_rank"].max().sort_values(ascending=False).index[:30].tolist()
    )
    plot = selected.loc[selected["tf_symbol"].isin(order)].copy()
    plot["region_label"] = plot["region_class"].map(CLASS_LABELS)
    matrix = plot.pivot_table(
        index="tf_symbol", columns="region_label", values="adjusted_log_odds", aggfunc="first"
    ).reindex(index=order, columns=[CLASS_LABELS["promoter"]] + [CLASS_LABELS[x] for x in PEAK_CLASSES])
    annot_source = plot.assign(
        marker=np.where(
            plot["analysis_fdr"].le(0.05) & plot["fraction_repeats_enriched"].ge(0.8),
            "*",
            "",
        )
    )
    annot = annot_source.pivot_table(
        index="tf_symbol", columns="region_label", values="marker", aggfunc="first", fill_value=""
    ).reindex(index=matrix.index, columns=matrix.columns, fill_value="")

    sns.set_theme(style="white")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, max(6.5, 0.27 * len(matrix) + 2.2)), gridspec_kw={"width_ratios": [1.0, 0.72]})
    vmax = max(1.0, float(np.nanquantile(np.abs(matrix.to_numpy(dtype=float)), 0.95)))
    sns.heatmap(
        matrix,
        ax=axes[0],
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        annot=annot,
        fmt="",
        cbar_kws={"label": "Adjusted motif log odds"},
    )
    axes[0].set_title("A. Covariate-adjusted JASPAR/HOCOMOCO consensus motif enrichment", loc="left", weight="bold")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("TF motif")
    axes[0].tick_params(axis="x", rotation=35)

    stability = plot.pivot_table(
        index="tf_symbol", columns="region_label", values="fraction_repeats_enriched", aggfunc="first"
    ).reindex(index=matrix.index, columns=matrix.columns)
    sns.heatmap(
        stability,
        ax=axes[1],
        cmap="Blues",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Fraction of 100 matches with OR > 1"},
    )
    axes[1].set_title("B. Repeated matching stability", loc="left", weight="bold")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].set_yticklabels([])
    axes[1].tick_params(axis="x", rotation=35)
    fig.text(
        0.01,
        0.01,
        "* analysis-wide FDR <= 0.05 and enrichment in >=80% of matched-control repeats. FDR is controlled across the three peak bins and separately across promoters.",
        fontsize=8,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output_dir / "robust_motif_enrichment_summary.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "robust_motif_enrichment_summary.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    output_dir: Path,
    peak_table: pd.DataFrame,
    promoter_table: pd.DataFrame,
    tf_symbols: list[str],
    results: pd.DataFrame,
    matching: pd.DataFrame,
    diagnostics: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    tested = results.loc[results["model_status"].eq("ok")].copy()
    supported = tested.loc[
        tested["analysis_fdr"].le(0.05)
        & tested["adjusted_odds_ratio"].gt(1)
        & tested["fraction_repeats_enriched"].ge(0.8)
    ].copy()
    peak_supported = supported.loc[supported["region_class"].isin(PEAK_CLASSES)]
    promoter_supported = supported.loc[supported["region_class"].eq("promoter")]
    stat_irf = supported.loc[supported["tf_symbol"].isin(STAT_IRF)].sort_values(
        ["analysis_fdr", "adjusted_odds_ratio"], ascending=[True, False]
    )
    max_smd = float(diagnostics["smd_after"].abs().max()) if not diagnostics.empty else np.nan
    control_label = "control" if args.controls_per_target == 1 else "controls"
    lines = [
        "# Robust motif enrichment validation",
        "",
        "## Design",
        "",
        f"- Linked-peak foreground: {int(peak_table['target'].sum())} target peak-by-cell-type rows.",
        f"- Linked-peak background: {int((1 - peak_table['target']).sum())} non-target linked peak-by-cell-type rows.",
        f"- Mutually exclusive linked-peak bins: <=2 kb, 2-50 kb, and >50 kb by absolute distance to the linked-gene TSS.",
        f"- Fixed-promoter foreground: {int(promoter_table['target'].sum())} unique target-gene TSS +/-2 kb sequences.",
        f"- Fixed-promoter background: {int((1 - promoter_table['target']).sum())} non-target reference promoters.",
        f"- TF symbols represented in both motif databases and screened: {len(tf_symbols)}.",
        "- Primary motif indicator: the same TF symbol had at least one motif hit from both JASPAR and HOCOMOCO in the region.",
        "- Motif scanning threshold inherited from the completed motifmatchr scan: p < 5e-5 per PWM/region scan.",
        "- Peak models adjusted for GC fraction, CpG density, peak width, baseline accessibility, shared-donor count, absolute peak-gene correlation, TSS distance, and major cell type.",
        "- Peak models used binomial GEE with genomic peak as the clustering unit because a peak can occur in more than one cell type.",
        "- Promoter models adjusted for GC fraction, CpG density, and median RNA abundance; promoter results are secondary and separate from observed ATAC-peak tests.",
        f"- Matched-control sensitivity analysis used {args.controls_per_target} optimally selected {control_label} per target over {args.match_repeats} repeats within cell type and distance bin.",
        "- FDR was controlled across all TF-by-bin tests in the three-bin linked-peak analysis; fixed-promoter tests were corrected separately because they use a distinct sequence universe and background.",
        "",
        "## Results",
        "",
        f"- Successfully fitted tests: {tested.shape[0]}.",
        f"- Supported linked-peak enrichments (peak-analysis FDR <= 0.05, adjusted OR > 1, OR > 1 in >=80% repeats): {peak_supported.shape[0]}.",
        f"- Supported fixed-promoter enrichments under the same rule: {promoter_supported.shape[0]}.",
        f"- Largest cell-type-specific absolute post-match standardized mean difference in the recorded first repeat: {max_smd:.3f}; matching is therefore treated as a direction-stability sensitivity analysis rather than a substitute for covariate adjustment.",
        "",
        "## STAT/IRF results passing the full rule",
        "",
    ]
    if stat_irf.empty:
        lines.append("No STAT/IRF motif passed the global-FDR and repeated-matching rule.")
    else:
        for row in stat_irf.itertuples(index=False):
            lines.append(
                f"- {row.tf_symbol}, {CLASS_LABELS[row.region_class]}: adjusted OR {row.adjusted_odds_ratio:.2f} "
                f"(95% CI {row.adjusted_ci_low:.2f}-{row.adjusted_ci_high:.2f}), analysis-wide FDR {row.analysis_fdr:.3g}, "
                f"enriched in {row.fraction_repeats_enriched:.0%} of matches."
            )
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "These are sequence-enrichment tests. A supported TF motif indicates compatibility with binding by that TF or a TF with a similar DNA-binding motif; it does not establish TF occupancy. Donor-aware motif-accessibility tests are reported separately as an orthogonal assay-level validation.",
            "",
        ]
    )
    (output_dir / "robust_motif_enrichment_summary.md").write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir = Path(args.fasta_dir)

    print("Loading sequence composition features...", flush=True)
    sequence_features = load_sequence_features(fasta_dir)
    print("Building mutually exclusive linked-peak and separate promoter tables...", flush=True)
    peak_table = build_peak_table(
        Path(args.target_regions),
        Path(args.background_peak_regions),
        Path(args.peak_links),
        sequence_features,
    )
    expression = build_expression_lookup(Path(args.reference_expression))
    promoter_table = build_promoter_table(
        Path(args.target_regions),
        Path(args.background_promoters),
        expression,
        sequence_features,
    )
    peak_table.to_csv(output_dir / "peak_analysis_regions.tsv.gz", sep="\t", index=False, compression="gzip")
    promoter_table.to_csv(output_dir / "promoter_analysis_regions.tsv.gz", sep="\t", index=False, compression="gzip")

    print("Building cross-database TF catalog...", flush=True)
    catalog, tf_symbols = build_motif_catalog(Path(args.target_hits))
    catalog.to_csv(output_dir / "shared_tf_motif_catalog.tsv", sep="\t", index=False)
    pd.Series(tf_symbols, name="tf_symbol").to_csv(output_dir / "tested_tf_symbols.tsv", sep="\t", index=False)

    print(f"Loading motif presence for {len(tf_symbols)} TF symbols across linked peaks...", flush=True)
    peak_matrices = load_presence_matrices(
        [Path(args.target_hits), Path(args.background_peak_hits)],
        peak_table["region_id"],
        catalog,
        tf_symbols,
        args.chunksize,
    )
    print("Loading motif presence across fixed promoters...", flush=True)
    promoter_matrices = load_presence_matrices(
        [Path(args.target_hits), Path(args.background_promoter_hits)],
        promoter_table["region_id"],
        catalog,
        tf_symbols,
        args.chunksize,
    )
    for database, matrix in peak_matrices.items():
        sparse.save_npz(output_dir / f"peak_motif_presence_{database.lower()}.npz", matrix)
    for database, matrix in promoter_matrices.items():
        sparse.save_npz(output_dir / f"promoter_motif_presence_{database.lower()}.npz", matrix)

    result_parts: list[pd.DataFrame] = []
    match_replicates: list[pd.DataFrame] = []
    match_summaries: list[pd.DataFrame] = []
    diagnostic_parts: list[pd.DataFrame] = []
    peak_covariates = [
        "gc_fraction",
        "cpg_density",
        "log_peak_width",
        "atac_median_logcpm",
        "shared_donors",
        "abs_pearson_r",
        "log_abs_distance",
    ]
    for class_index, region_class in enumerate(PEAK_CLASSES):
        print(f"Fitting adjusted models and repeated matches: {region_class}", flush=True)
        region_table = peak_table.loc[peak_table["region_class"].eq(region_class)].copy()
        result_parts.append(
            fit_adjusted_models(
                region_table,
                peak_matrices,
                tf_symbols,
                region_class,
                peak_covariates,
                ["major_label"],
                "peak_id",
                args.min_target_hits,
                args.min_background_hits,
            )
        )
        replicates, summary, diagnostics = run_matching_sensitivity(
            region_table,
            peak_matrices["consensus"],
            tf_symbols,
            region_class,
            peak_covariates,
            ["major_label"],
            args.controls_per_target,
            args.match_repeats,
            args.seed + class_index,
        )
        match_replicates.append(replicates)
        match_summaries.append(summary)
        diagnostics["region_class"] = region_class
        diagnostic_parts.append(diagnostics)

    promoter_covariates = ["gc_fraction", "cpg_density", "median_rna_logcpm"]
    print("Fitting separate fixed-promoter models and repeated matches...", flush=True)
    result_parts.append(
        fit_adjusted_models(
            promoter_table,
            promoter_matrices,
            tf_symbols,
            "promoter",
            promoter_covariates,
            [],
            None,
            args.min_target_hits,
            args.min_background_hits,
        )
    )
    replicates, summary, diagnostics = run_matching_sensitivity(
        promoter_table,
        promoter_matrices["consensus"],
        tf_symbols,
        "promoter",
        promoter_covariates,
        [],
        args.controls_per_target,
        args.match_repeats,
        args.seed + len(PEAK_CLASSES),
    )
    match_replicates.append(replicates)
    match_summaries.append(summary)
    diagnostics["region_class"] = "promoter"
    diagnostic_parts.append(diagnostics)

    results = pd.concat(result_parts, ignore_index=True)
    results["global_all_fdr"] = bh_fdr(results["adjusted_p"])
    peak_mask = results["region_class"].isin(PEAK_CLASSES)
    results["global_peak_fdr"] = np.nan
    results.loc[peak_mask, "global_peak_fdr"] = bh_fdr(results.loc[peak_mask, "adjusted_p"])
    promoter_mask = results["region_class"].eq("promoter")
    results["promoter_fdr"] = np.nan
    results.loc[promoter_mask, "promoter_fdr"] = bh_fdr(results.loc[promoter_mask, "adjusted_p"])
    results["analysis_fdr"] = np.where(
        peak_mask,
        results["global_peak_fdr"],
        results["promoter_fdr"],
    )

    replicate_table = pd.concat(match_replicates, ignore_index=True)
    matching = pd.concat(match_summaries, ignore_index=True)
    diagnostics = pd.concat(diagnostic_parts, ignore_index=True)
    results = results.merge(matching, on=["region_class", "tf_symbol"], how="left")
    results.to_csv(output_dir / "robust_motif_enrichment_results.tsv.gz", sep="\t", index=False, compression="gzip")
    replicate_table.to_csv(output_dir / "matched_enrichment_replicates.tsv.gz", sep="\t", index=False, compression="gzip")
    matching.to_csv(output_dir / "matched_enrichment_stability.tsv", sep="\t", index=False)
    diagnostics.to_csv(output_dir / "matching_balance_diagnostics.tsv", sep="\t", index=False)
    plot_results(results, matching, output_dir)
    write_summary(output_dir, peak_table, promoter_table, tf_symbols, results, matching, diagnostics, args)
    print(f"Wrote robust motif validation outputs to {output_dir}", flush=True)


if __name__ == "__main__":
    main()
