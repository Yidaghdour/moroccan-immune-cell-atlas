#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


TEXT = "#222222"
GRID = "#dddddd"
ROW_ORDER = [
    ("Monocytes", "F"),
    ("Naive CD4 T cells", "F"),
    ("Cytotoxic T/NK cells", "F"),
    ("B cells", "F"),
    ("Memory T cells", "F"),
    ("Memory T cells", "M"),
]
ROW_LABELS = {
    ("Monocytes", "F"): "Monocytes | Female",
    ("Naive CD4 T cells", "F"): "CD4 Naive | Female",
    ("Cytotoxic T/NK cells", "F"): "Cytotoxic T/NK | Female",
    ("B cells", "F"): "B cells | Female",
    ("Memory T cells", "F"): "Memory T | Female",
    ("Memory T cells", "M"): "Memory T | Male",
}
IRF_MOTIFS = ["IRF7", "IRF2", "IRF9", "IRF5", "IRF4", "IRF8", "IRF6"]
OTHER_ROBUST_MOTIFS = ["STAT3", "MAFB", "ZBED2", "HIC1", "HAND1"]
CONTEXT_MOTIFS = [
    "FOS",
    "JUN",
    "NFKB1",
    "REL",
    "BACH2",
    "RUNX1",
    "CTCF",
    "FOXA2",
    "HIF1A",
]
CONTEXT_FAMILY = {
    "FOS": "AP-1",
    "JUN": "AP-1",
    "NFKB1": "NF-kB",
    "REL": "NF-kB",
    "BACH2": "BACH",
    "RUNX1": "RUNX",
    "CTCF": "CTCF",
    "FOXA2": "FOX",
    "HIF1A": "HIF",
}
TF_ORDER = IRF_MOTIFS + OTHER_ROBUST_MOTIFS + CONTEXT_MOTIFS
N_SUPPORT_CONFIGURATIONS = 64


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add selected non-IRF family context to the threshold-robust "
            "Figure 2D TFBS summary."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=(
            "results/tfbs_gene_overlap_integrated/robust_motif_validation_report/"
            "sex_stratified_tfbs_enrichment"
        ),
    )
    return parser.parse_args()


def require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {sorted(missing)}")


def dot_size(fdr: pd.Series | np.ndarray | float) -> np.ndarray:
    values = np.asarray(fdr, dtype=float)
    evidence = np.clip(-np.log10(np.clip(values, 1e-6, 1)), 0, 6)
    return 20 + 34 * evidence


def build_plot_table(results: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    require_columns(
        results,
        {
            "major_label",
            "sex",
            "tf_symbol",
            "model_status",
            "adjusted_log2_odds",
            "adjusted_odds_ratio",
            "maximum_support_grid_fdr",
            "support_configurations_significant",
            "support_configurations_tested",
            "significant_enrichment",
            "target_regions",
        },
        "TFBS enrichment results",
    )
    require_columns(
        coverage,
        {"major_label", "sex", "submitted_genes"},
        "TFBS target coverage",
    )
    available = set(results["tf_symbol"].astype(str))
    missing = set(TF_ORDER).difference(available)
    if missing:
        raise ValueError(f"Requested context motifs were not tested: {sorted(missing)}")
    if not results["support_configurations_tested"].eq(N_SUPPORT_CONFIGURATIONS).all():
        raise ValueError("The extended figure requires the 64-configuration result table")

    row_map = {f"{cell}|{sex}": index for index, (cell, sex) in enumerate(ROW_ORDER)}
    tf_map = {tf: index for index, tf in enumerate(TF_ORDER)}
    plot = results.copy()
    plot["row_key"] = plot["major_label"].astype(str) + "|" + plot["sex"].astype(str)
    plot = plot.loc[
        plot["row_key"].isin(row_map)
        & plot["tf_symbol"].isin(TF_ORDER)
        & plot["model_status"].eq("ok")
    ].copy()
    expected_rows = len(ROW_ORDER) * len(TF_ORDER)
    if plot.shape[0] != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} extended-figure estimates, observed {plot.shape[0]}"
        )
    plot["y"] = plot["row_key"].map(row_map)
    plot["x"] = plot["tf_symbol"].map(tf_map)
    plot["display_fdr"] = plot["maximum_support_grid_fdr"]
    plot["display_group"] = np.select(
        [
            plot["tf_symbol"].isin(IRF_MOTIFS),
            plot["tf_symbol"].isin(OTHER_ROBUST_MOTIFS),
        ],
        ["IRF-like", "Other threshold-robust"],
        default="Selected context",
    )
    plot["context_family"] = plot["tf_symbol"].map(CONTEXT_FAMILY)
    plot = plot.merge(
        coverage[["major_label", "sex", "submitted_genes"]],
        on=["major_label", "sex"],
        how="left",
        validate="many_to_one",
    )
    return plot


def draw_group_header(
    axis: plt.Axes,
    start: int,
    end: int,
    label: str,
    color: str,
) -> None:
    axis.text(
        (start + end) / 2,
        -1.02,
        label,
        ha="center",
        va="bottom",
        fontsize=9.5,
        weight="bold",
        color=color,
        clip_on=False,
    )


def plot_extended(plot: pd.DataFrame, output_dir: Path) -> None:
    fig, axis = plt.subplots(figsize=(20.2, 8.7))
    plt.subplots_adjust(left=0.15, right=0.87, top=0.79, bottom=0.27)

    for row_index in range(len(ROW_ORDER)):
        if row_index % 2:
            axis.axhspan(row_index - 0.5, row_index + 0.5, color="#f8f8f8", zorder=0)
    irf_end = len(IRF_MOTIFS) - 0.5
    robust_end = len(IRF_MOTIFS) + len(OTHER_ROBUST_MOTIFS) - 0.5
    axis.axvspan(-0.5, irf_end, color="#f3f6fb", zorder=0)
    axis.axvspan(robust_end, len(TF_ORDER) - 0.5, color="#f1f1f1", zorder=0)
    axis.axvline(irf_end, color="#9ca9ba", linewidth=1.0, zorder=1)
    axis.axvline(robust_end, color="#999999", linewidth=1.2, zorder=1)

    color_limit = max(
        2.0,
        float(
            np.ceil(np.nanquantile(np.abs(plot["adjusted_log2_odds"]), 0.98) * 2)
            / 2
        ),
    )
    normalization = TwoSlopeNorm(vmin=-color_limit, vcenter=0, vmax=color_limit)
    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        ["#355ea9", "#f7f7f5", "#cd142b"],
    )
    scatter = axis.scatter(
        plot["x"],
        plot["y"],
        s=dot_size(plot["display_fdr"]),
        c=plot["adjusted_log2_odds"],
        cmap=cmap,
        norm=normalization,
        edgecolors="#888888",
        linewidths=0.45,
        zorder=3,
    )
    robust = plot.loc[plot["significant_enrichment"].astype(bool)]
    axis.scatter(
        robust["x"],
        robust["y"],
        s=dot_size(robust["display_fdr"]) + 24,
        facecolors="none",
        edgecolors="#111111",
        linewidths=1.45,
        zorder=4,
    )

    axis.set_xticks(np.arange(len(TF_ORDER)))
    axis.set_xticklabels(
        TF_ORDER,
        fontsize=9.2,
        rotation=42,
        ha="right",
        rotation_mode="anchor",
    )
    axis.set_yticks(np.arange(len(ROW_ORDER)))
    axis.set_yticklabels([ROW_LABELS[value] for value in ROW_ORDER], fontsize=10.5)
    axis.set_xlim(-0.55, len(TF_ORDER) + 1.55)
    axis.set_ylim(len(ROW_ORDER) - 0.5, -0.5)
    axis.grid(color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", length=0, colors=TEXT)
    for spine in axis.spines.values():
        spine.set_color("#bdbdbd")
        spine.set_linewidth(0.8)

    draw_group_header(axis, 0, len(IRF_MOTIFS) - 1, "IRF-like motifs", "#355ea9")
    draw_group_header(
        axis,
        len(IRF_MOTIFS),
        len(IRF_MOTIFS) + len(OTHER_ROBUST_MOTIFS) - 1,
        "Other threshold-robust motifs",
        "#555555",
    )
    draw_group_header(
        axis,
        len(IRF_MOTIFS) + len(OTHER_ROBUST_MOTIFS),
        len(TF_ORDER) - 1,
        "Selected family context",
        "#555555",
    )

    axis.text(
        len(TF_ORDER) + 0.45,
        -0.82,
        "Target genes /\nregulatory regions",
        ha="center",
        va="bottom",
        fontsize=8.5,
        weight="bold",
        clip_on=False,
    )
    coverage = (
        plot[["major_label", "sex", "submitted_genes", "target_regions"]]
        .drop_duplicates(["major_label", "sex"])
        .set_index(["major_label", "sex"])
    )
    for row_index, key in enumerate(ROW_ORDER):
        values = coverage.loc[key]
        axis.text(
            len(TF_ORDER) + 0.45,
            row_index,
            f"{int(values.submitted_genes)} / {int(values.target_regions)}",
            ha="center",
            va="center",
            fontsize=8.2,
            color="#555555",
        )

    colorbar_axis = fig.add_axes([0.895, 0.40, 0.012, 0.28])
    colorbar = fig.colorbar(scatter, cax=colorbar_axis)
    colorbar.set_label(
        "Adjusted log2 odds ratio\ntarget vs background",
        fontsize=9,
    )
    colorbar.ax.tick_params(labelsize=8)
    colorbar.outline.set_linewidth(0.7)

    fdr_values = [0.05, 0.01, 0.001]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#d8d8d8",
            markeredgecolor="#777777",
            markeredgewidth=0.45,
            markersize=np.sqrt(float(dot_size(value))),
            label=f"{value:g}",
        )
        for value in fdr_values
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor="#111111",
            markeredgewidth=1.45,
            markersize=8,
            label="FDR <= 0.05 and OR > 1\nin all 64 configurations",
        )
    )
    fig.legend(
        handles=handles,
        title="Dot size: worst-case\nglobal FDR",
        loc="upper left",
        bbox_to_anchor=(0.884, 0.35),
        frameon=False,
        fontsize=8,
        title_fontsize=8.5,
        borderaxespad=0,
    )

    fig.suptitle(
        "TFBS enrichment across concordant RNA-ATAC gene sets",
        x=0.15,
        y=0.955,
        ha="left",
        fontsize=20,
        weight="bold",
        color=TEXT,
    )
    fig.text(
        0.15,
        0.908,
        (
            "Threshold-robust individual motifs plus selected representatives "
            "of AP-1/NF-kB, BACH, RUNX, CTCF, FOX, and HIF programs"
        ),
        ha="left",
        va="center",
        fontsize=10.3,
        color="#4d4d4d",
    )
    fig.text(
        0.055,
        0.515,
        "Figure 2D gene set",
        rotation=90,
        ha="center",
        va="center",
        fontsize=10.5,
    )
    fig.text(0.515, 0.135, "Individual TF motif", ha="center", va="center", fontsize=10.5)
    fig.text(
        0.15,
        0.040,
        (
            "Black rings identify adjusted OR > 1 and global BH FDR <= 0.05 in all 64 support configurations; "
            "dot size uses the largest FDR across configurations. Motifs in the grey context block were selected "
            "and are shown even when not significant; they provide context, not positive evidence. Counts at right are "
            "target genes / foreground regulatory regions. Closely related IRF motifs cannot identify the exact bound IRF protein."
        ),
        ha="left",
        va="bottom",
        fontsize=8.3,
        color="#4d4d4d",
        wrap=True,
    )

    for extension in ["pdf", "png"]:
        fig.savefig(
            output_dir
            / f"figure2d_gene_set_tfbs_enrichment_extended_families.{extension}",
            dpi=320,
            facecolor="white",
        )
    plt.close(fig)


def write_notes(plot: pd.DataFrame, path: Path) -> None:
    context = plot.loc[plot["tf_symbol"].isin(CONTEXT_MOTIFS)].copy()
    context_robust = context.loc[context["significant_enrichment"].astype(bool)]
    robust = plot.loc[plot["significant_enrichment"].astype(bool)]
    lines = [
        "# Extended Figure 2D TFBS family context",
        "",
        "## Intended use",
        "",
        "This supplementary version retains every motif selected by the threshold-robust primary summary and adds manually selected representatives of additional TF families for context.",
        "",
        "## Statistical display",
        "",
        "- All estimates come from the same covariate-adjusted models as the primary Figure 2D TFBS summary.",
        "- Dot size represents the largest global BH FDR across the 64 support configurations.",
        "- Black rings require adjusted odds ratio >1 and global FDR <=0.05 in all 64 configurations.",
        f"- Threshold-robust combinations displayed: {robust.shape[0]}.",
        "",
        "## Selected context block",
        "",
    ]
    lines.extend(f"- {tf}: {CONTEXT_FAMILY[tf]}." for tf in CONTEXT_MOTIFS)
    lines.extend(
        [
            "",
            f"Threshold-robust context combinations: {context_robust.shape[0]}.",
            "A context estimate without a black ring did not satisfy the all-configuration enrichment rule.",
            "",
            "## Scope",
            "",
            "The grey context block is a descriptive comparison set and is not evidence that those TF families are enriched. Closely related motifs can be recognized by the same sequence and therefore do not establish which TF protein was bound.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = Path(args.input_dir)
    results = pd.read_csv(
        output_dir / "sex_stratified_tfbs_enrichment_results.tsv.gz",
        sep="\t",
    )
    coverage = pd.read_csv(
        output_dir / "sex_stratified_tfbs_target_coverage.tsv",
        sep="\t",
    )
    plot = build_plot_table(results, coverage)
    plot.to_csv(
        output_dir / "figure2d_gene_set_tfbs_enrichment_extended_families_values.tsv",
        sep="\t",
        index=False,
    )
    write_notes(
        plot,
        output_dir / "figure2d_gene_set_tfbs_enrichment_extended_families_notes.md",
    )
    plot_extended(plot, output_dir)
    print(f"Wrote extended TFBS family-context figure to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
