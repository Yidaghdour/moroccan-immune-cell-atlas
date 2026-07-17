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
IRF_ORDER = ["IRF7", "IRF2", "IRF9", "IRF5", "IRF4", "IRF8", "IRF6"]
OTHER_ORDER = ["STAT3", "MAFB", "ZBED2", "HIC1", "HAND1"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compact TFBS summary aligned to the Figure 2D gene sets."
    )
    parser.add_argument(
        "--input-dir",
        default=(
            "results/tfbs_gene_overlap_integrated/robust_motif_validation_report/"
            "sex_stratified_tfbs_enrichment"
        ),
    )
    return parser.parse_args()


def dot_size(fdr: pd.Series | np.ndarray | float) -> np.ndarray:
    values = np.asarray(fdr, dtype=float)
    evidence = np.clip(-np.log10(np.clip(values, 1e-6, 1)), 0, 6)
    return 20 + 34 * evidence


def motif_order(significant: pd.DataFrame) -> list[str]:
    observed = set(significant["tf_symbol"].astype(str))
    ordered = [tf for tf in IRF_ORDER + OTHER_ORDER if tf in observed]
    remaining = (
        significant.loc[~significant["tf_symbol"].isin(ordered)]
        .groupby("tf_symbol")["global_fdr"]
        .min()
        .sort_values()
        .index.tolist()
    )
    return ordered + remaining


def build_plot_table(
    results: pd.DataFrame,
    coverage: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    row_keys = {f"{cell}|{sex}" for cell, sex in ROW_ORDER}
    results = results.copy()
    results["row_key"] = results["major_label"].astype(str) + "|" + results["sex"].astype(str)
    significant = results.loc[
        results["significant_enrichment"].astype(bool) & results["row_key"].isin(row_keys)
    ].copy()
    tf_order = motif_order(significant)
    if not tf_order:
        raise ValueError("No significant individual TF motifs were available for Figure 2D strata")
    plot = results.loc[
        results["row_key"].isin(row_keys)
        & results["tf_symbol"].isin(tf_order)
        & results["model_status"].eq("ok")
    ].copy()
    row_map = {f"{cell}|{sex}": index for index, (cell, sex) in enumerate(ROW_ORDER)}
    tf_map = {tf: index for index, tf in enumerate(tf_order)}
    plot["y"] = plot["row_key"].map(row_map)
    plot["x"] = plot["tf_symbol"].map(tf_map)
    plot = plot.merge(
        coverage[["major_label", "sex", "submitted_genes"]],
        on=["major_label", "sex"],
        how="left",
        validate="many_to_one",
    )
    return plot, tf_order


def plot_summary(
    plot: pd.DataFrame,
    tf_order: list[str],
    output_dir: Path,
) -> None:
    fig, axis = plt.subplots(figsize=(9.2, 11.4))
    plt.subplots_adjust(left=0.23, right=0.77, top=0.76, bottom=0.14)

    for column_index in range(len(ROW_ORDER)):
        if column_index % 2:
            axis.axvspan(
                column_index - 0.5,
                column_index + 0.5,
                color="#f7f7f7",
                zorder=0,
            )
    irf_count = sum(tf in IRF_ORDER for tf in tf_order)
    if irf_count:
        axis.axhspan(-0.5, irf_count - 0.5, color="#f3f6fb", alpha=0.72, zorder=0)
        axis.axhline(irf_count - 0.5, color="#9ca9ba", linewidth=1.0, zorder=1)

    color_limit = max(
        2.0,
        float(np.ceil(np.nanquantile(np.abs(plot["adjusted_log2_odds"]), 0.98) * 2) / 2),
    )
    normalization = TwoSlopeNorm(vmin=-color_limit, vcenter=0, vmax=color_limit)
    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        ["#355ea9", "#f7f7f5", "#cd142b"],
    )
    scatter = axis.scatter(
        plot["y"],
        plot["x"],
        s=dot_size(plot["global_fdr"]),
        c=plot["adjusted_log2_odds"],
        cmap=cmap,
        norm=normalization,
        edgecolors="#888888",
        linewidths=0.45,
        zorder=3,
    )
    supported = plot.loc[plot["significant_enrichment"].astype(bool)]
    axis.scatter(
        supported["y"],
        supported["x"],
        s=dot_size(supported["global_fdr"]) + 24,
        facecolors="none",
        edgecolors="#111111",
        linewidths=1.45,
        zorder=4,
    )

    coverage = (
        plot[["major_label", "sex", "submitted_genes", "target_regions"]]
        .drop_duplicates(["major_label", "sex"])
        .set_index(["major_label", "sex"])
    )
    column_names = {
        ("Monocytes", "F"): "Monocytes\nFemale",
        ("Naive CD4 T cells", "F"): "CD4 Naive\nFemale",
        ("Cytotoxic T/NK cells", "F"): "Cytotoxic T/NK\nFemale",
        ("B cells", "F"): "B cells\nFemale",
        ("Memory T cells", "F"): "Memory T\nFemale",
        ("Memory T cells", "M"): "Memory T\nMale",
    }
    column_labels = []
    for key in ROW_ORDER:
        values = coverage.loc[key]
        column_labels.append(
            f"{column_names[key]}\n{int(values.submitted_genes)} / {int(values.target_regions)}"
        )

    axis.set_xticks(np.arange(len(ROW_ORDER)))
    axis.set_xticklabels(column_labels, fontsize=8.1, linespacing=1.25)
    axis.xaxis.tick_top()
    axis.tick_params(axis="x", labeltop=True, labelbottom=False, pad=8)
    axis.set_yticks(np.arange(len(tf_order)))
    axis.set_yticklabels(tf_order, fontsize=10)
    axis.set_xlim(-0.5, len(ROW_ORDER) - 0.5)
    axis.set_ylim(len(tf_order) - 0.5, -0.5)
    axis.grid(color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", length=0, colors=TEXT)
    for spine in axis.spines.values():
        spine.set_color("#bdbdbd")
        spine.set_linewidth(0.8)

    irf_positions = [index for index, tf in enumerate(tf_order) if tf in IRF_ORDER]
    other_positions = [index for index, tf in enumerate(tf_order) if tf not in IRF_ORDER]
    if irf_positions:
        axis.text(
            -1.30,
            np.mean(irf_positions),
            "IRF-like motifs",
            ha="center",
            va="center",
            rotation=90,
            fontsize=9.5,
            weight="bold",
            color="#355ea9",
            clip_on=False,
        )
    if other_positions:
        axis.text(
            -1.30,
            np.mean(other_positions),
            "Additional motifs",
            ha="center",
            va="center",
            rotation=90,
            fontsize=9.5,
            weight="bold",
            color="#555555",
            clip_on=False,
        )

    colorbar_axis = fig.add_axes([0.82, 0.56, 0.022, 0.20])
    colorbar = fig.colorbar(scatter, cax=colorbar_axis)
    colorbar.set_label("Adjusted log2 odds ratio\ntarget vs background", fontsize=9)
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
            label="FDR <= 0.05 in all 64\nsupport settings; OR > 1",
        )
    )
    fig.legend(
        handles=handles,
        title="Dot size: worst-case\nglobal FDR",
        loc="upper left",
        bbox_to_anchor=(0.80, 0.49),
        frameon=False,
        fontsize=8,
        title_fontsize=8.5,
        borderaxespad=0,
    )

    fig.suptitle(
        "TFBS enrichment across concordant RNA-ATAC gene sets",
        x=0.08,
        y=0.968,
        ha="left",
        fontsize=17.5,
        weight="bold",
        color=TEXT,
    )
    fig.text(
        0.08,
        0.925,
        "Gene sets correspond to the Unstimulated cell-type and sex strata in Figure 2D;\nindividual TF motifs are shown separately.",
        ha="left",
        va="center",
        fontsize=9.7,
        color="#4d4d4d",
    )
    fig.text(
        0.50,
        0.835,
        "Figure 2D gene set\nTarget genes / regulatory regions",
        ha="center",
        va="center",
        fontsize=9.5,
        weight="bold",
        linespacing=1.35,
    )
    fig.text(0.035, 0.465, "Individual TF motif", rotation=90, ha="center", va="center", fontsize=10.5)
    fig.text(
        0.08,
        0.035,
        "Dots show covariate-adjusted motif overrepresentation across each gene set's promoters and linked ATAC peaks. Dot size uses the largest global BH FDR across 64 support configurations; black rings mark FDR <= 0.05 with odds ratio > 1 in every configuration. "
        "Counts beneath column labels are target genes / foreground regulatory regions. Closely related IRF motifs cannot identify the exact bound IRF protein. "
        "Female and Male columns use different selected genes and are not a formal sex contrast.",
        ha="left",
        va="bottom",
        fontsize=8.2,
        color="#4d4d4d",
        wrap=True,
    )

    for extension in ["pdf", "png"]:
        fig.savefig(
            output_dir / f"figure2d_gene_set_tfbs_enrichment.{extension}",
            dpi=320,
            facecolor="white",
        )
    plt.close(fig)


def write_caption(plot: pd.DataFrame, tf_order: list[str], path: Path) -> None:
    significant = plot.loc[plot["significant_enrichment"].astype(bool)].copy()
    lines = [
        "# Figure 2D gene-set TFBS enrichment",
        "",
        "## Suggested caption",
        "",
        "TFBS enrichment across cell-type- and sex-stratified concordant RNA-ATAC gene sets. For each Unstimulated gene set shown in Figure 2D, reference promoters (TSS +/-2 kb) and linked ATAC peaks were compared with non-target background regions using JASPAR 2024 and HOCOMOCO v14 consensus TF-symbol annotations. Binomial models adjusted for regulatory-region class, GC fraction, CpG density, and sequence length. A 64-configuration support analysis varied target-set and motif-state count requirements while retaining a fixed multiple-testing universe. Dot color represents the adjusted log2 odds ratio, dot area represents the largest global BH FDR across the configurations, and black outlines identify combinations with FDR <=0.05 and odds ratio >1 in every configuration. Closely related motifs, particularly IRF motifs, cannot distinguish the exact TF family member or establish TF occupancy.",
        "",
        "## Main descriptive result",
        "",
        "Female lymphoid gene sets show a recurrent IRF-like TFBS signature. Female B-cell targets additionally support MAFB and STAT3 motifs, whereas the Male Memory T-cell set supports IRF7 together with ZBED2, HIC1, and HAND1 motifs. No motif passed global FDR in the Female Monocyte set. These are gene-set-specific enrichment patterns, not formal Female-versus-Male effects.",
        "",
        "## Displayed motifs",
        "",
        f"{', '.join(tf_order)}.",
        "",
        "## Significant combinations",
        "",
    ]
    for row in significant.sort_values(["y", "global_fdr", "tf_symbol"]).itertuples(index=False):
        lines.append(
            f"- {ROW_LABELS[(row.major_label, row.sex)]}: {row.tf_symbol}; "
            f"adjusted OR={row.adjusted_odds_ratio:.2f}, worst-case global FDR={row.global_fdr:.3g}."
        )
    lines.append("")
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
    plot, tf_order = build_plot_table(results, coverage)
    plot.to_csv(
        output_dir / "figure2d_gene_set_tfbs_enrichment_values.tsv",
        sep="\t",
        index=False,
    )
    write_caption(
        plot,
        tf_order,
        output_dir / "figure2d_gene_set_tfbs_enrichment_caption.md",
    )
    plot_summary(plot, tf_order, output_dir)
    print(f"Wrote Figure 2D-aligned TFBS summary to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
