#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
import numpy as np
import pandas as pd


COLORS = {
    "TF": "#FF3B3B",
    "strict": "#3699F6",
    "positive": "#A8DC5F",
    "negative": "#D8C4E8",
    "edge_positive": "#8A8A8A",
    "edge_negative": "#9B4A6D",
    "bar_sig": "#DDF5E7",
    "bar_trend": "#FFFFFF",
}

SELECTED_TERMS = [
    "Response to Molecule of Bacterial Origin (GO:0002237)",
    "IL-17 signaling pathway",
    "Defense Response to Bacterium (GO:0042742)",
    "Epithelium Development (GO:0060429)",
    "Epithelial Cell Differentiation (GO:0030855)",
]


def ellipse_width(label: str, strict: bool = False) -> float:
    return max(0.72, min(1.35, 0.35 + 0.090 * len(label))) + (0.08 if strict else 0)


def node_positions(edges: pd.DataFrame) -> dict[str, tuple[float, float]]:
    strict = edges.loc[edges.strict_filtered_link.astype(bool), "target"].tolist()
    other = edges.loc[~edges.strict_filtered_link.astype(bool), "target"].tolist()
    pos = {"KLF4": (0.0, 0.0)}
    strict_angles = np.deg2rad([48, 132, 228, 312][:len(strict)])
    for gene, angle in zip(strict, strict_angles):
        pos[gene] = (1.20 * np.cos(angle), 0.92 * np.sin(angle))
    rings = [(1.95, 1.38, 14), (2.78, 2.02, 16), (3.58, 2.70, 16)]
    offset = np.deg2rad(8)
    index = 0
    for ring_id, (rx, ry, count) in enumerate(rings):
        genes = other[index:index + count]
        angles = np.linspace(0, 2 * np.pi, max(len(genes), 1), endpoint=False)
        angles += offset + ring_id * np.deg2rad(7)
        for gene, angle in zip(genes, angles):
            pos[gene] = (rx * np.cos(angle), ry * np.sin(angle))
        index += count
    return pos


def draw_network(ax, edges: pd.DataFrame) -> None:
    positions = node_positions(edges)
    center = positions["KLF4"]
    for row in edges.itertuples(index=False):
        end = positions[row.target]
        positive = row.direction == "positive"
        patch = FancyArrowPatch(
            center,
            end,
            arrowstyle="-|>" if positive else "-[",
            mutation_scale=6.0,
            linewidth=0.55 if positive else 0.75,
            linestyle="-" if positive else (0, (3, 2)),
            color=COLORS["edge_positive"] if positive else COLORS["edge_negative"],
            shrinkA=18,
            shrinkB=10,
            zorder=1,
        )
        ax.add_patch(patch)

    for row in edges.itertuples(index=False):
        x, y = positions[row.target]
        if bool(row.strict_filtered_link):
            category, text_color = "strict", "black"
        elif row.direction == "positive":
            category, text_color = "positive", "black"
        else:
            category, text_color = "negative", "#5D2B54"
        node = Ellipse(
            (x, y),
            width=ellipse_width(row.target, bool(row.strict_filtered_link)),
            height=0.34,
            facecolor=COLORS[category],
            edgecolor="none",
            zorder=3,
        )
        ax.add_patch(node)
        ax.text(x, y, row.target, ha="center", va="center", fontsize=7,
                color=text_color, zorder=4)

    tf = Ellipse(center, width=1.22, height=0.46, facecolor=COLORS["TF"],
                 edgecolor="none", zorder=5)
    ax.add_patch(tf)
    ax.text(*center, "KLF4", ha="center", va="center", fontsize=10,
            fontweight="bold", zorder=6)
    ax.set_xlim(-4.35, 4.35)
    ax.set_ylim(-3.35, 3.35)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("KLF4-mediated Fib_K14 gene network", fontsize=11, pad=4)

    legend_y = -3.12
    entries = [
        ("Strict CellOracle target", COLORS["strict"]),
        ("Extended positive target", COLORS["positive"]),
        ("Extended negative target", COLORS["negative"]),
    ]
    start_x = -3.60
    for i, (label, color) in enumerate(entries):
        x = start_x + i * 2.45
        ax.add_patch(Ellipse((x, legend_y), 0.42, 0.20, facecolor=color,
                             edgecolor="none", zorder=3))
        ax.text(x + 0.30, legend_y, label, va="center", ha="left", fontsize=7)


def selected_enrichment(ora: pd.DataFrame) -> pd.DataFrame:
    order = {term: i for i, term in enumerate(SELECTED_TERMS)}
    out = ora[ora["Term"].isin(SELECTED_TERMS)].copy()
    out["figure_order"] = out["Term"].map(order)
    out = out.sort_values("figure_order")
    if len(out) != len(SELECTED_TERMS):
        missing = sorted(set(SELECTED_TERMS) - set(out.Term))
        raise RuntimeError(f"Selected enrichment terms missing: {missing}")
    out["minus_log10_FDR"] = -np.log10(out["Adjusted P-value"].clip(lower=1e-300))
    out["evidence_state"] = np.where(out["Adjusted P-value"] < 0.05,
                                     "FDR < 0.05", "epithelial trend; FDR >= 0.05")
    return out


def short_term(term: str) -> str:
    return term.split(" (GO:")[0]


def draw_enrichment(ax, selected: pd.DataFrame) -> None:
    x = np.arange(len(selected))
    colors = [COLORS["bar_sig"] if q < 0.05 else COLORS["bar_trend"]
              for q in selected["Adjusted P-value"]]
    bars = ax.bar(x, selected.minus_log10_FDR, width=0.84, color=colors,
                  edgecolor="#333333", linewidth=0.65)
    for bar, row in zip(bars, selected.to_dict("records")):
        epithelial = "Epithel" in row["Term"]
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            0.18,
            short_term(row["Term"]),
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=7,
            color="#D7191C" if epithelial else "black",
        )
        if row["Adjusted P-value"] >= 0.05:
            bar.set_hatch("///")
    ax.axhline(-np.log10(0.05), color="#777777", lw=0.8, ls=(0, (3, 2)))
    ax.text(len(selected) - 0.55, -np.log10(0.05) + 0.08, "FDR = 0.05",
            ha="right", va="bottom", fontsize=7, color="#555555")
    ax.set_ylabel("−log10(FDR)", fontsize=9)
    ax.set_xlabel("GO/KEGG terms", fontsize=10)
    ax.set_xticks([])
    ax.set_xlim(-0.55, len(selected) - 0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("KLF4 target enrichment", fontsize=11, pad=4)


def save(fig, stem: Path, tight: bool = True) -> None:
    extra = {"bbox_inches": "tight"} if tight else {}
    fig.savefig(stem.with_suffix(".pdf"), **extra)
    fig.savefig(stem.with_suffix(".png"), dpi=600, **extra)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"}, **extra)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", required=True, type=Path)
    ap.add_argument("--ora", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })
    edges = pd.read_csv(args.edges, sep="\t")
    ora = pd.read_csv(args.ora, sep="\t")
    selected = selected_enrichment(ora)
    selected.to_csv(args.output / "Fib_K14_KLF4_ORA_selected_for_figure.tsv",
                    sep="\t", index=False)

    fig_g, ax_g = plt.subplots(figsize=(8.2, 6.5))
    draw_network(ax_g, edges)
    fig_g.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.94)
    fig_g.text(0.015, 0.975, "G", fontsize=12, fontweight="bold", va="top")
    save(fig_g, args.output / "Panel_G_Fib_K14_KLF4_GRN")
    plt.close(fig_g)

    fig_h, ax_h = plt.subplots(figsize=(3.6, 6.5))
    draw_enrichment(ax_h, selected)
    fig_h.subplots_adjust(left=0.20, right=0.98, bottom=0.10, top=0.92)
    fig_h.text(0.025, 0.975, "H", fontsize=12, fontweight="bold", va="top")
    save(fig_h, args.output / "Panel_H_Fib_K14_KLF4_target_enrichment")
    plt.close(fig_h)

    fig = plt.figure(figsize=(11.69, 8.27))
    ax_g = fig.add_axes([0.035, 0.11, 0.69, 0.80])
    ax_h = fig.add_axes([0.755, 0.20, 0.225, 0.66])
    draw_network(ax_g, edges)
    draw_enrichment(ax_h, selected)
    fig.text(0.012, 0.955, "G", fontsize=12, fontweight="bold", va="top")
    fig.text(0.735, 0.955, "H", fontsize=12, fontweight="bold", va="top")
    save(fig, args.output / "Figure_GH_Fib_K14_KLF4_GRN_enrichment_A4_landscape", tight=False)
    plt.close(fig)

    (args.output / "figure_evidence_checkpoint.json").write_text(json.dumps({
        "status": "PASS",
        "analysis_version": "human_fibk14_klf4_grn_enrichment_v01",
        "visual_revision": "v02_csc_species_aware_faithful_style_A4_fontfix",
        "network_edges": int(len(edges)),
        "strict_network_edges": int(edges.strict_filtered_link.astype(bool).sum()),
        "positive_network_edges": int((edges.direction == "positive").sum()),
        "negative_network_edges": int((edges.direction == "negative").sum()),
        "enrichment_terms_displayed": int(len(selected)),
        "significant_terms_displayed": int((selected["Adjusted P-value"] < 0.05).sum()),
        "epithelial_trend_terms_displayed": int((selected["Adjusted P-value"] >= 0.05).sum()),
        "statistical_label": "ORA gene-set statistic; not donor-replicated clinical inference"
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
