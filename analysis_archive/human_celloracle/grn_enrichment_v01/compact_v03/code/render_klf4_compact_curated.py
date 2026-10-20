#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, Patch
import numpy as np
import pandas as pd


MM = 1 / 25.4
PANEL_SIZE = (52.5 * MM, 59.4 * MM)
COLORS = {
    "tf": "#FF3B3B",
    "strict": "#3699F6",
    "positive": "#A8DC5F",
    "negative": "#D8C4E8",
    "positive_edge": "#7A7A7A",
    "negative_edge": "#9B4A6D",
}

REQUESTED = {
    "TP63": "positive", "KRT5": "positive", "KRT14": "positive",
    "EPCAM": "positive", "TACSTD2": "positive",
    "VIM": "negative", "PDGFRA": "negative", "COL1A1": "negative",
}
STRICT = ["FABP5", "DMKN", "S100A8", "S100A9"]
POSITIONS = {
    "KLF4": (0.00, 0.00),
    "KRT14": (0.00, 1.38), "KRT5": (-1.05, 1.12), "TACSTD2": (1.05, 1.12),
    "TP63": (-1.38, 0.55), "EPCAM": (1.38, 0.55),
    "FABP5": (-0.63, 0.48), "DMKN": (0.63, 0.48),
    "S100A8": (-0.63, -0.52), "S100A9": (0.63, -0.52),
    "VIM": (-1.38, -0.43), "PDGFRA": (1.38, -0.43), "COL1A1": (0.00, -1.28),
}

SELECTED_TERMS = [
    "Response to Molecule of Bacterial Origin (GO:0002237)",
    "IL-17 signaling pathway",
    "Defense Response to Bacterium (GO:0042742)",
    "Epithelium Development (GO:0060429)",
    "Epithelial Cell Differentiation (GO:0030855)",
]


def build_evidence(edges: pd.DataFrame, universe: set[str]) -> pd.DataFrame:
    rows = []
    for gene, requested in REQUESTED.items():
        q = edges[edges.target.eq(gene)]
        if len(q):
            r = q.sort_values("evidence_score", ascending=False).iloc[0]
            matched = r.direction == requested
            significant = r.FDR < 0.05
            if matched and significant:
                status = "model_supported_direction_matched"
            elif matched:
                status = "hypothesis_direction_matched_not_FDR_significant"
            else:
                status = "hypothesis_direction_conflicts_with_model"
            rows.append({
                "gene": gene, "requested_direction": requested,
                "in_model_universe": gene in universe, "edge_present": True,
                "model_direction": r.direction, "coef_mean": r.coef_mean,
                "p": r.p, "FDR": r.FDR, "evidence_status": status,
                "display_as_model_edge": bool(matched and significant),
            })
        else:
            rows.append({
                "gene": gene, "requested_direction": requested,
                "in_model_universe": gene in universe, "edge_present": False,
                "model_direction": "not_available", "coef_mean": np.nan,
                "p": np.nan, "FDR": np.nan,
                "evidence_status": "hypothesis_no_CellOracle_edge",
                "display_as_model_edge": False,
            })
    return pd.DataFrame(rows)


def node_width(gene: str) -> float:
    return max(0.56, min(0.91, 0.29 + 0.070 * len(gene)))


def draw_edge(ax, gene: str, direction: str, supported: bool) -> None:
    end = POSITIONS[gene]
    if supported:
        color = COLORS["positive_edge"] if direction == "positive" else COLORS["negative_edge"]
        style = "-|>" if direction == "positive" else "-["
        ls = "-" if direction == "positive" else (0, (3, 2))
        lw = 0.55
    else:
        color = COLORS["positive"] if direction == "positive" else COLORS["negative_edge"]
        style = "-|>" if direction == "positive" else "-["
        ls = (0, (1.2, 1.8))
        lw = 0.75
    ax.add_patch(FancyArrowPatch(
        POSITIONS["KLF4"], end, arrowstyle=style, mutation_scale=5.2,
        linewidth=lw, linestyle=ls, color=color, shrinkA=12, shrinkB=8, zorder=1
    ))


def draw_network(ax, evidence: pd.DataFrame) -> None:
    strict_set = set(STRICT)
    evidence_map = evidence.set_index("gene").to_dict("index")
    for gene in STRICT:
        draw_edge(ax, gene, "positive", True)
    for gene, requested in REQUESTED.items():
        draw_edge(ax, gene, requested, bool(evidence_map[gene]["display_as_model_edge"]))

    for gene in STRICT + list(REQUESTED):
        x, y = POSITIONS[gene]
        if gene in strict_set:
            face, edge, ls, lw = COLORS["strict"], "none", "-", 0
        else:
            ev = evidence_map[gene]
            direction = ev["requested_direction"]
            if ev["display_as_model_edge"]:
                face = COLORS[direction]
                edge, ls, lw = "none", "-", 0
            else:
                face = "white"
                edge = COLORS["positive"] if direction == "positive" else COLORS["negative_edge"]
                ls, lw = (0, (2, 1.5)), 1.0
        ax.add_patch(Ellipse((x, y), node_width(gene), 0.29,
                             facecolor=face, edgecolor=edge, linestyle=ls,
                             linewidth=lw, zorder=3))
        ax.text(
            x, y, gene,
            ha="center", va="center", fontsize=7,
            fontstyle="italic", zorder=4,
        )

    ax.add_patch(Ellipse((0, 0), 0.86, 0.38, facecolor=COLORS["tf"],
                         edgecolor="none", zorder=5))
    ax.text(
        0, 0, "KLF4", fontsize=9,
        fontweight="bold", fontstyle="italic",
        ha="center", va="center", zorder=6,
    )
    ax.set_xlim(-1.83, 1.83); ax.set_ylim(-1.57, 1.67)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("KLF4 network in Fib_K14", fontsize=8.5, pad=1.5)


def add_network_legend(fig) -> None:
    handles = [
        Patch(facecolor=COLORS["strict"], edgecolor="none", label="Strict"),
        Patch(facecolor=COLORS["positive"], edgecolor="none", label="GRN +"),
        Patch(facecolor=COLORS["negative"], edgecolor="none", label="GRN −"),
        Patch(facecolor="white", edgecolor="#777777", linestyle=(0, (2, 1.5)), label="Hypothesis"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.50, 0.025),
               ncol=2, frameon=False, fontsize=7, handlelength=1.15,
               columnspacing=0.8, handletextpad=0.35)


def selected_enrichment(ora: pd.DataFrame) -> pd.DataFrame:
    order = {t: i for i, t in enumerate(SELECTED_TERMS)}
    q = ora[ora.Term.isin(SELECTED_TERMS)].copy()
    q["order"] = q.Term.map(order); q = q.sort_values("order")
    q["minus_log10_FDR"] = -np.log10(q["Adjusted P-value"].clip(lower=1e-300))
    return q


def short(term: str) -> str:
    mapping = {
        SELECTED_TERMS[0]: "Bacterial molecule response",
        SELECTED_TERMS[1]: "IL-17 signaling",
        SELECTED_TERMS[2]: "Defense response",
        SELECTED_TERMS[3]: "Epithelium development",
        SELECTED_TERMS[4]: "Epithelial differentiation",
    }
    return mapping[term]


def draw_enrichment(ax, q: pd.DataFrame) -> None:
    x = np.arange(len(q))
    sig = q["Adjusted P-value"].to_numpy() < 0.05
    bars = ax.bar(x, q.minus_log10_FDR, width=0.80,
                  color=np.where(sig, "#DDF5E7", "white"),
                  edgecolor="#333333", linewidth=0.6)
    for bar, is_sig, term in zip(bars, sig, q.Term):
        if not is_sig:
            bar.set_hatch("///")
        ax.text(bar.get_x() + bar.get_width()/2, 0.10, short(term),
                rotation=90, ha="center", va="bottom", fontsize=7,
                color="#D7191C" if "Epithel" in term else "black")
    cutoff = -np.log10(0.05)
    ax.axhline(cutoff, color="#777777", lw=0.7, ls=(0, (3, 2)))
    ax.set_ylabel("−log10(FDR)", fontsize=7.5)
    ax.set_xlabel("GO/KEGG terms", fontsize=7.5)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=7, length=2.5)
    ax.set_xlim(-0.52, 4.52)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_title("KLF4 target enrichment", fontsize=8.5, pad=2)


def save_exact(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-edges", required=True, type=Path)
    ap.add_argument("--universe", required=True, type=Path)
    ap.add_argument("--ora", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family": "Arial", "font.size": 7,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.linewidth": 0.7})
    edges = pd.read_csv(args.all_edges, sep="\t")
    universe = set(x.strip() for x in args.universe.read_text().splitlines() if x.strip())
    evidence = build_evidence(edges, universe)
    evidence.to_csv(args.output / "requested_target_evidence_audit.tsv", sep="\t", index=False)
    q = selected_enrichment(pd.read_csv(args.ora, sep="\t"))
    q.to_csv(args.output / "compact_enrichment_source.tsv", sep="\t", index=False)

    fig_g = plt.figure(figsize=PANEL_SIZE)
    ax_g = fig_g.add_axes([0.035, 0.18, 0.93, 0.72]); draw_network(ax_g, evidence)
    fig_g.text(0.015, 0.985, "G", fontsize=12, fontweight="bold", va="top")
    add_network_legend(fig_g)
    save_exact(fig_g, args.output / "Panel_G_KLF4_GRN_compact_52p5x59p4mm")
    plt.close(fig_g)

    fig_h = plt.figure(figsize=PANEL_SIZE)
    ax_h = fig_h.add_axes([0.20, 0.15, 0.76, 0.74]); draw_enrichment(ax_h, q)
    fig_h.text(0.015, 0.985, "H", fontsize=12, fontweight="bold", va="top")
    save_exact(fig_h, args.output / "Panel_H_KLF4_enrichment_compact_52p5x59p4mm")
    plt.close(fig_h)

    # Optional side-by-side assembly preview; each constituent panel retains
    # the requested 52.5 x 59.4 mm physical size.
    fig = plt.figure(figsize=(105.0 * MM, 59.4 * MM))
    ax1 = fig.add_axes([0.02, 0.18, 0.46, 0.72]); draw_network(ax1, evidence)
    ax2 = fig.add_axes([0.57, 0.15, 0.405, 0.74]); draw_enrichment(ax2, q)
    fig.text(0.008, 0.985, "G", fontsize=12, fontweight="bold", va="top")
    fig.text(0.505, 0.985, "H", fontsize=12, fontweight="bold", va="top")
    save_exact(fig, args.output / "Figure_GH_compact_preview_105x59p4mm")
    plt.close(fig)

    checkpoint = {
        "status": "PASS_WITH_HYPOTHESIS_OVERLAY",
        "analysis_version": "human_fibk14_klf4_grn_curated_overlay_v02",
        "visual_revision": "v03_compact_52.5x59.4mm",
        "panel_width_mm": 52.5, "panel_height_mm": 59.4,
        "font_family": "Arial", "font_range_pt": [7, 12],
        "strict_targets": STRICT,
        "requested_targets": REQUESTED,
        "formal_ORA_excludes_hypothesis_nodes": True,
        "evidence_counts": evidence.evidence_status.value_counts().to_dict(),
    }
    (args.output / "compact_figure_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
