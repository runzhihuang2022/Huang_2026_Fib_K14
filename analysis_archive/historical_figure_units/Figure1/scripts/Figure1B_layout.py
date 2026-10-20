from __future__ import annotations

import csv

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from figure_svg_utils import ROOT
from figure1_data_plot_utils import CELL_ORDER, COLORS, DISPLAY_NAMES, SCRNA_H5, read_h5_panel, save_panel, style


WIDTH_MM, HEIGHT_MM = 61, 92
STEM = "Figure1B_v02"


def broad_group(state: str) -> str:
    if state.startswith("KC_"): return "KC"
    if state.startswith("SAC_") or state == "Melanocyte": return "SAC"
    if state.startswith("Fib_"): return "Fib"
    if state.startswith("Endo_") or state in {"Pericyte", "Schwann"}: return "Endo"
    return "Immune"


def main() -> None:
    style()
    umap, labels = read_h5_panel(SCRNA_H5)
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")
    gs = GridSpec(2, 1, figure=fig, height_ratios=[0.63, 0.37], hspace=0.03, left=0.035, right=0.985, top=0.94, bottom=0.025)
    ax = fig.add_subplot(gs[0])
    for state in CELL_ORDER:
        mask = labels == state
        if mask.any():
            ax.scatter(umap[mask, 0], umap[mask, 1], s=0.055, c=COLORS[state], linewidths=0, rasterized=True)
    for idx, state in enumerate(CELL_ORDER, 1):
        mask = labels == state
        if mask.any():
            x, y = np.median(umap[mask], axis=0)
            ax.text(x, y, str(idx), ha="center", va="center", fontsize=5.2, fontweight="bold",
                    bbox={"boxstyle": "circle,pad=0.16", "facecolor": "white", "edgecolor": "#333333", "linewidth": 0.45, "alpha": 0.92})
    for group in ["KC", "SAC", "Fib", "Endo", "Immune"]:
        states = [s for s in CELL_ORDER if broad_group(s) == group]
        mask = np.isin(labels, states)
        x, y = np.median(umap[mask], axis=0)
        offsets = {"KC": (-2.5, 0.0), "SAC": (1.2, 0.3), "Fib": (1.8, -0.6), "Endo": (1.0, 0.8), "Immune": (-2.4, -0.5)}
        dx, dy = offsets[group]
        ax.text(x + dx, y + dy, group, fontsize=7.2, fontweight="bold", ha="center")
    ax.set_aspect("equal", adjustable="datalim")
    ax.axis("off")
    ax.set_title("Human single-cell atlas", fontsize=8.5, fontweight="bold", pad=1)
    x0, y0 = 0.055, 0.06
    ax.annotate("", xy=(x0 + .13, y0), xytext=(x0, y0), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": .7, "color": "black"})
    ax.annotate("", xy=(x0, y0 + .13), xytext=(x0, y0), xycoords="axes fraction", arrowprops={"arrowstyle": "-|>", "lw": .7, "color": "black"})
    ax.text(x0 + .065, y0 - .035, "UMAP 1", transform=ax.transAxes, ha="center", fontsize=5.5)
    ax.text(x0 - .04, y0 + .065, "UMAP 2", transform=ax.transAxes, va="center", rotation=90, fontsize=5.5)

    key = fig.add_subplot(gs[1]); key.set_xlim(0, 1); key.set_ylim(0, 1); key.axis("off")
    key.text(0.0, 0.985, "47 cell states", fontsize=6.3, fontweight="bold", va="top")
    n_rows = 16
    for i, (state, name) in enumerate(zip(CELL_ORDER, DISPLAY_NAMES)):
        col, row = divmod(i, n_rows)
        x = [0.0, 0.335, 0.67][col]; y = 0.90 - row * 0.057
        key.add_patch(Rectangle((x, y - 0.018), 0.020, 0.025, facecolor=COLORS[state], edgecolor="#555555", linewidth=0.25))
        key.text(x + 0.026, y, f"{i+1} {name}", fontsize=5.1, va="center")
    fig.text(0.012, 0.982, "B", fontsize=12, fontweight="bold", va="top")
    save_panel(fig, STEM, WIDTH_MM, HEIGHT_MM)

    with (ROOT / "source_data" / "Figure1B_sources.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["element", "source", "processing"])
        writer.writerow(["UMAP coordinates and labels", SCRNA_H5, "read directly from X_umap and sub_labels; no PPT extraction"])
        writer.writerow(["colors", "original annotation color registry", "same colors used in UMAP and key swatches"])
    (ROOT / "QC" / f"{STEM}_content_proof.txt").write_text(
        "Figure 1B v02\n- Recomputed directly from pbmc_final.h5ad (279,305 cells).\n"
        "- Original cell-state color registry retained and added to the numbered key.\n"
        "- Dense point layer is code-rasterized; labels and key remain vector.\n- No PPT-derived bioinformatics image is used.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
