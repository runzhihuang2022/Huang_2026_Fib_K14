#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d


PALETTE = {
    "Fib_K14": "#4C91D9", "KC_Basal": "#1B9E77", "KC_Basal_Mig": "#66C2A5",
    "KC_Basal_Prolif": "#E6AB02", "KC_Spinous": "#7570B3", "KC_Spinous_Mig": "#E7298A",
    "KC_Spinous_Mat": "#D95F02", "KC_Granular": "#A6761D",
}


def setup_style() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
        "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.8, "savefig.facecolor": "white",
    })


def smooth_xy(x, y, bins=80):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    edges = np.linspace(x.min(), x.max(), bins + 1)
    ids = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    bx = np.array([np.median(x[ids == i]) if np.any(ids == i) else np.nan for i in range(bins)])
    by = np.array([np.mean(y[ids == i]) if np.any(ids == i) else np.nan for i in range(bins)])
    keep = np.isfinite(bx) & np.isfinite(by)
    bx, by = bx[keep], by[keep]
    if len(bx) < 5:
        return bx, by
    grid = np.linspace(bx.min(), bx.max(), 300)
    # Smooth binned means without cubic overshoot below the non-negative scale.
    by_smooth = gaussian_filter1d(by, sigma=2.0, mode="nearest")
    return grid, np.clip(np.interp(grid, bx, by_smooth), 0, None)


def save_all(fig, stem: Path, width, height):
    fig.set_size_inches(width, height)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-dir", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=20260811)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_style()
    cells = pd.read_csv(args.result_dir / "human_FibK14_allKC_monocle3_cells.tsv.gz", sep="\t")
    expr = pd.read_csv(args.result_dir / "human_FibK14_allKC_focus_expression.tsv.gz", sep="\t")
    df = expr.merge(cells[["cell","pseudotime","sub_labels"]], on="cell", how="left")
    finite = np.isfinite(cells.pseudotime)
    rng = np.random.default_rng(args.seed)
    display_idx = np.flatnonzero(finite)
    if len(display_idx) > 70000:
        display_idx = np.sort(rng.choice(display_idx, 70000, replace=False))

    fig, axes = plt.subplots(1, 2)
    c = cells.iloc[display_idx]
    for lab in sorted(c.sub_labels.unique()):
        q = c[c.sub_labels == lab]
        axes[0].scatter(q.UMAP_1, q.UMAP_2, s=0.8, c=PALETTE.get(lab, "#999999"),
                        label=lab, linewidths=0, alpha=.75, rasterized=True)
    axes[0].set_title("Fib_K14 + keratinocyte states")
    axes[0].legend(frameon=False, markerscale=4, bbox_to_anchor=(.5, -.12),
                   loc="upper center", ncol=4, columnspacing=1.0, handletextpad=.35)
    sc = axes[1].scatter(c.UMAP_1, c.UMAP_2, s=.8, c=c.pseudotime, cmap="viridis",
                         linewidths=0, alpha=.8, rasterized=True)
    axes[1].set_title("Monocle3 pseudotime")
    axes[1].text(.02, .98, "Root: Fib_K14", transform=axes[1].transAxes,
                 ha="left", va="top", fontsize=8)
    fig.colorbar(sc, ax=axes[1], label="Pseudotime", fraction=.046, pad=.04)
    for ax in axes:
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top","right","bottom","left"]].set_visible(False)
    fig.subplots_adjust(wspace=.28, bottom=.24)
    save_all(fig, args.output_dir / "Human_FibK14_allKC_Monocle3_UMAP", 7.0, 3.5)
    plt.close(fig)

    selected = [g for g in ["KLF4", "KRT5"] if g in set(df.gene)]
    fig, axes = plt.subplots(2, 2, sharex=True)
    for row, gene in enumerate(selected):
        g = df[(df.gene == gene) & np.isfinite(df.pseudotime)].copy()
        if len(g) > 80000:
            g = g.iloc[np.sort(rng.choice(len(g), 80000, replace=False))]
        for lab in sorted(g.sub_labels.unique()):
            q = g[g.sub_labels == lab]
            axes[row, 0].scatter(q.pseudotime, q.expression, s=.7, alpha=.55,
                                 c=PALETTE.get(lab,"#999999"), label=lab,
                                 linewidths=0, rasterized=True)
        sx, sy = smooth_xy(g.pseudotime, g.expression)
        axes[row, 0].plot(sx, sy, color="black", lw=1.4)
        axes[row, 0].set_title(f"{gene} by cell type", fontweight="bold")
        m = axes[row, 1].scatter(g.pseudotime, g.expression, c=g.expression, cmap="magma",
                                 s=.7, alpha=.55, linewidths=0, rasterized=True)
        axes[row, 1].plot(sx, sy, color="black", lw=1.4)
        axes[row, 1].set_title(f"{gene} expression trend", fontweight="bold")
        fig.colorbar(m, ax=axes[row, 1], label="Expression", fraction=.046, pad=.04)
        for col in range(2):
            axes[row, col].set_ylabel("log-normalized expression")
            axes[row, col].set_xlabel("Pseudotime")
            axes[row, col].spines[["top","right"]].set_visible(False)
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="center right", bbox_to_anchor=(1.14,.5), markerscale=4)
    fig.subplots_adjust(hspace=.42, wspace=.35, right=.86)
    save_all(fig, args.output_dir / "Human_FibK14_allKC_KLF4_KRT5_pseudotime", 7.0, 5.5)
    plt.close(fig)

    df.to_csv(args.output_dir / "Human_FibK14_allKC_focus_gene_source_data.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    main()
