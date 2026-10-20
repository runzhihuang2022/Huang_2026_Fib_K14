#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression


MM = 1 / 25.4
PANEL_SIZE = (52.5 * MM, 49.5 * MM)
GENES = ["KLF4", "KRT5", "TP63", "KRT14", "EPCAM", "TACSTD2", "VIM", "PDGFRA", "COL1A1"]

# Original atlas category palette. Fib_K14 is restored to the red used in the
# source atlas rather than the blue proxy used in the first trajectory draft.
PALETTE = {
    "Fib_K14": "#E41A1C",
    "KC_Basal": "#008941",
    "KC_Basal_Mig": "#FF3030",
    "KC_Basal_Prolif": "#FFD700",
    "KC_Spinous": "#00D93A",
    "KC_Spinous_Mig": "#1A0099",
    "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17",
}
ORDER = [
    "Fib_K14", "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif",
    "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
]
LEGEND_LABELS = {
    "Fib_K14": "Fib_K14 (root state)",
    "KC_Basal": "KC_Basal",
    "KC_Basal_Mig": "KC_Basal_Mig",
    "KC_Basal_Prolif": "KC_Basal_Prolif",
    "KC_Spinous": "KC_Spinous",
    "KC_Spinous_Mig": "KC_Spinous_Mig",
    "KC_Spinous_Mat": "KC_Spinous_Mat",
    "KC_Granular": "KC_Granular",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    })


def fit_monotonic(q: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, str, float]:
    bins = np.linspace(q.pseudotime.min(), q.pseudotime.max(), 121)
    bid = np.clip(np.digitize(q.pseudotime, bins) - 1, 0, 119)
    b = (
        q.assign(bin=bid)
        .groupby("bin", observed=True)
        .agg(pseudotime=("pseudotime", "median"),
             expression=("expression", "mean"),
             n=("expression", "size"))
        .query("n >= 30")
        .reset_index(drop=True)
    )
    increasing = bool(b.expression.iloc[-10:].mean() >= b.expression.iloc[:10].mean())
    model = IsotonicRegression(increasing=increasing, out_of_bounds="clip")
    yhat = model.fit_transform(b.pseudotime, b.expression, sample_weight=b.n)
    mean = np.average(b.expression, weights=b.n)
    ss_res = float(np.sum(b.n * (b.expression - yhat) ** 2))
    ss_tot = float(np.sum(b.n * (b.expression - mean) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return b, yhat, "increasing" if increasing else "decreasing", r2


def draw_gene(ax, plot_q: pd.DataFrame, gene: str):
    for label in ORDER:
        z = plot_q[plot_q.sub_labels == label]
        ax.scatter(
            z.pseudotime, z.expression,
            s=0.45, alpha=0.34, color=PALETTE[label],
            linewidths=0, rasterized=True,
        )
    b, yhat, direction, r2 = fit_monotonic(plot_q)
    ax.plot(b.pseudotime, yhat, color="black", lw=1.4, zorder=8)
    ax.set_title(rf"$\it{{{gene}}}$ expression", fontweight="bold", pad=2)
    ax.set_xlabel("Monocle3 pseudotime", labelpad=2)
    ax.set_ylabel("Log-normalized expression", labelpad=2)
    ax.text(0.025, 0.955, "Fib_K14 root", transform=ax.transAxes,
            ha="left", va="top", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=0.8, pad=2)
    ax.margins(x=0.01, y=0.04)
    return {
        "gene": gene,
        "direction": direction,
        "weighted_R2": r2,
        "n_cells_fit": int(len(plot_q)),
        "n_points_displayed": int(len(plot_q)),
    }


def save_exact(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})


def legend_handles():
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=4.2,
               markerfacecolor=PALETTE[label], markeredgecolor="none",
               label=LEGEND_LABELS[label])
        for label in ORDER
    ]
    handles.append(Line2D([0], [0], color="black", lw=1.4,
                          label="Monotonic fit"))
    return handles


def render_shared_legend(output: Path) -> None:
    fig = plt.figure(figsize=(157.5 * MM, 17 * MM))
    fig.legend(handles=legend_handles(), frameon=False, ncol=5,
               loc="center", bbox_to_anchor=(0.5, 0.52),
               handlelength=1.25, handletextpad=0.35,
               columnspacing=1.0, labelspacing=0.55)
    save_exact(fig, output / "Shared_celltype_legend_157p5x17mm")
    plt.close(fig)


def render_montage(sampled: dict[str, pd.DataFrame], output: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(157.5 * MM, 165.5 * MM))
    for ax, gene in zip(axes.flat, GENES):
        draw_gene(ax, sampled[gene], gene)
    fig.legend(handles=legend_handles(), frameon=False, ncol=5,
               loc="lower center", bbox_to_anchor=(0.5, 0.012),
               handlelength=1.25, handletextpad=0.35,
               columnspacing=0.9, labelspacing=0.45)
    fig.subplots_adjust(left=0.08, right=0.985, top=0.975, bottom=0.13,
                        wspace=0.34, hspace=0.43)
    save_exact(fig, output / "Human_FibK14_allKC_9gene_monotonic_montage")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260811)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    configure()

    data = pd.read_csv(args.source, sep="\t", compression="gzip")
    data = data[data.gene.isin(GENES) & np.isfinite(data.pseudotime)].copy()
    missing_genes = sorted(set(GENES) - set(data.gene.unique()))
    missing_labels = sorted(set(ORDER) - set(data.sub_labels.unique()))
    if missing_genes or missing_labels:
        raise RuntimeError(f"Missing genes={missing_genes}; missing labels={missing_labels}")

    # Match the v01 KLF4/KRT5 renderer exactly: one seeded generator, genes in
    # declared order, and both display and isotonic fit use the same 90,000-cell
    # deterministic subset. New genes inherit this unchanged rule.
    rng = np.random.default_rng(args.seed)
    sampled = {}
    source_n = {}
    for gene in GENES:
        q = data[data.gene == gene].copy()
        source_n[gene] = int(len(q))
        if len(q) > 90000:
            take = np.sort(rng.choice(len(q), 90000, replace=False))
            q = q.iloc[take].copy()
        sampled[gene] = q

    rows = []
    for gene in GENES:
        fig = plt.figure(figsize=PANEL_SIZE)
        ax = fig.add_axes([0.20, 0.20, 0.76, 0.68])
        row = draw_gene(ax, sampled[gene], gene)
        row["n_cells_source"] = source_n[gene]
        rows.append(row)
        save_exact(fig, args.output / f"Trajectory_{gene}_52p5x49p5mm")
        plt.close(fig)

    render_shared_legend(args.output)
    render_montage(sampled, args.output)

    pd.DataFrame(rows).to_csv(
        args.output / "nine_gene_monotonic_fit_statistics.tsv", sep="\t", index=False
    )
    pd.DataFrame({
        "sub_labels": ORDER,
        "legend_label": [LEGEND_LABELS[x] for x in ORDER],
        "hex_color": [PALETTE[x] for x in ORDER],
        "source": ["original atlas category palette"] * len(ORDER),
    }).to_csv(args.output / "celltype_palette_and_legend.tsv", sep="\t", index=False)
    checkpoint = {
        "status": "PASS_PENDING_PHYSICAL_QA",
        "analysis_version": "human_fibk14_allkc_monocle3_trajectory9gene_v02",
        "parent_monocle3_analysis": "human_fibk14_allkc_monocle3_20260811_v01",
        "visual_revision": "trajectory_compact_v02_52.5x49.5mm",
        "source_sha256": sha256(args.source),
        "genes": GENES,
        "cell_types": ORDER,
        "panel_width_mm": 52.5,
        "panel_height_mm": 49.5,
        "font_family": "Arial",
        "font_range_pt": [8, 12],
        "points_are_display_downsampled": True,
        "trend_fit_uses_deterministic_90000_cell_subset": True,
        "fit_rule_matches_original_KLF4_KRT5_renderer": True,
        "scatter_seed": args.seed,
    }
    (args.output / "trajectory_compact_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
