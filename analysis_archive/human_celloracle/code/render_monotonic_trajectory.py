#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression


PALETTE = {
    "KC_Basal": "#008941", "KC_Basal_Mig": "#FF3030", "KC_Basal_Prolif": "#FFD700",
    "KC_Spinous": "#00FF00", "KC_Spinous_Mig": "#1A0099", "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17", "Fib_K14": "#4682B4",
}
ORDER = ["Fib_K14", "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
         "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular"]


def save(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--seed", type=int, default=20260811)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family":"Arial", "font.size":8, "axes.titlesize":10,
                         "axes.labelsize":9, "legend.fontsize":7, "pdf.fonttype":42,
                         "ps.fonttype":42, "axes.linewidth":0.8})
    d = pd.read_csv(args.source, sep="\t", compression="gzip")
    d = d[d.gene.isin(["KLF4", "KRT5"]) & np.isfinite(d.pseudotime)].copy()
    rng = np.random.default_rng(args.seed)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2), sharex=True)
    rows = []
    for ax, gene in zip(axes, ["KLF4", "KRT5"]):
        q = d[d.gene == gene].copy()
        if len(q) > 90000:
            q = q.iloc[np.sort(rng.choice(len(q), 90000, replace=False))]
        for lab in ORDER:
            z = q[q.sub_labels == lab]
            ax.scatter(z.pseudotime, z.expression, s=.65, alpha=.42,
                       color=PALETTE[lab], linewidths=0, rasterized=True, label=lab)
        bins = np.linspace(q.pseudotime.min(), q.pseudotime.max(), 121)
        bid = np.clip(np.digitize(q.pseudotime, bins)-1, 0, 119)
        b = q.assign(bin=bid).groupby("bin", observed=True).agg(
            pseudotime=("pseudotime","median"), expression=("expression","mean"), n=("expression","size")
        ).query("n >= 30")
        direction = bool(b.expression.iloc[-10:].mean() >= b.expression.iloc[:10].mean())
        iso = IsotonicRegression(increasing=direction, out_of_bounds="clip")
        yhat = iso.fit_transform(b.pseudotime, b.expression, sample_weight=b.n)
        ax.plot(b.pseudotime, yhat, color="black", lw=2.0,
                label=f"Isotonic fit ({'increasing' if direction else 'decreasing'})")
        ss_res = float(np.sum(b.n * (b.expression-yhat)**2))
        ss_tot = float(np.sum(b.n * (b.expression-np.average(b.expression,weights=b.n))**2))
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
        rows.append({"gene":gene,"direction":"increasing" if direction else "decreasing",
                     "weighted_R2":r2,"n_cells":len(q)})
        ax.set_title(f"{gene} expression", fontweight="bold")
        ax.set_xlabel("Monocle3 pseudotime")
        ax.set_ylabel("Log-normalized expression")
        ax.spines[["top","right"]].set_visible(False)
        ax.text(.02,.96,"Fib_K14 root",transform=ax.transAxes,ha="left",va="top")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:8], labels[:8], frameon=False, ncol=4, loc="lower center",
               bbox_to_anchor=(.5,-.05), markerscale=4, columnspacing=1.0)
    fig.subplots_adjust(bottom=.25, wspace=.30)
    save(fig, args.output / "Human_FibK14_allKC_KLF4_KRT5_monotonic_trajectory")
    pd.DataFrame(rows).to_csv(args.output / "monotonic_fit_statistics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
