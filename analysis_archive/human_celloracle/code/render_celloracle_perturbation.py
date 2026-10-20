#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull


TARGET_STATES = ["KC_Basal", "KC_Basal_Mig", "KC_Spinous_Mig"]
CMAP = LinearSegmentedColormap.from_list("perturb", ["#A5005A", "#F7F7F7", "#2E7D16"])


def orient_grid(a):
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError(f"Expected 2D grid/flow array, got {a.shape}")
    return a.T if a.shape[0] == 2 else a


def cell_scores(d: pd.DataFrame):
    target = d[d.sub_labels.isin(TARGET_STATES)][["UMAP_1","UMAP_2"]].median().to_numpy()
    xy = d[["UMAP_1","UMAP_2"]].to_numpy()
    vec = d[["delta_UMAP_1","delta_UMAP_2"]].to_numpy()
    toward = target[None,:] - xy
    toward /= np.maximum(np.linalg.norm(toward, axis=1, keepdims=True), 1e-8)
    score = np.sum(vec * toward, axis=1)
    scale = np.quantile(np.abs(score), .99)
    return score, max(float(scale), 1e-12), target


def grid_scores(grid, flow, target, scale):
    toward = target[None,:] - grid
    toward /= np.maximum(np.linalg.norm(toward, axis=1, keepdims=True), 1e-8)
    return np.clip(np.sum(flow * toward, axis=1) / scale, -1, 1)


def add_hull(ax, d):
    pts = d[d.sub_labels == "Fib_K14"][["UMAP_1","UMAP_2"]].to_numpy()
    # A full convex hull is distorted by sparse boundary cells.  Use the
    # central 65% of Fib_K14 cells to mark the core state, while retaining all
    # cells in the perturbation scatter and source table.
    center = np.median(pts, axis=0)
    distance = np.linalg.norm(pts - center, axis=1)
    pts = pts[distance <= np.quantile(distance, 0.65)]
    if len(pts) > 2000:
        pts = pts[np.linspace(0, len(pts)-1, 2000, dtype=int)]
    hull = ConvexHull(pts)
    poly = pts[hull.vertices]
    poly = np.vstack([poly, poly[0]])
    ax.plot(poly[:,0], poly[:,1], color="#D7191C", lw=1.0, ls=(0,(3,2)))
    c = np.median(pts, axis=0)
    ax.text(c[0], c[1], "Fib_K14", color="#B30000", fontsize=8, fontweight="bold",
            ha="center", va="center")


def save(fig, stem):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight",
                pil_kwargs={"compression":"tiff_lzw"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, type=Path)
    ap.add_argument("--analysis", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":10,
                         "pdf.fonttype":42,"ps.fonttype":42})
    atlas = pd.read_csv(args.atlas, sep="\t", compression="gzip")
    rng = np.random.default_rng(20260811)
    if len(atlas) > 180000:
        atlas_plot = atlas.iloc[np.sort(rng.choice(len(atlas),180000,replace=False))]
    else:
        atlas_plot = atlas

    fig, axes = plt.subplots(1,2,figsize=(7.0,3.5))
    summaries=[]
    for ax, condition, title in zip(
        axes, ["KO","OE"], ["In silico KLF4 knockout", "In silico KLF4 overexpression"]
    ):
        d = pd.read_csv(args.analysis/f"KLF4_{condition}_cell_vectors.tsv.gz",sep="\t")
        npz = np.load(args.analysis/f"KLF4_{condition}_celloracle_grid.npz")
        grid, flow = orient_grid(npz["flow_grid"]), orient_grid(npz["flow"])
        # CellOracle defines mass_filter=True for low-mass grid positions that
        # must be excluded (not positions that pass the filter).
        mass_excluded = np.asarray(npz["mass_filter"]).ravel().astype(bool)
        mass_pass = ~mass_excluded
        score, scale, target = cell_scores(d)
        d["KC_proxy_perturbation_score"] = np.clip(score/scale,-1,1)
        d.to_csv(args.output/f"KLF4_{condition}_cell_vectors_with_score.tsv.gz",sep="\t",index=False)
        gs = grid_scores(grid,flow,target,scale)
        pd.DataFrame({"UMAP_1":grid[:,0],"UMAP_2":grid[:,1],"delta_1":flow[:,0],
                      "delta_2":flow[:,1],"score":gs,"mass_pass":mass_pass}).to_csv(
            args.output/f"KLF4_{condition}_grid_vectors.tsv",sep="\t",index=False)

        ax.scatter(atlas_plot.UMAP_1,atlas_plot.UMAP_2,s=.08,c="#D9D9D9",alpha=.35,
                   linewidths=0,rasterized=True)
        norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1)
        ax.scatter(d.UMAP_1,d.UMAP_2,s=.65,c=d.KC_proxy_perturbation_score,cmap=CMAP,norm=norm,
                   alpha=.72,linewidths=0,rasterized=True)
        magnitude = np.linalg.norm(flow, axis=1)
        threshold = np.quantile(magnitude[mass_pass], 0.65)
        valid=mass_pass & np.isfinite(gs) & (magnitude >= threshold)
        spacing = np.median(np.r_[np.diff(np.unique(grid[:,0])),
                                  np.diff(np.unique(grid[:,1]))])
        direction = flow[valid] / np.maximum(magnitude[valid,None], 1e-12)
        display_flow = direction * (0.62 * spacing)
        ax.quiver(grid[valid,0],grid[valid,1],display_flow[:,0],display_flow[:,1],
                  color="black",angles="xy",scale_units="xy",scale=1.0,width=.0026,
                  headwidth=3.2,headlength=4.0,headaxislength=3.5)
        add_hull(ax,d)
        ax.set_title(title)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_aspect("equal",adjustable="datalim")
        summaries.append({"condition":condition,"median_score":float(np.median(score)),
                          "mean_score":float(np.mean(score)),"positive_fraction":float(np.mean(score>0)),
                          "robust_scale_q99":scale,"n_cells":len(d)})
    sm=mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-1,vcenter=0,vmax=1),cmap=CMAP)
    cbar=fig.colorbar(sm,ax=axes,orientation="horizontal",fraction=.08,pad=.08,aspect=35)
    cbar.set_label("Perturbation score toward KC proxy")
    cbar.set_ticks([-1,0,1]); cbar.set_ticklabels(["Away / inhibition","0","Toward / promotion"])
    fig.subplots_adjust(bottom=.18,wspace=.08)
    save(fig,args.output/"Human_FibK14_KLF4_CellOracle_KO_OE")
    pd.DataFrame(summaries).to_csv(args.output/"KLF4_KO_OE_summary.tsv",sep="\t",index=False)
    (args.output/"perturbation_score_definition.json").write_text(json.dumps({
        "engine":"CellOracle 0.18.0","score_definition":"dot product of CellOracle simulated embedding shift and unit vector from each cell toward median original-UMAP position of KC_Basal, KC_Basal_Mig and KC_Spinous_Mig; divided by absolute q99",
        "positive":"predicted motion toward KC proxy","negative":"predicted motion away from KC proxy",
        "embedding":"original parent X_umap","background":"full original atlas",
        "arrow_display":"CellOracle grid-flow directions; arrows are direction-normalized to 0.62 grid spacing and restricted to mass-pass vectors above the 65th magnitude percentile for legibility; raw magnitudes are retained in the source table"
    },indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
