#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree


MM = 1 / 25.4
FIG_SIZE = (52.5 * MM, 52.5 * MM)
CMAP = LinearSegmentedColormap.from_list(
    "perturbation", ["#A5005A", "#F7F7F7", "#2E7D16"]
)
TARGET_KC = ["KC_Basal", "KC_Basal_Mig", "KC_Spinous_Mig"]


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
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    })


def central_hull(points: np.ndarray, retain: float = 0.68) -> np.ndarray:
    center = np.median(points, axis=0)
    dist = np.linalg.norm(points - center, axis=1)
    core = points[dist <= np.quantile(dist, retain)]
    if len(core) > 3000:
        core = core[np.linspace(0, len(core) - 1, 3000, dtype=int)]
    hull = ConvexHull(core)
    poly = core[hull.vertices]
    return np.vstack([poly, poly[0]])


def state_centers(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    fib = data[data.sub_labels == "Fib_K14"][["UMAP_1", "UMAP_2"]].median().to_numpy()
    # KC_Basal is the epithelial core directly connected to Fib_K14 and avoids
    # the upper/right SAC-adjacent portion of the full KC manifold.
    kc = data[data.sub_labels == "KC_Basal"][["UMAP_1", "UMAP_2"]].median().to_numpy()
    return fib, kc


def prepare_algorithm_arrows(grid: pd.DataFrame, atlas: pd.DataFrame,
                             condition: str) -> pd.DataFrame:
    """Mask raw CellOracle grid vectors by atlas-defined Fib_K14/KC regions.

    Region membership is determined from 25 nearest full-atlas cells. SAC-rich
    positions are explicitly excluded even if a nearby KC point is present.
    """
    xy_atlas = atlas[["UMAP_1", "UMAP_2"]].to_numpy()
    xy_grid = grid[["UMAP_1", "UMAP_2"]].to_numpy()
    _, nn = cKDTree(xy_atlas).query(xy_grid, k=25)
    bulk = atlas.bulk_labels.astype(str).to_numpy()[nn]
    sub = atlas.sub_labels.astype(str).to_numpy()[nn]
    kc_fraction = np.mean(bulk == "Keratinocyte", axis=1)
    fib_fraction = np.mean(sub == "Fib_K14", axis=1)
    sac_fraction = np.mean(np.char.startswith(bulk.astype(str), "SAC"), axis=1)
    region = np.where(kc_fraction >= fib_fraction, "KC", "Fib_K14")
    eligible_region = ((kc_fraction >= 0.56) | (fib_fraction >= 0.16)) & (sac_fraction < 0.20)
    magnitude = np.hypot(grid.delta_1.to_numpy(), grid.delta_2.to_numpy())
    base = eligible_region & grid.mass_pass.astype(bool).to_numpy() & np.isfinite(magnitude)
    threshold = float(np.quantile(magnitude[base], 0.25)) if np.any(base) else np.inf
    display = base & (magnitude >= threshold)
    out = grid.copy()
    out["condition"] = condition
    out["nearest25_KC_fraction"] = kc_fraction
    out["nearest25_Fib_K14_fraction"] = fib_fraction
    out["nearest25_SAC_fraction"] = sac_fraction
    out["assigned_display_region"] = region
    out["magnitude"] = magnitude
    out["magnitude_threshold_q25"] = threshold
    out["display_arrow"] = display
    return out


def add_algorithm_arrows(ax, arrow_data: pd.DataFrame) -> None:
    q = arrow_data[arrow_data.display_arrow].copy()
    magnitude = q.magnitude.to_numpy()
    unit_x = q.delta_1.to_numpy() / np.maximum(magnitude, 1e-12)
    unit_y = q.delta_2.to_numpy() / np.maximum(magnitude, 1e-12)
    ux = np.unique(arrow_data.UMAP_1.to_numpy())
    uy = np.unique(arrow_data.UMAP_2.to_numpy())
    spacing = float(np.median(np.r_[np.diff(ux), np.diff(uy)]))
    ax.quiver(
        q.UMAP_1, q.UMAP_2,
        unit_x * (0.50 * spacing), unit_y * (0.50 * spacing),
        color="#111111", alpha=0.78,
        angles="xy", scale_units="xy", scale=1.0,
        width=0.0022, headwidth=3.0, headlength=3.8,
        headaxislength=3.3, zorder=6,
    )


def aggregate_local_direction(arrow_data: pd.DataFrame, atlas: pd.DataFrame,
                              condition: str, seed: int) -> dict:
    """Project Fib_K14-local CellOracle arrows onto the prespecified Fib-to-KC axis."""
    fib, kc = state_centers(atlas)
    axis = kc - fib
    axis = axis / np.linalg.norm(axis)
    q = arrow_data[
        arrow_data.display_arrow
        & arrow_data.assigned_display_region.eq("Fib_K14")
    ].copy()
    magnitude = np.hypot(q.delta_1.to_numpy(), q.delta_2.to_numpy())
    projected = (
        q.delta_1.to_numpy() * axis[0] + q.delta_2.to_numpy() * axis[1]
    ) / np.maximum(magnitude, 1e-12)
    rng = np.random.default_rng(seed + (0 if condition == "KO" else 1))
    boot = np.array([
        rng.choice(projected, len(projected), replace=True).mean()
        for _ in range(20000)
    ])
    mean_projection = float(np.mean(projected))
    return {
        "condition": condition,
        "arrow_subset": "display arrows assigned to Fib_K14 region",
        "n_arrows": int(len(projected)),
        "axis_definition": "Fib_K14 median to KC_Basal median",
        "mean_unit_projection": mean_projection,
        "median_unit_projection": float(np.median(projected)),
        "fraction_toward_KC": float(np.mean(projected > 0)),
        "bootstrap_ci_low": float(np.quantile(boot, 0.025)),
        "bootstrap_ci_high": float(np.quantile(boot, 0.975)),
        "display_direction": "toward_KC" if mean_projection > 0 else "toward_Fib_K14",
        "inference_status": "descriptive; bootstrap CI includes zero"
            if np.quantile(boot, 0.025) <= 0 <= np.quantile(boot, 0.975)
            else "descriptive; bootstrap CI excludes zero",
    }


def add_aggregate_arrow(ax, atlas: pd.DataFrame, stats: dict) -> None:
    fib, kc = state_centers(atlas)
    if stats["mean_unit_projection"] > 0:
        start = fib + 0.08 * (kc - fib)
        end = fib + 0.62 * (kc - fib)
        color = "#2E7D16"
    else:
        start = fib + 0.62 * (kc - fib)
        end = fib + 0.08 * (kc - fib)
        color = "#A5005A"
    ax.add_patch(FancyArrowPatch(
        posA=start, posB=end, arrowstyle="-|>", mutation_scale=12,
        connectionstyle="arc3,rad=0", lw=2.0, color=color, zorder=10,
    ))
    mid = (start + end) / 2
    ax.text(
        mid[0] + 3.5, mid[1], "Mean projected\ndirection",
        color=color, fontsize=7, fontweight="bold", ha="left", va="center",
        bbox={"boxstyle": "round,pad=0.10", "fc": "white", "ec": "none", "alpha": 0.84},
        zorder=11,
    )


def draw_panel(ax, atlas_plot: pd.DataFrame, data: pd.DataFrame,
               arrow_data: pd.DataFrame,
               condition: str, title: str,
               atlas: pd.DataFrame, summary_stats: dict | None = None) -> None:
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    ax.scatter(atlas_plot.UMAP_1, atlas_plot.UMAP_2,
               s=0.05, c="#D9D9D9", alpha=0.22,
               linewidths=0, rasterized=True, zorder=1)
    ax.scatter(data.UMAP_1, data.UMAP_2,
               s=0.34, c=data.KC_proxy_perturbation_score,
               cmap=CMAP, norm=norm, alpha=0.64,
               linewidths=0, rasterized=True, zorder=3)
    add_algorithm_arrows(ax, arrow_data)
    if summary_stats is not None:
        add_aggregate_arrow(ax, atlas, summary_stats)

    ax.set_title(title, pad=1.5)
    ax.set_xlim(-18.2, 18.5)
    ax.set_ylim(-19.3, 13.7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def add_color_key(fig, rect) -> None:
    cax = fig.add_axes(rect)
    sm = mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), cmap=CMAP)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_ticks([-1, 0, 1])
    cb.set_ticklabels(["Inhibition", "0", "Promotion"])
    cb.ax.tick_params(labelsize=8, length=2, pad=1)
    cb.outline.set_linewidth(0.6)
    fig.text(rect[0] + rect[2] / 2, rect[1] + rect[3] + 0.022,
             "CellOracle perturbation score", ha="center", va="bottom", fontsize=8)
    cb.ax.get_xticklabels()[0].set_color("#A5005A")
    cb.ax.get_xticklabels()[-1].set_color("#2E7D16")


def save_exact(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=600)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260811)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    configure()

    atlas = pd.read_csv(args.atlas, sep="\t", compression="gzip")
    rng = np.random.default_rng(args.seed)
    if len(atlas) > 180000:
        atlas_plot = atlas.iloc[np.sort(rng.choice(len(atlas), 180000, replace=False))].copy()
    else:
        atlas_plot = atlas.copy()
    data = {
        condition: pd.read_csv(
            args.analysis / f"KLF4_{condition}_cell_vectors_with_score.tsv.gz", sep="\t"
        )
        for condition in ["KO", "OE"]
    }
    grid = {
        condition: pd.read_csv(args.analysis / f"KLF4_{condition}_grid_vectors.tsv", sep="\t")
        for condition in ["KO", "OE"]
    }
    arrows = {
        condition: prepare_algorithm_arrows(grid[condition], atlas, condition)
        for condition in ["KO", "OE"]
    }

    summary = {
        condition: aggregate_local_direction(arrows[condition], atlas, condition, args.seed)
        for condition in ["KO", "OE"]
    }
    for condition, title in [
        ("KO", r"$\it{In\ silico}$ $\it{KLF4}$ knockout"),
        ("OE", r"$\it{In\ silico}$ $\it{KLF4}$ overexpression"),
    ]:
        for version, stats in [("local_arrows", None), ("aggregate_direction", summary[condition])]:
            f = plt.figure(figsize=FIG_SIZE)
            ax = f.add_axes([0.055, 0.205, 0.89, 0.70])
            draw_panel(
                ax, atlas_plot, data[condition], arrows[condition],
                condition, title, atlas, stats,
            )
            add_color_key(f, [0.205, 0.065, 0.59, 0.026])
            save_exact(f, args.output / f"KLF4_{condition}_{version}_52p5mm_square")
            plt.close(f)

    pd.DataFrame(summary.values()).to_csv(
        args.output / "aggregate_direction_statistics.tsv", sep="\t", index=False
    )
    pd.concat(arrows.values(), ignore_index=True).to_csv(
        args.output / "algorithm_grid_arrow_mask_audit.tsv", sep="\t", index=False
    )
    checkpoint = {
        "status": "PASS_PENDING_PHYSICAL_QA",
        "analysis_version": "human_fibk14_klf4_celloracle_v01",
        "visual_revision": "perturbation_square_v04_52.5mm",
        "analysis_policy": "visual_revision_only",
        "atlas_sha256": sha256(args.atlas),
        "KO_source_sha256": sha256(args.analysis / "KLF4_KO_cell_vectors_with_score.tsv.gz"),
        "OE_source_sha256": sha256(args.analysis / "KLF4_OE_cell_vectors_with_score.tsv.gz"),
        "width_mm": 52.5,
        "height_mm": 52.5,
        "font_family": "Arial",
        "font_range_pt": [7, 12],
        "raw_grid_arrows_removed": False,
        "algorithm_grid_arrows": "raw CellOracle directions, normalized only for display",
        "algorithm_arrow_region_mask": "25-nearest full-atlas labels; KC or Fib_K14 region; SAC fraction < 0.20",
        "algorithm_arrow_magnitude_filter": "mass-pass and condition-specific magnitude >= q25",
        "contours": "removed at user request",
        "Fib_K14_text_label": "removed at user request",
        "direction_arrows": "descriptive signed mean projection of Fib_K14-local CellOracle arrows onto Fib_K14-to-KC_Basal axis",
        "direction_uncertainty": "bootstrap confidence intervals recorded; do not interpret as significance when CI spans zero",
        "SAC_arrows": False,
        "score_background": "inherited CellOracle KC-proxy perturbation score",
    }
    (args.output / "perturbation_square_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
