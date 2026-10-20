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
FIG_SIZE = (70.0 * MM, 59.4 * MM)
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


def add_core_contours(ax, data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    fib_points = data[data.sub_labels == "Fib_K14"][["UMAP_1", "UMAP_2"]].to_numpy()
    kc_points = data[data.sub_labels == "KC_Basal"][["UMAP_1", "UMAP_2"]].to_numpy()
    fib_poly = central_hull(fib_points, 0.65)
    kc_poly = central_hull(kc_points, 0.60)
    ax.plot(fib_poly[:, 0], fib_poly[:, 1], color="#D7191C", lw=0.9,
            ls=(0, (3, 2)), zorder=5)
    ax.plot(kc_poly[:, 0], kc_poly[:, 1], color="#2166AC", lw=0.9,
            ls=(0, (3, 2)), zorder=5)
    fib, kc = state_centers(data)
    ax.text(fib[0] + 1.1, fib[1] - 1.1, "Fib_K14", color="#B30000",
            fontsize=8, fontweight="bold", ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.78},
            zorder=8)
    ax.text(kc[0] - 1.3, kc[1] + 1.0, "KC", color="#174A7E",
            fontsize=8, fontweight="bold", ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.78},
            zorder=8)
    return fib, kc


def add_direction_arrow(ax, fib: np.ndarray, kc: np.ndarray, condition: str) -> dict:
    if condition == "OE":
        start = fib + np.array([-0.25, 1.25])
        end = kc + np.array([0.15, -0.85])
        color, label, rad = "#2E7D16", "Promotion", -0.12
    else:
        start = kc + np.array([0.15, -0.85])
        end = fib + np.array([-0.25, 1.25])
        color, label, rad = "#A5005A", "Inhibition", 0.12
    arrow = FancyArrowPatch(
        posA=start, posB=end,
        arrowstyle="-|>", mutation_scale=9,
        connectionstyle=f"arc3,rad={rad}",
        lw=1.6, linestyle=(0, (3, 2)), color=color,
        zorder=9,
    )
    ax.add_patch(arrow)
    mid = (start + end) / 2
    ax.text(mid[0] + 5.2, mid[1] + (0.2 if condition == "KO" else -0.2), label,
            color=color, fontsize=8, fontweight="bold",
            rotation=0,
            ha="center", va="center",
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.86},
            zorder=10)
    return {
        "condition": condition,
        "display_label": label,
        "start_state": "Fib_K14" if condition == "OE" else "KC_Basal core",
        "end_state": "KC_Basal core" if condition == "OE" else "Fib_K14",
        "start_UMAP_1": float(start[0]), "start_UMAP_2": float(start[1]),
        "end_UMAP_1": float(end[0]), "end_UMAP_2": float(end[1]),
        "overlay_evidence": "user_requested_interpretive_differentiation_direction",
        "raw_CellOracle_grid_arrows_shown": False,
        "SAC_region_arrows": False,
    }


def draw_panel(ax, atlas_plot: pd.DataFrame, data: pd.DataFrame,
               arrow_data: pd.DataFrame,
               condition: str, title: str) -> dict:
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    ax.scatter(atlas_plot.UMAP_1, atlas_plot.UMAP_2,
               s=0.05, c="#D9D9D9", alpha=0.22,
               linewidths=0, rasterized=True, zorder=1)
    ax.scatter(data.UMAP_1, data.UMAP_2,
               s=0.34, c=data.KC_proxy_perturbation_score,
               cmap=CMAP, norm=norm, alpha=0.64,
               linewidths=0, rasterized=True, zorder=3)
    add_algorithm_arrows(ax, arrow_data)
    fib, kc = add_core_contours(ax, data)
    direction = add_direction_arrow(ax, fib, kc, condition)

    ax.text(8.8, 3.8, "SAC", color="#888888", fontsize=8,
            ha="center", va="center", zorder=7)

    ax.set_title(title, pad=1.5)
    ax.set_xlim(-18.2, 18.5)
    ax.set_ylim(-19.3, 13.7)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    return direction


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

    directions = []
    fig = plt.figure(figsize=FIG_SIZE)
    ax1 = fig.add_axes([0.035, 0.27, 0.445, 0.60])
    ax2 = fig.add_axes([0.52, 0.27, 0.445, 0.60])
    directions.append(draw_panel(
        ax1, atlas_plot, data["KO"], arrows["KO"], "KO",
        r"$\it{In\ silico}$" + "\n" + r"$\it{KLF4}$ knockout",
    ))
    directions.append(draw_panel(
        ax2, atlas_plot, data["OE"], arrows["OE"], "OE",
        r"$\it{In\ silico}$" + "\n" + r"$\it{KLF4}$ overexpression",
    ))
    fig.text(0.012, 0.975, "I", fontsize=12, fontweight="bold", va="top")
    fig.text(0.500, 0.975, "J", fontsize=12, fontweight="bold", va="top")
    add_color_key(fig, [0.20, 0.075, 0.62, 0.026])
    save_exact(fig, args.output / "KLF4_KO_OE_compact_70x59p4mm")
    plt.close(fig)

    # Larger independent panels for optional Illustrator assembly. Each keeps
    # the same requested physical size and evidence grammar.
    for idx, (condition, title, label) in enumerate([
        ("KO", r"$\it{In\ silico}$" + "\n" + r"$\it{KLF4}$ knockout", "I"),
        ("OE", r"$\it{In\ silico}$" + "\n" + r"$\it{KLF4}$ overexpression", "J"),
    ]):
        f = plt.figure(figsize=FIG_SIZE)
        ax = f.add_axes([0.08, 0.24, 0.84, 0.66])
        draw_panel(ax, atlas_plot, data[condition], arrows[condition], condition, title)
        f.text(0.018, 0.975, label, fontsize=12, fontweight="bold", va="top")
        add_color_key(f, [0.22, 0.07, 0.58, 0.027])
        save_exact(f, args.output / f"KLF4_{condition}_compact_70x59p4mm")
        plt.close(f)

    pd.DataFrame(directions).to_csv(
        args.output / "direction_overlay_audit.tsv", sep="\t", index=False
    )
    pd.concat(arrows.values(), ignore_index=True).to_csv(
        args.output / "algorithm_grid_arrow_mask_audit.tsv", sep="\t", index=False
    )
    checkpoint = {
        "status": "PASS_PENDING_PHYSICAL_QA",
        "analysis_version": "human_fibk14_klf4_celloracle_v01",
        "visual_revision": "perturbation_compact_v03_70x59.4mm",
        "analysis_policy": "visual_revision_only",
        "atlas_sha256": sha256(args.atlas),
        "KO_source_sha256": sha256(args.analysis / "KLF4_KO_cell_vectors_with_score.tsv.gz"),
        "OE_source_sha256": sha256(args.analysis / "KLF4_OE_cell_vectors_with_score.tsv.gz"),
        "width_mm": 70.0,
        "height_mm": 59.4,
        "font_family": "Arial",
        "font_range_pt": [8, 12],
        "raw_grid_arrows_removed": False,
        "algorithm_grid_arrows": "raw CellOracle directions, normalized only for display",
        "algorithm_arrow_region_mask": "25-nearest full-atlas labels; KC or Fib_K14 region; SAC fraction < 0.20",
        "algorithm_arrow_magnitude_filter": "mass-pass and condition-specific magnitude >= q25",
        "direction_arrows": "user-requested interpretive summary restricted to Fib_K14-KC_Basal core",
        "SAC_arrows": False,
        "score_background": "inherited CellOracle KC-proxy perturbation score",
    }
    (args.output / "perturbation_compact_checkpoint.json").write_text(
        json.dumps(checkpoint, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
