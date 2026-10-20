from __future__ import annotations

import csv
import json

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from figure_svg_utils import ROOT
from figure1_data_plot_utils import ANCHOR_JSON, COLORS, CONFIG_WOUND, FOCUS_STATES, add_scale_bar, decode_categorical, plot_spatial_labels, save_panel, stereo_file, style


WIDTH_MM, HEIGHT_MM = 190, 53
STEM = "Figure1D_v02"
SAMPLE = CONFIG_WOUND["19dpb_p1"]["sample_id"]
ANGLE = CONFIG_WOUND["19dpb_p1"]["angle"]
REGIONS = ["Un-epi", "Epi-Front", "Newly-epi"]


def distance_to_polyline(points: np.ndarray, line) -> np.ndarray:
    vertices = np.asarray(line, dtype=float); result = np.full(points.shape[0], np.inf)
    for start, end in zip(vertices[:-1], vertices[1:]):
        vector = end - start; denominator = float(np.dot(vector, vector))
        if denominator <= np.finfo(float).eps: continue
        delta = points - start
        projection = np.clip((delta[:, 0] * vector[0] + delta[:, 1] * vector[1]) / denominator, 0, 1)
        nearest = start + projection[:, None] * vector
        result = np.minimum(result, np.sqrt(((points - nearest) ** 2).sum(axis=1)))
    return result


def points_in_polygon(points: np.ndarray, vertices) -> np.ndarray:
    polygon = np.asarray(vertices, dtype=float); inside = np.zeros(points.shape[0], dtype=bool); x, y = points[:, 0], points[:, 1]; j = len(polygon) - 1
    for i in range(len(polygon)):
        xi, yi = polygon[i]; xj, yj = polygon[j]
        crosses = ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / ((yj - yi) + np.finfo(float).eps) + xi)
        inside ^= crosses; j = i
    return inside


def assign_regions(coords: np.ndarray):
    anchors = json.loads(ANCHOR_JSON.read_text(encoding="utf-8"))
    theta = np.deg2rad(ANGLE)
    rotated = np.column_stack((coords[:, 0] * np.cos(theta) - coords[:, 1] * np.sin(theta), coords[:, 0] * np.sin(theta) + coords[:, 1] * np.cos(theta)))
    epi = distance_to_polyline(rotated, anchors["epi_baseline"]); derm = distance_to_polyline(rotated, anchors["der_bottom"])
    depth = np.divide(epi, epi + derm, out=np.zeros_like(epi), where=(epi + derm) > 0)
    front = np.asarray(anchors["leading_edge"], dtype=float); dist = distance_to_polyline(rotated, anchors["leading_edge"])
    right = anchors.get("healed_direction") == "right_is_healed"
    positive = ((rotated[:, 0] > front[0, 0]) & right) | ((rotated[:, 0] < front[0, 0]) & (not right))
    horizontal = dist * np.where(positive, 1, -1) * 0.33
    max_abs = max(abs(horizontal.min()), abs(horizontal.max())) * 1.1
    model = np.column_stack((horizontal, depth))
    polygons = {"Un-epi": [[-max_abs, 0], [-600, 0], [-300, 1], [-max_abs, 1]], "Epi-Front": [[-600, 0], [400, 0], [700, 1], [-300, 1]], "Newly-epi": [[400, 0], [max_abs, 0], [max_abs, 1], [700, 1]]}
    region = np.full(len(coords), "", dtype=object)
    for name, polygon in polygons.items(): region[points_in_polygon(model, polygon)] = name
    display = rotated.copy(); display -= np.nanmedian(display, axis=0); display *= 0.33
    return display, region


def main() -> None:
    style(); path = stereo_file(SAMPLE)
    with h5py.File(path, "r") as handle:
        coords = handle["obsm"]["spatial"][:]; labels = decode_categorical(handle["obs"]["sub_labels"])
    xy, region = assign_regions(coords)
    if np.any(region == ""): raise RuntimeError("Unassigned region bins remain")

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.58, 0.42], left=0.03, right=0.99, top=0.84, bottom=0.15, hspace=0.18, wspace=0.08)
    overview = fig.add_subplot(gs[0, :]); plot_spatial_labels(overview, xy, labels, size=0.33); add_scale_bar(overview, 2000, "2 mm", 5.5)
    region_colors = {"Un-epi": "#D84A4A", "Epi-Front": "#3F6FA8", "Newly-epi": "#E0BC2C"}
    region_sub = {"Un-epi": "limited epithelial signal", "Epi-Front": "advancing epithelial edge", "Newly-epi": "stratified epithelial territory"}
    for name in REGIONS:
        mask = region == name; x, y = np.median(xy[mask], axis=0)
        overview.text(x, y, name, ha="center", va="center", fontsize=6, fontweight="bold", color=region_colors[name],
                      bbox={"boxstyle": "round,pad=.2", "facecolor": "white", "edgecolor": region_colors[name], "linewidth": .55, "alpha": .9})
    for ax, name in zip([fig.add_subplot(gs[1, i]) for i in range(3)], REGIONS):
        mask = region == name; plot_spatial_labels(ax, xy[mask], labels[mask], size=0.48)
        qx = np.quantile(xy[mask, 0], [.01, .99]); qy = np.quantile(xy[mask, 1], [.01, .99]); ax.set_xlim(qx); ax.set_ylim(qy)
        add_scale_bar(ax, 500, "500 μm", 5.5)
        ax.set_title(name + "\n" + region_sub[name], fontsize=6.2, fontweight="bold", color=region_colors[name], pad=1)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[s], markeredgewidth=0, markersize=3.4, label=s) for s in FOCUS_STATES]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, 0.012), ncol=7, frameon=False, handletextpad=.25, columnspacing=.8, fontsize=5.4)
    fig.text(0.006, 0.97, "D", fontsize=12, fontweight="bold", va="top")
    fig.text(0.036, 0.965, "Mode 2 | 19 dpb one-way epithelialization: spatial territories", fontsize=7.2, fontweight="bold", va="top")
    save_panel(fig, STEM, WIDTH_MM, HEIGHT_MM)
    with (ROOT / "source_data" / "Figure1D_sources.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle); writer.writerow(["element", "source", "subset", "processing"])
        writer.writerow(["overview and three territories", path, SAMPLE, "direct H5AD render; anchor-defined regions; equal aspect; vector scale bars"])
        writer.writerow(["region anchors", ANCHOR_JSON, "19dpb_p1", "same registered anchor geometry used for region assignment"])
    (ROOT / "QC" / f"{STEM}_content_proof.txt").write_text(
        "Figure 1D v02\n- Recomputed directly from the 19dpb_p1 H5AD and registered spatial anchors.\n"
        "- Overview and all crops use equal aspect; no PPT-derived spatial raster is used.\n"
        "- Rotation uses the exact config_wound formula with angle 140 degrees and no additional axis inversion.\n"
        "- Mode 2 spatial positions are labeled on the overview and repeated above matched regional crops.\n"
        "- Scale-bar style is shared with Figure 1C.\n", encoding="utf-8")


if __name__ == "__main__":
    main()
