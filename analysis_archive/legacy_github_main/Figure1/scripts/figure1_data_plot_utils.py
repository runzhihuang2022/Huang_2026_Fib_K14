from __future__ import annotations

import os
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from figure_svg_utils import ROOT, ensure_dirs, render_svg


SCRNA_H5 = Path(os.environ.get("FIGURE1_HUMAN_SCRNA_H5", r"data/external/local_decd995edda5"))
STEREO_ROOT = Path(os.environ.get(
    "FIGURE1_HUMAN_STEREO_ROOT",
    r"data/external/local_be81937a46ed",
))
ANCHOR_JSON = Path(os.environ.get(
    "FIGURE1_19DPB_ANCHORS",
    r"data/external/local_f6f72cdda076",
))

CELL_ORDER = [
    "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
    "SAC_SG_Progenitor", "SAC_SG_Clear", "SAC_SG_Dark", "SAC_SG_Ductal", "SAC_HF_IRS", "SAC_HF_ORS", "SAC_HF_HFSC",
    "SAC_HF_DP_DS", "SAC_HF_Matrix", "Melanocyte", "Fib_Papi", "Fib_SFRP2", "Fib_Myo", "Fib_Inflama", "Fib_Fasci",
    "Fib_EN1", "Fib_Prolif", "Fib_K14", "Fib_CD45", "Endo_Capillary", "Endo_Arterial", "Endo_Venous", "Endo_Lymphatic",
    "Endo_Prolif", "Pericyte", "Schwann", "M1_Macrophage", "M2_Macrophage", "LAM", "Neutrophil", "Mast", "cDC1", "cDC2",
    "pDC", "LGHS", "B_Mature", "B_Plasma", "CD4_T", "CD8_T", "NK",
]

DISPLAY_NAMES = [
    "Basal", "Basal-Mig", "Basal-Prolif", "Spinous", "Spinous-Mig", "Spinous-Mat", "Granular", "SG-Prog", "SG-Clear", "SG-Dark",
    "SG-Ductal", "HF-IRS", "HF-ORS", "HF-HFSC", "HF-DP/DS", "HF-Matrix", "Melanocyte", "Fib-Papi", "Fib-SFRP2", "Fib-Myo",
    "Fib-Inflama", "Fib-Fasci", "Fib-EN1", "Fib-Prolif", "Fib-K14", "Fib-CD45", "Endo-Cap", "Endo-Art", "Endo-Ven", "Endo-Lymph",
    "Endo-Prolif", "Pericyte", "Schwann", "M1-Macro", "M2-Macro", "LAM", "Neutrophil", "Mast", "cDC1", "cDC2", "pDC", "LGHS",
    "B_Mature", "B_Plasma", "CD4_T", "CD8_T", "NK",
]

COLORS = dict(zip(CELL_ORDER, [
    "#008B45", "#FF8C00", "#FFD700", "#00C853", "#00008B", "#FF69B4", "#8B4513",
    "#0099CC", "#AFEEEE", "#0A2B00", "#CD853F", "#DA70D6", "#D8BFD8", "#20B2AA", "#B8860B", "#48D1CC", "#A0522D",
    "#4682B4", "#87CEFA", "#DDA0DD", "#BA55D3", "#9370DB", "#F0E68C", "#FFB6C1", "#E31A1C", "#9932CC",
    "#9ACD32", "#6B8E23", "#556B2F", "#8FBC8F", "#BDB76B", "#DAA520", "#D2B48C", "#1E90FF", "#4169E1", "#0000CD",
    "#483D8B", "#800080", "#6495ED", "#7B68EE", "#8470FF", "#00CED1", "#000080", "#4B0082", "#5F9EA0", "#6A5ACD", "#008080",
]))

FOCUS_STATES = ["KC_Basal", "KC_Spinous_Mig", "KC_Spinous", "KC_Spinous_Mat", "KC_Granular", "SAC_SG_Progenitor", "SAC_SG_Dark"]

CONFIG_WOUND = {
    "5dpb": {"sample_id": "BW32_A01597A3_SDSDB_5dpb", "angle": 135},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb", "angle": 45},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb", "angle": 135},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb", "angle": 135},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb", "angle": 135},
    "19dpb": {"sample_id": "BW15D_C6_SDSDB_19dpb", "angle": -225},
    "19dpb_p1": {"sample_id": "BW15D_1_D1_SDSDB_19dpb", "angle": 140},
    "26dpb_p1": {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1", "angle": 0},
    "26dpb_p2": {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2", "angle": 0},
    "2mph": {"sample_id": "2mph_A03699G6.SCT", "angle": 45},
    "Normal": {"sample_id": "NS_C02847B1", "angle": 180},
}


def style() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 6, "axes.titlesize": 7, "axes.labelsize": 6,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.55,
        "svg.fonttype": "none",
    })


def decode_categorical(group: h5py.Group) -> np.ndarray:
    categories = np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in group["categories"][:]], dtype=object)
    codes = np.asarray(group["codes"], dtype=int)
    out = np.full(codes.shape, "NA", dtype=object)
    valid = codes >= 0
    out[valid] = categories[codes[valid]]
    return out


def read_h5_panel(path: Path, obs_label: str = "sub_labels") -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        coords = np.asarray(handle["obsm"]["X_umap" if "X_umap" in handle["obsm"] else "spatial"], dtype=float)
        labels = decode_categorical(handle["obs"][obs_label])
    return coords, labels


def stereo_file(sample_id: str) -> Path:
    path = STEREO_ROOT / f"{sample_id}_Annotated.h5ad"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def rotate_xy(coords: np.ndarray, angle: float, pixel_um: float = 0.33) -> np.ndarray:
    xy = np.asarray(coords, dtype=float)
    center = np.nanmedian(xy, axis=0)
    xy = xy - center
    theta = np.deg2rad(angle)
    rotated = np.column_stack((xy[:, 0] * np.cos(theta) - xy[:, 1] * np.sin(theta), xy[:, 0] * np.sin(theta) + xy[:, 1] * np.cos(theta)))
    return rotated * pixel_um


def add_scale_bar(ax: plt.Axes, length_um: float, label: str, fontsize: float = 6) -> None:
    xmin, xmax = ax.get_xlim(); ymin, ymax = ax.get_ylim()
    x0 = xmin + 0.055 * (xmax - xmin)
    y0 = ymin + 0.075 * (ymax - ymin)
    ax.plot([x0, x0 + length_um], [y0, y0], color="black", lw=0.9, solid_capstyle="butt", zorder=20)
    ax.text(x0 + length_um / 2, y0 + 0.035 * (ymax - ymin), label, ha="center", va="bottom", fontsize=fontsize, fontweight="bold", zorder=20,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.4})


def plot_spatial_labels(ax: plt.Axes, xy_um: np.ndarray, labels: np.ndarray, focus_states: list[str] | None = None, size: float = 0.35) -> None:
    focus = focus_states or FOCUS_STATES
    ax.scatter(xy_um[:, 0], xy_um[:, 1], s=size * 0.55, c="#E3E3E3", linewidths=0, rasterized=True, zorder=1)
    for state in focus:
        mask = labels == state
        if mask.any():
            ax.scatter(xy_um[mask, 0], xy_um[mask, 1], s=size, c=COLORS[state], linewidths=0, rasterized=True, zorder=2)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def save_panel(fig: plt.Figure, stem: str, width_mm: float, height_mm: float) -> None:
    ensure_dirs()
    svg = ROOT / "outputs" / "panels" / f"{stem}.svg"
    fig.savefig(svg, format="svg", transparent=False)
    plt.close(fig)
    render_svg(svg, width_mm, height_mm, stem)
