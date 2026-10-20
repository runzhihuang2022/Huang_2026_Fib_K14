from __future__ import annotations

import json
from pathlib import Path

import h5py
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
VECTOR = OUT / "vector"
RASTER = OUT / "raster"
ASSETS = ROOT / "assets" / "linked_rasters"
QC = ROOT / "QC"

STEREO = Path(
    r"data/external/local_91ac4b6618ea"
    r"\时空组学\stereoseq\pdf_output\Wound_Healing_Annotation_Output_White"
)
STAIN = Path(
    r"data/external/local_3413b1de7846"
    r"\修回最终版文件\主图和附件\张伟染色结果\20260723染色数据给师兄"
    r"\SH  20230915  张伟 VIM+krt10 krt14 X30  DQ3046"
    r"\VIM KRT14X15\截图"
)

SECTIONS = [
    ("5 dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12 dpb DPTDI-1", "BW13_1_B1_DSDB_12dpb", 45),
    ("12 dpb SPTDI-1", "BW14_1_C1_SSDB_12dpb", 135),
    ("12 dpb DPTDI-2", "BW13_A3_DSDB_12dpb", 135),
    ("12 dpb SPTDI-2", "BW14_B3_SSDB_12dpb", 135),
    ("19 dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("19 dpb p1", "BW15D_1_D1_SDSDB_19dpb", 140),
    ("26 dpb p1", "BW81_C02846B6_SDSDB_26dpb_part1", 0),
    ("26 dpb p2", "BW81_C02846B6_SDSDB_26dpb_part2", 0),
    ("2 mph", "2mph_A03699G6.SCT", 45),
    ("Normal", "NS_C02847B1", 180),
]
TEMPORAL_SECTIONS = [
    ("Normal", "NS_C02847B1", 180),
    ("5 dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12 dpb SPTDI", "BW14_B3_SSDB_12dpb", 135),
    ("19 dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("2 mph", "2mph_A03699G6.SCT", 45),
]
KC_STATES = [
    "KC_Basal", "KC_Basal_Mig", "KC_Spinous",
    "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
]
KC_COLORS = {
    "KC_Basal": "#008B45",
    "KC_Basal_Mig": "#FF8C00",
    "KC_Spinous": "#00D62F",
    "KC_Spinous_Mig": "#00008B",
    "KC_Spinous_Mat": "#FF69B4",
    "KC_Granular": "#8B4513",
}

PAGE_W_MM, PAGE_H_MM = 196.0, 139.0
GAP_MM = 3.0
E_HEIGHT_MM = 71.0
FG_HEIGHT_MM = PAGE_H_MM - E_HEIGHT_MM - GAP_MM

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.titlesize": 8,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.55,
})


def rect_mm(x: float, top: float, width: float, height: float) -> list[float]:
    return [x / PAGE_W_MM, 1 - (top + height) / PAGE_H_MM, width / PAGE_W_MM, height / PAGE_H_MM]


def panel_header(fig: plt.Figure, rect: list[float], letter: str, title: str) -> None:
    x, y, w, h = rect
    fig.text(x, y + h - 0.010, letter, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(x + 0.025, y + h - 0.010, title, ha="left", va="top", fontsize=8, fontweight="bold")


def decode_strings(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values)
    if arr.dtype.kind == "S" or (arr.dtype.kind == "O" and arr.size and isinstance(arr.flat[0], (bytes, np.bytes_))):
        return np.array([v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in arr], dtype=object)
    return arr.astype(object)


def decode_categorical(node: h5py.Group) -> np.ndarray:
    categories = decode_strings(np.asarray(node["categories"]))
    codes = np.asarray(node["codes"], dtype=int)
    result = np.empty(codes.shape, dtype=object)
    result[:] = ""
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    return result


def read_obs_field(handle: h5py.File, field: str) -> np.ndarray:
    node = handle["obs"][field]
    return decode_categorical(node) if isinstance(node, h5py.Group) else decode_strings(np.asarray(node))


def rotate(coords: np.ndarray, angle: float) -> np.ndarray:
    theta = np.deg2rad(angle)
    centered = coords - np.nanmedian(coords, axis=0)
    matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    result = centered @ matrix.T
    if np.nanpercentile(result[:, 1], 95) < np.nanpercentile(result[:, 1], 5):
        result[:, 1] *= -1
    return result


def read_spatial(sample_id: str, fields: list[str]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    path = STEREO / f"{sample_id}_Annotated.h5ad"
    with h5py.File(path, "r") as handle:
        coords = np.asarray(handle["obsm"]["spatial"], dtype=float)
        obs = {field: read_obs_field(handle, field) for field in fields}
    return coords, obs


def add_scale_bar(
    ax: plt.Axes, length_units: float = 2000.0, label: str = "1 mm", position: str = "left"
) -> None:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x0 = xmin + (0.70 if position == "right" else 0.045) * (xmax - xmin)
    x0 = min(x0, xmax - length_units - 0.02 * (xmax - xmin))
    y0 = ymin + 0.075 * (ymax - ymin)
    ax.plot([x0, x0 + length_units], [y0, y0], color="black", lw=1.0, solid_capstyle="butt")
    ax.text(x0 + length_units / 2, y0 + 0.045 * (ymax - ymin), label, ha="center", va="bottom", fontsize=6)


def spatial_cache() -> list[tuple[str, str, float, np.ndarray, dict[str, np.ndarray]]]:
    result = []
    for label, sample_id, angle in SECTIONS:
        coords, obs = read_spatial(sample_id, ["nCount_Spatial", "nFeature_Spatial"])
        result.append((label, sample_id, angle, rotate(coords, angle), obs))
    return result


def panel_e(fig: plt.Figure, rect: list[float], cached: list[tuple[str, str, float, np.ndarray, dict[str, np.ndarray]]]) -> None:
    panel_header(fig, rect, "E", "Spatial QC across all 11 registered Stereo-seq sections")
    x, y, w, h = rect
    metrics = ("nCount_Spatial", "nFeature_Spatial")
    logs = {
        metric: np.concatenate([np.log10(np.asarray(item[4][metric], dtype=float) + 1) for item in cached])
        for metric in metrics
    }
    limits = {metric: np.nanquantile(values, [0.01, 0.99]) for metric, values in logs.items()}
    cmaps = {"nCount_Spatial": "viridis", "nFeature_Spatial": "plasma"}
    cols, rows = 4, 3
    left, right, top, bottom = 0.018, 0.012, 0.15, 0.14
    gx, gy = 0.012, 0.035
    tile_w = (w * (1 - left - right) - w * gx * (cols - 1)) / cols
    tile_h = (h * (1 - top - bottom) - h * gy * (rows - 1)) / rows
    for index, (label, _sample_id, _angle, coords, obs) in enumerate(cached):
        row, col = divmod(index, cols)
        tile_x = x + w * left + col * (tile_w + w * gx)
        tile_y = y + h * (1 - top) - (row + 1) * tile_h - row * h * gy
        fig.text(tile_x + tile_w / 2, tile_y + tile_h + 0.0015, label, ha="center", va="bottom", fontsize=7, fontweight="bold")
        map_gap = tile_h * 0.06
        map_h = (tile_h - map_gap) / 2
        for metric_index, metric in enumerate(metrics):
            map_y = tile_y + (1 - metric_index) * (map_h + map_gap)
            ax = fig.add_axes([tile_x, map_y, tile_w, map_h])
            values = np.log10(np.asarray(obs[metric], dtype=float) + 1)
            ax.scatter(
                coords[:, 0], coords[:, 1], c=values, s=0.17,
                cmap=cmaps[metric], vmin=limits[metric][0], vmax=limits[metric][1],
                linewidths=0, rasterized=True,
            )
            ax.set_aspect("equal")
            ax.axis("off")
            if col == 0:
                badge = "UMI" if metric_index == 0 else "Gene"
                ax.text(-0.02, 0.5, badge, transform=ax.transAxes, ha="right", va="center", fontsize=6, fontweight="bold")
            if metric_index == 0:
                add_scale_bar(ax, position="right")
    legend_specs = [
        ("nCount_Spatial", "log10(Number of UMI + 1)", 0.18),
        ("nFeature_Spatial", "log10(Number of Gene + 1)", 0.58),
    ]
    for metric, label, xpos in legend_specs:
        cax = fig.add_axes([x + xpos * w, y + 0.026 * h, 0.24 * w, 0.012 * h])
        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(*limits[metric]), cmap=cmaps[metric])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=6, length=1, pad=1)
        fig.text(x + (xpos + 0.12) * w, y + 0.066 * h, label, ha="center", va="bottom", fontsize=6)


def panel_f(fig: plt.Figure, rect: list[float]) -> list[str]:
    panel_header(fig, rect, "F", "Clinical VIM/KRT14 immunofluorescence")
    x, y, w, h = rect
    source_files = {
        "Merge": STAIN / "OCT-3-1 VIM KRT14_20.0x.jpg",
        "DAPI": STAIN / "1.jpg",
        "VIM": STAIN / "2.jpg",
        "KRT14": STAIN / "3.jpg",
    }
    for src in source_files.values():
        if not src.exists():
            raise FileNotFoundError(src)
    merge_ax = fig.add_axes([x + 0.025 * w, y + 0.37 * h, 0.95 * w, 0.46 * h])
    merge_ax.imshow(Image.open(source_files["Merge"]))
    merge_ax.axis("off")
    merge_ax.set_title("OCT-3-1 (20x), merged channels", fontsize=7, pad=1)
    labels = [("DAPI", "#2047D7"), ("VIM", "#00A651"), ("KRT14", "#E64B35")]
    thumb_y, thumb_h = y + 0.055 * h, 0.23 * h
    gap = 0.018 * w
    thumb_w = (0.95 * w - 2 * gap) / 3
    for index, (label, color) in enumerate(labels):
        ax = fig.add_axes([x + 0.025 * w + index * (thumb_w + gap), thumb_y, thumb_w, thumb_h])
        ax.imshow(Image.open(source_files[label]))
        ax.axis("off")
        ax.set_title(label, fontsize=7, color=color, fontweight="bold", pad=1)
    return [str(path) for path in source_files.values()]


def panel_g(fig: plt.Figure, rect: list[float]) -> None:
    panel_header(fig, rect, "G", "Keratinocyte states across healing")
    x, y, w, h = rect
    left, right, top, bottom, gap_x, gap_y = 0.025, 0.018, 0.16, 0.25, 0.022, 0.035
    map_w = (w * (1 - left - right) - w * gap_x * 2) / 3
    map_h = h * (1 - top - bottom - gap_y) / 2
    for index, (label, sample_id, angle) in enumerate(TEMPORAL_SECTIONS):
        coords, obs = read_spatial(sample_id, ["sub_labels"])
        coords = rotate(coords, angle)
        labels = obs["sub_labels"]
        row, col = divmod(index, 3)
        row_count = 3 if row == 0 else 2
        row_offset = 0 if row == 0 else (3 - row_count) * (map_w + w * gap_x) / 2
        map_x = x + w * left + row_offset + col * (map_w + w * gap_x)
        map_y = y + h * bottom + (1 - row) * (map_h + h * gap_y)
        ax = fig.add_axes([map_x, map_y, map_w, map_h])
        ax.scatter(coords[:, 0], coords[:, 1], c="#E9EDF2", s=0.09, linewidths=0, rasterized=True)
        for state in KC_STATES:
            mask = labels == state
            ax.scatter(coords[mask, 0], coords[mask, 1], c=KC_COLORS[state], s=1.6, linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(label, fontsize=8, fontweight="bold", pad=1)
        add_scale_bar(ax)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=KC_COLORS[state],
               markeredgecolor="none", markersize=3.5, label=state)
        for state in KC_STATES
    ]
    fig.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(x + w / 2, y + 0.012 * h),
        ncol=3, frameon=False, fontsize=7, handletextpad=0.25, columnspacing=0.75,
    )


def main() -> None:
    for directory in (VECTOR, RASTER, ASSETS, QC):
        directory.mkdir(parents=True, exist_ok=True)
    cached = spatial_cache()
    fig = plt.figure(figsize=(PAGE_W_MM / 25.4, PAGE_H_MM / 25.4), facecolor="none")
    fig.patch.set_alpha(0)
    e_rect = rect_mm(0, 0, PAGE_W_MM, E_HEIGHT_MM)
    f_width = 77.0
    g_width = PAGE_W_MM - f_width - GAP_MM
    f_rect = rect_mm(0, E_HEIGHT_MM + GAP_MM, f_width, FG_HEIGHT_MM)
    g_rect = rect_mm(f_width + GAP_MM, E_HEIGHT_MM + GAP_MM, g_width, FG_HEIGHT_MM)
    panel_e(fig, e_rect, cached)
    sources = panel_f(fig, f_rect)
    panel_g(fig, g_rect)
    stem = "FigureS1_EFG_revised_v07"
    svg = VECTOR / f"{stem}_editable.svg"
    pdf = VECTOR / f"{stem}.pdf"
    png = RASTER / f"{stem}_600dpi.png"
    fig.savefig(svg, dpi=600, transparent=True)
    fig.savefig(pdf, dpi=600, transparent=True)
    fig.savefig(png, dpi=600, transparent=True)
    plt.close(fig)
    report = {
        "canvas_mm": [PAGE_W_MM, PAGE_H_MM],
        "font": "Arial",
        "font_range_pt": [6, 12],
        "stereo_sections": len(SECTIONS),
        "temporal_sections": [item[0] for item in TEMPORAL_SECTIONS],
        "panel_F_source_files": sources,
        "panel_F_interpretation": "Representative staining only; no clinical stage inferred from filenames.",
        "background": "transparent composition on white final artboard; black editorial text",
        "rasterization_dpi_for_dense_spatial_points": 600,
    }
    (QC / "FigureS1_v07_EFG_content.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(svg)


if __name__ == "__main__":
    main()
