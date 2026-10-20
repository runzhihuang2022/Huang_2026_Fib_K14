from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import build_figure_s1_efg_v07 as base


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
VECTOR = OUT / "vector"
RASTER = OUT / "raster"
QC = ROOT / "QC"

base.STEREO = Path(os.environ.get("FIGS1_STEREO_ROOT", str(base.STEREO)))

PAGE_W_MM, PAGE_H_MM = 196.0, 100.0
E_HEIGHT_MM = 70.0
PANEL_GAP_MM = 3.0
F_HEIGHT_MM = PAGE_H_MM - E_HEIGHT_MM - PANEL_GAP_MM

LEFT_MM = 7.0
RIGHT_MM = 3.0
COL_GAP_MM = 2.0
COLS = 5
COL_W_MM = (PAGE_W_MM - LEFT_MM - RIGHT_MM - COL_GAP_MM * (COLS - 1)) / COLS
MAP_H_MM = 6.0

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


def rect_mm(left: float, top: float, width: float, height: float) -> list[float]:
    return [
        left / PAGE_W_MM,
        1.0 - (top + height) / PAGE_H_MM,
        width / PAGE_W_MM,
        height / PAGE_H_MM,
    ]


def add_scale_bar(ax: plt.Axes) -> None:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    length_units = 2000.0
    x0 = xmax - length_units - 0.05 * (xmax - xmin)
    y0 = ymin + 0.10 * (ymax - ymin)
    ax.plot([x0, x0 + length_units], [y0, y0], color="black", lw=0.9, solid_capstyle="butt")
    ax.text(
        x0 + length_units / 2,
        y0 + 0.05 * (ymax - ymin),
        "1 mm",
        ha="center",
        va="bottom",
        fontsize=6,
        color="black",
    )


def panel_header(fig: plt.Figure, letter: str, title: str, top_mm: float) -> None:
    fig.text(
        0.0,
        1.0 - (top_mm + 1.0) / PAGE_H_MM,
        letter,
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        9.0 / PAGE_W_MM,
        1.0 - (top_mm + 1.0) / PAGE_H_MM,
        title,
        ha="left",
        va="top",
        fontsize=8,
        fontweight="bold",
    )


def panel_e(
    fig: plt.Figure,
    cached: list[tuple[str, str, float, np.ndarray, dict[str, np.ndarray]]],
) -> None:
    panel_header(fig, "E", "Spatial QC across all 11 registered Stereo-seq sections", 0.0)
    metrics = ("nCount_Spatial", "nFeature_Spatial")
    logs = {
        metric: np.concatenate(
            [np.log10(np.asarray(item[4][metric], dtype=float) + 1) for item in cached]
        )
        for metric in metrics
    }
    limits = {metric: np.nanquantile(values, [0.01, 0.99]) for metric, values in logs.items()}
    cmaps = {"nCount_Spatial": "viridis", "nFeature_Spatial": "plasma"}

    row_start_mm = 8.5
    label_h_mm = 2.5
    within_gap_mm = 0.7
    row_gap_mm = 0.8
    row_h_mm = label_h_mm + 2 * MAP_H_MM + within_gap_mm

    for index, (label, _sample_id, _angle, coords, obs) in enumerate(cached):
        row, col = divmod(index, COLS)
        col_count = COLS if row < 2 else 1
        row_offset = 0.0 if col_count == COLS else (PAGE_W_MM - COL_W_MM) / 2 - LEFT_MM
        left_mm = LEFT_MM + row_offset + col * (COL_W_MM + COL_GAP_MM)
        top_mm = row_start_mm + row * (row_h_mm + row_gap_mm)
        fig.text(
            (left_mm + COL_W_MM / 2) / PAGE_W_MM,
            1.0 - top_mm / PAGE_H_MM,
            label,
            ha="center",
            va="top",
            fontsize=6,
            fontweight="bold",
        )
        for metric_index, metric in enumerate(metrics):
            map_top = top_mm + label_h_mm + metric_index * (MAP_H_MM + within_gap_mm)
            ax = fig.add_axes(rect_mm(left_mm, map_top, COL_W_MM, MAP_H_MM))
            values = np.log10(np.asarray(obs[metric], dtype=float) + 1)
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=0.18,
                cmap=cmaps[metric],
                vmin=limits[metric][0],
                vmax=limits[metric][1],
                linewidths=0,
                rasterized=True,
            )
            ax.set_aspect("equal")
            ax.axis("off")
            add_scale_bar(ax)
            if col == 0 and row < 2:
                fig.text(
                    6.2 / PAGE_W_MM,
                    1.0 - (map_top + MAP_H_MM / 2) / PAGE_H_MM,
                    "UMI" if metric_index == 0 else "Gene",
                    ha="right",
                    va="center",
                    fontsize=6,
                    fontweight="bold",
                )

    legend_top = 61.0
    legend_specs = [
        ("nCount_Spatial", "log10(Number of UMI + 1)", 24.0),
        ("nFeature_Spatial", "log10(Number of Gene + 1)", 112.0),
    ]
    for metric, label, left_mm in legend_specs:
        cax = fig.add_axes(rect_mm(left_mm, legend_top, 60.0, 1.8))
        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(*limits[metric]), cmap=cmaps[metric])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=6, length=1, pad=1)
        fig.text(
            (left_mm + 30.0) / PAGE_W_MM,
            1.0 - (legend_top - 1.2) / PAGE_H_MM,
            label,
            ha="center",
            va="bottom",
            fontsize=6,
        )


def panel_f(fig: plt.Figure) -> None:
    top_panel = E_HEIGHT_MM + PANEL_GAP_MM
    panel_header(fig, "F", "Keratinocyte states across healing", top_panel)
    map_top = top_panel + 8.0

    for index, (label, sample_id, angle) in enumerate(base.TEMPORAL_SECTIONS):
        coords, obs = base.read_spatial(sample_id, ["sub_labels"])
        coords = base.rotate(coords, angle)
        labels = obs["sub_labels"]
        left_mm = LEFT_MM + index * (COL_W_MM + COL_GAP_MM)
        ax = fig.add_axes(rect_mm(left_mm, map_top, COL_W_MM, MAP_H_MM))
        ax.scatter(coords[:, 0], coords[:, 1], c="#E9EDF2", s=0.10, linewidths=0, rasterized=True)
        for state in base.KC_STATES:
            mask = labels == state
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=base.KC_COLORS[state],
                s=1.6,
                linewidths=0,
                rasterized=True,
            )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(label, fontsize=7, fontweight="bold", pad=1)
        add_scale_bar(ax)

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=base.KC_COLORS[state],
            markeredgecolor="none",
            markersize=3.5,
            label=state,
        )
        for state in base.KC_STATES
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.008),
        ncol=6,
        frameon=False,
        fontsize=6,
        handletextpad=0.25,
        columnspacing=0.8,
    )


def main() -> None:
    for directory in (VECTOR, RASTER, QC):
        directory.mkdir(parents=True, exist_ok=True)
    cached = base.spatial_cache()
    fig = plt.figure(figsize=(PAGE_W_MM / 25.4, PAGE_H_MM / 25.4), facecolor="none")
    fig.patch.set_alpha(0)
    panel_e(fig, cached)
    panel_f(fig)

    stem = "FigureS1_EF_revised_v08"
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
        "panel_E_sections": len(base.SECTIONS),
        "panel_F_temporal_sections": [item[0] for item in base.TEMPORAL_SECTIONS],
        "common_column_width_mm": COL_W_MM,
        "common_map_height_mm": MAP_H_MM,
        "common_scale_bar": "1 mm on every E and F map",
        "removed_panel": "former clinical VIM/KRT14 immunofluorescence panel F",
        "renumbering": "former panel G is panel F",
        "rasterization_dpi_for_dense_spatial_points": 600,
    }
    (QC / "FigureS1_v08_EF_content.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    print(svg)


if __name__ == "__main__":
    main()
