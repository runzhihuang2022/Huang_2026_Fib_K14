from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(
    r"data/external/local_3413b1de7846"
    r"\20260813找曾老师指导\免疫荧光结果统计\Figure 3G overlap"
)
SOURCE = BASE / "09_BW9_single_slide_spatial_statistics"
OUTPUT = BASE / "10_BW9_square_triptych_A4"

A4_WIDTH_MM = 210.0
PANEL_MM = A4_WIDTH_MM / 5.0
FIGURE_WIDTH_MM = PANEL_MM * 3.0
PT_PER_INCH = 72.0
MM_PER_INCH = 25.4


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.titleweight": "bold",
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def format_p_value(p: float) -> str:
    return "P < 0.001" if p < 0.001 else f"P = {p:.3f}"


def plot_panel(ax: mpl.axes.Axes, panel: str, by_radius: pd.DataFrame, null: pd.DataFrame) -> None:
    radii = by_radius["radius_um"].to_numpy(float)
    observed = by_radius["symmetric_observed_pct"].to_numpy(float)
    null_mean = by_radius["symmetric_null_mean_pct"].to_numpy(float)
    null_low = by_radius["symmetric_null_95ci_low_pct"].to_numpy(float)
    null_high = by_radius["symmetric_null_95ci_high_pct"].to_numpy(float)

    if panel == "a":
        values = null["symmetric_proximity_10um_pct"].to_numpy(float)
        ax.hist(values, bins=35, color="#BDBDBD", edgecolor="white", linewidth=0.25)
        ax.axvline(observed[0], color="#6A3D9A", linewidth=1.3)
        ax.axvline(values.mean(), color="#272727", linestyle="--", linewidth=0.9)
        ax.text(0.98, 0.95, f"Observed {observed[0]:.2f}%", transform=ax.transAxes, ha="right", va="top", color="#6A3D9A", fontsize=8)
        ax.text(0.98, 0.86, f"Null {values.mean():.2f}%", transform=ax.transAxes, ha="right", va="top", color="#272727", fontsize=8)
        p = float(by_radius.iloc[0]["permutation_p_one_sided"])
        ax.text(0.98, 0.77, format_p_value(p), transform=ax.transAxes, ha="right", va="top", fontsize=8)
        ax.set_xlabel("Proximity within 10 µm (%)")
        ax.set_ylabel("Permutations")
        ax.set_title("10-µm randomization")
    elif panel == "b":
        ax.fill_between(radii, null_low, null_high, color="#BDBDBD", alpha=0.45, linewidth=0)
        ax.plot(radii, null_mean, "--o", color="#272727", linewidth=1.0, markersize=2.7, label="Null")
        ax.plot(radii, observed, "-o", color="#6A3D9A", linewidth=1.2, markersize=2.7, label="Observed")
        ax.set_xscale("log")
        ax.set_xticks(radii, [str(int(x)) for x in radii])
        ax.set_xlabel("Radius (µm)")
        ax.set_ylabel("Symmetric proximity (%)")
        ax.set_title("Radius sensitivity")
        ax.legend(loc="upper left", handlelength=1.8, borderaxespad=0.2)
    elif panel == "c":
        effect = observed - null_mean
        ax.axhline(0, color="#272727", linewidth=0.7)
        ax.plot(radii, effect, "-o", color="#1B9E77", linewidth=1.2, markersize=2.7)
        ax.set_xscale("log")
        ax.set_xticks(radii, [str(int(x)) for x in radii])
        ax.set_xlabel("Radius (µm)")
        ax.set_ylabel("Observed − null (pp)")
        ax.set_title("Spatial enrichment")
    else:
        raise ValueError(panel)

    ax.tick_params(width=0.8, length=3)


def make_panel(panel: str, by_radius: pd.DataFrame, null: pd.DataFrame) -> mpl.figure.Figure:
    size = PANEL_MM / MM_PER_INCH
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes([0.22, 0.22, 0.74, 0.68])
    plot_panel(ax, panel, by_radius, null)
    return fig


def save_all(fig: mpl.figure.Figure, stem: Path, dpi: int = 600) -> None:
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=300)
    # Uncompressed TIFF is lossless and avoids optional imagecodecs/LZW runtime failures.
    fig.savefig(stem.with_suffix(".tiff"), dpi=dpi)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    set_style()
    by_radius = pd.read_csv(SOURCE / "BW9_spatial_randomization_by_radius.csv")
    null = pd.read_csv(SOURCE / "BW9_null_distribution_source_data.csv")

    for panel in ("a", "b", "c"):
        fig = make_panel(panel, by_radius, null)
        save_all(fig, OUTPUT / f"BW9_panel_{panel}_42mm_square")
        plt.close(fig)

    width = FIGURE_WIDTH_MM / MM_PER_INCH
    height = PANEL_MM / MM_PER_INCH
    fig = plt.figure(figsize=(width, height))
    for index, panel in enumerate(("a", "b", "c")):
        x0 = index / 3.0
        # Each allocated panel cell is exactly 42 x 42 mm.
        ax = fig.add_axes([x0 + 0.0733, 0.22, 0.2467, 0.68])
        plot_panel(ax, panel, by_radius, null)
    save_all(fig, OUTPUT / "BW9_three_square_panels_horizontal_A4_3of5")
    plt.close(fig)

    audit = f"""BW9 square triptych layout audit
A4 width: {A4_WIDTH_MM:.1f} mm
Each panel allocation: {PANEL_MM:.1f} x {PANEL_MM:.1f} mm (A4 width / 5)
Combined figure: {FIGURE_WIDTH_MM:.1f} x {PANEL_MM:.1f} mm
Arrangement: 1 row x 3 columns
Font: Arial with Helvetica/DejaVu Sans fallback
Font sizes: tick/axis/legend 8 pt; titles 9 pt bold
Panel labels: none
P-value display: P < 0.001, otherwise at most 3 decimals
Italics: none; no conventional biological nomenclature in these plot labels
Primary editable export: SVG; PDF uses TrueType text
Raster exports: PNG 300 dpi; TIFF 600 dpi uncompressed (lossless)
No image crop, microscopy ZoomIn selection, or statistical recomputation was performed.
"""
    (OUTPUT / "BW9_square_triptych_layout_audit.txt").write_text(audit, encoding="utf-8")


if __name__ == "__main__":
    main()
