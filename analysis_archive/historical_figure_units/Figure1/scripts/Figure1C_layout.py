from __future__ import annotations

import csv

import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from figure_svg_utils import ROOT
from figure1_data_plot_utils import COLORS, CONFIG_WOUND, FOCUS_STATES, add_scale_bar, decode_categorical, plot_spatial_labels, rotate_xy, save_panel, stereo_file, style


WIDTH_MM, HEIGHT_MM = 190, 39
STEM = "Figure1C_v02"
SECTIONS = [
    ("Mode 1 | 5 dpb SPTDI & DPTDI", "Superficial  ↔  Deep", CONFIG_WOUND["5dpb"]),
    ("Mode 2 | 19 dpb SPTDI & DPTDI", "Newly-epi  →  Epi-Front  →  Un-epi", CONFIG_WOUND["19dpb_p1"]),
    ("Mode 3 | 12 dpb SPTDI", "Healed  →  Focal unhealed", CONFIG_WOUND["12dpb_SPTDI2"]),
    ("Mode 4 | 2 mph SPTDI & DPTDI", "Scar / hyperplasia  →  Scarless healed", CONFIG_WOUND["2mph"]),
]


def read_section(sample_id: str, angle: float):
    path = stereo_file(sample_id)
    with h5py.File(path, "r") as handle:
        xy = handle["obsm"]["spatial"][:]
        labels = decode_categorical(handle["obs"]["sub_labels"])
    return path, rotate_xy(xy, angle), labels


def main() -> None:
    style()
    fig, axes = plt.subplots(1, 4, figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")
    fig.subplots_adjust(left=0.025, right=0.99, top=0.77, bottom=0.22, wspace=0.07)
    source_rows = []
    for ax, (title, spatial_note, config) in zip(axes, SECTIONS):
        sid, angle = config["sample_id"], config["angle"]
        path, xy, labels = read_section(sid, angle)
        plot_spatial_labels(ax, xy, labels, size=0.30)
        ax.set_title(title + "\n" + spatial_note, fontsize=6.4, fontweight="bold", pad=1.5)
        add_scale_bar(ax, 1000, "1 mm", fontsize=5.5)
        source_rows.append([title, sid, str(path), angle, len(labels), "direct H5AD render; equal aspect; 1 mm vector scale bar"])
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[s], markeredgewidth=0, markersize=3.5, label=s) for s in FOCUS_STATES]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, 0.015), ncol=7, frameon=False, handletextpad=0.25, columnspacing=0.8, fontsize=5.5)
    fig.text(0.005, 0.97, "C", fontsize=12, fontweight="bold", va="top")
    fig.text(0.035, 0.965, "Single-cell label projection onto representative Stereo-seq sections", fontsize=7, fontweight="bold", va="top")
    save_panel(fig, STEM, WIDTH_MM, HEIGHT_MM)
    with (ROOT / "source_data" / "Figure1C_sources.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle); writer.writerow(["display", "sample_id", "source", "rotation_deg", "n_bins", "processing"]); writer.writerows(source_rows)
    (ROOT / "QC" / f"{STEM}_content_proof.txt").write_text(
        "Figure 1C v02\n- All four maps were regenerated from registered per-section H5AD objects.\n"
        "- No PPT crop is used; every axes uses equal aspect and preserve-by-code sizing.\n"
        "- Rotation uses the exact config_wound formula x'=x cos(theta)-y sin(theta), y'=x sin(theta)+y cos(theta), with no additional axis inversion.\n"
        "- Scale bars share one vector style and 5.5 pt text.\n- Mode-derived spatial positions are stated above each section.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
