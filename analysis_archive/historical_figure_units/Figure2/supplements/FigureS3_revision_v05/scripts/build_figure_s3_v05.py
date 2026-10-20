from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import kruskal

from scripts.build_figure_s3_v04 import (
    CELL_COLORS,
    DISPLAY_LABELS,
    TIME_ORDER,
    add_letter,
    clean,
    configure_style,
    panel_a,
)


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
VEC = ROOT / "outputs" / "vector"
RAS = ROOT / "outputs" / "raster"
QC = ROOT / "QC"
for directory in (VEC, RAS, QC):
    directory.mkdir(parents=True, exist_ok=True)


# Exact palette used by the source 26dpb_p2 co-occurrence panel.
SOURCE_COLORS = {
    "KC_Basal": "#008941",
    "KC_Spinous": "#00FF00",
    "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17",
    "KC_Spinous_Mig": "#1A0099",
    "KC_Basal_Mig": "#FF3030",
    "KC_Basal_Prolif": "#FFD700",
    "SAC_SG_Progenitor": "#0099CC",
    "SAC_SG_Clear": "#B7E4F9",
    "SAC_SG_Dark": "#193006",
    "SAC_SG_Ductal": "#8F7700",
    "SAC_HF_IRS": "#CC9900",
    "SAC_HF_ORS": "#A30059",
    "SAC_HF_HFSC": "#00F5FF",
    "SAC_HF_DP_DS": "#CDBE70",
    "SAC_HF_Matrix": "#7FFFD4",
    "Melanocyte": "#FF4500",
    "Fib_Papi": "#3DA873",
    "Fib_SFRP2": "#CCFF00",
    "Fib_Myo": "#B79762",
    "Fib_Inflama": "#800000",
    "Fib_Fasci": "#9370DB",
    "Fib_EN1": "#BDB76B",
    "Fib_CD45": "#DA70D6",
    "Fib_Prolif": "#EED2EE",
    "Fib_K14": "#4682B4",
    "Endo_Capillary": "#CF4E9C",
    "Endo_Arterial": "#FF00FF",
    "Endo_Venous": "#4169E1",
    "Endo_Lymphatic": "#ADFF2F",
    "Endo_Prolif": "#F08080",
    "Pericyte": "#FFFF00",
    "Schwann": "#7CFC00",
    "M1_Macrophage": "#FF4040",
    "M2_Macrophage": "#54FF9F",
    "LAM": "#8B4513",
    "Neutrophil": "#FF7F50",
    "Mast": "#FF1493",
    "cDC1": "#9B30FF",
    "cDC2": "#EE9A00",
    "pDC": "#D1EEEE",
    "LGHS": "#FFB5C5",
    "B_Mature": "#0000CD",
    "B_Plasma": "#FAF0E6",
    "CD4_T": "#00BFFF",
    "CD8_T": "#1E90FF",
    "NK": "#00E5EE",
}


def time_group_statistics(section_stats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for state in ["KC_Spinous_Mig", "SAC_SG_Progenitor", "Fib_K14"]:
        subset = section_stats[section_stats["cell_state"].eq(state)]
        groups = [
            subset.loc[subset["timepoint"].eq(timepoint), "cooccurrence_score"].to_numpy()
            for timepoint in TIME_ORDER
        ]
        statistic, p_value = kruskal(*groups)
        rows.append(
            {
                "cell_state": state,
                "comparison": "co-occurrence score across six healing time groups",
                "test": "Kruskal-Wallis",
                "n_sections": len(subset),
                "statistic": statistic,
                "p_value": p_value,
                "displayed_in_figure": bool(p_value < 0.05),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(SRC / "FigureS3_v05_time_group_statistics.csv", index=False)
    return result


def plot_time_series(
    ax: mpl.axes.Axes,
    section_stats: pd.DataFrame,
    states: list[str],
    title: str,
    show_legend: bool,
) -> None:
    x = np.arange(len(TIME_ORDER))
    for state in states:
        raw = section_stats[section_stats["cell_state"].eq(state)]
        means = raw.groupby("timepoint")["cooccurrence_score"].mean().reindex(TIME_ORDER)
        color = CELL_COLORS[state]
        ax.plot(
            x,
            means.to_numpy(float),
            color=color,
            lw=1.5,
            marker="o",
            markersize=4,
            label=DISPLAY_LABELS[state],
            zorder=3,
        )
        for x_index, timepoint in enumerate(TIME_ORDER):
            values = raw.loc[
                raw["timepoint"].eq(timepoint), "cooccurrence_score"
            ].to_numpy(float)
            jitter = np.linspace(-0.055, 0.055, len(values)) if len(values) > 1 else np.zeros(len(values))
            ax.scatter(
                x_index + jitter,
                values,
                s=12,
                color=color,
                edgecolor="white",
                linewidth=0.3,
                zorder=4,
            )
    ax.axhline(1, color="#666666", lw=0.7, ls="--", zorder=1)
    ax.set_xticks(x, TIME_ORDER, rotation=25, ha="right")
    ax.set_ylabel("Fib_K14-conditioned\nco-occurrence score")
    ax.set_title(title, fontsize=8, fontweight="bold", pad=3)
    if show_legend:
        ax.legend(frameon=False, loc="upper left", ncol=2, fontsize=6.5)
    clean(ax)


def panel_b(fig: mpl.figure.Figure, section_stats: pd.DataFrame) -> None:
    add_letter(fig, 0.515, 0.975, "B")

    cross_ax = fig.add_axes([0.585, 0.735, 0.375, 0.185])
    plot_time_series(
        cross_ax,
        section_stats,
        ["KC_Spinous_Mig", "SAC_SG_Progenitor"],
        "Temporal Fib_K14–cell co-occurrence",
        show_legend=True,
    )

    self_ax = fig.add_axes([0.585, 0.555, 0.375, 0.115])
    plot_time_series(
        self_ax,
        section_stats,
        ["Fib_K14"],
        "Fib_K14 self-spatial aggregation",
        show_legend=False,
    )
    self_ax.set_ylabel("Self-co-occurrence\nscore")


def panel_c(
    fig: mpl.figure.Figure,
    curves: pd.DataFrame,
    closest: pd.DataFrame,
) -> None:
    """Reproduce the author-specified 26dpb_p2 source panel without redesign."""
    add_letter(fig, 0.025, 0.505, "C")

    closest_without_self = closest[~closest["cell_state"].eq("Fib_K14")]
    top10 = closest_without_self.head(10).copy()
    top10_states = top10["cell_state"].tolist()
    available_states = set(curves["cell_state"])
    all_states = [state for state in SOURCE_COLORS if state in available_states]
    all_states.extend(
        state
        for state in curves["cell_state"].drop_duplicates()
        if state not in all_states
    )

    curve_ax = fig.add_axes([0.075, 0.145, 0.445, 0.295])
    handles: list[Line2D] = []
    for state in all_states:
        state_data = curves[curves["cell_state"].eq(state)].sort_values(
            "distance_coordinate_units"
        )
        color = SOURCE_COLORS.get(state, "#888888")
        if state == "Fib_K14":
            linewidth, alpha, zorder = 2.6, 1.0, 15
        elif state in top10_states:
            linewidth, alpha, zorder = 1.9, 1.0, 10
        else:
            linewidth, alpha, zorder = 0.6, 0.15, 1
        curve_ax.plot(
            state_data["distance_coordinate_units"],
            state_data["cooccurrence_probability_ratio"],
            color=color,
            lw=linewidth,
            alpha=alpha,
            zorder=zorder,
        )
        handles.append(Line2D([0], [0], color=color, lw=max(linewidth, 1.0), label=state))

    curve_ax.set_xlim(0, 10000)
    curve_ax.set_ylim(bottom=0)
    curve_ax.set_xlabel("Distance (µm)")
    curve_ax.set_ylabel("Co-occurrence Probability")
    curve_ax.set_title(
        "[26dpb_p2] Co-occurrence Probability Conditioned on Fib_K14",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    curve_ax.tick_params(labelsize=7)
    for spine in curve_ax.spines.values():
        spine.set_linewidth(0.8)
    curve_ax.spines["top"].set_visible(True)
    curve_ax.spines["right"].set_visible(True)

    bar_ax = fig.add_axes([0.615, 0.145, 0.345, 0.295])
    plot_df = top10.iloc[::-1]
    bars = bar_ax.barh(
        plot_df["cell_state"],
        plot_df["closest_distance_score"],
        color=[SOURCE_COLORS.get(state, "#888888") for state in plot_df["cell_state"]],
        edgecolor="none",
        height=0.60,
    )
    maximum = float(plot_df["closest_distance_score"].max())
    for bar in bars:
        width = bar.get_width()
        bar_ax.text(
            width + maximum * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{width:.3f}",
            va="center",
            ha="left",
            fontsize=7,
            fontweight="bold",
        )
    bar_ax.set_xlim(0, maximum * 1.15)
    bar_ax.set_xlabel("Co-occurrence Score at Closest Distance")
    bar_ax.set_title(
        "[26dpb_p2] Top 10 Co-occurring Cells (Score)",
        fontsize=10,
        fontweight="bold",
        pad=8,
    )
    bar_ax.tick_params(labelsize=7)
    for spine in bar_ax.spines.values():
        spine.set_linewidth(0.8)
    bar_ax.spines["top"].set_visible(True)
    bar_ax.spines["right"].set_visible(True)

    legend_ax = fig.add_axes([0.095, 0.045, 0.81, 0.075])
    legend_ax.axis("off")
    legend = legend_ax.legend(
        handles=handles,
        labels=all_states,
        loc="center",
        ncol=7,
        frameon=False,
        fontsize=6,
        handletextpad=0.4,
        columnspacing=1.0,
    )
    emphasized = set(top10_states + ["Fib_K14"])
    for text_item in legend.get_texts():
        if text_item.get_text() in emphasized:
            text_item.set_color("black")
            text_item.set_fontweight("bold")
        else:
            text_item.set_color("#BBBBBB")


def build() -> None:
    configure_style()
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )
    top10 = pd.read_csv(SRC / "FigureS3_FibK14_KC_SAC_timepoint_TOP10.csv")
    section_stats = pd.read_csv(SRC / "FigureS3_v04_prespecified_section_statistics.csv")
    curves = pd.read_csv(SRC / "FigureS3_v04_26dpb_p2_distance_curves.csv")
    closest = pd.read_csv(SRC / "FigureS3_v04_26dpb_p2_closest_scores.csv")
    group_stats = time_group_statistics(section_stats)

    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4), facecolor="white")
    panel_a(fig, top10)
    panel_b(fig, section_stats)
    panel_c(fig, curves, closest)

    svg = VEC / "FigureS3_revised_v05_editable.svg"
    pdf = VEC / "FigureS3_revised_v05.pdf"
    png = RAS / "FigureS3_revised_v05_600dpi.png"
    tiff = RAS / "FigureS3_revised_v05_600dpi.tiff"
    fig.savefig(svg)
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    fig.savefig(tiff, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    with (QC / "FigureS3_v05_output_sha256.txt").open("w", encoding="utf-8") as handle:
        for path in (svg, pdf, png, tiff):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")
    (QC / "FigureS3_v05_group_statistics.txt").write_text(
        group_stats.to_string(index=False), encoding="utf-8"
    )


if __name__ == "__main__":
    build()
