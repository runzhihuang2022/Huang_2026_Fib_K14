from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
VEC = ROOT / "outputs" / "vector"
RAS = ROOT / "outputs" / "raster"
QC = ROOT / "QC"
for directory in (VEC, RAS, QC):
    directory.mkdir(parents=True, exist_ok=True)

TIME_ORDER = ["Normal", "5 dpb", "12 dpb", "19 dpb", "26 dpb", "2 mph"]
CELL_ORDER = [
    "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
    "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
    "SAC_SG_Progenitor", "SAC_SG_Clear", "SAC_SG_Dark", "SAC_SG_Ductal",
    "SAC_HF_IRS", "SAC_HF_ORS", "SAC_HF_HFSC", "SAC_HF_DP_DS",
    "SAC_HF_Matrix",
]
CELL_COLORS = {
    "Fib_K14": "#E31A1C",
    "KC_Basal": "#008941", "KC_Basal_Mig": "#FF8C00",
    "KC_Basal_Prolif": "#FFD700", "KC_Spinous": "#00D600",
    "KC_Spinous_Mig": "#1A0099", "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17", "SAC_SG_Progenitor": "#0099CC",
    "SAC_SG_Clear": "#A8E6E6", "SAC_SG_Dark": "#193006",
    "SAC_SG_Ductal": "#CD853F", "SAC_HF_IRS": "#CC9900",
    "SAC_HF_ORS": "#A30059", "SAC_HF_HFSC": "#20B2AA",
    "SAC_HF_DP_DS": "#CDBE70", "SAC_HF_Matrix": "#48C9C9",
}
DISPLAY_STATES = ["Fib_K14", "KC_Spinous_Mig", "SAC_SG_Progenitor"]
DISPLAY_LABELS = {
    "Fib_K14": "Fib_K14 (self)",
    "KC_Spinous_Mig": "KC_Spinous_Mig",
    "SAC_SG_Progenitor": "SAC_SG_Progenitor",
}


def add_letter(fig: mpl.figure.Figure, x: float, y: float, letter: str) -> None:
    fig.text(x, y, letter, fontsize=14, fontweight="bold", va="top", ha="left")


def clean(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=2.5)


def star(q_value: float) -> str:
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return "ns"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def vector_horizontal_colorbar(
    fig: mpl.figure.Figure,
    position: list[float],
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
    maximum: float,
) -> None:
    ax = fig.add_axes(position)
    edges = np.linspace(0, maximum, 41)
    for low, high in zip(edges[:-1], edges[1:]):
        midpoint = (low + high) / 2
        ax.add_patch(
            mpl.patches.Rectangle(
                (low, 0),
                high - low,
                1,
                facecolor=cmap(norm(midpoint)),
                edgecolor="none",
            )
        )
    ax.set_xlim(0, maximum)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([0, 1, 2])
    ax.tick_params(axis="x", labelsize=6, length=2, pad=1)
    ax.set_xlabel("Mean co-occurrence score", fontsize=7, labelpad=2)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)


def panel_a(fig: mpl.figure.Figure, top10: pd.DataFrame) -> None:
    add_letter(fig, 0.025, 0.975, "A")
    ax = fig.add_axes([0.155, 0.565, 0.315, 0.355])
    score = pd.DataFrame(np.nan, index=CELL_ORDER, columns=TIME_ORDER)
    rank = pd.DataFrame(np.nan, index=CELL_ORDER, columns=TIME_ORDER)
    for _, row in top10.iterrows():
        cell = row["cell_state"]
        timepoint = row["timepoint"]
        if cell in score.index and timepoint in score.columns:
            score.loc[cell, timepoint] = float(row["mean_score"])
            rank.loc[cell, timepoint] = int(row["rank_within_time"])
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#EFEFEF")
    maximum = max(2.0, float(np.nanmax(score.to_numpy(float))))
    mesh = ax.pcolormesh(
        np.arange(len(TIME_ORDER) + 1),
        np.arange(len(CELL_ORDER) + 1),
        score.to_numpy(float),
        cmap=cmap,
        norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=maximum),
        shading="flat",
        edgecolors="white",
        linewidth=0.35,
    )
    ax.set_xlim(0, len(TIME_ORDER))
    ax.set_ylim(len(CELL_ORDER), 0)
    ax.set_xticks(np.arange(len(TIME_ORDER)) + 0.5, TIME_ORDER)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(CELL_ORDER)) + 0.5, CELL_ORDER)
    ax.tick_params(length=0, pad=2)
    ax.set_title(
        "Time-resolved TOP10 Fib_K14–KC/SAC\nco-occurrence scores",
        fontweight="bold",
        pad=6,
    )
    for tick in ax.get_yticklabels():
        tick.set_color(CELL_COLORS.get(tick.get_text(), "#222222"))
    for row_index in range(len(CELL_ORDER)):
        for column_index in range(len(TIME_ORDER)):
            value = score.iat[row_index, column_index]
            rank_value = rank.iat[row_index, column_index]
            if np.isfinite(value):
                text_color = "white" if value > max(1.45, maximum * 0.56) else "black"
                ax.text(
                    column_index + 0.5,
                    row_index + 0.5,
                    f"{value:.2f}\n#{int(rank_value)}",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color=text_color,
                    linespacing=0.85,
                )
    vector_horizontal_colorbar(
        fig,
        [0.155, 0.535, 0.315, 0.012],
        cmap,
        mesh.norm,
        maximum,
    )


def panel_b(
    fig: mpl.figure.Figure,
    section_stats: pd.DataFrame,
    time_stats: pd.DataFrame,
) -> None:
    add_letter(fig, 0.515, 0.975, "B")
    ax = fig.add_axes([0.585, 0.705, 0.375, 0.215])
    x = np.arange(len(TIME_ORDER))
    for state in DISPLAY_STATES:
        summary = (
            time_stats[time_stats["cell_state"] == state]
            .set_index("timepoint")
            .reindex(TIME_ORDER)
        )
        color = CELL_COLORS[state]
        ax.plot(
            x,
            summary["mean_score"],
            marker="o",
            markersize=4.2,
            lw=1.6,
            color=color,
            label=DISPLAY_LABELS[state],
        )
        raw = section_stats[section_stats["cell_state"] == state]
        for x_index, timepoint in enumerate(TIME_ORDER):
            values = raw[raw["timepoint"] == timepoint]["cooccurrence_score"].to_numpy()
            jitter = (
                np.linspace(-0.06, 0.06, len(values))
                if len(values) > 1 else np.zeros(len(values))
            )
            ax.scatter(
                x_index + jitter,
                values,
                s=14,
                color=color,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
    ax.axhline(1, color="#555555", lw=0.75, ls="--")
    ax.set_xticks(x, TIME_ORDER, rotation=25, ha="right")
    ax.set_ylabel("Fib_K14-conditioned\nco-occurrence score")
    ax.set_title(
        "Temporal Fib_K14-centered co-occurrence",
        fontweight="bold",
        pad=5,
    )
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=1)
    clean(ax)

    q_ax = fig.add_axes([0.585, 0.565, 0.375, 0.095])
    display_time = time_stats[time_stats["cell_state"].isin(DISPLAY_STATES)].copy()
    q_matrix = (
        display_time.pivot(
            index="cell_state",
            columns="timepoint",
            values="fdr_across_102_timepoint_state_tests",
        )
        .reindex(index=DISPLAY_STATES, columns=TIME_ORDER)
    )
    for row_index, state in enumerate(DISPLAY_STATES):
        for column_index, timepoint in enumerate(TIME_ORDER):
            q_value = float(q_matrix.loc[state, timepoint])
            q_ax.add_patch(
                mpl.patches.Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1,
                    1,
                    facecolor="#FADBD8" if q_value < 0.05 else "#F3F3F3",
                    edgecolor="white",
                    linewidth=0.5,
                )
            )
    q_ax.set_xlim(-0.5, len(TIME_ORDER) - 0.5)
    q_ax.set_ylim(len(DISPLAY_STATES) - 0.5, -0.5)
    q_ax.set_xticks(np.arange(len(TIME_ORDER)), TIME_ORDER, fontsize=6)
    q_ax.set_yticks(
        np.arange(len(DISPLAY_STATES)),
        [DISPLAY_LABELS[state] for state in DISPLAY_STATES],
        fontsize=6,
    )
    q_ax.tick_params(length=0, pad=1)
    for row_index, state in enumerate(DISPLAY_STATES):
        q_ax.get_yticklabels()[row_index].set_color(CELL_COLORS[state])
        for column_index, timepoint in enumerate(TIME_ORDER):
            q_value = float(q_matrix.loc[state, timepoint])
            label = f"q={q_value:.3f}" if q_value < 0.1 else f"q={q_value:.2f}"
            if q_value < 0.05:
                label += f" {star(q_value)}"
            q_ax.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=6,
                fontweight="bold" if q_value < 0.05 else "normal",
            )
    q_ax.set_title(
        "Time-point spatial-label permutation statistics",
        fontsize=7,
        fontweight="bold",
        pad=3,
    )
    for spine in q_ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#777777")
    fig.text(
        0.585,
        0.525,
        "999 permutations; BH correction across 17 states × 6 time points.",
        fontsize=6,
        ha="left",
    )


def panel_c(
    fig: mpl.figure.Figure,
    curves: pd.DataFrame,
    closest: pd.DataFrame,
    section_stats: pd.DataFrame,
) -> None:
    add_letter(fig, 0.025, 0.485, "C")
    fig.text(
        0.50,
        0.478,
        "26 dpb_p2: Fib_K14-centered spatial co-occurrence",
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )
    top_states = closest.head(10)["cell_state"].tolist()
    curve_ax = fig.add_axes([0.075, 0.165, 0.40, 0.27])
    for state, state_data in curves.groupby("cell_state", sort=False):
        state_data = state_data.sort_values("distance_um")
        if state == "Fib_K14":
            color, width, alpha, zorder = CELL_COLORS[state], 2.4, 1.0, 5
        elif state == "KC_Spinous_Mig":
            color, width, alpha, zorder = CELL_COLORS[state], 2.4, 1.0, 6
        elif state in top_states:
            color, width, alpha, zorder = "#9E9E9E", 1.0, 0.55, 3
        else:
            color, width, alpha, zorder = "#D9D9D9", 0.55, 0.22, 1
        curve_ax.plot(
            state_data["distance_um"],
            state_data["cooccurrence_probability_ratio"],
            color=color,
            lw=width,
            alpha=alpha,
            zorder=zorder,
        )
    curve_ax.axhline(1, color="#555555", lw=0.75, ls="--")
    curve_ax.set_xlim(left=0)
    curve_ax.set_ylim(bottom=0)
    curve_ax.set_xlabel("Distance (µm)")
    curve_ax.set_ylabel("Co-occurrence probability ratio")
    curve_ax.set_title(
        "Distance-dependent co-occurrence conditioned on Fib_K14",
        fontweight="bold",
        pad=5,
    )
    curve_ax.legend(
        handles=[
            Line2D([0], [0], color=CELL_COLORS["Fib_K14"], lw=2.4, label="Fib_K14 (self)"),
            Line2D([0], [0], color=CELL_COLORS["KC_Spinous_Mig"], lw=2.4, label="KC_Spinous_Mig"),
            Line2D([0], [0], color="#9E9E9E", lw=1.2, label="Other top co-occurring states"),
        ],
        frameon=False,
        loc="upper right",
        fontsize=7,
    )
    clean(curve_ax)

    bar_ax = fig.add_axes([0.635, 0.165, 0.31, 0.27])
    top10 = closest.head(10).sort_values("closest_distance_score", ascending=True)
    bar_colors = [
        CELL_COLORS.get(state, "#BDBDBD")
        if state in ("Fib_K14", "KC_Spinous_Mig")
        else "#BDBDBD"
        for state in top10["cell_state"]
    ]
    bars = bar_ax.barh(
        top10["cell_state"],
        top10["closest_distance_score"],
        color=bar_colors,
        edgecolor="none",
        height=0.65,
    )
    maximum = float(top10["closest_distance_score"].max())
    for bar, (_, row) in zip(bars, top10.iterrows()):
        bar_ax.text(
            bar.get_width() + maximum * 0.018,
            bar.get_y() + bar.get_height() / 2,
            f"{float(row['closest_distance_score']):.3f}",
            ha="left",
            va="center",
            fontsize=6,
            fontweight=(
                "bold"
                if row["cell_state"] in ("Fib_K14", "KC_Spinous_Mig")
                else "normal"
            ),
        )
    for tick in bar_ax.get_yticklabels():
        if tick.get_text() in ("Fib_K14", "KC_Spinous_Mig"):
            tick.set_color(CELL_COLORS[tick.get_text()])
            tick.set_fontweight("bold")
    bar_ax.tick_params(axis="y", labelsize=6)
    bar_ax.set_xlim(0, maximum * 1.18)
    bar_ax.set_xlabel("Closest-distance co-occurrence score")
    bar_ax.set_title(
        "Top 10 states including Fib_K14 self-reference",
        fontweight="bold",
        pad=5,
    )
    clean(bar_ax)

    stats_26 = section_stats[
        (section_stats["section"] == "26dpb_p2")
        & (section_stats["cell_state"].isin(["Fib_K14", "KC_Spinous_Mig"]))
    ].set_index("cell_state")
    fib = stats_26.loc["Fib_K14"]
    kc = stats_26.loc["KC_Spinous_Mig"]
    fig.text(
        0.075,
        0.105,
        (
            "Fixed 50-µm permutation test (999 permutations; BH across Fib_K14 + all KC/SAC states):\n"
            f"Fib_K14 self: score={float(fib['cooccurrence_score']):.3f}, "
            f"P={float(fib['p_enrichment']):.3f}, "
            f"q={float(fib['fdr_within_17_FibK14_KC_SAC_states']):.3f}; "
            f"KC_Spinous_Mig: score={float(kc['cooccurrence_score']):.3f}, "
            f"P={float(kc['p_enrichment']):.3f}, "
            f"q={float(kc['fdr_within_17_FibK14_KC_SAC_states']):.3f} (ns)."
        ),
        fontsize=7,
        ha="left",
        va="top",
        linespacing=1.35,
    )
    fig.text(
        0.075,
        0.055,
        (
            "Curves reproduce the cumulative-radius probability-ratio algorithm; "
            "distance was converted using 0.5 µm per coordinate unit."
        ),
        fontsize=6,
        ha="left",
        va="top",
    )


def build() -> None:
    configure_style()
    top10 = pd.read_csv(SRC / "FigureS3_FibK14_KC_SAC_timepoint_TOP10.csv")
    section_stats = pd.read_csv(SRC / "FigureS3_v04_prespecified_section_statistics.csv")
    time_stats = pd.read_csv(SRC / "FigureS3_v04_prespecified_timepoint_statistics.csv")
    curves = pd.read_csv(SRC / "FigureS3_v04_26dpb_p2_distance_curves.csv")
    closest = pd.read_csv(SRC / "FigureS3_v04_26dpb_p2_closest_scores.csv")

    fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4), facecolor="white")
    panel_a(fig, top10)
    panel_b(fig, section_stats, time_stats)
    panel_c(fig, curves, closest, section_stats)

    svg = VEC / "FigureS3_revised_v04_editable.svg"
    pdf = VEC / "FigureS3_revised_v04.pdf"
    png = RAS / "FigureS3_revised_v04_600dpi.png"
    tiff = RAS / "FigureS3_revised_v04_600dpi.tiff"
    fig.savefig(svg)
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    fig.savefig(tiff, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    with (QC / "FigureS3_v04_output_sha256.txt").open("w", encoding="utf-8") as handle:
        for path in (svg, pdf, png, tiff):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")


if __name__ == "__main__":
    build()
