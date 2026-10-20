from __future__ import annotations

from pathlib import Path
import hashlib

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
VEC = ROOT / "outputs" / "vector"
RAS = ROOT / "outputs" / "raster"
QC = ROOT / "QC"
for directory in (VEC, RAS, QC):
    directory.mkdir(parents=True, exist_ok=True)

SECTION_ORDER = [
    "Normal", "5dpb", "12dpb_DPTDI1", "12dpb_SPTDI1", "12dpb_DPTDI2",
    "12dpb_SPTDI2", "19dpb", "19dpb_p1", "26dpb_p1", "26dpb_p2", "2mph",
]
SECTION_LABELS = [
    "Normal", "5 dpb", "12D1", "12S1", "12D2", "12S2",
    "19-1", "19-2", "26-1", "26-2", "2 mph",
]
TIME_ORDER = ["Normal", "5 dpb", "12 dpb", "19 dpb", "26 dpb", "2 mph"]
CELL_ORDER = [
    "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
    "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
    "SAC_SG_Progenitor", "SAC_SG_Clear", "SAC_SG_Dark", "SAC_SG_Ductal",
    "SAC_HF_IRS", "SAC_HF_ORS", "SAC_HF_HFSC", "SAC_HF_DP_DS",
    "SAC_HF_Matrix",
]
CELL_COLORS = {
    "KC_Basal": "#008941", "KC_Basal_Mig": "#FF8C00",
    "KC_Basal_Prolif": "#FFD700", "KC_Spinous": "#00D600",
    "KC_Spinous_Mig": "#1A0099", "KC_Spinous_Mat": "#EE4C97",
    "KC_Granular": "#882D17", "SAC_SG_Progenitor": "#0099CC",
    "SAC_SG_Clear": "#A8E6E6", "SAC_SG_Dark": "#193006",
    "SAC_SG_Ductal": "#CD853F", "SAC_HF_IRS": "#CC9900",
    "SAC_HF_ORS": "#A30059", "SAC_HF_HFSC": "#20B2AA",
    "SAC_HF_DP_DS": "#CDBE70", "SAC_HF_Matrix": "#48C9C9",
}


def add_letter(fig: mpl.figure.Figure, x: float, y: float, letter: str) -> None:
    fig.text(x, y, letter, fontsize=14, fontweight="bold", va="top")


def clean(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.7, length=2.5)


def fdr_label(q: float) -> str:
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def vector_color_scale(
    fig: mpl.figure.Figure,
    position: list[float],
    cmap: mpl.colors.Colormap,
    norm: mpl.colors.Normalize,
    vmin: float,
    vmax: float,
    ticks: list[float],
    label: str,
) -> mpl.axes.Axes:
    ax = fig.add_axes(position)
    edges = np.linspace(vmin, vmax, 41)
    for low, high in zip(edges[:-1], edges[1:]):
        mid = (low + high) / 2
        ax.add_patch(
            mpl.patches.Rectangle(
                (0, low), 1, high - low,
                facecolor=cmap(norm(mid)), edgecolor="none",
            )
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(vmin, vmax)
    ax.set_xticks([])
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_yticks(ticks)
    ax.set_ylabel(label, rotation=90, labelpad=3)
    ax.tick_params(width=0.6, length=2, pad=2)
    return ax


def build() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    enrichment = pd.read_csv(SRC / "Figure2B_S2AB_19dpb_p1_region_enrichment.csv")
    section = pd.read_csv(SRC / "FigureS3_FibK14_KC_SAC_section_scores.csv")
    time_stats = pd.read_csv(SRC / "FigureS3_FibK14_KC_SAC_timepoint_summary.csv")
    top10 = pd.read_csv(SRC / "FigureS3_FibK14_KC_SAC_timepoint_TOP10.csv")

    fig = plt.figure(figsize=(8.2677, 11.6929))

    # A. Rank cell states by positive Epi-Front enrichment and show corrected statistics.
    add_letter(fig, 0.025, 0.975, "A")
    ax_a = fig.add_axes([0.17, 0.665, 0.31, 0.275])
    epi = (
        enrichment[(enrichment["region"] == "Epi-Front") & (enrichment["log2_enrichment"] > 0)]
        .nlargest(12, "log2_enrichment")
        .sort_values("log2_enrichment")
    )
    y = np.arange(len(epi))
    significant = epi["fisher_fdr"].to_numpy(float) < 0.05
    bar_colors = np.where(significant, "#D7191C", "#C9C9C9")
    ax_a.barh(
        y, epi["log2_enrichment"], color=bar_colors,
        edgecolor="#333333", linewidth=0.45,
    )
    ax_a.set_yticks(y, epi["cell_type"])
    ax_a.set_xlabel("log2 enrichment score")
    ax_a.set_title("Cell states enriched in Epi-Front (19 dpb_p1)", fontweight="bold")
    ax_a.axvline(0, color="#333333", lw=0.7)
    xmax = float(epi["log2_enrichment"].max()) * 1.25
    ax_a.set_xlim(0, xmax)
    for yi, (_, row) in enumerate(epi.iterrows()):
        q = float(row["fisher_fdr"])
        label = f"{fdr_label(q)}" if q < 0.05 else ""
        if label:
            ax_a.text(
                float(row["log2_enrichment"]) + xmax * 0.018, yi, label,
                va="center", ha="left", fontsize=8, fontweight="bold",
            )
    for tick, is_sig in zip(ax_a.get_yticklabels(), significant):
        if is_sig:
            tick.set_color("#D7191C")
            tick.set_fontweight("bold")
    clean(ax_a)
    # B. Section-resolved Fib_K14-conditioned co-occurrence.
    add_letter(fig, 0.515, 0.975, "B")
    ax_b = fig.add_axes([0.665, 0.67, 0.275, 0.27])
    score_matrix = (
        section.pivot(index="cell_state", columns="section", values="cooccurrence_score")
        .reindex(index=CELL_ORDER, columns=SECTION_ORDER)
    )
    mesh = ax_b.pcolormesh(
        np.arange(len(SECTION_ORDER) + 1),
        np.arange(len(CELL_ORDER) + 1),
        score_matrix.to_numpy(float),
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=3),
        shading="flat",
        edgecolors="white",
        linewidth=0.2,
    )
    ax_b.set_xlim(0, len(SECTION_ORDER))
    ax_b.set_ylim(len(CELL_ORDER), 0)
    ax_b.set_xticks(
        np.arange(len(SECTION_ORDER)) + 0.5,
        SECTION_LABELS, rotation=50, ha="left",
    )
    ax_b.xaxis.tick_top()
    ax_b.set_yticks(np.arange(len(CELL_ORDER)) + 0.5, CELL_ORDER)
    ax_b.tick_params(length=0, pad=1)
    ax_b.set_title("Fib_K14-conditioned co-occurrence by section", fontweight="bold", pad=5)
    for tick in ax_b.get_yticklabels():
        tick.set_color(CELL_COLORS.get(tick.get_text(), "#222222"))
    vector_color_scale(
        fig, [0.948, 0.735, 0.012, 0.14],
        mpl.colormaps["RdBu_r"], mesh.norm, 0, 3, [0, 1, 2, 3], "Score",
    )
    fig.text(
        0.665, 0.645,
        "Fixed 50-µm radius; score >1 denotes enrichment",
        fontsize=8,
    )

    # C. Time-resolved TOP10 displayed as one compact rank/score matrix.
    add_letter(fig, 0.025, 0.60, "C")
    ax_c = fig.add_axes([0.20, 0.315, 0.70, 0.25])
    score_top = pd.DataFrame(np.nan, index=CELL_ORDER, columns=TIME_ORDER)
    rank_top = pd.DataFrame(np.nan, index=CELL_ORDER, columns=TIME_ORDER)
    for _, row in top10.iterrows():
        cell = row["cell_state"]
        timepoint = row["timepoint"]
        if cell in score_top.index and timepoint in score_top.columns:
            score_top.loc[cell, timepoint] = float(row["mean_score"])
            rank_top.loc[cell, timepoint] = int(row["rank_within_time"])
    cmap_c = mpl.colormaps["RdBu_r"].copy()
    cmap_c.set_bad("#F0F0F0")
    max_top = max(2.0, float(np.nanmax(score_top.to_numpy(float))))
    heat_c = ax_c.pcolormesh(
        np.arange(len(TIME_ORDER) + 1),
        np.arange(len(CELL_ORDER) + 1),
        score_top.to_numpy(float),
        cmap=cmap_c,
        norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=max_top),
        shading="flat",
        edgecolors="white",
        linewidth=0.35,
    )
    ax_c.set_xlim(0, len(TIME_ORDER))
    ax_c.set_ylim(len(CELL_ORDER), 0)
    ax_c.set_xticks(np.arange(len(TIME_ORDER)) + 0.5, TIME_ORDER)
    ax_c.set_yticks(np.arange(len(CELL_ORDER)) + 0.5, CELL_ORDER)
    ax_c.xaxis.tick_top()
    ax_c.tick_params(length=0, pad=3)
    ax_c.set_title(
        "Time-resolved TOP10 Fib_K14–KC/SAC co-occurrence scores",
        fontweight="bold", pad=7,
    )
    for tick in ax_c.get_yticklabels():
        tick.set_color(CELL_COLORS.get(tick.get_text(), "#222222"))
    for row_i in range(len(CELL_ORDER)):
        for col_i in range(len(TIME_ORDER)):
            score = score_top.iat[row_i, col_i]
            rank = rank_top.iat[row_i, col_i]
            if np.isfinite(score):
                color = "white" if score > max(1.45, max_top * 0.56) else "black"
                ax_c.text(
                    col_i + 0.5, row_i + 0.5,
                    f"{score:.2f}\n#{int(rank)}",
                    ha="center", va="center", fontsize=7, color=color,
                )
    vector_color_scale(
        fig, [0.915, 0.355, 0.014, 0.16],
        cmap_c, heat_c.norm, 0, max_top, [0, 1, max_top], "Mean score",
    )
    fig.text(
        0.20, 0.292,
        "Values are time-point means; rank is shown after #. Gray cells fall outside that time point's TOP10.",
        fontsize=8,
    )

    # D. Prespecified target trajectories with individual-section points.
    add_letter(fig, 0.025, 0.245, "D")
    ax_d = fig.add_axes([0.13, 0.075, 0.80, 0.135])
    target_states = ["KC_Spinous_Mig", "SAC_SG_Progenitor"]
    x = np.arange(len(TIME_ORDER))
    for state in target_states:
        summary = (
            time_stats[time_stats["cell_state"] == state]
            .set_index("timepoint")
            .reindex(TIME_ORDER)
        )
        color = CELL_COLORS[state]
        ax_d.plot(
            x, summary["mean_score"], marker="o", markersize=4.5,
            lw=1.7, color=color, label=state,
        )
        raw = section[section["cell_state"] == state]
        for xi, timepoint in enumerate(TIME_ORDER):
            vals = raw[raw["timepoint"] == timepoint]["cooccurrence_score"].to_numpy()
            jitter = np.linspace(-0.055, 0.055, len(vals)) if len(vals) > 1 else np.zeros(len(vals))
            ax_d.scatter(
                xi + jitter, vals, s=18, color=color,
                edgecolor="white", linewidth=0.4, zorder=3,
            )
    ax_d.axhline(1, color="#555555", lw=0.8, ls="--")
    ax_d.set_xticks(x, TIME_ORDER)
    ax_d.set_ylabel("Fib_K14-conditioned\nco-occurrence score")
    ax_d.set_title(
        "Temporal co-occurrence of prespecified epithelial/appendage states",
        fontweight="bold",
    )
    ax_d.legend(frameon=False, ncol=2, loc="upper right")
    clean(ax_d)
    fig.text(
        0.50, 0.035,
        "Statistics: 999 spatial-label permutations per section; BH correction across all KC/SAC states. "
        "No positive comparison passed FDR < 0.05.",
        ha="center", fontsize=8, color="#B2182B", fontweight="bold",
    )

    svg = VEC / "FigureS3_revised_v03_editable.svg"
    pdf = VEC / "FigureS3_revised_v03.pdf"
    png = RAS / "FigureS3_revised_v03_600dpi.png"
    tif = RAS / "FigureS3_revised_v03_600dpi.tiff"
    fig.savefig(svg)
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    fig.savefig(tif, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    with (QC / "output_sha256.txt").open("w", encoding="utf-8") as handle:
        for path in (svg, pdf, png, tif):
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{digest}  {path.name}\n")


if __name__ == "__main__":
    build()
