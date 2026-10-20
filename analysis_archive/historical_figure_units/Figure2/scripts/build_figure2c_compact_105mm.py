from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source_data"
OUT = ROOT / "panels" / "Figure2C_compact_105mm"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "KC_Spinous_Mig": "#25378C",
    "SAC_SG_Progenitor": "#13A2B3",
    "Fib_K14": "#D94B45",
}

MM_TO_IN = 1 / 25.4
FIGURE_WIDTH_MM = 105
FIGURE_HEIGHT_MM = 48

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 7,
        "axes.titlesize": 10,
        "axes.titleweight": "normal",
        "axes.labelsize": 8,
        "axes.labelweight": "normal",
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.75,
        "lines.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", width=0.75, length=2.6, pad=2.0)


def build() -> None:
    data_path = DATA / "Figure2C_trajectory_centroids.csv"
    frame = pd.read_csv(data_path)
    frame = frame[frame["Cell_Type"].isin(COLORS)].copy()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(FIGURE_WIDTH_MM * MM_TO_IN, FIGURE_HEIGHT_MM * MM_TO_IN),
        gridspec_kw={"wspace": 0.28},
    )
    fig.subplots_adjust(left=0.085, right=0.992, bottom=0.22, top=0.84)

    for state, color in COLORS.items():
        data = frame[frame["Cell_Type"].eq(state)].sort_values("Day")
        axes[0].plot(
            data["Day"],
            data["Median_X"] / 1000,
            marker="o",
            markersize=3.6,
            markeredgewidth=0,
            color=color,
            label=state,
        )
        axes[1].plot(
            data["Day"],
            data["Median_Y"],
            marker="o",
            markersize=3.6,
            markeredgewidth=0,
            color=color,
            label=state,
        )

    axes[0].axhline(
        0,
        color="#707070",
        linewidth=0.7,
        linestyle=(0, (3, 2)),
        zorder=0,
    )
    axes[0].set_title("Horizontal trajectory", pad=5)
    axes[0].set_xlabel("Healing time (days)", labelpad=3)
    axes[0].set_ylabel("Distance to front (mm)", labelpad=3)

    axes[1].set_title("Vertical trajectory", pad=5)
    axes[1].set_xlabel("Healing time (days)", labelpad=3)
    axes[1].set_ylabel("Relative depth", labelpad=3)
    axes[1].invert_yaxis()

    for ax in axes:
        ax.set_xticks([5, 12, 19, 26, 60])
        ax.margins(x=0.05)
        clean_axis(ax)

    axes[1].legend(
        frameon=False,
        loc="lower right",
        handlelength=1.7,
        handletextpad=0.5,
        borderaxespad=0.1,
        labelspacing=0.25,
        columnspacing=0.6,
    )

    base = OUT / "Figure2C_compact_105mm"
    fig.savefig(base.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches=None)
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches=None,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    readme = "\n".join(
        [
            "# Figure 2C compact 105-mm panel",
            "",
            "- Final width: 105 mm (one-half of A4 width).",
            "- Final height: 48 mm.",
            "- Typeface: Arial.",
            "- Typography: titles 10 pt; axis labels 8 pt; ticks and legend 7 pt.",
            "- Data and axis direction are unchanged from the 170-mm source panel.",
            "- No statistical annotations were added because "
            "`Figure2C_statistical_tests.csv` contains no verified tests.",
        ]
    )
    (OUT / "README.md").write_text(readme + "\n", encoding="utf-8")


if __name__ == "__main__":
    build()
