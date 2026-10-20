from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source_data"
OUT = ROOT / "panels" / "Figure2C_optimized"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "KC_Spinous_Mig": "#25378C",
    "SAC_SG_Progenitor": "#13A2B3",
    "Fib_K14": "#D94B45",
}

MM_TO_IN = 1 / 25.4
FIGURE_WIDTH_MM = 170  # Matches the target width occupied by Figure 2A.
FIGURE_HEIGHT_MM = 66

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 8,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def p_to_symbol(p: float) -> str:
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def clean_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", width=0.8, length=3, pad=2.5)


def add_statistics(axes: list[plt.Axes], stats_path: Path) -> None:
    """Draw only verified, explicitly supplied adjusted P values.

    Expected columns: projection, day, cell_type_a, cell_type_b, adjusted_p,
    test, correction, biological_n. Rows lacking a finite adjusted P value are
    ignored. Annotations are placed above the corresponding time point; the
    exact test metadata remains in the companion CSV for the legend/methods.
    """
    if not stats_path.exists():
        return
    stats = pd.read_csv(stats_path)
    required = {
        "projection",
        "day",
        "cell_type_a",
        "cell_type_b",
        "adjusted_p",
        "test",
        "correction",
        "biological_n",
    }
    if not required.issubset(stats.columns):
        raise ValueError(f"Statistics table is missing: {sorted(required - set(stats.columns))}")
    stats = stats[pd.to_numeric(stats["adjusted_p"], errors="coerce").notna()].copy()
    if stats.empty:
        return
    stats["adjusted_p"] = stats["adjusted_p"].astype(float)

    projection_to_axis = {"horizontal": axes[0], "vertical": axes[1]}
    for projection, group in stats.groupby("projection", sort=False):
        ax = projection_to_axis.get(str(projection).lower())
        if ax is None:
            continue
        y0, y1 = ax.get_ylim()
        span = abs(y1 - y0)
        # Vertical depth is inverted; axes coordinates keep placement stable.
        for level, (_, row) in enumerate(group.sort_values("day").iterrows()):
            x = float(row["day"])
            y_axes = 0.94 - 0.075 * (level % 3)
            ax.text(
                x,
                y_axes,
                p_to_symbol(float(row["adjusted_p"])),
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                clip_on=False,
            )


def build() -> None:
    data_path = DATA / "Figure2C_trajectory_centroids.csv"
    stats_path = DATA / "Figure2C_statistical_tests.csv"
    df = pd.read_csv(data_path)
    df = df[df["Cell_Type"].isin(COLORS)].copy()

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(FIGURE_WIDTH_MM * MM_TO_IN, FIGURE_HEIGHT_MM * MM_TO_IN),
        gridspec_kw={"wspace": 0.25},
    )
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.20, top=0.86)

    for state, color in COLORS.items():
        d = df[df["Cell_Type"].eq(state)].sort_values("Day")
        axes[0].plot(
            d["Day"],
            d["Median_X"] / 1000,
            marker="o",
            markersize=4.5,
            markeredgewidth=0,
            color=color,
            label=state,
        )
        axes[1].plot(
            d["Day"],
            d["Median_Y"],
            marker="o",
            markersize=4.5,
            markeredgewidth=0,
            color=color,
            label=state,
        )

    axes[0].axhline(0, color="#707070", linewidth=0.8, linestyle=(0, (3, 2)), zorder=0)
    axes[0].set_title("Horizontal trajectory", pad=8)
    axes[0].set_xlabel("Healing time (days)")
    axes[0].set_ylabel("Distance to front (mm)")

    axes[1].set_title("Vertical trajectory", pad=8)
    axes[1].set_xlabel("Healing time (days)")
    axes[1].set_ylabel("Relative depth")
    axes[1].invert_yaxis()

    for ax in axes:
        ax.set_xticks([5, 12, 19, 26, 60])
        ax.margins(x=0.05)
        clean_axis(ax)

    axes[1].legend(
        frameon=False,
        loc="lower right",
        handlelength=1.8,
        handletextpad=0.6,
        borderaxespad=0.2,
        labelspacing=0.4,
    )
    add_statistics(list(axes), stats_path)

    fig.text(0.012, 0.97, "C", fontsize=12, fontweight="bold", ha="left", va="top")
    base = OUT / "Figure2C_optimized_170mm"
    fig.savefig(base.with_suffix(".svg"), bbox_inches=None)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches=None)
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches=None)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches=None, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    build()
