"""Render compact publication panels from the Monocle2 DDRTree coordinates."""

from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


root = Path(sys.argv[1])
df = pd.read_csv(root / "tables" / "S5E_monocle2_coordinates_pseudotime.csv")
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
    "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 5.5, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.linewidth": 0.7,
})
states = [
    "Fib_K14", "SAC_Progenitor", "KC_Basal", "KC_Basal_Mig",
    "KC_Spinous_Mig", "KC_Spinous",
]
state_colors = {
    "Fib_K14": "#E31A1C", "SAC_Progenitor": "#0099CC",
    "KC_Basal": "#008B45", "KC_Basal_Mig": "#FF8C00",
    "KC_Spinous_Mig": "#00008B", "KC_Spinous": "#00A000",
}
display = {x: ("Fib_K14-like" if x == "Fib_K14" else x) for x in states}
time_order = [
    "Normal", "0dpi", "0.5dpi", "1dpi", "2dpi", "3dpi",
    "5dpi", "7dpi", "10dpi", "11dpi", "15dpi", "19dpi",
]
time_colors = dict(zip(time_order, mpl.colormaps["Spectral"](np.linspace(0.02, 0.98, len(time_order)))))

fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.15))
x, y = df["Component_1"].to_numpy(), df["Component_2"].to_numpy()
order = np.argsort(df["Pseudotime"].to_numpy())
scatter = axes[0].scatter(
    x[order], y[order], c=df["Pseudotime"].to_numpy()[order],
    s=1.4, cmap="plasma", linewidths=0, rasterized=True,
)
axes[0].set_title("Monocle2 pseudotime")
cb = fig.colorbar(scatter, ax=axes[0], fraction=0.045, pad=0.02)
cb.set_label("Pseudotime", fontsize=6); cb.ax.tick_params(labelsize=5.5)

for state in states:
    mask = df["sub_labels"].astype(str).eq(state).to_numpy()
    axes[1].scatter(
        x[mask], y[mask], s=1.4, color=state_colors[state],
        label=display[state], linewidths=0, rasterized=True,
    )
axes[1].set_title("Cell-state mapping")
axes[1].legend(
    frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.16),
    markerscale=2.2, columnspacing=0.8, handletextpad=0.25,
)

for time in time_order:
    mask = df["time_point"].astype(str).eq(time).to_numpy()
    if mask.any():
        axes[2].scatter(
            x[mask], y[mask], s=1.4, color=time_colors[time],
            label=time, linewidths=0, rasterized=True,
        )
axes[2].set_title("Observed wound time")
axes[2].legend(
    frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.16),
    markerscale=2.2, columnspacing=0.7, handletextpad=0.25,
)

for ax in axes:
    ax.set_xlabel("Component 1"); ax.set_ylabel("Component 2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
fig.subplots_adjust(left=0.06, right=0.99, top=0.90, bottom=0.30, wspace=0.30)
stem = root / "figures" / "panels" / "FigureS5_E_rat_monocle2_trajectory"
fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
fig.savefig(stem.with_suffix(".svg"), facecolor="white")
plt.close(fig)
