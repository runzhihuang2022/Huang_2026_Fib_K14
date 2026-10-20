"""Render the precomputed 11-section density grids with a stable backend."""

from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


LABELS = [
    "Normal", "5 dpb", "12 dpb DPTDI-1", "12 dpb SPTDI-1",
    "12 dpb DPTDI-2", "12 dpb SPTDI-2", "19 dpb", "19 dpb p1",
    "26 dpb p1", "26 dpb p2", "2 mph",
]


def main() -> None:
    root = Path(sys.argv[1])
    data = np.load(root / "tables" / "S4F_spatial_density_grids.npz")
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 7,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })
    fig, axes = plt.subplots(3, 4, figsize=(7.08, 4.95))
    axes = axes.ravel()
    last = None
    for i, (ax, label) in enumerate(zip(axes, LABELS)):
        grid = data[f"grid_{i}"]
        extent = data[f"extent_{i}"]
        vmax = max(np.nanquantile(grid, .995), 1e-6)
        last = ax.imshow(
            grid, origin="lower", extent=extent, cmap="viridis",
            interpolation="bilinear", aspect="equal", norm=Normalize(0, vmax),
            rasterized=True,
        )
        ax.set_title(label, pad=2)
        ax.axis("off")
    for ax in axes[len(LABELS):]:
        ax.axis("off")
    cax = fig.add_axes([.76, .07, .18, .018])
    cb = fig.colorbar(last, cax=cax, orientation="horizontal")
    cb.set_label("PDGFRA + KRT14 joint spatial density", fontsize=6)
    cb.set_ticks([])
    fig.subplots_adjust(left=.02, right=.98, top=.96, bottom=.12,
                        wspace=.08, hspace=.15)
    stem = root / "figures" / "panels" / "FigureS4_F_all_human_spatial_joint_density"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
