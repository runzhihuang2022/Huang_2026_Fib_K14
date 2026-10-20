"""Render the approved PAGA content with non-overlapping labels."""

from pathlib import Path
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


root = Path(sys.argv[1])
matrix = pd.read_csv(root / "tables" / "S5D_rat_PAGA_connectivity.csv", index_col=0)
colors = {
    "Fib_K14": "#E31A1C", "SAC_Progenitor": "#0099CC",
    "KC_Basal": "#008B45", "KC_Basal_Mig": "#FF8C00",
    "KC_Basal_Prolif": "#F1C40F", "KC_Spinous": "#00C000",
    "KC_Spinous_Mig": "#00008B", "KC_Spinous_Mat": "#FF69B4",
    "KC_Granular": "#8B4513",
}
pos = {
    "Fib_K14": (-1.55, 0.15), "SAC_Progenitor": (-0.75, 0.78),
    "KC_Basal": (0.00, -0.40), "KC_Basal_Mig": (1.48, 0.38),
    "KC_Basal_Prolif": (0.82, 0.72), "KC_Spinous": (0.25, 0.12),
    "KC_Spinous_Mig": (1.05, -0.62), "KC_Spinous_Mat": (-0.05, 0.62),
    "KC_Granular": (0.78, -0.18),
}
label_pos = {
    "Fib_K14": (-1.55, -0.08), "SAC_Progenitor": (-0.75, 1.00),
    "KC_Basal": (-0.25, -0.67), "KC_Basal_Mig": (1.48, 0.10),
    "KC_Basal_Prolif": (0.98, 0.94), "KC_Spinous": (0.20, -0.12),
    "KC_Spinous_Mig": (1.28, -0.82), "KC_Spinous_Mat": (-0.20, 0.88),
    "KC_Granular": (0.75, -0.45),
}
display = {x: ("Fib_K14-like" if x == "Fib_K14" else x) for x in matrix.index}
mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "pdf.fonttype": 42,
    "ps.fonttype": 42, "svg.fonttype": "none",
})
G = nx.Graph()
for state in matrix.index:
    G.add_node(state)
for i, left in enumerate(matrix.index):
    for right in matrix.index[i + 1:]:
        value = float(matrix.loc[left, right])
        if value >= 0.02:
            G.add_edge(left, right, weight=value)
fig, ax = plt.subplots(figsize=(7.08, 1.85))
nx.draw_networkx_edges(
    G, pos, width=[0.35 + 4.2 * G.edges[e]["weight"] for e in G.edges],
    edge_color="#666666", alpha=0.72, ax=ax,
)
nx.draw_networkx_nodes(
    G, pos, node_size=220, node_color=[colors[x] for x in G.nodes],
    edgecolors="black", linewidths=0.55, ax=ax,
)
for node, xy in label_pos.items():
    ax.text(
        *xy, display[node], ha="center", va="center", fontsize=6.2,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.4},
    )
ax.set_title("Rat Fib_K14-like–keratinocyte PAGA connectivity", fontsize=8)
ax.set_xlim(-1.9, 1.85); ax.set_ylim(-1.02, 1.15); ax.axis("off")
stem = root / "figures" / "panels" / "FigureS5_D_rat_PAGA"
fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
plt.close(fig)
