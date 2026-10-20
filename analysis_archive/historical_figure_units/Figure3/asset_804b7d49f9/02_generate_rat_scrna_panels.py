"""Generate rat Figure S5 panels A, C and D and export Monocle2 input.

Krt14 is absent from the registered rat gene universe. The script therefore
renders an explicit unavailable tile and computes the sixth tile from the
prespecified surrogate program Pdgfra/Vim/Krt5/Tacstd2. It never substitutes
the surrogate for a strict Pdgfra+Krt14 double-positive definition.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.io import mmwrite
from scipy.ndimage import gaussian_filter


FIB_ORDER = [
    "Fib_Papi", "Fib_Sfrp2", "Fib_Myo", "Fib_Inflama", "Fib_Fasci",
    "Fib_En1", "Fib_Prolif", "Fib_K14", "Fib_Cd45",
]
FIB_COLORS = {
    "Fib_Papi": "#4682B4", "Fib_Sfrp2": "#87CEFA", "Fib_Myo": "#DDA0DD",
    "Fib_Inflama": "#BA55D3", "Fib_Fasci": "#9370DB", "Fib_En1": "#F0E68C",
    "Fib_Prolif": "#FFB6C1", "Fib_K14": "#E31A1C", "Fib_Cd45": "#9932CC",
}
DISPLAY = {x: ("Fib_K14-like" if x == "Fib_K14" else x) for x in FIB_ORDER}
MARKERS = ["Pdgfra", "Vim", "Krt14", "Krt5", "Tacstd2"]
SURROGATE = ["Pdgfra", "Vim", "Krt5", "Tacstd2"]
TRAJECTORY_STATES = [
    "Fib_K14", "SAC_Progenitor", "KC_Basal", "KC_Basal_Mig",
    "KC_Spinous_Mig", "KC_Spinous",
]
STATE_COLORS = {
    "Fib_K14": "#E31A1C", "SAC_Progenitor": "#0099CC",
    "KC_Basal": "#008B45", "KC_Basal_Mig": "#FF8C00",
    "KC_Basal_Prolif": "#F1C40F", "KC_Spinous": "#00C000",
    "KC_Spinous_Mig": "#00008B", "KC_Spinous_Mat": "#FF69B4",
    "KC_Granular": "#8B4513",
}


def setup() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.linewidth": 0.7,
    })
    sns.set_style("white")


def save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def dense(x) -> np.ndarray:
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def local_density(coords: np.ndarray, weights: np.ndarray, bins=240, sigma=3.2) -> np.ndarray:
    weights = np.nan_to_num(weights, nan=0.0)
    lo = coords.min(0); hi = coords.max(0)
    ij = np.clip(
        ((coords - lo) / np.maximum(hi - lo, 1e-9) * (bins - 1)).astype(int),
        0, bins - 1,
    )
    weighted = np.zeros((bins, bins)); occupied = np.zeros((bins, bins))
    np.add.at(weighted, (ij[:, 1], ij[:, 0]), weights)
    np.add.at(occupied, (ij[:, 1], ij[:, 0]), 1.0)
    sm_w = gaussian_filter(weighted, sigma=sigma)
    sm_n = gaussian_filter(occupied, sigma=sigma)
    den = np.divide(sm_w, sm_n, out=np.zeros_like(sm_w), where=sm_n > 1e-6)
    return den[ij[:, 1], ij[:, 0]]


def panel_a_c(adata: ad.AnnData, figures: Path, tables: Path) -> None:
    labels = adata.obs["sub_labels"].astype(str)
    fib_idx = np.flatnonzero(labels.str.startswith("Fib_").to_numpy())
    harmony = np.asarray(adata.obsm["X_harmony"])[fib_idx, :40]
    obs = adata.obs.iloc[fib_idx].copy()
    obs["sub_labels"] = pd.Categorical(obs["sub_labels"].astype(str), FIB_ORDER)
    fib = ad.AnnData(X=harmony, obs=obs)
    sc.pp.neighbors(fib, n_neighbors=20, use_rep="X", random_state=17)
    sc.tl.umap(fib, min_dist=0.35, spread=1.0, random_state=17)
    sc.tl.leiden(fib, resolution=0.8, key_added="fib_leiden", random_state=17)
    coords = np.asarray(fib.obsm["X_umap"])
    fib.obs.assign(UMAP1=coords[:, 0], UMAP2=coords[:, 1]).to_csv(
        tables / "S5A_rat_fibroblast_reclustered_coordinates.csv.gz"
    )

    fig, ax = plt.subplots(figsize=(3.8, 2.75))
    for state in FIB_ORDER:
        m = fib.obs["sub_labels"].astype(str).eq(state).to_numpy()
        if m.any():
            ax.scatter(
                coords[m, 0], coords[m, 1], s=0.45, linewidths=0,
                color=FIB_COLORS[state], label=DISPLAY[state], rasterized=True,
            )
    ax.set(title="Rat fibroblast-only re-clustering", xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(frameon=False, ncol=3, markerscale=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.12))
    sns.despine(ax=ax)
    save(fig, figures / "FigureS5_A_rat_fibroblast_reclustered_umap")

    available = [g for g in MARKERS if g in adata.var_names]
    raw = dense(adata[fib_idx, available].X)
    counts = adata.obs.iloc[fib_idx]["nCount_RNA.x"].to_numpy(float)
    norm = np.log1p(raw / np.maximum(counts[:, None], 1.0) * 1e4)
    expr = {g: norm[:, j] for j, g in enumerate(available)}
    out = pd.DataFrame({"cell_id": fib.obs_names.astype(str)})
    for gene in MARKERS:
        out[gene] = expr.get(gene, np.nan)
    out.to_csv(tables / "S5C_rat_marker_log_normalized_expression.csv.gz", index=False)

    fig, axes = plt.subplots(1, 6, figsize=(7.08, 1.42))
    for j, (ax, gene) in enumerate(zip(axes[:5], MARKERS)):
        if gene not in expr:
            ax.set_facecolor("#F2F2F2")
            ax.text(0.5, 0.58, "Not available", ha="center", va="center",
                    transform=ax.transAxes, fontsize=7, fontweight="bold")
            ax.text(0.5, 0.42, "absent from rat\ngene universe", ha="center",
                    va="center", transform=ax.transAxes, fontsize=6)
            ax.set_title(gene, fontstyle="italic")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#999999")
            continue
        values = expr[gene]
        positive = values[values > 0]
        vmax = np.quantile(positive, 0.99) if len(positive) else 1.0
        order = np.argsort(values)
        scatter = ax.scatter(
            coords[order, 0], coords[order, 1], c=values[order],
            s=0.32, cmap="viridis", vmin=0, vmax=vmax, linewidths=0,
            rasterized=True,
        )
        ax.set_title(gene, fontstyle="italic")
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if j == 4:
            cb = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.01)
            cb.ax.set_title("Expr.", fontsize=6, pad=2); cb.ax.tick_params(labelsize=6)

    scaled = []
    for gene in SURROGATE:
        value = expr[gene]
        positive = value[value > 0]
        q = np.quantile(positive, 0.99) if len(positive) else 1.0
        scaled.append(np.clip(value / max(q, 1e-9), 0, 1))
    weight = np.prod(np.column_stack(scaled), axis=1) ** (1 / len(scaled))
    density = local_density(coords, weight)
    order = np.argsort(density)
    scatter = axes[5].scatter(
        coords[order, 0], coords[order, 1], c=density[order],
        s=0.32, cmap="viridis", linewidths=0, rasterized=True,
    )
    axes[5].set_title("Fib_K14-like\nsurrogate density")
    axes[5].set_aspect("equal"); axes[5].set_xticks([]); axes[5].set_yticks([])
    cb = fig.colorbar(scatter, ax=axes[5], fraction=0.046, pad=0.01)
    cb.ax.set_title("Density", fontsize=6, pad=2); cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(wspace=0.13)
    save(fig, figures / "FigureS5_C_rat_marker_featureplots")


def panel_d(adata: ad.AnnData, figures: Path, tables: Path) -> None:
    states = [
        "Fib_K14", "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif",
        "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
        "SAC_Progenitor",
    ]
    labels = adata.obs["sub_labels"].astype(str)
    rng = np.random.default_rng(17)
    selected: list[int] = []
    for state in states:
        idx = np.flatnonzero(labels.eq(state).to_numpy())
        if len(idx):
            selected.extend(rng.choice(idx, min(4000, len(idx)), replace=False))
    selected = np.asarray(selected)
    obs = adata.obs.iloc[selected][["sub_labels"]].copy()
    obs["sub_labels"] = pd.Categorical(obs["sub_labels"].astype(str), states)
    rep = np.asarray(adata.obsm["X_harmony"])[selected, :40]
    graph = ad.AnnData(X=rep, obs=obs)
    sc.pp.neighbors(graph, n_neighbors=20, use_rep="X", random_state=17)
    sc.tl.paga(graph, groups="sub_labels")
    connectivity = graph.uns["paga"]["connectivities"].toarray()
    pd.DataFrame(connectivity, index=states, columns=states).to_csv(
        tables / "S5D_rat_PAGA_connectivity.csv"
    )

    G = nx.Graph()
    for state in states:
        G.add_node(state)
    for i in range(len(states)):
        for j in range(i + 1, len(states)):
            if connectivity[i, j] >= 0.02:
                G.add_edge(states[i], states[j], weight=float(connectivity[i, j]))
    pos = nx.spring_layout(G, seed=17, weight="weight", k=0.85)
    fig, ax = plt.subplots(figsize=(7.08, 1.85))
    widths = [0.4 + 4.0 * G.edges[e]["weight"] for e in G.edges]
    nx.draw_networkx_edges(G, pos, width=widths, edge_color="#555555", alpha=0.75, ax=ax)
    nx.draw_networkx_nodes(
        G, pos, node_size=260,
        node_color=[STATE_COLORS.get(x, "#888888") for x in G.nodes],
        edgecolors="black", linewidths=0.6, ax=ax,
    )
    labels_display = {x: ("Fib_K14-like" if x == "Fib_K14" else x) for x in G.nodes}
    nx.draw_networkx_labels(G, pos, labels=labels_display, font_family="Arial",
                            font_size=7, ax=ax)
    ax.set_title("Rat Fib_K14-like–keratinocyte PAGA connectivity")
    ax.axis("off")
    save(fig, figures / "FigureS5_D_rat_PAGA")


def export_monocle(adata: ad.AnnData, tables: Path) -> None:
    labels = adata.obs["sub_labels"].astype(str)
    rng = np.random.default_rng(17)
    selected: list[int] = []
    represented = []
    for state in TRAJECTORY_STATES:
        idx = np.flatnonzero(labels.eq(state).to_numpy())
        if len(idx):
            represented.append(state)
        selected.extend(rng.choice(idx, min(1000, len(idx)), replace=False))
    selected = np.asarray(selected)
    traj = adata[selected].to_memory()
    sc.pp.filter_genes(traj, min_cells=20)
    counts = traj.X.tocsr() if sparse.issparse(traj.X) else sparse.csr_matrix(traj.X)
    mmwrite(tables / "S5E_monocle2_counts_cells_by_genes.mtx", counts)
    pd.DataFrame({"gene": traj.var_names.astype(str)}).to_csv(
        tables / "S5E_monocle2_genes.csv", index=False
    )
    pheno = traj.obs.copy()
    pheno["cell_id"] = traj.obs_names.astype(str)
    pheno.to_csv(tables / "S5E_monocle2_phenodata.csv", index=False)
    pd.DataFrame({"represented_state": represented}).to_csv(
        tables / "S5E_represented_states.csv", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scrna", required=True)
    parser.add_argument("--package", required=True)
    args = parser.parse_args()
    package = Path(args.package)
    figures = package / "figures" / "panels"
    tables = package / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    setup()
    adata = ad.read_h5ad(args.scrna, backed="r")
    panel_a_c(adata, figures, tables)
    panel_d(adata, figures, tables)
    export_monocle(adata, tables)


if __name__ == "__main__":
    main()
