"""Generate human-only Figure S4 source panels A, C, D and F.

Panel A is a fibroblast-only re-embedding. Panel C uses that same embedding for
five expression feature plots and a PDGFRA+KRT14 joint-density plot. Panel D
is regenerated without changing the previously approved statistic or layout.
Panel F uses the registered 11-section human Stereo-seq subset only.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.ndimage import gaussian_filter
from scipy.io import mmwrite


FIB_ORDER = [
    "Fib_Papi", "Fib_SFRP2", "Fib_Myo", "Fib_Inflama", "Fib_Fasci",
    "Fib_EN1", "Fib_Prolif", "Fib_K14", "Fib_CD45",
]
FIB_COLORS = {
    "Fib_Papi": "#4682B4", "Fib_SFRP2": "#87CEFA", "Fib_Myo": "#DDA0DD",
    "Fib_Inflama": "#BA55D3", "Fib_Fasci": "#9370DB", "Fib_EN1": "#F0E68C",
    "Fib_Prolif": "#FFB6C1", "Fib_K14": "#E31A1C", "Fib_CD45": "#9932CC",
}
GENES = ["PDGFRA", "VIM", "KRT14", "KRT5", "TACSTD2"]
SAMPLES = [
    ("Normal", "NS_C02847B1", 180),
    ("5 dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12 dpb DPTDI-1", "BW13_1_B1_DSDB_12dpb", 45),
    ("12 dpb SPTDI-1", "BW14_1_C1_SSDB_12dpb", 135),
    ("12 dpb DPTDI-2", "BW13_A3_DSDB_12dpb", 135),
    ("12 dpb SPTDI-2", "BW14_B3_SSDB_12dpb", 135),
    ("19 dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("19 dpb p1", "BW15D_1_D1_SDSDB_19dpb", 140),
    ("26 dpb p1", "BW81_C02846B6_SDSDB_26dpb_part1", 0),
    ("26 dpb p2", "BW81_C02846B6_SDSDB_26dpb_part2", 0),
    ("2 mph", "2mph_A03699G6.SCT", 45),
]
TRAJECTORY_STATES = [
    "Fib_K14", "SAC_SG_Progenitor", "KC_Basal", "KC_Basal_Mig",
    "KC_Spinous_Mig", "KC_Spinous",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--scrna", required=True)
    p.add_argument("--stereo", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--old-d-table", required=True)
    p.add_argument("--skip-f", action="store_true")
    p.add_argument("--only-f", action="store_true")
    p.add_argument("--skip-trajectory-export", action="store_true")
    return p.parse_args()


def style() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.linewidth": 0.7,
    })
    sns.set_style("white")


def save(fig: plt.Figure, path: Path, dpi: int = 600) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".png"), dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def dense(x) -> np.ndarray:
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def weighted_density(xy: np.ndarray, weights: np.ndarray, bins: int = 220,
                     sigma: float = 3.0) -> np.ndarray:
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    if weights.max() > 0:
        weights = weights / weights.max()
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    xi = np.clip(((xy[:, 0] - xmin) / max(xmax - xmin, 1e-9) * (bins - 1)).astype(int), 0, bins - 1)
    yi = np.clip(((xy[:, 1] - ymin) / max(ymax - ymin, 1e-9) * (bins - 1)).astype(int), 0, bins - 1)
    weighted = np.zeros((bins, bins), dtype=float)
    occupancy = np.zeros((bins, bins), dtype=float)
    np.add.at(weighted, (yi, xi), weights)
    np.add.at(occupancy, (yi, xi), 1.0)
    sm_w = gaussian_filter(weighted, sigma=sigma)
    sm_n = gaussian_filter(occupancy, sigma=sigma)
    den = np.divide(sm_w, sm_n, out=np.zeros_like(sm_w), where=sm_n > 1e-5)
    return den[yi, xi]


def rotate(xy: np.ndarray, angle: float) -> np.ndarray:
    theta = np.deg2rad(angle)
    center = np.median(xy, axis=0)
    r = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    return (xy - center) @ r.T


def panel_a_c_and_export_trajectory(adata: ad.AnnData, figures: Path,
                                    tables: Path, skip_trajectory: bool = False) -> None:
    labels = adata.obs["sub_labels"].astype(str)
    fib_mask = labels.str.startswith("Fib_").to_numpy()
    fib_idx = np.flatnonzero(fib_mask)
    # Rebuild the graph and UMAP using only fibroblasts in the registered
    # Harmony space. This is a fibroblast-only re-clustering while avoiding a
    # second dense copy of the 52,399-gene count matrix.
    harmony = np.asarray(adata.obsm["X_harmony"])[fib_idx, :40]
    fib_obs = adata.obs.iloc[fib_idx].copy()
    fib_obs["sub_labels"] = pd.Categorical(
        fib_obs["sub_labels"].astype(str), categories=FIB_ORDER,
    )
    fib = ad.AnnData(X=harmony, obs=fib_obs)
    sc.pp.neighbors(fib, n_neighbors=20, use_rep="X", random_state=17)
    sc.tl.umap(fib, min_dist=0.35, spread=1.0, random_state=17)
    sc.tl.leiden(fib, resolution=0.8, key_added="fib_leiden", random_state=17)
    coords = np.asarray(fib.obsm["X_umap"])
    expr = dense(adata[fib_idx, GENES].X)
    libsize = adata.obs.iloc[fib_idx]["nCount_RNA"].to_numpy(dtype=float)
    expr = np.log1p(expr / np.maximum(libsize[:, None], 1.0) * 1e4)
    fib.obs[["sub_labels", "sample_id", "time_point", "fib_leiden"]].assign(
        UMAP1=coords[:, 0], UMAP2=coords[:, 1],
    ).to_csv(tables / "S4A_fibroblast_reclustered_coordinates.csv.gz")

    fig, ax = plt.subplots(figsize=(3.85, 2.75))
    for state in FIB_ORDER:
        m = fib.obs["sub_labels"].astype(str).eq(state).to_numpy()
        if m.any():
            ax.scatter(coords[m, 0], coords[m, 1], s=0.55,
                       c=FIB_COLORS[state], label=state, linewidths=0)
    ax.set(title="Fibroblast-only re-clustering", xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(frameon=False, ncol=3, markerscale=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.12))
    sns.despine(ax=ax)
    save(fig, figures / "FigureS4_A_fibroblast_reclustered_umap")

    # Feature plots use the log-normalized values captured before HVG scaling.
    fig, axes = plt.subplots(1, 6, figsize=(7.08, 1.38))
    for j, (ax, gene) in enumerate(zip(axes[:5], GENES)):
        values = expr[:, j]
        vmax = np.quantile(values[values > 0], 0.99) if np.any(values > 0) else 1
        order = np.argsort(values)
        sca = ax.scatter(coords[order, 0], coords[order, 1], c=values[order],
                         s=0.38, cmap="viridis", vmin=0, vmax=vmax, linewidths=0)
        ax.set_title(gene, fontstyle="italic")
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        if j == 4:
            cb = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.01)
            cb.ax.set_title("Expr.", fontsize=6, pad=2)
            cb.ax.tick_params(labelsize=6)
    pdg = expr[:, GENES.index("PDGFRA")]
    k14 = expr[:, GENES.index("KRT14")]
    pdg_n = pdg / max(np.quantile(pdg[pdg > 0], 0.99), 1e-9) if np.any(pdg > 0) else pdg
    k14_n = k14 / max(np.quantile(k14[k14 > 0], 0.99), 1e-9) if np.any(k14 > 0) else k14
    joint_weight = np.sqrt(np.clip(pdg_n, 0, 1) * np.clip(k14_n, 0, 1))
    joint_density = weighted_density(coords, joint_weight, bins=240, sigma=3.2)
    order = np.argsort(joint_density)
    sca = axes[5].scatter(coords[order, 0], coords[order, 1],
                          c=joint_density[order], s=0.38, cmap="viridis",
                          linewidths=0)
    axes[5].set_title("PDGFRA + KRT14\ndensity")
    axes[5].set_aspect("equal")
    axes[5].set_xticks([]); axes[5].set_yticks([])
    cb = fig.colorbar(sca, ax=axes[5], fraction=0.045, pad=0.01)
    cb.ax.set_title("Density", fontsize=6, pad=2)
    cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(wspace=0.12)
    save(fig, figures / "FigureS4_C_marker_featureplots")

    # Export a balanced cell set for Monocle2. All six Figure 3B states are
    # represented in one model; per-state sampling avoids domination by KC_Basal.
    if skip_trajectory:
        return
    all_labels = adata.obs["sub_labels"].astype(str)
    rng = np.random.default_rng(17)
    selected: list[int] = []
    for state in TRAJECTORY_STATES:
        idx = np.flatnonzero(all_labels.eq(state).to_numpy())
        selected.extend(rng.choice(idx, min(2500, len(idx)), replace=False))
    selected = np.asarray(selected)
    traj = adata[selected].to_memory() if adata.isbacked else adata[selected].copy()
    sc.pp.filter_genes(traj, min_cells=20)
    counts = traj.X.tocsr() if sparse.issparse(traj.X) else sparse.csr_matrix(traj.X)
    sparse.save_npz(tables / "S4E_monocle2_counts_cells_by_genes.npz", counts)
    mmwrite(tables / "S4E_monocle2_counts_cells_by_genes.mtx", counts)
    pd.DataFrame({"gene": traj.var_names.astype(str)}).to_csv(
        tables / "S4E_monocle2_genes.csv", index=False,
    )
    traj.obs.assign(cell_id=traj.obs_names.astype(str)).to_csv(
        tables / "S4E_monocle2_phenodata.csv", index=False,
    )


def panel_d(old_table: Path, figures: Path, tables: Path) -> None:
    df = pd.read_csv(old_table)
    df.to_csv(tables / "S4D_sample_level_proxy.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.08, 1.85))
    order = [x for x in FIB_ORDER if x in set(df["subcluster"])]
    sns.boxplot(data=df, x="subcluster", y="KC_proxy", order=order, ax=ax,
                showfliers=False, color="white", linewidth=0.7)
    sns.stripplot(data=df, x="subcluster", y="KC_proxy", order=order, ax=ax,
                  size=2, palette=[FIB_COLORS[x] for x in order], jitter=0.18)
    ax.set(title="Sample-level transcriptional proximity to keratinocyte states",
           ylabel="KC proxy score", xlabel="")
    ax.tick_params(axis="x", rotation=30)
    sns.despine(ax=ax)
    save(fig, figures / "FigureS4_D_proxy_statistics")


def panel_f(stereo_path: str, figures: Path, tables: Path) -> None:
    stereo = ad.read_h5ad(stereo_path, backed="r")
    available = set(stereo.obs["sample_batch_new"].astype(str).unique())
    rows = []
    fig, axes = plt.subplots(3, 4, figsize=(7.08, 4.95))
    axes = axes.ravel()
    for ax, (label, sample_id, angle) in zip(axes, SAMPLES):
        if sample_id not in available:
            raise KeyError(f"Registered section missing: {sample_id}")
        mask = stereo.obs["sample_batch_new"].astype(str).eq(sample_id).to_numpy()
        idx = np.flatnonzero(mask)
        xy = np.asarray(stereo.obsm["spatial"])[idx]
        xy = rotate(xy, angle)
        # The backed Raw view in anndata 0.12 ignores the gene slice for this
        # legacy file; use the registered log-normalized X matrix explicitly.
        values = dense(stereo[idx, ["PDGFRA", "KRT14"]].X)
        pdg, k14 = values[:, 0], values[:, 1]
        pdg_n = pdg / max(np.quantile(pdg[pdg > 0], 0.99), 1e-9) if np.any(pdg > 0) else pdg
        k14_n = k14 / max(np.quantile(k14[k14 > 0], 0.99), 1e-9) if np.any(k14 > 0) else k14
        joint = np.sqrt(np.clip(pdg_n, 0, 1) * np.clip(k14_n, 0, 1))
        density = weighted_density(xy, joint, bins=260, sigma=2.8)
        order = np.argsort(density)
        sca = ax.scatter(xy[order, 0], xy[order, 1], c=density[order],
                         s=0.10, cmap="viridis", linewidths=0, rasterized=True,
                         norm=Normalize(0, max(np.quantile(density, 0.995), 1e-6)))
        ax.set_title(label, fontsize=7, pad=2)
        ax.set_aspect("equal")
        ax.axis("off")
        rows.append(pd.DataFrame({
            "section": label, "sample_id": sample_id,
            "x_rotated": xy[:, 0], "y_rotated": xy[:, 1],
            "PDGFRA": pdg, "KRT14": k14, "joint_density": density,
        }))
    for ax in axes[len(SAMPLES):]:
        ax.axis("off")
    cax = fig.add_axes([0.76, 0.075, 0.18, 0.018])
    cb = fig.colorbar(sca, cax=cax, orientation="horizontal")
    cb.set_label("PDGFRA + KRT14 joint spatial density", fontsize=6)
    cb.set_ticks([])
    fig.subplots_adjust(left=0.02, right=0.98, top=0.96, bottom=0.12,
                        wspace=0.08, hspace=0.15)
    pd.concat(rows, ignore_index=True).to_csv(
        tables / "S4F_human_spatial_joint_density.csv.gz", index=False,
    )
    save(fig, figures / "FigureS4_F_all_human_spatial_joint_density")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    figures = outdir / "figures" / "panels"
    tables = outdir / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    style()
    if not args.only_f:
        scrna = ad.read_h5ad(args.scrna, backed="r")
        panel_a_c_and_export_trajectory(
            scrna, figures, tables, skip_trajectory=args.skip_trajectory_export
        )
        panel_d(Path(args.old_d_table), figures, tables)
        del scrna
    if not args.skip_f:
        panel_f(args.stereo, figures, tables)


if __name__ == "__main__":
    main()
