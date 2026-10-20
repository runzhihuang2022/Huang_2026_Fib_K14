"""Re-render Figure S4/S5 with the matched A4 layout requested on 2026-07-25.

The script reuses audited derived tables for doublets and Monocle2 coordinates,
but reads the registered H5AD objects for the added COL1A1/Col1a1 feature plot
and for sample/section-level Fib_K14 statistics. Biological samples or spatial
sections, never cells, are the statistical replicates.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import anndata as ad
import h5py
try:
    import fitz
except ImportError:
    fitz = None
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.ndimage import gaussian_filter
from PIL import Image
import openslide
from skimage import measure, morphology
from skimage.filters import threshold_otsu


A4_WIDTH_PT = 210 / 25.4 * 72
A4_HEIGHT_PT = 297 / 25.4 * 72
FIB_COLORS = {
    "Fib_Papi": "#1B9E77",
    "Fib_SFRP2": "#377EB8",
    "Fib_Sfrp2": "#377EB8",
    "Fib_Myo": "#984EA3",
    "Fib_Inflama": "#FF7F00",
    "Fib_Fasci": "#A65628",
    "Fib_EN1": "#F4D03F",
    "Fib_En1": "#F4D03F",
    "Fib_Prolif": "#7F8C8D",
    "Fib_CD45": "#00A6A6",
    "Fib_Cd45": "#00A6A6",
    "Fib_K14": "#E31A1C",
}
STATE_COLORS = {
    "Fib_K14": "#E31A1C",
    "SAC_SG_Progenitor": "#00A6D6",
    "SAC_Progenitor": "#00A6D6",
    "KC_Basal": "#009E73",
    "KC_Basal_Mig": "#E69F00",
    "KC_Spinous_Mig": "#1F1A9E",
    "KC_Spinous": "#18C83E",
}
DOUBLE_POSITIVE_COLOR = "#CC3377"
MATCHED_CATEGORY_COLORS = {
    "Normal": "#4D4D4D",
    "SPTDI": "#E6AB02",
    "DPTDI": "#D62728",
    "Embryonic": "#4D4D4D",
    "P5": "#E6AB02",
    "Adult": "#D62728",
}
STAINING_MARKER_COLORS = {
    "VIM": "#00A651",
    "KRT14": "#E64B35",
    "KRT5": "#E64B35",
    "DAPI": "#2F55D4",
}
HUMAN_TIME_ORDER = ["Normal", "0–7 dpb", "8–14 dpb", "15–28 dpb", "1–2 mph"]
RAT_TIME_ORDER = ["Normal", "0–7 dpi", "8–14 dpi", "15–28 dpi"]
TIME_COLORS = {
    "Normal": "#4D4D4D",
    "0–7 dpb": "#377EB8",
    "8–14 dpb": "#4DAF4A",
    "15–28 dpb": "#984EA3",
    "1–2 mph": "#FF7F00",
    "0–7 dpi": "#377EB8",
    "8–14 dpi": "#4DAF4A",
    "15–28 dpi": "#984EA3",
}
HUMAN_SECTIONS = [
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", choices=["human", "rat"], required=True)
    parser.add_argument("--package", required=True)
    parser.add_argument("--scrna", required=True)
    parser.add_argument("--spatial", required=True)
    parser.add_argument("--skip-assemble", action="store_true")
    parser.add_argument(
        "--panels",
        default="ABCDEF",
        help="Subset of panels to rebuild, for example CD.",
    )
    return parser.parse_args()


def setup() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.7,
        }
    )
    sns.set_style("white")


def save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".png"),
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def dense(matrix) -> np.ndarray:
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)


def h5ad_csr_columns(path: str, genes: list[str], chunk_rows: int = 20000) -> np.ndarray:
    """Read selected columns from an on-disk CSR matrix without materializing X."""
    with h5py.File(path, "r") as handle:
        var_names = [
            value.decode("utf-8") if isinstance(value, bytes) else str(value)
            for value in handle["var"]["_index"][:]
        ]
        gene_positions = []
        for gene in genes:
            if gene not in var_names:
                raise KeyError(f"{gene} is absent from {path}")
            gene_positions.append(var_names.index(gene))
        matrix = handle["X"]
        shape = tuple(int(value) for value in matrix.attrs["shape"])
        output = np.zeros((shape[0], len(genes)), dtype=np.float32)
        indptr = matrix["indptr"][:]
        for row_start in range(0, shape[0], chunk_rows):
            row_end = min(row_start + chunk_rows, shape[0])
            data_start = int(indptr[row_start])
            data_end = int(indptr[row_end])
            indices = matrix["indices"][data_start:data_end]
            data = matrix["data"][data_start:data_end]
            counts = np.diff(indptr[row_start : row_end + 1])
            rows = np.repeat(np.arange(row_start, row_end), counts)
            for column, gene_position in enumerate(gene_positions):
                selected = indices == gene_position
                output[rows[selected], column] = data[selected]
    return output


def h5ad_categorical_obs(path: str, column: str) -> pd.Series:
    with h5py.File(path, "r") as handle:
        group = handle["obs"][column]
        categories = np.asarray(
            [
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in group["categories"][:]
            ],
            dtype=object,
        )
        codes = group["codes"][:]
    return pd.Series(categories[codes], dtype="object")


def h5ad_obsm_array(path: str, key: str) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        return handle["obsm"][key][:]


def rotate(xy: np.ndarray, angle: float) -> np.ndarray:
    theta = np.deg2rad(angle)
    center = np.median(xy, axis=0)
    rotation = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )
    return (xy - center) @ rotation.T


def weighted_density(
    xy: np.ndarray, weights: np.ndarray, bins: int = 220, sigma: float = 3.0
) -> np.ndarray:
    weights = np.nan_to_num(weights, nan=0.0, posinf=0.0, neginf=0.0)
    if np.max(weights) > 0:
        weights = weights / np.max(weights)
    low = np.min(xy, axis=0)
    high = np.max(xy, axis=0)
    scaled = (xy - low) / np.maximum(high - low, 1e-9)
    ij = np.clip((scaled * (bins - 1)).astype(int), 0, bins - 1)
    weighted = np.zeros((bins, bins), dtype=float)
    occupied = np.zeros((bins, bins), dtype=float)
    np.add.at(weighted, (ij[:, 1], ij[:, 0]), weights)
    np.add.at(occupied, (ij[:, 1], ij[:, 0]), 1.0)
    smooth_weight = gaussian_filter(weighted, sigma=sigma)
    smooth_count = gaussian_filter(occupied, sigma=sigma)
    grid = np.divide(
        smooth_weight,
        smooth_count,
        out=np.zeros_like(smooth_weight),
        where=smooth_count > 1e-5,
    )
    return grid[ij[:, 1], ij[:, 0]]


def cell_ids_and_coordinates(package: Path, species: str) -> pd.DataFrame:
    filename = (
        "S4A_fibroblast_reclustered_coordinates.csv.gz"
        if species == "human"
        else "S5A_rat_fibroblast_reclustered_coordinates.csv.gz"
    )
    frame = pd.read_csv(package / "tables" / filename)
    frame = frame.rename(columns={"Unnamed: 0": "cell_id"})
    return frame


def render_a(package: Path, species: str) -> None:
    frame = cell_ids_and_coordinates(package, species)
    order = [state for state in FIB_COLORS if state in set(frame["sub_labels"])]
    fig, ax = plt.subplots(figsize=(3.38, 2.35))
    for state in order:
        current = frame["sub_labels"].astype(str).eq(state).to_numpy()
        ax.scatter(
            frame.loc[current, "UMAP1"],
            frame.loc[current, "UMAP2"],
            s=0.45,
            color=FIB_COLORS[state],
            linewidths=0,
            rasterized=True,
            label="Fib_K14-like" if species == "rat" and state == "Fib_K14" else state,
        )
    ax.set(title="Fibroblast-only re-clustering", xlabel="UMAP1", ylabel="UMAP2")
    ax.legend(
        frameon=False,
        ncol=3,
        markerscale=5,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        columnspacing=0.8,
        handletextpad=0.3,
    )
    sns.despine(ax=ax)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.91, bottom=0.36)
    prefix = "FigureS4" if species == "human" else "FigureS5"
    save(fig, package / "figures" / "panels" / f"{prefix}_A_fibroblast_umap_v02")


def render_b(package: Path, species: str) -> None:
    filename = (
        "S4B_QC_doublet_plot_data.csv.gz"
        if species == "human"
        else "S5B_QC_doublet_plot_data.csv.gz"
    )
    frame = pd.read_csv(package / "tables" / filename)
    if species == "human":
        aliases = {
            "Keratinocyte": "Keratinocytes",
            "Other Fib": "Other fibroblasts",
            "Fib_K14": "Fib_K14",
        }
    else:
        aliases = {
            "Keratinocyte": "Keratinocytes",
            "Other Fib": "Other fibroblasts",
            "Fib_K14-like": "Fib_K14-like",
        }
    frame["Display group"] = frame["Group"].map(aliases).fillna(frame["Group"])
    order = list(aliases.values())
    palette = {
        order[0]: "#BDBDBD",
        order[1]: "#377EB8",
        order[2]: "#E31A1C",
    }
    rng = np.random.default_rng(17)
    selected: list[int] = []
    for group in order:
        indices = np.flatnonzero(frame["Display group"].eq(group).to_numpy())
        selected.extend(rng.choice(indices, min(3000, len(indices)), replace=False))
    plot = frame.iloc[selected]
    fig, axes = plt.subplots(2, 2, figsize=(3.38, 2.35))
    metrics = [
        "Number of Gene",
        "Number of UMI",
        "Mitochondrial (%)",
        "Doublet score",
    ]
    for ax, metric in zip(axes.ravel(), metrics):
        sns.violinplot(
            data=plot,
            x="Display group",
            y=metric,
            order=order,
            hue="Display group",
            palette=palette,
            legend=False,
            cut=0,
            inner="quartile",
            linewidth=0.55,
            ax=ax,
        )
        ax.set_title(metric, pad=2)
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=18, pad=1)
        for label in ax.get_xticklabels():
            label.set_ha("right")
        sns.despine(ax=ax)
    fractions = (
        frame.groupby("Display group")["Predicted doublet"].mean().reindex(order)
    )
    axes[1, 1].text(
        0.98,
        0.98,
        "Predicted doublets\n"
        + "\n".join(
            f"{name}: {value:.1%}" for name, value in fractions.items()
        ),
        transform=axes[1, 1].transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1},
    )
    fig.subplots_adjust(
        left=0.15, right=0.98, top=0.92, bottom=0.25, wspace=0.48, hspace=0.70
    )
    prefix = "FigureS4" if species == "human" else "FigureS5"
    save(fig, package / "figures" / "panels" / f"{prefix}_B_qc_doublets_v02")


def normalized_expression(
    adata: ad.AnnData,
    cell_ids: list[str],
    genes: list[str],
    total_column: str,
) -> np.ndarray:
    positions = adata.obs_names.get_indexer(cell_ids)
    if np.any(positions < 0):
        raise KeyError(f"{np.sum(positions < 0)} fibroblast cells missing from source")
    present = [gene for gene in genes if gene in adata.var_names]
    expression = np.zeros((len(cell_ids), len(genes)), dtype=np.float32)
    if present:
        raw = dense(adata[positions, present].X).astype(np.float32)
        library = adata.obs.iloc[positions][total_column].to_numpy(float)
        raw = np.log1p(raw / np.maximum(library[:, None], 1.0) * 1e4)
        for column, gene in enumerate(present):
            expression[:, genes.index(gene)] = raw[:, column]
    return expression


def render_c(package: Path, species: str, scrna_path: str) -> None:
    coordinates = cell_ids_and_coordinates(package, species)
    if species == "human":
        genes = ["PDGFRA", "VIM", "KRT14", "KRT5", "TACSTD2", "COL1A1"]
        total_column = "nCount_RNA"
    else:
        genes = ["Pdgfra", "Vim", "Krt5", "Tacstd2", "Col1a1"]
        total_column = "nCount_RNA.x"
    adata = ad.read_h5ad(scrna_path, backed="r")
    expression = normalized_expression(
        adata,
        coordinates["cell_id"].astype(str).tolist(),
        genes,
        total_column,
    )
    adata.file.close()
    output = coordinates[["cell_id", "sub_labels", "UMAP1", "UMAP2"]].copy()
    for index, gene in enumerate(genes):
        output[gene] = expression[:, index]
    table_name = (
        "S4C_marker_log_normalized_expression_v02.csv.gz"
        if species == "human"
        else "S5C_marker_log_normalized_expression_v02.csv.gz"
    )
    output.to_csv(package / "tables" / table_name, index=False)

    xy = coordinates[["UMAP1", "UMAP2"]].to_numpy(float)
    fig, axes = plt.subplots(1, len(genes) + 1, figsize=(7.08, 1.35))
    last_scatter = None
    measurable_columns = []
    for index, (ax, gene) in enumerate(zip(axes[: len(genes)], genes)):
        values = expression[:, index]
        measurable_columns.append(index)
        positive = values[values > 0]
        vmax = np.quantile(positive, 0.99) if len(positive) else 1.0
        sort_order = np.argsort(values)
        last_scatter = ax.scatter(
            xy[sort_order, 0],
            xy[sort_order, 1],
            c=values[sort_order],
            s=0.22,
            cmap="viridis",
            vmin=0,
            vmax=max(vmax, 1e-8),
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(gene, fontstyle="italic", pad=2)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    scaled = []
    for index in measurable_columns:
        values = expression[:, index]
        positive = values[values > 0]
        scale = np.quantile(positive, 0.99) if len(positive) else 1.0
        scaled.append(np.clip(values / max(scale, 1e-9), 0, 1))
    merged_weight = np.exp(
        np.mean(np.log(np.maximum(np.vstack(scaled), 1e-4)), axis=0)
    )
    merged_density = weighted_density(xy, merged_weight, bins=240, sigma=3.0)
    sort_order = np.argsort(merged_density)
    merged_axis = axes[-1]
    merged_scatter = merged_axis.scatter(
        xy[sort_order, 0],
        xy[sort_order, 1],
        c=merged_density[sort_order],
        s=0.22,
        cmap="viridis",
        linewidths=0,
        rasterized=True,
    )
    merged_label = "Six-gene\nmerged density" if species == "human" else "Five-gene\nmerged density"
    merged_axis.set_title(merged_label, pad=2)
    merged_axis.set_aspect("equal")
    merged_axis.set_xticks([])
    merged_axis.set_yticks([])
    for spine in merged_axis.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(
        merged_scatter,
        ax=axes,
        orientation="horizontal",
        fraction=0.05,
        pad=0.12,
        aspect=45,
    )
    colorbar.set_label("Expression / merged density", fontsize=8.0)
    colorbar.ax.tick_params(labelsize=8.0)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.84, bottom=0.23, wspace=0.12)
    prefix = "FigureS4" if species == "human" else "FigureS5"
    save(fig, package / "figures" / "panels" / f"{prefix}_C_markers_merged_v02")


def render_d(package: Path, species: str) -> None:
    prefix = "FigureS4" if species == "human" else "FigureS5"
    if species == "human":
        frame = pd.read_csv(package / "tables" / "S4D_sample_level_proxy.csv")
        order = [
            state
            for state in FIB_COLORS
            if state in set(frame["subcluster"].astype(str))
        ]
        fig, ax = plt.subplots(figsize=(2.25, 2.15))
        sns.boxplot(
            data=frame,
            y="subcluster",
            x="KC_proxy",
            order=order,
            showfliers=False,
            color="white",
            linewidth=0.65,
            ax=ax,
        )
        sns.stripplot(
            data=frame,
            y="subcluster",
            x="KC_proxy",
            order=order,
            palette=[FIB_COLORS[state] for state in order],
            size=1.8,
            jitter=0.16,
            ax=ax,
        )
        ax.set(title="Keratinocyte proxy", xlabel="Sample-level proxy score", ylabel="")
        sns.despine(ax=ax)
        fig.subplots_adjust(left=0.35, right=0.98, top=0.90, bottom=0.18)
    else:
        matrix = pd.read_csv(
            package / "tables" / "S5D_rat_PAGA_connectivity.csv", index_col=0
        )
        states = list(matrix.index)
        preset = {
            "Fib_K14": (0.0, 0.0),
            "KC_Basal": (-0.72, 0.70),
            "KC_Basal_Mig": (0.00, 0.98),
            "KC_Basal_Prolif": (-0.98, 0.12),
            "KC_Spinous": (-0.68, -0.66),
            "KC_Spinous_Mig": (0.02, -0.98),
            "KC_Spinous_Mat": (0.72, -0.68),
            "KC_Granular": (0.98, 0.10),
            "SAC_Progenitor": (0.68, 0.72),
        }
        positions = {
            state: np.asarray(preset.get(state, (0.0, 0.0)), dtype=float)
            for state in states
        }
        fig, ax = plt.subplots(figsize=(2.25, 2.15))
        threshold = np.nanquantile(matrix.to_numpy()[np.triu_indices(len(states), 1)], 0.55)
        for left_index, left in enumerate(states):
            for right_index in range(left_index + 1, len(states)):
                right = states[right_index]
                weight = float(matrix.loc[left, right])
                if weight >= threshold:
                    ax.plot(
                        [positions[left][0], positions[right][0]],
                        [positions[left][1], positions[right][1]],
                        color="#555555",
                        lw=0.4 + 2.2 * weight,
                        alpha=0.72,
                        zorder=1,
                    )
        for state in states:
            color = (
                FIB_COLORS["Fib_K14"]
                if state == "Fib_K14"
                else STATE_COLORS.get(state, "#999999")
            )
            ax.scatter(
                *positions[state],
                s=48,
                color=color,
                edgecolor="black",
                linewidth=0.45,
                zorder=2,
            )
            label = "Fib_K14-like" if state == "Fib_K14" else state
            vector = positions[state]
            if state == "Fib_K14":
                text_position = vector + np.array([0.06, 0.08])
                horizontal = "left"
            else:
                text_position = vector * 1.17
                horizontal = "left" if vector[0] >= 0 else "right"
            ax.text(
                text_position[0],
                text_position[1],
                label,
                ha=horizontal,
                va="center",
                fontsize=7.9,
            )
        ax.set_title("Fib–KC PAGA connectivity", pad=2)
        ax.set_aspect("equal")
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.25, 1.30)
        ax.axis("off")
        fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.02)
    save(fig, package / "figures" / "panels" / f"{prefix}_D_compact_v02")


def normalized_time(value: str, species: str) -> str:
    value = str(value).replace("dpi", "dpb" if species == "human" else "dpi")
    mapping = {
        "0-7dpb": "0–7 dpb",
        "8-14dpb": "8–14 dpb",
        "15-28dpb": "15–28 dpb",
        "0-7dpi": "0–7 dpi",
        "8-14dpi": "8–14 dpi",
        "15-28dpi": "15–28 dpi",
        "1-2mph": "1–2 mph",
        "Normal": "Normal",
    }
    return mapping.get(value, value)


def add_trajectory_arrow(ax: plt.Axes, frame: pd.DataFrame) -> None:
    ordered = frame.sort_values("Pseudotime")
    start = ordered.iloc[int(0.18 * (len(ordered) - 1))]
    end = ordered.iloc[int(0.58 * (len(ordered) - 1))]
    arrow = FancyArrowPatch(
        (start["Component_1"], start["Component_2"]),
        (end["Component_1"], end["Component_2"]),
        arrowstyle="-|>",
        mutation_scale=9,
        lw=1.0,
        color="black",
        zorder=5,
    )
    ax.add_patch(arrow)


def render_e(package: Path, species: str) -> None:
    prefix = "S4" if species == "human" else "S5"
    frame = pd.read_csv(
        package / "tables" / f"{prefix}E_monocle2_coordinates_pseudotime.csv"
    )
    frame["Observed healing time"] = frame["Time_category"].map(
        lambda value: normalized_time(value, species)
    )
    state_order = [
        state for state in STATE_COLORS if state in set(frame["sub_labels"].astype(str))
    ]
    time_order = HUMAN_TIME_ORDER if species == "human" else RAT_TIME_ORDER
    fig, axes = plt.subplots(1, 3, figsize=(4.65, 2.15))
    rng = np.random.default_rng(17)
    if len(frame) > 15000:
        draw = frame.iloc[rng.choice(len(frame), 15000, replace=False)]
    else:
        draw = frame
    for state in state_order:
        current = draw["sub_labels"].astype(str).eq(state)
        axes[0].scatter(
            draw.loc[current, "Component_1"],
            draw.loc[current, "Component_2"],
            s=0.7,
            color=STATE_COLORS[state],
            linewidths=0,
            rasterized=True,
        )
    pseudotime = axes[1].scatter(
        draw["Component_1"],
        draw["Component_2"],
        s=0.7,
        c=draw["Pseudotime"],
        cmap="viridis",
        linewidths=0,
        rasterized=True,
    )
    for time in time_order:
        current = draw["Observed healing time"].eq(time)
        if current.any():
            axes[2].scatter(
                draw.loc[current, "Component_1"],
                draw.loc[current, "Component_2"],
                s=0.7,
                color=TIME_COLORS[time],
                linewidths=0,
                rasterized=True,
            )
    titles = ["Cell state", "Pseudotime", "Observed healing time"]
    for ax, title in zip(axes, titles):
        ax.set_title(title, pad=2)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        add_trajectory_arrow(ax, frame)
    cbar = fig.colorbar(
        pseudotime,
        ax=axes[1],
        orientation="horizontal",
        fraction=0.045,
        pad=0.04,
    )
    cbar.ax.tick_params(labelsize=6.5)
    state_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=STATE_COLORS[state],
            markeredgecolor="none",
            markersize=3.2,
            label="Fib_K14-like" if species == "rat" and state == "Fib_K14" else state,
        )
        for state in state_order
    ]
    time_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=TIME_COLORS[time],
            markeredgecolor="none",
            markersize=3.2,
            label=time,
        )
        for time in time_order
        if time in set(frame["Observed healing time"])
    ]
    fig.legend(
        handles=state_handles + time_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=5,
        frameon=False,
        columnspacing=0.7,
        handletextpad=0.25,
        fontsize=6.5,
    )
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.25, wspace=0.08)
    panel_prefix = "FigureS4" if species == "human" else "FigureS5"
    save(
        fig,
        package / "figures" / "panels" / f"{panel_prefix}_E_monocle_triplet_v02",
    )


def add_scale_bar(
    ax: plt.Axes, xy: np.ndarray, length: float = 1000, label: str = "1 mm"
) -> None:
    xmin, ymin = np.min(xy, axis=0)
    xmax, ymax = np.max(xy, axis=0)
    x0 = xmin + 0.06 * (xmax - xmin)
    y0 = ymin + 0.08 * (ymax - ymin)
    available = xmax - x0
    actual = min(length, 0.28 * (xmax - xmin), available)
    ax.plot([x0, x0 + actual], [y0, y0], color="black", lw=1.1)
    ax.text(
        x0 + actual / 2,
        y0 + 0.035 * (ymax - ymin),
        label,
        ha="center",
        va="bottom",
        fontsize=7.0,
    )


def sample_fib_k14_fraction(scrna_path: str, species: str) -> pd.DataFrame:
    adata = ad.read_h5ad(scrna_path, backed="r")
    obs = adata.obs[["sub_labels", "sample_id", "Time_category", "Grade"]].copy()
    labels = obs["sub_labels"].astype(str)
    obs["is_fib"] = labels.str.startswith("Fib_")
    obs["is_fib_k14"] = labels.eq("Fib_K14")
    rows = []
    for sample_id, current in obs.groupby("sample_id", observed=True):
        fibroblasts = int(current["is_fib"].sum())
        if fibroblasts == 0:
            continue
        rows.append(
            {
                "sample_id": sample_id,
                "Time": normalized_time(current["Time_category"].iloc[0], species),
                "Group": str(current["Grade"].iloc[0]),
                "Fib_K14_fraction": float(current["is_fib_k14"].sum() / fibroblasts),
                "n_fibroblasts": fibroblasts,
            }
        )
    adata.file.close()
    return pd.DataFrame(rows)


def human_f(package: Path, scrna_path: str, spatial_path: str) -> None:
    print("human F: read spatial metadata", flush=True)
    stereo = ad.read_h5ad(spatial_path, backed="r")
    sample_column = (
        "sample_batch_new"
        if "sample_batch_new" in stereo.obs
        else "sample_id"
    )
    available = set(stereo.obs[sample_column].astype(str))
    fig = plt.figure(figsize=(7.08, 3.35))
    grid = fig.add_gridspec(
        3, 6, height_ratios=[1, 1, 1.42], hspace=0.18, wspace=0.08
    )
    map_summary = []
    for index, (title, sample_id, angle) in enumerate(HUMAN_SECTIONS):
        row, column = divmod(index, 6)
        ax = fig.add_subplot(grid[row, column])
        if sample_id not in available:
            ax.text(0.5, 0.5, "Section unavailable", ha="center", va="center")
            ax.axis("off")
            continue
        selected = np.flatnonzero(
            stereo.obs[sample_column].astype(str).eq(sample_id).to_numpy()
        )
        xy = rotate(np.asarray(stereo.obsm["spatial"])[selected], angle)
        values = dense(stereo[selected, ["PDGFRA", "KRT14"]].X)
        normalized = []
        for column_values in values.T:
            positive = column_values[column_values > 0]
            scale = np.quantile(positive, 0.99) if len(positive) else 1.0
            normalized.append(np.clip(column_values / max(scale, 1e-9), 0, 1))
        joint = np.sqrt(normalized[0] * normalized[1])
        density = weighted_density(xy, joint, bins=240, sigma=2.8)
        draw_order = np.argsort(density)
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=0.025,
            color="#E5E8EB",
            linewidths=0,
            rasterized=True,
        )
        scatter = ax.scatter(
            xy[draw_order, 0],
            xy[draw_order, 1],
            s=0.055,
            c=density[draw_order],
            cmap="viridis",
            norm=Normalize(0, max(np.quantile(density, 0.995), 1e-6)),
            linewidths=0,
            rasterized=True,
        )
        add_scale_bar(ax, xy)
        ax.set_title(title, fontsize=7.0, pad=1)
        ax.set_aspect("equal")
        ax.axis("off")
        map_summary.append(
            {
                "section": title,
                "sample_id": sample_id,
                "n_bins": len(selected),
                "joint_density_mean": float(np.mean(density)),
            }
        )
    empty_ax = fig.add_subplot(grid[1, 5])
    empty_ax.axis("off")
    color_ax = empty_ax.inset_axes([0.18, 0.30, 0.68, 0.20])
    colorbar = fig.colorbar(scatter, cax=color_ax, orientation="horizontal")
    colorbar.set_ticks([])
    colorbar.set_label("Merged spatial density", fontsize=7.0)

    stats_ax = fig.add_subplot(grid[2, :])
    fractions = sample_fib_k14_fraction(scrna_path, "human")
    fractions.to_csv(
        package / "tables" / "S4F_human_sample_level_FibK14_fraction_v02.csv",
        index=False,
    )
    fractions = fractions[fractions["Time"].isin(HUMAN_TIME_ORDER)]
    group_colors = {"Normal": "#222222", "SPTDI": "#D8AA00", "DPTDI": "#E31A1C"}
    x_positions = np.arange(len(HUMAN_TIME_ORDER))
    rng = np.random.default_rng(17)
    for group in ["Normal", "SPTDI", "DPTDI"]:
        means, errors = [], []
        for time in HUMAN_TIME_ORDER:
            values = fractions.loc[
                fractions["Group"].eq(group) & fractions["Time"].eq(time),
                "Fib_K14_fraction",
            ]
            means.append(values.mean() if len(values) else np.nan)
            errors.append(
                values.sem() if len(values) > 1 else (0.0 if len(values) else np.nan)
            )
            if len(values):
                stats_ax.scatter(
                    np.full(len(values), x_positions[HUMAN_TIME_ORDER.index(time)])
                    + rng.uniform(-0.06, 0.06, len(values)),
                    values * 100,
                    s=8,
                    color=group_colors[group],
                    alpha=0.60,
                    linewidths=0,
                    zorder=2,
                )
        stats_ax.errorbar(
            x_positions,
            np.asarray(means) * 100,
            yerr=np.asarray(errors) * 100,
            color=group_colors[group],
            marker="o",
            lw=1.2,
            ms=3.4,
            capsize=2,
            label=group,
            zorder=3,
        )
    stats_ax.set(
        ylabel="Fib_K14 among fibroblasts (%)",
        xlabel="Observed healing time",
        xticks=x_positions,
        xticklabels=HUMAN_TIME_ORDER,
    )
    stats_ax.legend(frameon=False, ncol=3, loc="upper left")
    stats_ax.grid(axis="y", color="#E6E6E6", lw=0.5)
    sns.despine(ax=stats_ax)
    pd.DataFrame(map_summary).to_csv(
        package / "tables" / "S4F_human_spatial_map_summary_v02.csv", index=False
    )
    stereo.file.close()
    fig.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.12)
    save(fig, package / "figures" / "panels" / "FigureS4_F_spatial_time_v02")


def rat_stage_from_filename(name: str) -> str:
    if name.startswith(("E16.5", "E17.5", "E18.5")):
        return "Embryonic"
    if name.startswith("P5"):
        return "P5"
    return "Adult"


def rat_phase_from_filename(name: str) -> str:
    lower = name.lower()
    if "nor" in lower or "control" in lower:
        return "Normal"
    match = re.search(r"_(\d+(?:\.\d+)?)(h|d)", lower)
    if not match:
        return "Unresolved"
    value = float(match.group(1))
    if match.group(2) == "h":
        value /= 24
    if value <= 1:
        return "Inflammation"
    if value <= 3:
        return "Proliferation"
    return "Remodeling"


def rat_f(package: Path, spatial_root: str) -> None:
    paths = sorted(
        glob.glob(
            str(Path(spatial_root) / "**" / "*_Annotated_Spatial.h5ad"),
            recursive=True,
        )
    )
    metadata = pd.DataFrame(
        {
            "path": paths,
            "filename": [Path(path).name for path in paths],
        }
    )
    metadata["Stage"] = metadata["filename"].map(rat_stage_from_filename)
    metadata["Phase"] = metadata["filename"].map(rat_phase_from_filename)
    representatives = []
    desired = [
        ("Embryonic", "Normal"),
        ("Embryonic", "Proliferation"),
        ("Embryonic", "Remodeling"),
        ("P5", "Normal"),
        ("P5", "Proliferation"),
        ("P5", "Remodeling"),
        ("Adult", "Normal"),
        ("Adult", "Proliferation"),
        ("Adult", "Remodeling"),
    ]
    for stage, phase in desired:
        current = metadata[
            metadata["Stage"].eq(stage) & metadata["Phase"].eq(phase)
        ]
        if len(current):
            representatives.append(current.iloc[0])
    fig = plt.figure(figsize=(7.08, 3.35))
    grid = fig.add_gridspec(
        4, 3, height_ratios=[1, 1, 1, 1.42], hspace=0.18, wspace=0.08
    )
    rows = []
    for index, record in enumerate(representatives):
        ax = fig.add_subplot(grid[index // 3, index % 3])
        adata = ad.read_h5ad(record["path"], backed="r")
        xy = np.asarray(adata.obsm["X_rot_spatial"])
        labels = adata.obs["sub_labels_gated"].astype(str).to_numpy()
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=0.025,
            color="#E5E8EB",
            linewidths=0,
            rasterized=True,
        )
        current = labels == "Fib_K14"
        ax.scatter(
            xy[current, 0],
            xy[current, 1],
            s=0.065,
            color="#E31A1C",
            linewidths=0,
            rasterized=True,
        )
        add_scale_bar(ax, xy)
        ax.set_title(f'{record["Stage"]} · {record["Phase"]}', fontsize=7.0, pad=1)
        ax.set_aspect("equal")
        ax.axis("off")
        adata.file.close()

    for path in paths:
        adata = ad.read_h5ad(path, backed="r")
        labels = adata.obs["sub_labels_gated"].astype(str)
        fibroblasts = labels.str.startswith("Fib_")
        rows.append(
            {
                "section": Path(path).name,
                "Stage": rat_stage_from_filename(Path(path).name),
                "Phase": rat_phase_from_filename(Path(path).name),
                "Fib_K14_like_fraction": float(
                    labels.eq("Fib_K14").sum() / max(fibroblasts.sum(), 1)
                ),
                "n_fibroblasts": int(fibroblasts.sum()),
            }
        )
        adata.file.close()
    fractions = pd.DataFrame(rows)
    fractions.to_csv(
        package / "tables" / "S5F_rat_section_level_FibK14_fraction_v02.csv",
        index=False,
    )
    stats_ax = fig.add_subplot(grid[3, :])
    phase_order = ["Normal", "Inflammation", "Proliferation", "Remodeling"]
    stage_colors = {"Embryonic": "#2CA02C", "P5": "#FF8C00", "Adult": "#D62728"}
    x_positions = np.arange(len(phase_order))
    rng = np.random.default_rng(17)
    for stage in ["Embryonic", "P5", "Adult"]:
        means, errors = [], []
        for phase in phase_order:
            values = fractions.loc[
                fractions["Stage"].eq(stage) & fractions["Phase"].eq(phase),
                "Fib_K14_like_fraction",
            ]
            means.append(values.mean() if len(values) else np.nan)
            errors.append(
                values.sem() if len(values) > 1 else (0.0 if len(values) else np.nan)
            )
            if len(values):
                stats_ax.scatter(
                    np.full(len(values), x_positions[phase_order.index(phase)])
                    + rng.uniform(-0.06, 0.06, len(values)),
                    values * 100,
                    s=8,
                    color=stage_colors[stage],
                    alpha=0.60,
                    linewidths=0,
                    zorder=2,
                )
        stats_ax.errorbar(
            x_positions,
            np.asarray(means) * 100,
            yerr=np.asarray(errors) * 100,
            color=stage_colors[stage],
            marker="o",
            lw=1.2,
            ms=3.4,
            capsize=2,
            label=stage,
            zorder=3,
        )
    stats_ax.set(
        ylabel="Fib_K14-like among fibroblasts (%)",
        xlabel="Repair phase",
        xticks=x_positions,
        xticklabels=phase_order,
    )
    stats_ax.legend(frameon=False, ncol=3, loc="upper left")
    stats_ax.grid(axis="y", color="#E6E6E6", lw=0.5)
    sns.despine(ax=stats_ax)
    fig.text(
        0.995,
        0.005,
        "Rat Krt14 is absent; maps use the annotated Fib_K14-like state.",
        ha="right",
        va="bottom",
        fontsize=7.0,
        color="#555555",
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.12)
    save(fig, package / "figures" / "panels" / "FigureS5_F_spatial_time_v02")


def normalized_time_clean(value: str, species: str) -> str:
    text = str(value).strip().replace("–", "-").replace("—", "-")
    if species == "human":
        text = text.replace("dpi", "dpb")
    mapping = {
        "Normal": "Normal",
        "0-7dpb": "0-7 dpb",
        "8-14dpb": "8-14 dpb",
        "15-28dpb": "15-28 dpb",
        "1-2mph": "1-2 mph",
        "0-7dpi": "0-7 dpi",
        "8-14dpi": "8-14 dpi",
        "15-28dpi": "15-28 dpi",
    }
    return mapping.get(text, text)


def clinical_group(value: str) -> str:
    text = str(value).upper()
    if "NORMAL" in text or text in {"HNS", "NS"}:
        return "Normal"
    if "SPTDI" in text or "SSDB" in text:
        return "SPTDI"
    if "DPTDI" in text or "DSDB" in text:
        return "DPTDI"
    return str(value)


def human_sample_double_positive_fraction(scrna_path: str) -> pd.DataFrame:
    adata = ad.read_h5ad(scrna_path, backed="r")
    required = ["VIM", "KRT14"]
    missing = [gene for gene in required if gene not in adata.var_names]
    if missing:
        raise KeyError(f"Missing human genes: {missing}")
    values = dense(adata[:, required].X).astype(np.float32)
    library = adata.obs["nCount_RNA"].to_numpy(float)
    expression = np.log1p(values / np.maximum(library[:, None], 1.0) * 1e4)
    obs = adata.obs[["sample_id", "Time_category", "Grade"]].copy()
    obs["double_positive"] = (expression[:, 0] > 0.5) & (expression[:, 1] > 0.5)
    rows = []
    for sample_id, current in obs.groupby("sample_id", observed=True):
        rows.append(
            {
                "sample_id": str(sample_id),
                "Time": normalized_time_clean(current["Time_category"].iloc[0], "human"),
                "Group": clinical_group(current["Grade"].iloc[0]),
                "double_positive_fraction": float(current["double_positive"].mean()),
                "n_cells": int(len(current)),
                "n_double_positive": int(current["double_positive"].sum()),
            }
        )
    adata.file.close()
    return pd.DataFrame(rows)


def human_f(package: Path, scrna_path: str, spatial_path: str) -> None:
    sample_column = "sample_batch_new"
    sample_values = h5ad_categorical_obs(spatial_path, sample_column).astype(str)
    available = set(sample_values)
    all_spatial = h5ad_obsm_array(spatial_path, "spatial")
    print("human F: read VIM/KRT14 columns", flush=True)
    all_marker_values = h5ad_csr_columns(spatial_path, ["VIM", "KRT14"])
    print("human F: render spatial maps", flush=True)
    fig = plt.figure(figsize=(7.08, 2.78))
    grid = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.45], hspace=0.22, wspace=0.10)
    map_summary = []
    for index, (title, sample_id, angle) in enumerate(HUMAN_SECTIONS):
        ax = fig.add_subplot(grid[index // 6, index % 6])
        if sample_id not in available:
            ax.text(0.5, 0.5, "Unavailable", ha="center", va="center")
            ax.axis("off")
            continue
        selected = np.flatnonzero(sample_values.eq(sample_id).to_numpy())
        xy = rotate(all_spatial[selected], angle)
        values = all_marker_values[selected]
        double_positive = (values[:, 0] > 0.5) & (values[:, 1] > 0.5)
        ax.scatter(
            xy[:, 0], xy[:, 1], s=0.028, color="#E2E6EA",
            linewidths=0, rasterized=True,
        )
        ax.scatter(
            xy[double_positive, 0], xy[double_positive, 1], s=0.050,
            color=DOUBLE_POSITIVE_COLOR, linewidths=0, rasterized=True,
        )
        add_scale_bar(ax, xy)
        ax.set_title(title, fontsize=7.0, pad=1)
        ax.set_aspect("equal")
        ax.axis("off")
        map_summary.append(
            {
                "section": title,
                "sample_id": sample_id,
                "n_bins": int(len(selected)),
                "n_double_positive": int(double_positive.sum()),
                "double_positive_fraction": float(double_positive.mean()),
                "threshold": "VIM > 0.5 and KRT14 > 0.5 in stored spatial expression",
            }
        )
    legend_ax = fig.add_subplot(grid[1, 5])
    legend_ax.axis("off")
    legend_ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#E2E6EA", markeredgecolor="none", label="Other bins"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=DOUBLE_POSITIVE_COLOR, markeredgecolor="none", label="VIM+KRT14 double-positive"),
        ],
        frameon=False, loc="center", fontsize=6.2,
    )
    stats_ax = fig.add_subplot(grid[2, :])
    print("human F: compute scRNA sample fractions", flush=True)
    fractions = human_sample_double_positive_fraction(scrna_path)
    fractions.to_csv(package / "tables" / "S4F_human_VIM_KRT14_double_positive_samples_v03.csv", index=False)
    pd.DataFrame(map_summary).to_csv(
        package / "tables" / "S4F_human_VIM_KRT14_double_positive_spatial_v03.csv", index=False
    )
    time_order = ["Normal", "0-7 dpb", "8-14 dpb", "15-28 dpb", "1-2 mph"]
    fractions = fractions[fractions["Time"].isin(time_order)]
    x_positions = np.arange(len(time_order))
    rng = np.random.default_rng(17)
    for group in ["Normal", "SPTDI", "DPTDI"]:
        means, errors = [], []
        for time in time_order:
            current = fractions.loc[
                fractions["Group"].eq(group) & fractions["Time"].eq(time),
                "double_positive_fraction",
            ]
            means.append(current.mean() if len(current) else np.nan)
            errors.append(current.sem() if len(current) > 1 else (0.0 if len(current) else np.nan))
            if len(current):
                stats_ax.scatter(
                    np.full(len(current), x_positions[time_order.index(time)]) + rng.uniform(-0.055, 0.055, len(current)),
                    current * 100, s=8, color=MATCHED_CATEGORY_COLORS[group],
                    alpha=0.62, linewidths=0, zorder=2,
                )
        stats_ax.errorbar(
            x_positions, np.asarray(means) * 100, yerr=np.asarray(errors) * 100,
            color=MATCHED_CATEGORY_COLORS[group], marker="o", lw=1.2, ms=3.3,
            capsize=2, label=group, zorder=3,
        )
    stats_ax.set(
        ylabel="VIM+KRT14 double-positive cells (%)",
        xlabel="Observed healing time",
        xticks=x_positions,
        xticklabels=time_order,
    )
    stats_ax.legend(frameon=False, ncol=3, loc="upper left")
    stats_ax.grid(axis="y", color="#E6E6E6", lw=0.5)
    sns.despine(ax=stats_ax)
    fig.subplots_adjust(left=0.09, right=0.995, top=0.96, bottom=0.15)
    print("human F: save", flush=True)
    save(fig, package / "figures" / "panels" / "FigureS4_F_spatial_time_v03")


def rat_f(package: Path, spatial_root: str) -> None:
    paths = sorted(glob.glob(str(Path(spatial_root) / "**" / "*_Annotated_Spatial.h5ad"), recursive=True))
    metadata = pd.DataFrame({"path": paths, "filename": [Path(path).name for path in paths]})
    metadata["Stage"] = metadata["filename"].map(rat_stage_from_filename)
    metadata["Phase"] = metadata["filename"].map(rat_phase_from_filename)
    desired = [
        ("Embryonic", "Normal"), ("Embryonic", "Proliferation"), ("Embryonic", "Remodeling"),
        ("P5", "Normal"), ("P5", "Proliferation"), ("P5", "Remodeling"),
        ("Adult", "Normal"), ("Adult", "Proliferation"), ("Adult", "Remodeling"),
    ]
    representatives = []
    for stage, phase in desired:
        current = metadata[metadata["Stage"].eq(stage) & metadata["Phase"].eq(phase)]
        if len(current):
            representatives.append(current.iloc[0])
    fig = plt.figure(figsize=(7.08, 3.25))
    grid = fig.add_gridspec(4, 3, height_ratios=[1, 1, 1, 1.38], hspace=0.20, wspace=0.10)
    rows = []
    for index, record in enumerate(representatives):
        ax = fig.add_subplot(grid[index // 3, index % 3])
        adata = ad.read_h5ad(record["path"], backed="r")
        xy = np.asarray(adata.obsm["X_rot_spatial"])
        genes = ["Vim", "Krt5"]
        missing = [gene for gene in genes if gene not in adata.var_names]
        if missing:
            raise KeyError(f"{Path(record['path']).name}: missing {missing}")
        expression = dense(adata[:, genes].X).astype(float)
        double_positive = (expression[:, 0] > 0.5) & (expression[:, 1] > 0.5)
        ax.scatter(xy[:, 0], xy[:, 1], s=0.025, color="#E2E6EA", linewidths=0, rasterized=True)
        ax.scatter(
            xy[double_positive, 0], xy[double_positive, 1], s=0.055,
            color=DOUBLE_POSITIVE_COLOR, linewidths=0, rasterized=True,
        )
        add_scale_bar(ax, xy)
        ax.set_title(f"{record['Stage']} | {record['Phase']}", fontsize=7.0, pad=1)
        ax.set_aspect("equal")
        ax.axis("off")
        adata.file.close()
    for path in paths:
        adata = ad.read_h5ad(path, backed="r")
        genes = ["Vim", "Krt5"]
        missing = [gene for gene in genes if gene not in adata.var_names]
        if missing:
            adata.file.close()
            continue
        expression = dense(adata[:, genes].X).astype(float)
        double_positive = (expression[:, 0] > 0.5) & (expression[:, 1] > 0.5)
        rows.append(
            {
                "section": Path(path).name,
                "Stage": rat_stage_from_filename(Path(path).name),
                "Phase": rat_phase_from_filename(Path(path).name),
                "double_positive_fraction": float(double_positive.mean()),
                "n_bins": int(adata.n_obs),
                "n_double_positive": int(double_positive.sum()),
                "threshold": "Vim > 0.5 and Krt5 > 0.5 in stored spatial expression",
            }
        )
        adata.file.close()
    fractions = pd.DataFrame(rows)
    fractions.to_csv(package / "tables" / "S5F_rat_Vim_Krt5_double_positive_sections_v03.csv", index=False)
    stats_ax = fig.add_subplot(grid[3, :])
    phase_order = ["Normal", "Inflammation", "Proliferation", "Remodeling"]
    x_positions = np.arange(len(phase_order))
    rng = np.random.default_rng(17)
    for stage in ["Embryonic", "P5", "Adult"]:
        means, errors = [], []
        for phase in phase_order:
            current = fractions.loc[
                fractions["Stage"].eq(stage) & fractions["Phase"].eq(phase),
                "double_positive_fraction",
            ]
            means.append(current.mean() if len(current) else np.nan)
            errors.append(current.sem() if len(current) > 1 else (0.0 if len(current) else np.nan))
            if len(current):
                stats_ax.scatter(
                    np.full(len(current), x_positions[phase_order.index(phase)]) + rng.uniform(-0.055, 0.055, len(current)),
                    current * 100, s=8, color=MATCHED_CATEGORY_COLORS[stage],
                    alpha=0.62, linewidths=0, zorder=2,
                )
        stats_ax.errorbar(
            x_positions, np.asarray(means) * 100, yerr=np.asarray(errors) * 100,
            color=MATCHED_CATEGORY_COLORS[stage], marker="o", lw=1.2, ms=3.3,
            capsize=2, label=stage, zorder=3,
        )
    stats_ax.set(
        ylabel="Vim+Krt5 double-positive bins (%)",
        xlabel="Repair phase",
        xticks=x_positions,
        xticklabels=phase_order,
    )
    stats_ax.legend(frameon=False, ncol=3, loc="upper left")
    stats_ax.grid(axis="y", color="#E6E6E6", lw=0.5)
    sns.despine(ax=stats_ax)
    fig.subplots_adjust(left=0.09, right=0.995, top=0.96, bottom=0.13)
    save(fig, package / "figures" / "panels" / "FigureS5_F_spatial_time_v03")


def place_pdf(page: fitz.Page, source: Path, rectangle: fitz.Rect) -> None:
    document = fitz.open(source)
    page.show_pdf_page(rectangle, document, 0, keep_proportion=True, overlay=True)
    document.close()


def assemble(package: Path, species: str) -> Path:
    if fitz is None:
        raise RuntimeError("PyMuPDF is required for PDF assembly")
    prefix = "FigureS4" if species == "human" else "FigureS5"
    panels = package / "figures" / "panels"
    names = {
        "A": f"{prefix}_A_fibroblast_umap_v02.pdf",
        "B": f"{prefix}_B_qc_doublets_v02.pdf",
        "C": f"{prefix}_C_markers_merged_v02.pdf",
        "D": f"{prefix}_D_compact_v02.pdf",
        "E": f"{prefix}_E_monocle_triplet_v02.pdf",
        "F": f"{prefix}_F_spatial_time_v02.pdf",
    }
    document = fitz.open()
    page = document.new_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
    margin = 18
    gap = 10
    usable_width = A4_WIDTH_PT - 2 * margin
    cursor = 20
    row_heights = {"AB": 164, "C": 100, "DE": 155, "F": 318}

    half_width = (usable_width - gap) / 2
    place_pdf(
        page,
        panels / names["A"],
        fitz.Rect(margin, cursor, margin + half_width, cursor + row_heights["AB"]),
    )
    place_pdf(
        page,
        panels / names["B"],
        fitz.Rect(
            margin + half_width + gap,
            cursor,
            A4_WIDTH_PT - margin,
            cursor + row_heights["AB"],
        ),
    )
    page.insert_text((7, cursor + 12), "A", fontsize=14, fontname="helv")
    page.insert_text(
        (margin + half_width + 1, cursor + 12), "B", fontsize=14, fontname="helv"
    )
    cursor += row_heights["AB"] + gap

    place_pdf(
        page,
        panels / names["C"],
        fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, cursor + row_heights["C"]),
    )
    page.insert_text((7, cursor + 12), "C", fontsize=14, fontname="helv")
    cursor += row_heights["C"] + gap

    third_width = (usable_width - gap) / 3
    place_pdf(
        page,
        panels / names["D"],
        fitz.Rect(margin, cursor, margin + third_width, cursor + row_heights["DE"]),
    )
    place_pdf(
        page,
        panels / names["E"],
        fitz.Rect(
            margin + third_width + gap,
            cursor,
            A4_WIDTH_PT - margin,
            cursor + row_heights["DE"],
        ),
    )
    page.insert_text((7, cursor + 12), "D", fontsize=14, fontname="helv")
    page.insert_text(
        (margin + third_width + 1, cursor + 12), "E", fontsize=14, fontname="helv"
    )
    cursor += row_heights["DE"] + gap

    place_pdf(
        page,
        panels / names["F"],
        fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, A4_HEIGHT_PT - margin),
    )
    page.insert_text((7, cursor + 12), "F", fontsize=14, fontname="helv")

    output = package / "figures" / f"{prefix}_revised_layout_v02.pdf"
    document.save(output, garbage=4, deflate=True)
    document.close()
    rendered = fitz.open(output)
    pixmap = rendered[0].get_pixmap(
        matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False
    )
    pixmap.save(
        package / "figures" / f"{prefix}_revised_layout_v02_600dpi.png"
    )
    rendered.close()
    return output


def write_legend(package: Path, species: str) -> None:
    if species == "human":
        legend = """Figure S4. Human single-cell and spatial validation of the Fib_K14 state.
(A) Fibroblast-only UMAP using a high-contrast state palette; Fib_K14 is shown in red.
(B) Single-cell quality-control distributions and independently recomputed Scrublet doublet scores for keratinocytes, other fibroblasts and Fib_K14.
(C) Fibroblast-UMAP feature plots for PDGFRA, VIM, KRT14, KRT5, TACSTD2 and COL1A1, followed by the merged six-gene density.
(D) Sample-level keratinocyte-proxy scores across human fibroblast states.
(E) Equal-sized Monocle2 maps colored by cell state, pseudotime and observed healing-time group. Arrows denote increasing pseudotime; legends are shared below the plots.
(F) PDGFRA/KRT14 merged spatial-density maps across the registered 11 human Stereo-seq sections, with 1-mm scale bars, and sample-level Fib_K14 fractions across Normal, SPTDI and DPTDI groups. Points are biological samples; lines show group means and s.e.m.
"""
    else:
        legend = """Figure S5. Rat single-cell and spatial validation of the Fib_K14-like state.
(A) Rat fibroblast-only UMAP using a high-contrast state palette; the annotated Fib_K14-like state is shown in red.
(B) Single-cell quality-control distributions and independently recomputed Scrublet doublet scores for keratinocytes, other fibroblasts and Fib_K14-like cells.
(C) Rat fibroblast-UMAP feature plots for Pdgfra, Vim, Krt14, Krt5, Tacstd2 and Col1a1, followed by the merged measurable five-gene density. Krt14 is absent from the registered rat gene universe and is not imputed.
(D) Rat Fib–keratinocyte PAGA connectivity.
(E) Equal-sized Monocle2 maps colored by cell state, pseudotime and observed healing-time group. Arrows denote increasing pseudotime; legends are shared below the plots.
(F) Representative rat spatial maps of the annotated Fib_K14-like state across embryonic, P5 and adult repair, with 1-mm scale bars, and section-level temporal fractions. Points are spatial sections; lines show stage means and s.e.m. Because Krt14 is absent, this panel does not claim direct Pdgfra/Krt14 co-expression.
"""
    (package / f"{'FigureS4' if species == 'human' else 'FigureS5'}_legend_v02.txt").write_text(
        legend, encoding="utf-8"
    )


def write_legend_v03(package: Path, species: str) -> None:
    if species == "human":
        legend = """Figure S4. Human single-cell, spatial and staining validation of the epithelial-like fibroblast phenotype.
(A) Fibroblast-only UMAP using a high-contrast state palette.
(B) Single-cell quality-control distributions and independently recomputed Scrublet doublet scores for keratinocytes and fibroblast states.
(C) Fibroblast-UMAP feature plots for PDGFRA, VIM, KRT14, KRT5, TACSTD2 and COL1A1, followed by the merged six-gene density.
(D) Sample-level keratinocyte-proxy scores across human fibroblast states.
(E) Equal-sized Monocle2 maps colored by cell state, pseudotime and observed healing-time group. Arrows denote increasing pseudotime; legends are shared below the plots.
(F) VIM/KRT14 double-positive spatial bins across the registered 11 human Stereo-seq sections. A bin was defined as double positive only when both stored VIM and KRT14 expression exceeded 0.5. The line plot independently summarizes scRNA-seq cells satisfying the same threshold after library-size normalization and log1p transformation, stratified by Normal, SPTDI and DPTDI. Points are biological samples; lines show group means and s.e.m.
(G) Representative VIM/KRT14/DAPI multiplex immunofluorescence fields and slide-level temporal quantification. Nuclei were detected from DAPI; pooled Otsu thresholds were applied unchanged to per-nucleus VIM and KRT14 signals. Three fields were averaged per slide, and one slide/sample was the statistical unit. Scale bars, 200 µm.
"""
    else:
        legend = """Figure S5. Rat single-cell and spatial validation of the Vim/Krt5 double-positive phenotype.
(A) Rat fibroblast-only UMAP using the same high-contrast palette system as Figure S4.
(B) Single-cell quality-control distributions and independently recomputed Scrublet doublet scores for keratinocytes and fibroblast states.
(C) Rat fibroblast-UMAP feature plots for Pdgfra, Vim, Krt5, Tacstd2 and Col1a1, followed by the merged five-gene density. Krt14 is absent and is not displayed or imputed.
(D) Rat fibroblast-keratinocyte PAGA connectivity.
(E) Equal-sized Monocle2 maps colored by cell state, pseudotime and observed healing-time group. Arrows denote increasing pseudotime; legends are shared below the plots.
(F) Rat spatial bins satisfying Vim > 0.5 and Krt5 > 0.5 across embryonic, P5 and adult repair, with 1-mm scale bars, and section-level temporal fractions. Points are spatial sections; lines show stage means and s.e.m. No Fib_K14 or Fib_K14-like label is used.
"""
    prefix = "FigureS4" if species == "human" else "FigureS5"
    (package / f"{prefix}_legend_v03.txt").write_text(legend, encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup()
    package = Path(args.package)
    requested = set(args.panels.upper())
    if "A" in requested:
        print("render A", flush=True)
        render_a(package, args.species)
    if "B" in requested:
        print("render B", flush=True)
        render_b(package, args.species)
    if "C" in requested:
        print("render C", flush=True)
        render_c(package, args.species, args.scrna)
    if "D" in requested:
        print("render D", flush=True)
        render_d(package, args.species)
    if "E" in requested:
        print("render E", flush=True)
        render_e(package, args.species)
    if "F" in requested:
        print("render F", flush=True)
        if args.species == "human":
            human_f(package, args.scrna, args.spatial)
        else:
            rat_f(package, args.spatial)
    write_legend_v03(package, args.species)
    print("legend written", flush=True)
    if not args.skip_assemble:
        output = assemble(package, args.species)
        print(output)


if __name__ == "__main__":
    main()
