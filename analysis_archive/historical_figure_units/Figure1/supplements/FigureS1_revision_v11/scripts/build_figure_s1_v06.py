from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[3]
SOURCE_QC = PROJECT / "Figure1_Supplement_QC_20260721"
FIG1 = PROJECT / "Final_Figure_Units" / "Figure1"
OUT = ROOT / "outputs"
PANELS = OUT / "panels"
VECTOR = OUT / "vector"
RASTER = OUT / "raster"
ASSETS = ROOT / "assets" / "linked_rasters"
QC = ROOT / "QC"

SCRNA = Path(os.environ.get("FIGS1_SCRNA_H5AD", r"data/external/local_decd995edda5"))
STEREO = Path(os.environ.get(
    "FIGS1_STEREO_ROOT",
    r"data/external/local_d4d64fe45cab"
    r"\stereoseq\pdf_output\Wound_Healing_Annotation_Output_White",
))
VISIUM = Path(os.environ.get(
    "FIGS1_VISIUM_H5AD",
    r"data/external/local_d4d64fe45cab"
    r"\10X visium\h5ad_output\Combined_Aligned_800k_20260403.h5ad",
))
OLD_FIG3 = Path(r"data/external/local_343563411aff")
OLD_S31 = ROOT / "inputs" / "old_supplement_docx" / "word" / "media" / "image31.png"
FIG1A_COMPOSITE = FIG1 / "inputs" / "raw_panels" / "ppt_slide2_media" / "image3.png"

PAGE_W_MM, PAGE_H_MM = 210.0, 297.0
MARGIN_X, MARGIN_Y, GAP = 7.0, 6.0, 3.0
CONTENT_W = PAGE_W_MM - 2 * MARGIN_X
AB_PANEL_HEIGHT = 70.0
LOWER_PANEL_HEIGHTS = [42.0, 62.0, 44.0]
BOTTOM_PANEL_HEIGHT = 38.0

SECTIONS = [
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
    ("Normal", "NS_C02847B1", 180),
]
TEMPORAL_SECTIONS = [
    ("Normal", "NS_C02847B1", 180),
    ("5 dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12 dpb SPTDI", "BW14_B3_SSDB_12dpb", 135),
    ("19 dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("2 mph", "2mph_A03699G6.SCT", 45),
]
KC_STATES = ["KC_Basal", "KC_Basal_Mig", "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular"]
KC_COLORS = {
    "KC_Basal": "#008B45",
    "KC_Basal_Mig": "#FF8C00",
    "KC_Spinous": "#00FF00",
    "KC_Spinous_Mig": "#00008B",
    "KC_Spinous_Mat": "#FF69B4",
    "KC_Granular": "#8B4513",
}

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 6,
    "axes.titlesize": 8,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.55,
})


def decode_strings(values: np.ndarray) -> np.ndarray:
    return np.asarray([value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values], dtype=object)


def decode_categorical(node: h5py.Group) -> np.ndarray:
    categories = decode_strings(node["categories"][:])
    codes = np.asarray(node["codes"], dtype=int)
    result = np.full(codes.shape, "NA", dtype=object)
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    return result


def read_obs_field(handle: h5py.File, field: str) -> np.ndarray:
    node = handle["obs"][field]
    return decode_categorical(node) if isinstance(node, h5py.Group) else np.asarray(node)


def rotate(coords: np.ndarray, angle: float) -> np.ndarray:
    values = np.asarray(coords, dtype=float).copy()
    values -= np.nanmedian(values, axis=0)
    theta = np.deg2rad(angle)
    return np.column_stack((
        values[:, 0] * np.cos(theta) - values[:, 1] * np.sin(theta),
        values[:, 0] * np.sin(theta) + values[:, 1] * np.cos(theta),
    ))


def rect_mm(x: float, top: float, width: float, height: float) -> list[float]:
    return [x / PAGE_W_MM, (PAGE_H_MM - top - height) / PAGE_H_MM, width / PAGE_W_MM, height / PAGE_H_MM]


def panel_header(fig: plt.Figure, rect: list[float], letter: str, title: str) -> None:
    x, y, _w, h = rect
    fig.text(x, y + h, letter, ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(x + 0.024, y + h - 0.001, title, ha="left", va="top", fontsize=9, fontweight="bold")


def add_scale_bar(ax: plt.Axes, length_units: float = 2000.0, label: str = "1 mm", color: str = "black") -> None:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x0 = xmin + 0.06 * (xmax - xmin)
    y0 = ymin + 0.08 * (ymax - ymin)
    ax.plot([x0, x0 + length_units], [y0, y0], color=color, linewidth=1.2, solid_capstyle="butt")
    ax.text(x0 + length_units / 2, y0 + 0.025 * (ymax - ymin), label, ha="center", va="bottom", fontsize=8, color=color)


def crop_legacy_assets() -> dict[str, Path]:
    ASSETS.mkdir(parents=True, exist_ok=True)
    outputs = {
        "clinical": ASSETS / "FigureS1A_clinical_modes_legacy_native_crop.png",
        "he": ASSETS / "FigureS1C_HE_DPTDI_SPTDI_native_crop.tiff",
        "stats": ASSETS / "FigureS1C_appendage_statistics_native_crop.png",
        "schematic": ASSETS / "FigureS1C_wound_outcome_schematic_imagegen_v01.png",
        "burn_schematic": ASSETS / "FigureS1C_burn_outcome_illustration_base_v01.png",
    }
    with Image.open(FIG1A_COMPOSITE) as image:
        image.crop((28, 190, 1024, 379)).save(outputs["clinical"], dpi=image.info.get("dpi", (96, 96)))
    with Image.open(OLD_FIG3) as image:
        image.crop((135, 0, 1735, 475)).save(outputs["he"], dpi=image.info.get("dpi", (600, 600)), compression="tiff_lzw")
    with Image.open(OLD_S31) as image:
        image.crop((610, 15, 900, 365)).save(outputs["stats"], dpi=image.info.get("dpi", (168, 168)))
    return outputs


def panel_a(fig: plt.Figure, rect: list[float], assets: dict[str, Path]) -> None:
    panel_header(fig, rect, "A", "Clinical wound-healing modes")
    x, y, w, h = rect
    image = np.asarray(Image.open(assets["clinical"]).convert("RGB"))
    crop_width = image.shape[1] / 4
    labels = [
        "Mode 1\nUnhealed gradient",
        "Mode 2\nOne-way epithelialization",
        "Mode 3\nAlmost healed",
        "Mode 4\nScar / hyperplasia",
    ]
    row_top = y + 0.87 * h
    row_h = 0.205 * h
    for index, label in enumerate(labels):
        row_y = row_top - (index + 1) * row_h
        ax = fig.add_axes([x + 0.04 * w, row_y + 0.03 * row_h, 0.48 * w, 0.88 * row_h])
        left = int(index * crop_width)
        right = int((index + 1) * crop_width)
        ax.imshow(image[:, left:right])
        ax.axis("off")
        # Approximate scale only: inferred from the visible adult lower-leg anatomy.
        bar_x0, bar_x1, bar_y = 0.08, 0.25, 0.08
        ax.plot([bar_x0, bar_x1], [bar_y, bar_y], transform=ax.transAxes, color="white", linewidth=2.4, solid_capstyle="butt")
        ax.plot([bar_x0, bar_x1], [bar_y, bar_y], transform=ax.transAxes, color="black", linewidth=1.1, solid_capstyle="butt")
        ax.text((bar_x0 + bar_x1) / 2, bar_y + 0.04, "≈2 cm", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=6, color="black", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.80, "pad": 0.3})
        fig.text(
            x + 0.57 * w,
            row_y + 0.48 * row_h,
            label,
            ha="left",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            linespacing=1.05,
        )


def sample_qc_tables() -> dict[str, pd.DataFrame]:
    tables = {}
    with h5py.File(SCRNA, "r") as handle:
        sample = read_obs_field(handle, "sample_id")
        stage = read_obs_field(handle, "Time_category")
        values = pd.DataFrame({
            "sample": sample,
            "stage": stage,
            "umi": np.asarray(read_obs_field(handle, "nCount_RNA"), dtype=float),
            "gene": np.asarray(read_obs_field(handle, "nFeature_RNA"), dtype=float),
        })
        order = {"Normal": 0, "0-7dpi": 1, "0-7dpb": 1, "8-14dpi": 2, "8-14dpb": 2, "15-28dpi": 3, "15-28dpb": 3, "1-2mph": 4}
        sample_order = (
            values[["sample", "stage"]].drop_duplicates()
            .assign(rank=lambda frame: frame["stage"].map(order).fillna(9))
            .sort_values(["rank", "sample"])["sample"].tolist()
        )
        values["sample"] = pd.Categorical(values["sample"], categories=sample_order, ordered=True)
        tables["scRNA-seq"] = values

    rows = []
    for label, sample_id, _angle in SECTIONS:
        with h5py.File(STEREO / f"{sample_id}_Annotated.h5ad", "r") as handle:
            rows.append(pd.DataFrame({
                "sample": label,
                "stage": label,
                "umi": np.asarray(read_obs_field(handle, "nCount_Spatial"), dtype=float),
                "gene": np.asarray(read_obs_field(handle, "nFeature_Spatial"), dtype=float),
            }))
    tables["Stereo-seq"] = pd.concat(rows, ignore_index=True)

    with h5py.File(VISIUM, "r") as handle:
        sample = read_obs_field(handle, "sample_batch_new")
        tables["10x Visium"] = pd.DataFrame({
            "sample": sample,
            "stage": sample,
            "umi": np.asarray(read_obs_field(handle, "total_counts"), dtype=float),
            "gene": np.asarray(read_obs_field(handle, "n_genes_by_counts"), dtype=float),
        })
    return tables


def short_sample(value: str) -> str:
    result = str(value).split("-Skin-")[0]
    return result.replace("_SDSDB_", "\n").replace("_DSDB_", "\n").replace("_SSDB_", "\n")


def draw_boxplot(ax: plt.Axes, table: pd.DataFrame, metric: str, color: str, prefix: str) -> None:
    groups = [(str(name), group[metric].dropna().to_numpy()) for name, group in table.groupby("sample", observed=True, sort=False)]
    positions = np.arange(len(groups))
    ax.boxplot(
        [np.log10(values + 1) for _, values in groups],
        positions=positions,
        widths=0.66,
        showfliers=False,
        patch_artist=False,
        medianprops={"color": "#222222", "linewidth": 0.55},
        whiskerprops={"color": "#444444", "linewidth": 0.45},
        capprops={"color": "#444444", "linewidth": 0.45},
        boxprops={"color": color, "linewidth": 0.65},
    )
    ax.set_xticks(positions, [f"{prefix}{index + 1:02d}" for index in range(len(groups))], rotation=90, ha="center")
    ax.tick_params(axis="x", labelsize=6, length=1, pad=0.5)
    ax.tick_params(axis="y", labelsize=6, length=2, pad=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.4)


def panel_b(fig: plt.Figure, rect: list[float], tables: dict[str, pd.DataFrame]) -> None:
    panel_header(fig, rect, "B", "Per-sample sequencing depth and detected genes across all human datasets")
    x, y, w, h = rect
    colors = {"scRNA-seq": "#4C78A8", "Stereo-seq": "#59A14F", "10x Visium": "#F28E2B"}
    prefixes = {"scRNA-seq": "R", "Stereo-seq": "ST", "10x Visium": "V"}
    top_gap, bottom_gap, row_gap = 0.11 * h, 0.02 * h, 0.025 * h
    row_h = (h - top_gap - bottom_gap - 2 * row_gap) / 3
    left = x + 0.058 * w
    plot_w = 0.43 * w
    right_x = x + 0.55 * w
    for row, modality in enumerate(("scRNA-seq", "Stereo-seq", "10x Visium")):
        bottom = y + h - top_gap - (row + 1) * row_h - row * row_gap
        for col, metric in enumerate(("umi", "gene")):
            ax = fig.add_axes([left if col == 0 else right_x, bottom + 0.31 * row_h, plot_w, 0.69 * row_h])
            draw_boxplot(ax, tables[modality], metric, colors[modality], prefixes[modality])
            if row == 0:
                ax.set_title("Number of UMI" if metric == "umi" else "Number of Gene", fontsize=8, pad=1.5)
        fig.text(x + 0.005 * w, bottom + row_h / 2, modality, rotation=90, ha="left", va="center", fontsize=7, fontweight="bold", color=colors[modality])


def draw_skin_cartoon(ax: plt.Axes, deep: bool, title: str, outcome: str, color: str) -> None:
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_facecolor("#FFF7ED")
    x = np.linspace(0.5, 9.5, 200)
    surface = 3.45 + 0.12 * np.sin(x * 1.8)
    ax.fill_between(x, 0.7, surface, color="#F6D6CF")
    ax.fill_between(x, surface, surface + 0.35, color="#D98A8A")
    center = 4.9
    width = 2.6 if deep else 1.6
    depth = 2.6 if deep else 1.45
    wound = depth * np.exp(-((x - center) / width) ** 4)
    ax.fill_between(x, surface - wound, surface + 0.36, where=wound > 0.08, color="white")
    for follicle_x in ([1.8, 7.7] if deep else [1.4, 3.0, 6.9, 8.5]):
        ax.plot([follicle_x, follicle_x - 0.15], [3.1, 1.4], color="#7C3F2C", linewidth=1.3)
        ax.scatter([follicle_x - 0.15], [1.25], s=10, c="#7C3F2C", linewidths=0)
    ax.text(0.5, 4.7, title, fontsize=6.5, fontweight="bold", color=color, va="top")
    ax.text(9.4, 0.55, outcome, fontsize=6.2, fontweight="bold", color=color, ha="right", va="bottom")


def panel_c(fig: plt.Figure, rect: list[float], assets: dict[str, Path]) -> None:
    panel_header(fig, rect, "C", "Burn-depth-dependent healing outcomes and residual appendages")
    x, y, w, h = rect
    base = Image.open(assets["burn_schematic"]).convert("RGB")
    outcome_crops = [
        base.crop((1166, 15, 1683, 376)),
        base.crop((1166, 450, 1683, 824)),
    ]
    captions = [
        ("SPTDI: ≈2 weeks\nScarless healing\nMany residual appendages", "#007F87"),
        ("DPTDI: ≈4 weeks\nScar-forming healing\nFew residual appendages", "#62358A"),
    ]
    for index, (crop, (caption, color)) in enumerate(zip(outcome_crops, captions)):
        tile_y = y + (0.48 if index == 0 else 0.08) * h
        ax = fig.add_axes([x + 0.025 * w, tile_y, 0.245 * w, 0.31 * h])
        ax.imshow(crop)
        ax.axis("off")
        fig.text(
            x + 0.285 * w,
            tile_y + 0.155 * h,
            caption,
            ha="left",
            va="center",
            fontsize=6.5,
            fontweight="bold",
            color=color,
            linespacing=1.05,
        )

    he = np.asarray(Image.open(assets["he"]).convert("RGB"))
    he_ax = fig.add_axes([x + 0.46 * w, y + 0.10 * h, 0.36 * w, 0.69 * h])
    he_ax.imshow(he)
    he_ax.axis("off")
    he_ax.set_title("H&E at 5 dpb", fontsize=7, pad=1, fontweight="bold")

    stats = np.asarray(Image.open(assets["stats"]).convert("RGB"))
    stats_ax = fig.add_axes([x + 0.835 * w, y + 0.10 * h, 0.155 * w, 0.69 * h])
    stats_ax.imshow(stats)
    stats_ax.axis("off")
    stats_ax.set_title("Residual appendages", fontsize=7, pad=1, fontweight="bold")
    fig.text(x + 0.42 * w, y + 0.015 * h, "Outcome illustrations are explanatory and not to scale.", ha="right", va="bottom", fontsize=6, color="#555555")


def marker_matrix() -> tuple[list[str], list[str], np.ndarray]:
    markers = pd.read_csv(SOURCE_QC / "cell_state_markers.csv")
    states = markers["cell_state"].astype(str).tolist()
    chosen, used = [], set()
    with h5py.File(SCRNA, "r") as handle:
        var_names = decode_strings(handle["var"]["_index"][:])
        lookup = {gene: index for index, gene in enumerate(var_names)}
        for _, row in markers.iterrows():
            candidates = [str(row[f"marker_{index}"]) for index in range(1, 5)]
            gene = next((candidate for candidate in candidates if candidate in lookup and candidate not in used), None)
            gene = gene or next((candidate for candidate in candidates if candidate in lookup), None)
            if gene is None:
                raise RuntimeError(f"No detected marker for {row['cell_state']}")
            chosen.append(gene)
            used.add(gene)
        gene_indices = [lookup[gene] for gene in chosen]
        labels = decode_categorical(handle["obs"]["sub_labels"])
        state_index = {state: index for index, state in enumerate(states)}
        sums = np.zeros((len(states), len(chosen)), dtype=float)
        counts = np.asarray([(labels == state).sum() for state in states], dtype=float)
        group = handle["X"]
        indptr, indices, data = group["indptr"], group["indices"], group["data"]
        for start in range(0, len(labels), 5000):
            stop = min(start + 5000, len(labels))
            pointers = np.asarray(indptr[start:stop + 1], dtype=np.int64)
            low, high = int(pointers[0]), int(pointers[-1])
            block = sparse.csr_matrix(
                (np.asarray(data[low:high]), np.asarray(indices[low:high]), pointers - low),
                shape=(stop - start, len(var_names)),
            )[:, gene_indices]
            block_labels = labels[start:stop]
            for state in np.unique(block_labels):
                if state in state_index:
                    sums[state_index[state]] += np.asarray(block[block_labels == state].sum(axis=0)).ravel()
    means = sums / np.maximum(counts[:, None], 1)
    values = (means - means.mean(axis=0)) / (means.std(axis=0) + 1e-8)
    return states, chosen, np.clip(values, -2.5, 2.5)


def draw_heatmap_block(fig: plt.Figure, box: list[float], states: list[str], genes: list[str], values: np.ndarray, show_genes: bool, title: str) -> mpl.image.AxesImage:
    ax = fig.add_axes(box)
    image = ax.imshow(values, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto", interpolation="nearest", rasterized=True)
    ax.set_yticks(np.arange(len(states)), states)
    ax.tick_params(axis="y", labelsize=6, length=0, pad=0.8)
    if show_genes:
        ax.set_xticks(np.arange(len(genes)), genes, rotation=90, ha="center", va="top")
        ax.tick_params(axis="x", labelsize=6, length=0, pad=0.8)
        for label in ax.get_xticklabels():
            label.set_fontstyle("italic")
    else:
        ax.set_xticks([])
    ax.set_title(title, fontsize=7, fontweight="bold", pad=1.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return image


def panel_d(fig: plt.Figure, rect: list[float]) -> None:
    panel_header(fig, rect, "D", "Canonical marker expression across 47 annotated cell states")
    x, y, w, h = rect
    states, genes, values = marker_matrix()
    boundaries = [(0, 16), (16, 32), (32, len(states))]
    titles = ["Epithelial and appendage states", "Stromal and vascular states", "Immune states"]
    boxes = [
        [x + 0.105 * w, y + 0.27 * h, 0.190 * w, 0.61 * h],
        [x + 0.435 * w, y + 0.27 * h, 0.190 * w, 0.61 * h],
        [x + 0.765 * w, y + 0.27 * h, 0.190 * w, 0.61 * h],
    ]
    image = None
    for box, (start, stop), title in zip(boxes, boundaries, titles):
        image = draw_heatmap_block(
            fig,
            box,
            states[start:stop],
            genes[start:stop],
            values[start:stop, start:stop],
            True,
            title,
        )
    cax = fig.add_axes([x + 0.976 * w, y + 0.39 * h, 0.008 * w, 0.34 * h])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label("Mean expression z-score", fontsize=6)
    colorbar.ax.tick_params(labelsize=6, length=2, width=0.4)


def read_spatial(sample_id: str, fields: list[str]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with h5py.File(STEREO / f"{sample_id}_Annotated.h5ad", "r") as handle:
        coords = np.asarray(handle["obsm"]["spatial"], dtype=float)
        obs = {field: read_obs_field(handle, field) for field in fields}
    return coords, obs


def spatial_cache() -> list[tuple[str, str, float, np.ndarray, dict[str, np.ndarray]]]:
    result = []
    for label, sample_id, angle in SECTIONS:
        coords, obs = read_spatial(sample_id, ["nCount_Spatial", "nFeature_Spatial"])
        result.append((label, sample_id, angle, rotate(coords, angle), obs))
    return result


def panel_e(fig: plt.Figure, rect: list[float], cached: list[tuple[str, str, float, np.ndarray, dict[str, np.ndarray]]]) -> None:
    panel_header(fig, rect, "E", "Spatial QC across all 11 registered Stereo-seq sections")
    x, y, w, h = rect
    logs = {
        metric: np.concatenate([np.log10(np.asarray(item[4][metric], dtype=float) + 1) for item in cached])
        for metric in ("nCount_Spatial", "nFeature_Spatial")
    }
    limits = {metric: np.nanquantile(values, [0.01, 0.99]) for metric, values in logs.items()}
    cols, rows = 4, 3
    left, right, bottom, top = 0.018, 0.015, 0.21, 0.22
    gx, gy = 0.012, 0.055
    tile_w = (w * (1 - left - right) - w * gx * (cols - 1)) / cols
    tile_h = (h * (1 - bottom - top) - h * gy * (rows - 1)) / rows
    for index, (label, _sample_id, _angle, coords, obs) in enumerate(cached):
        row, col = divmod(index, cols)
        tile_x = x + w * left + col * (tile_w + w * gx)
        tile_y = y + h * (1 - top) - (row + 1) * tile_h - row * h * gy
        for metric_index, metric in enumerate(("nCount_Spatial", "nFeature_Spatial")):
            ax = fig.add_axes([tile_x + metric_index * tile_w * 0.51, tile_y, tile_w * 0.48, tile_h])
            values = np.log10(np.asarray(obs[metric], dtype=float) + 1)
            ax.scatter(coords[:, 0], coords[:, 1], c=values, s=0.05, cmap="viridis", vmin=limits[metric][0], vmax=limits[metric][1], linewidths=0, rasterized=True)
            ax.set_aspect("equal")
            ax.axis("off")
            if metric_index == 0:
                add_scale_bar(ax)
        fig.text(tile_x + tile_w / 2, tile_y + tile_h + 0.002, label, ha="center", va="bottom", fontsize=8, fontweight="bold")
    fig.text(
        x + 0.99 * w,
        y + 0.885 * h,
        "Each section pair:  Number of UMI  |  Number of Gene",
        ha="right",
        va="center",
        fontsize=8,
        fontweight="bold",
    )
    for index, metric in enumerate(("nCount_Spatial", "nFeature_Spatial")):
        cax = fig.add_axes([x + (0.28 + index * 0.42) * w, y + 0.075 * h, 0.22 * w, 0.018 * h])
        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(*limits[metric]), cmap="viridis")
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=8, length=1, pad=1)
        fig.text(
            x + (0.39 + index * 0.42) * w,
            y + 0.145 * h,
            "log10(Number of UMI + 1)" if index == 0 else "log10(Number of Gene + 1)",
            ha="center",
            va="bottom",
            fontsize=8,
        )


def panel_f(fig: plt.Figure, rect: list[float]) -> None:
    panel_header(fig, rect, "F", "Reserved for temporal clinical immunofluorescence")
    x, y, w, h = rect
    ax = fig.add_axes([x + 0.025 * w, y + 0.10 * h, 0.95 * w, 0.70 * h])
    ax.set_facecolor("#F9FAFB")
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color="#9CA3AF", linewidth=0.8, linestyle=(0, (3, 2)))
    ax.text(0.5, 0.5, "Placeholder: normal skin + wound / healed skin / scar;\nKRT14 or KRT5, KRT1 and DAPI",
            ha="center", va="center", fontsize=8, color="#6B7280")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def panel_g(fig: plt.Figure, rect: list[float]) -> None:
    panel_header(fig, rect, "G", "Temporal spatial distribution of keratinocyte states")
    x, y, w, h = rect
    left, right, top, bottom, gap_x, gap_y = 0.04, 0.03, 0.18, 0.22, 0.035, 0.05
    map_w = (w * (1 - left - right) - w * gap_x * 2) / 3
    map_h = h * (1 - top - bottom - gap_y) / 2
    for index, (label, sample_id, angle) in enumerate(TEMPORAL_SECTIONS):
        coords, obs = read_spatial(sample_id, ["sub_labels"])
        coords = rotate(coords, angle)
        labels = obs["sub_labels"]
        row, col = divmod(index, 3)
        row_count = 3 if row == 0 else 2
        row_offset = 0 if row == 0 else (3 - row_count) * (map_w + w * gap_x) / 2
        map_x = x + w * left + row_offset + col * (map_w + w * gap_x)
        map_y = y + h * bottom + (1 - row) * (map_h + h * gap_y)
        ax = fig.add_axes([map_x, map_y, map_w, map_h])
        ax.scatter(coords[:, 0], coords[:, 1], c="#E5E7EB", s=0.04, linewidths=0, rasterized=True)
        for state in KC_STATES:
            mask = labels == state
            ax.scatter(coords[mask, 0], coords[mask, 1], c=KC_COLORS[state], s=0.16, linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(label, fontsize=8, fontweight="bold", pad=1)
        add_scale_bar(ax)
    handles = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=KC_COLORS[state], markeredgecolor="none", markersize=3, label=state) for state in KC_STATES]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(x + w / 2, y + 0.015 * h), ncol=3,
               frameon=False, fontsize=8, handletextpad=0.3, columnspacing=0.8)


def main() -> None:
    for directory in (PANELS, VECTOR, RASTER, ASSETS, QC):
        directory.mkdir(parents=True, exist_ok=True)
    assets = crop_legacy_assets()
    qc_tables = sample_qc_tables()
    cached_spatial = spatial_cache()
    data_dir = ROOT / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    prefixes = {"scRNA-seq": "R", "Stereo-seq": "ST", "10x Visium": "V"}
    mapping_rows = []
    for modality, table in qc_tables.items():
        names = [str(name) for name, _group in table.groupby("sample", observed=True, sort=False)]
        mapping_rows.extend(
            {"modality": modality, "plot_code": f"{prefixes[modality]}{index + 1:02d}", "sample_id": name}
            for index, name in enumerate(names)
        )
    pd.DataFrame(mapping_rows).to_csv(data_dir / "FigureS1B_sample_code_mapping.csv", index=False)

    top = MARGIN_Y
    a_width = CONTENT_W * 0.25
    b_width = CONTENT_W - a_width - GAP
    rects = {
        "A": rect_mm(MARGIN_X, top, a_width, AB_PANEL_HEIGHT),
        "B": rect_mm(MARGIN_X + a_width + GAP, top, b_width, AB_PANEL_HEIGHT),
    }
    top += AB_PANEL_HEIGHT + GAP
    for letter, height in zip("CDE", LOWER_PANEL_HEIGHTS):
        rects[letter] = rect_mm(MARGIN_X, top, CONTENT_W, height)
        top += height + GAP
    half_width = (CONTENT_W - GAP) / 2
    rects["F"] = rect_mm(MARGIN_X, top, half_width, BOTTOM_PANEL_HEIGHT)
    rects["G"] = rect_mm(MARGIN_X + half_width + GAP, top, half_width, BOTTOM_PANEL_HEIGHT)

    fig = plt.figure(figsize=(PAGE_W_MM / 25.4, PAGE_H_MM / 25.4), facecolor="white")
    panel_a(fig, rects["A"], assets)
    panel_b(fig, rects["B"], qc_tables)
    panel_c(fig, rects["C"], assets)
    panel_d(fig, rects["D"])
    panel_e(fig, rects["E"], cached_spatial)
    panel_f(fig, rects["F"])
    panel_g(fig, rects["G"])

    stem = "FigureS1_revised_v06"
    pdf = VECTOR / f"{stem}.pdf"
    svg = VECTOR / f"{stem}_editable.svg"
    png = RASTER / f"{stem}_600dpi.png"
    tiff = RASTER / f"{stem}_600dpi.tiff"
    fig.savefig(svg)
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(tiff, dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(pdf)
    plt.close(fig)

    layout = {
        "page_mm": [PAGE_W_MM, PAGE_H_MM],
        "page_count": 1,
        "font": "Arial",
        "font_range_pt": [6.0, 12.0],
        "panel_gap_mm": GAP,
        "panels": rects,
        "human_scrna_samples": int(qc_tables["scRNA-seq"]["sample"].nunique()),
        "human_stereo_sections": int(qc_tables["Stereo-seq"]["sample"].nunique()),
        "human_visium_sections": int(qc_tables["10x Visium"]["sample"].nunique()),
        "spatial_scale": "1 mm = 2000 Stereo-seq coordinate units (0.5 um DNB pitch)",
        "clinical_photo_scale_status": "Approximate 2-cm bars inferred from adult lower-leg anatomy; not acquisition-calibrated.",
        "panel_F_status": "Intentionally reserved placeholder at author request.",
        "panel_C_interpretation": "Second author statement treated as DPTDI ~4-week healing with scar, consistent with supplied H&E/appendage evidence.",
    }
    (QC / "FigureS1_v06_layout_and_content.json").write_text(json.dumps(layout, indent=2), encoding="utf-8")
    print(pdf)


if __name__ == "__main__":
    main()
