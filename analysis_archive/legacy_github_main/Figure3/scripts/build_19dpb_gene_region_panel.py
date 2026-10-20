from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from scipy import stats
from shapely import distance, points
from shapely.geometry import LineString
from statsmodels.stats.multitest import multipletests


PROJECT = Path(r"data/external/local_91ac4b6618ea")
H5AD = PROJECT / "时空组学" / "stereoseq" / "pdf_output" / "BW15D_C6_SDSDB_19dpb_Annotated.h5ad"
ANCHORS = PROJECT / "时空组学" / "stereoseq" / "pdf_output" / "19dpb_spatial_anchors.json"
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "panels" / "19dpb_gene_region_expression"
DATA = ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

PHENOTYPES = {
    "iPSC Reprogramming": ["KLF4", "SOX2", "POU5F1", "MYC"],
    "Epithelialization": ["GRHL3", "TP63", "TACSTD2", "IGFL1"],
    "Scar Formation": ["COL1A1", "ACTA2", "FN1", "TGFB1"],
}
GENES = [gene for genes in PHENOTYPES.values() for gene in genes]
REGIONS = ["Un-epi", "Epi-Front", "Newly-epi"]
REGION_COLORS = {"Un-epi": "#E53935", "Epi-Front": "#FDD835", "Newly-epi": "#3A77B7"}
ZISSOU = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
CMAP = LinearSegmentedColormap.from_list("zissou", ZISSOU)

# Within-section spatial blocking. These are analytical blocks, not patients.
X_BLOCK_UM = 500
Y_BLOCK = 0.10
MIN_BINS_PER_BLOCK = 20

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial"],
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "text.color": "white",
        "axes.labelcolor": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "axes.edgecolor": "white",
        "axes.facecolor": "black",
        "figure.facecolor": "black",
        "savefig.facecolor": "black",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def load_data() -> tuple[pd.DataFrame, np.ndarray, dict]:
    with open(ANCHORS, encoding="utf-8") as handle:
        anchor = json.load(handle)

    a = ad.read_h5ad(H5AD, backed="r")
    upper = np.asarray([str(v).upper() for v in a.var_names])
    missing = [g for g in GENES if not np.any(upper == g)]
    if missing:
        raise ValueError(f"Genes absent from the 19dpb object: {missing}")
    gene_idx = [int(np.flatnonzero(upper == g)[0]) for g in GENES]
    expr = a[:, gene_idx].X
    expr = expr.toarray() if sp.issparse(expr) else np.asarray(expr)
    raw_xy = np.asarray(a.obsm["spatial"])
    obs_names = np.asarray(a.obs_names.astype(str))
    a.file.close()

    angle = np.radians(float(anchor.get("rotation_applied", -225)))
    rot_x = raw_xy[:, 0] * np.cos(angle) - raw_xy[:, 1] * np.sin(angle)
    rot_y = raw_xy[:, 0] * np.sin(angle) + raw_xy[:, 1] * np.cos(angle)
    pts = points(rot_x, rot_y)
    d_epi = np.asarray(distance(pts, LineString(anchor["epi_baseline"])))
    d_der = np.asarray(distance(pts, LineString(anchor["der_bottom"])))
    rel_depth = np.divide(d_epi, d_epi + d_der, out=np.zeros_like(d_epi), where=(d_epi + d_der) > 0)
    front = LineString(anchor["leading_edge"])
    front_dist = np.asarray(distance(pts, front))
    front_x = float(anchor["leading_edge"][0][0])
    if anchor.get("healed_direction") == "right_is_healed":
        sign = np.where(rot_x > front_x, 1.0, -1.0)
    else:
        sign = np.where(rot_x < front_x, 1.0, -1.0)
    distance_um = front_dist * sign * 0.33

    b1 = -600 + 300 * rel_depth
    b2 = 400 + 300 * rel_depth
    region = np.where(distance_um < b1, "Un-epi", np.where(distance_um <= b2, "Epi-Front", "Newly-epi"))
    frame = pd.DataFrame(
        {
            "bin_id": obs_names,
            "distance_to_front_um": distance_um,
            "relative_depth": rel_depth,
            "region": region,
        }
    )
    return frame, expr.astype(np.float32, copy=False), anchor


def module_scores(expr: np.ndarray) -> dict[str, np.ndarray]:
    scores: dict[str, np.ndarray] = {}
    gene_to_col = {g: i for i, g in enumerate(GENES)}
    for name, genes in PHENOTYPES.items():
        scaled = []
        for gene in genes:
            values = expr[:, gene_to_col[gene]]
            positive = values[values > 0]
            cap = float(np.quantile(positive, 0.99)) if len(positive) else 1.0
            scaled.append(np.clip(values / max(cap, np.finfo(float).eps), 0, 1))
        scores[name] = np.mean(np.column_stack(scaled), axis=1)
    return scores


def make_block_table(frame: pd.DataFrame, expr: np.ndarray) -> pd.DataFrame:
    block_frame = frame.copy()
    block_frame["x_block"] = np.floor(block_frame["distance_to_front_um"] / X_BLOCK_UM).astype(int)
    block_frame["y_block"] = np.floor(block_frame["relative_depth"] / Y_BLOCK).astype(int)
    rows = []
    for gene_idx, gene in enumerate(GENES):
        temp = block_frame[["region", "x_block", "y_block"]].copy()
        temp["expression"] = expr[:, gene_idx]
        grouped = temp.groupby(["region", "x_block", "y_block"], observed=True)["expression"].agg(["mean", "size"]).reset_index()
        grouped = grouped[grouped["size"] >= MIN_BINS_PER_BLOCK].copy()
        grouped["gene"] = gene
        grouped.rename(columns={"mean": "block_mean_expression", "size": "n_bins"}, inplace=True)
        rows.append(grouped)
    return pd.concat(rows, ignore_index=True)


def compute_statistics(blocks: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    global_rows = []
    pair_rows = []
    for gene in GENES:
        subset = blocks[blocks["gene"].eq(gene)]
        groups = {r: subset.loc[subset["region"].eq(r), "block_mean_expression"].to_numpy() for r in REGIONS}
        if all(len(v) > 0 for v in groups.values()):
            h, p = stats.kruskal(*(groups[r] for r in REGIONS))
        else:
            h, p = np.nan, np.nan
        global_rows.append({"gene": gene, "test": "Kruskal-Wallis", "statistic": h, "p_raw": p, **{f"n_blocks_{r}": len(groups[r]) for r in REGIONS}})
        for r1, r2 in combinations(REGIONS, 2):
            if len(groups[r1]) and len(groups[r2]):
                u, pair_p = stats.mannwhitneyu(groups[r1], groups[r2], alternative="two-sided")
            else:
                u, pair_p = np.nan, np.nan
            pair_rows.append({"gene": gene, "region_1": r1, "region_2": r2, "test": "two-sided Mann-Whitney U", "statistic": u, "p_raw": pair_p, "n_blocks_1": len(groups[r1]), "n_blocks_2": len(groups[r2])})

    global_df = pd.DataFrame(global_rows)
    pair_df = pd.DataFrame(pair_rows)
    for table in (global_df, pair_df):
        valid = table["p_raw"].notna()
        table["p_adjusted_bh"] = np.nan
        table.loc[valid, "p_adjusted_bh"] = multipletests(table.loc[valid, "p_raw"], method="fdr_bh")[1]
    return global_df, pair_df


def p_symbol(p: float) -> str:
    if not np.isfinite(p):
        return ""
    if p <= 0.0001:
        return "****"
    if p <= 0.001:
        return "***"
    if p <= 0.01:
        return "**"
    if p <= 0.05:
        return "*"
    return "ns"


def add_region_background(ax: plt.Axes, x_limit: float) -> None:
    ax.add_patch(Polygon([[-x_limit, 0], [-600, 0], [-300, 1], [-x_limit, 1]], color=REGION_COLORS["Un-epi"], alpha=0.15, lw=0))
    ax.add_patch(Polygon([[-600, 0], [400, 0], [700, 1], [-300, 1]], color=REGION_COLORS["Epi-Front"], alpha=0.15, lw=0))
    ax.add_patch(Polygon([[400, 0], [x_limit, 0], [x_limit, 1], [700, 1]], color=REGION_COLORS["Newly-epi"], alpha=0.15, lw=0))
    ax.plot([-600, -300], [0, 1], color="white", ls=(0, (2, 2)), lw=0.6, alpha=0.8)
    ax.plot([400, 700], [0, 1], color="white", ls=(0, (2, 2)), lw=0.6, alpha=0.8)


def draw_density_contours(ax: plt.Axes, frame: pd.DataFrame, x_limit: float) -> None:
    hist, xedges, yedges = np.histogram2d(
        frame["distance_to_front_um"], frame["relative_depth"], bins=(80, 45), range=[[-x_limit, x_limit], [0, 1]]
    )
    xc = (xedges[:-1] + xedges[1:]) / 2
    yc = (yedges[:-1] + yedges[1:]) / 2
    positive = hist[hist > 0]
    if len(positive):
        levels = np.unique(np.quantile(positive, [0.40, 0.58, 0.72, 0.84, 0.93]))
        if len(levels) >= 2:
            ax.contour(xc, yc, hist.T, levels=levels, colors="white", linewidths=0.35, alpha=0.42, zorder=2)


def draw_spatial(ax: plt.Axes, frame: pd.DataFrame, values: np.ndarray, title: str, phenotype: str, show_y: bool, show_regions: bool, show_colorbar: bool) -> None:
    x = frame["distance_to_front_um"].to_numpy()
    y = frame["relative_depth"].to_numpy()
    x_limit = float(np.quantile(np.abs(x), 0.995) * 1.05)
    add_region_background(ax, x_limit)
    draw_density_contours(ax, frame, x_limit)
    order = np.argsort(values)
    positive = values[values > 0]
    vmax = float(np.quantile(positive, 0.995)) if len(positive) else 1.0
    ax.scatter(x[order], y[order], c=np.where(values[order] > 0, values[order], np.nan), cmap=CMAP, vmin=0, vmax=vmax, s=1.15, linewidths=0, rasterized=True, zorder=3)
    zero = values <= 0
    ax.scatter(x[zero], y[zero], color="#343434", s=0.7, alpha=0.45, linewidths=0, rasterized=True, zorder=1)
    if phenotype in {"iPSC Reprogramming", "Epithelialization"}:
        ax.annotate("", xy=(1800, 0.20), xytext=(-500, 0.75), arrowprops={"arrowstyle": "-|>", "color": "#00A6D6", "lw": 1.0})
    if phenotype in {"iPSC Reprogramming", "Scar Formation"}:
        ax.annotate("", xy=(-3200, 0.97), xytext=(-500, 0.75), arrowprops={"arrowstyle": "-|>", "color": "#FFF200", "lw": 1.0})
    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(1.03, -0.03)
    ax.set_box_aspect(1)
    ax.set_title(
        title if title in PHENOTYPES else rf"$\it{{{title}}}$",
        color="white",
        y=1.115 if show_regions else 1.02,
        pad=0,
        fontweight="bold" if title in PHENOTYPES else "normal",
    )
    ax.set_xticks([-3000, 0, 3000])
    ax.set_yticks([0, 0.5, 1.0])
    ax.tick_params(length=2, width=0.6, pad=1)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    if show_y:
        ax.set_ylabel("Epidermis-Dermis axis")
    else:
        ax.set_yticklabels([])
    if show_regions:
        trans = ax.get_xaxis_transform()
        ax.text((-x_limit - 600) / 2, 1.005, "Un-epi", color=REGION_COLORS["Un-epi"], transform=trans, ha="center", va="bottom", fontsize=7)
        ax.text(-100, 1.005, "Epi-Front", color=REGION_COLORS["Epi-Front"], transform=trans, ha="center", va="bottom", fontsize=7)
        ax.text((400 + x_limit) / 2, 1.005, "Newly-epi", color=REGION_COLORS["Newly-epi"], transform=trans, ha="center", va="bottom", fontsize=7)
    if show_colorbar:
        sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax), cmap=CMAP)
        cax = ax.inset_axes([1.05, 0.12, 0.035, 0.76])
        cb = plt.colorbar(sm, cax=cax)
        cb.set_ticks([0, vmax])
        cb.set_ticklabels(["Min", "Max"])
        cb.ax.tick_params(length=0, labelsize=7, pad=1)
        cb.set_label("Expression", fontsize=8, labelpad=2)


def draw_bar(ax: plt.Axes, gene: str, blocks: pd.DataFrame, global_df: pd.DataFrame, pair_df: pd.DataFrame, show_y: bool) -> None:
    sub = blocks[blocks["gene"].eq(gene)]
    rng = np.random.default_rng(abs(hash(gene)) % (2**32))
    tops = []
    for i, region in enumerate(REGIONS):
        vals = sub.loc[sub["region"].eq(region), "block_mean_expression"].to_numpy()
        mean = float(np.mean(vals)) if len(vals) else 0
        sem = float(stats.sem(vals)) if len(vals) > 1 else 0
        ax.bar(i, mean, yerr=sem, width=0.68, color=REGION_COLORS[region], edgecolor="white", linewidth=0.5, capsize=2, error_kw={"lw": 0.7, "capthick": 0.7})
        jitter = rng.uniform(-0.16, 0.16, len(vals))
        ax.scatter(i + jitter, vals, s=5, color="white", alpha=0.42, linewidths=0, zorder=3)
        tops.append(max(vals.max() if len(vals) else 0, mean + sem))
    ax.set_xticks(range(3), ["Un-epi", "Epi-Front", "Newly-epi"], rotation=28, ha="right")
    ax.set_ylabel("Block mean expression" if show_y else "")
    ax.tick_params(length=2, width=0.6, pad=1)
    ax.spines[["top", "right"]].set_visible(False)
    gq = global_df.loc[global_df["gene"].eq(gene), "p_adjusted_bh"].iloc[0]
    ax.set_title(rf"$\it{{{gene}}}$" + f"  (KW q={gq:.2g})", fontsize=9, pad=3)

    ymax = max(tops) if tops else 1
    step = max(ymax * 0.12, 0.003)
    significant_pairs = []
    for r1, r2 in combinations(REGIONS, 2):
        row = pair_df[(pair_df["gene"].eq(gene)) & (pair_df["region_1"].eq(r1)) & (pair_df["region_2"].eq(r2))]
        if not row.empty and float(row["p_adjusted_bh"].iloc[0]) <= 0.05:
            significant_pairs.append((r1, r2, float(row["p_adjusted_bh"].iloc[0])))
    for level, (r1, r2, q) in enumerate(significant_pairs):
        y = ymax + step * (1.2 + level * 1.25)
        x1, x2 = REGIONS.index(r1), REGIONS.index(r2)
        ax.plot([x1, x1, x2, x2], [y - step * 0.15, y, y, y - step * 0.15], color="white", lw=0.65, clip_on=False)
        ax.text((x1 + x2) / 2, y + step * 0.05, p_symbol(q), ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.set_ylim(0, ymax + step * (2.8 + 1.25 * len(significant_pairs)))


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), facecolor="black")
    fig.savefig(stem.with_suffix(".svg"), facecolor="black")
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="black")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, facecolor="black", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def main() -> None:
    frame, expr, _ = load_data()
    scores = module_scores(expr)
    blocks = make_block_table(frame, expr)
    global_df, pair_df = compute_statistics(blocks)
    blocks.to_csv(DATA / "Figure3_19dpb_gene_expression_spatial_blocks.csv", index=False)
    global_df.to_csv(DATA / "Figure3_19dpb_gene_expression_Kruskal_Wallis.csv", index=False)
    pair_df.to_csv(DATA / "Figure3_19dpb_gene_expression_pairwise_Mann_Whitney.csv", index=False)
    frame.groupby("region", observed=True).size().rename("n_bins").reset_index().to_csv(DATA / "Figure3_19dpb_region_bin_counts.csv", index=False)

    # Spatial maps: same 3 x 5 organization as the supplied reference.
    fig, axes = plt.subplots(3, 5, figsize=(11.3, 7.2), gridspec_kw={"hspace": 0.38, "wspace": 0.28})
    gene_to_col = {g: i for i, g in enumerate(GENES)}
    for row, (phenotype, genes) in enumerate(PHENOTYPES.items()):
        features = genes + [phenotype]
        for col, feature in enumerate(features):
            values = scores[phenotype] if feature == phenotype else expr[:, gene_to_col[feature]]
            draw_spatial(axes[row, col], frame, values, feature, phenotype, col == 0, row == 0, col == 4)
        axes[row, 0].text(-0.35, 0.5, phenotype, transform=axes[row, 0].transAxes, rotation=90, ha="center", va="center", fontsize=10, fontweight="bold")
    fig.supxlabel("Distance to leading edge (µm)", y=0.02, fontsize=10)
    fig.subplots_adjust(left=0.075, right=0.945, bottom=0.075, top=0.945)
    save_figure(fig, OUT / "Figure3_19dpb_spatial_gene_programs_black")

    # Region-level statistical bar plots for all 12 genes.
    fig, axes = plt.subplots(3, 4, figsize=(9.2, 6.6), gridspec_kw={"hspace": 0.62, "wspace": 0.34})
    for row, (phenotype, genes) in enumerate(PHENOTYPES.items()):
        for col, gene in enumerate(genes):
            draw_bar(axes[row, col], gene, blocks, global_df, pair_df, col == 0)
        axes[row, 0].text(-0.37, 0.5, phenotype, transform=axes[row, 0].transAxes, rotation=90, ha="center", va="center", fontsize=10, fontweight="bold")
    fig.text(0.5, 0.015, "Bars: mean ± SEM of spatial-block means; points: spatial blocks", ha="center", fontsize=8)
    fig.subplots_adjust(left=0.115, right=0.99, bottom=0.11, top=0.96)
    save_figure(fig, OUT / "Figure3_19dpb_regional_expression_statistics_black")


if __name__ == "__main__":
    main()
