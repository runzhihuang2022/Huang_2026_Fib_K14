"""Build rat Figure S5 panel B/G, copy human comparator F, and assemble.

Statistics use biological samples as replicates. Rat Krt14 is absent, so the
rat state is consistently labelled Fib_K14-like and cross-species conservation
is limited to the measurable Pdgfra/Vim/Krt5/Tacstd2 program.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
import shutil

import anndata as ad
import fitz
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import kruskal, mannwhitneyu


STATE_COLORS = {
    "KC_Basal": "#008B45", "KC_Basal_Mig": "#FF8C00",
    "KC_Basal_Prolif": "#F1C40F", "KC_Spinous": "#00D000",
    "KC_Spinous_Mig": "#00008B", "KC_Spinous_Mat": "#FF69B4",
    "KC_Granular": "#8B4513", "SAC_Progenitor": "#0099CC",
    "SAC_SG_Progenitor": "#0099CC", "SAC_SG_Dark": "#004225",
    "SAC_SG_Ductal": "#CD853F", "SAC_HF_IRS": "#DDA0DD",
    "SAC_HF_ORS": "#C71585", "Fib_K14": "#E31A1C",
}
G_MAPS = [
    ("Embryonic E18.5 5 dpw", "E18.5_5D_C02926D4_Annotated_Spatial.h5ad"),
    ("Postnatal P5 3 dpw", "P5_3D_B02621C2_Annotated_Spatial.h5ad"),
    ("Adult 19 dpw", "Adult_19d_B03424F4_Annotated_Spatial.h5ad"),
    ("Adult normal", "Adult_nor_A02988C1_Annotated_Spatial.h5ad"),
]
ORTHOLOGS = [
    ("PDGFRA", "Pdgfra"), ("VIM", "Vim"), ("KRT5", "Krt5"),
    ("TACSTD2", "Tacstd2"), ("KRT17", "Krt17"), ("COL1A1", "Col1a1"),
    ("COL3A1", "Col3a1"), ("DCN", "Dcn"), ("LUM", "Lum"),
]


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


def build_b(scrna_path: str, package: Path) -> None:
    figures = package / "figures" / "panels"
    tables = package / "tables"
    adata = ad.read_h5ad(scrna_path, backed="r")
    obs = adata.obs.copy()
    scores = pd.read_csv(tables / "S5B_scrublet_cell_scores.csv.gz")
    scores = scores.set_index("cell_id").reindex(obs.index.astype(str))
    labels = obs["sub_labels"].astype(str)
    keep = labels.str.startswith("KC_") | labels.str.startswith("Fib_")
    idx = np.flatnonzero(keep.to_numpy())
    group = np.where(
        labels.eq("Fib_K14"), "Fib_K14-like",
        np.where(labels.str.startswith("Fib_"), "Other Fib", "Keratinocyte"),
    )
    mito_genes = [g for g in adata.var_names if str(g).lower().startswith("mt-")]
    mt = adata[idx, mito_genes].X
    mt_counts = np.asarray(mt.sum(axis=1)).ravel()
    total_counts = obs.iloc[idx]["nCount_RNA.x"].to_numpy(float)
    df = pd.DataFrame({
        "Group": group[idx],
        "Number of Gene": obs.iloc[idx]["nFeature_RNA.x"].to_numpy(float),
        "Number of UMI": total_counts,
        "Mitochondrial (%)": mt_counts / np.maximum(total_counts, 1) * 100,
        "Doublet score": scores.iloc[idx]["doublet_score"].to_numpy(float),
        "Predicted doublet": scores.iloc[idx]["predicted_doublet"].to_numpy(bool),
    })
    df.to_csv(tables / "S5B_QC_doublet_plot_data.csv.gz", index=False)
    order = ["Keratinocyte", "Other Fib", "Fib_K14-like"]
    colors = ["#BDBDBD", "#4682B4", "#E31A1C"]
    rng = np.random.default_rng(17)
    plot_idx: list[int] = []
    for name in order:
        current = np.flatnonzero(df["Group"].eq(name).to_numpy())
        plot_idx.extend(rng.choice(current, min(3500, len(current)), replace=False))
    plot = df.iloc[plot_idx]
    fig, axes = plt.subplots(2, 2, figsize=(3.05, 2.75))
    for ax, metric in zip(axes.ravel(), [
        "Number of Gene", "Number of UMI", "Mitochondrial (%)", "Doublet score",
    ]):
        sns.violinplot(
            data=plot, x="Group", y=metric, order=order, ax=ax,
            hue="Group", palette=dict(zip(order, colors)), legend=False,
            cut=0, inner="quartile", linewidth=0.55,
        )
        ax.set_title(metric, pad=2); ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=28); sns.despine(ax=ax)
    frac = df.groupby("Group")["Predicted doublet"].mean().reindex(order)
    axes[1, 1].text(
        0.02, 0.98,
        f"Pred. doublets: KC {frac.iloc[0]:.1%}; other Fib {frac.iloc[1]:.1%}; "
        f"Fib_K14-like {frac.iloc[2]:.1%}",
        transform=axes[1, 1].transAxes, ha="left", va="top", fontsize=4.5,
    )
    fig.subplots_adjust(left=0.16, right=0.98, top=0.94, bottom=0.20,
                        wspace=0.48, hspace=0.62)
    save(fig, figures / "FigureS5_B_rat_QC_and_doublet_reanalysis")


def stage_from_name(name: str) -> str:
    if name.startswith(("E16.5", "E17.5", "E18.5")):
        return "Embryonic"
    if name.startswith("P5"):
        return "Postnatal"
    return "Adult"


def bh_adjust(pvalues: list[float]) -> np.ndarray:
    p = np.asarray(pvalues, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked); out[order] = np.minimum(ranked, 1)
    return out


def sample_state_expression(
    path: str, human: bool, label: str, genes: list[str]
) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    labels = adata.obs["sub_labels"].astype(str)
    idx_all = np.flatnonzero(labels.eq(label).to_numpy())
    sample_ids = adata.obs.iloc[idx_all]["sample_id"].astype(str)
    present = [g for g in genes if g in adata.var_names]
    matrix = adata[idx_all, present].X
    matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    total_col = "nCount_RNA" if human else "nCount_RNA.x"
    lib = adata.obs.iloc[idx_all][total_col].to_numpy(float)
    matrix = matrix.multiply(1e4 / np.maximum(lib, 1)[:, None]).tocsr()
    matrix.data = np.log1p(matrix.data)
    rows = []
    for sample in sorted(sample_ids.unique()):
        local = np.flatnonzero(sample_ids.eq(sample).to_numpy())
        mean = np.asarray(matrix[local].mean(axis=0)).ravel()
        detected = np.asarray((matrix[local] > 0).mean(axis=0)).ravel()
        rows.extend({
            "species": "Human" if human else "Rat",
            "sample_id": sample, "gene": gene, "mean_log1p_CPM": value,
            "fraction_expressing": fraction,
        } for gene, value, fraction in zip(present, mean, detected))
    return pd.DataFrame(rows)


def add_scale_bar(ax: plt.Axes, xy: np.ndarray, label: bool = True) -> None:
    xmin, xmax = np.min(xy[:, 0]), np.max(xy[:, 0])
    ymin, ymax = np.min(xy[:, 1]), np.max(xy[:, 1])
    x0 = xmin + 0.06 * (xmax - xmin)
    y0 = ymin + 0.08 * (ymax - ymin)
    ax.plot([x0, x0 + 1000], [y0, y0], color="black", lw=1.5, solid_capstyle="butt")
    if label:
        ax.text(x0 + 500, y0 + 0.025 * (ymax - ymin), "1 mm",
                ha="center", va="bottom", fontsize=6)


def build_g(
    spatial_root: str, human_scrna: str, rat_scrna: str, package: Path
) -> None:
    figures = package / "figures" / "panels"
    tables = package / "tables"
    all_paths = sorted(glob.glob(str(Path(spatial_root) / "**" / "*_Annotated_Spatial.h5ad"),
                                 recursive=True))
    lookup = {Path(p).name: p for p in all_paths}
    rows = []
    for path in all_paths:
        adata = ad.read_h5ad(path, backed="r")
        labels = adata.obs["sub_labels_gated"].astype(str)
        fib = labels.str.startswith("Fib_")
        rows.append({
            "sample": Path(path).name,
            "stage": stage_from_name(Path(path).name),
            "condition": "Normal" if "nor" in Path(path).name.lower() else "Wound",
            "fib_k14_like": int(labels.eq("Fib_K14").sum()),
            "fibroblasts": int(fib.sum()),
            "fraction_among_fibroblasts": float(labels.eq("Fib_K14").sum() / fib.sum()),
        })
        adata.file.close()
    fractions = pd.DataFrame(rows)
    fractions.to_csv(tables / "S5G_rat_spatial_FibK14_like_sample_fractions.csv", index=False)
    groups = {s: fractions.loc[fractions.stage.eq(s), "fraction_among_fibroblasts"]
              for s in ["Embryonic", "Postnatal", "Adult"]}
    kw = kruskal(*groups.values())
    comparisons = []
    raw_p = []
    for left, right in [
        ("Embryonic", "Postnatal"), ("Embryonic", "Adult"), ("Postnatal", "Adult")
    ]:
        test = mannwhitneyu(groups[left], groups[right], alternative="two-sided")
        comparisons.append({
            "comparison": f"{left} vs {right}",
            "n_left": len(groups[left]), "n_right": len(groups[right]),
            "median_left": groups[left].median(), "median_right": groups[right].median(),
            "U": test.statistic, "p_value": test.pvalue,
        })
        raw_p.append(test.pvalue)
    stats = pd.DataFrame(comparisons)
    stats["q_value_BH"] = bh_adjust(raw_p)
    stats["kruskal_wallis_H"] = kw.statistic
    stats["kruskal_wallis_p"] = kw.pvalue
    stats.to_csv(tables / "S5G_rat_stage_statistics.csv", index=False)

    human_genes = [x[0] for x in ORTHOLOGS]
    rat_genes = [x[1] for x in ORTHOLOGS]
    h = sample_state_expression(human_scrna, True, "Fib_K14", human_genes)
    r = sample_state_expression(rat_scrna, False, "Fib_K14", rat_genes)
    rat_to_human = {r: h for h, r in ORTHOLOGS}
    r["gene"] = r["gene"].map(rat_to_human)
    conservation = pd.concat([h, r], ignore_index=True)
    conservation.to_csv(tables / "S5G_cross_species_state_pseudobulk.csv", index=False)
    matrix = conservation.groupby(["gene", "species"])["fraction_expressing"].mean().unstack()
    matrix = matrix.reindex(human_genes)

    fig = plt.figure(figsize=(7.08, 3.55))
    grid = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.15], hspace=0.30, wspace=0.12)
    for col, (title, filename) in enumerate(G_MAPS):
        ax = fig.add_subplot(grid[0, col])
        adata = ad.read_h5ad(lookup[filename], backed="r")
        xy = np.asarray(adata.obsm["X_rot_spatial"])
        labels = adata.obs["sub_labels_gated"].astype(str).to_numpy()
        ax.scatter(xy[:, 0], xy[:, 1], s=0.06, color="#E4E8EB",
                   linewidths=0, rasterized=True)
        selected_states = [s for s in STATE_COLORS if s != "Fib_K14"]
        for state in selected_states:
            mask = labels == state
            if np.any(mask):
                ax.scatter(xy[mask, 0], xy[mask, 1], s=0.13,
                           color=STATE_COLORS[state], linewidths=0, rasterized=True)
        mask = labels == "Fib_K14"
        ax.scatter(xy[mask, 0], xy[mask, 1], s=0.20,
                   color=STATE_COLORS["Fib_K14"], linewidths=0, rasterized=True)
        add_scale_bar(ax, xy, label=True)
        ax.set_title(title, fontsize=7, pad=2)
        ax.set_aspect("equal"); ax.axis("off")
        adata.file.close()

    ax_hm = fig.add_subplot(grid[1, :2])
    sns.heatmap(
        matrix, cmap="YlGnBu", vmin=0, vmax=1, linewidths=0.4, linecolor="white",
        cbar=False, ax=ax_hm,
    )
    ax_hm.set_title("Conserved measurable Fib_K14/Fib_K14-like program")
    ax_hm.set_xlabel(""); ax_hm.set_ylabel("")
    ax_hm.tick_params(axis="x", rotation=0)
    ax_hm.tick_params(axis="y", rotation=0)
    ax_hm.text(
        0.0, -0.18,
        "Color scale: 0–1 mean fraction expressing. Rat Krt14 unavailable; KRT14 excluded.",
        transform=ax_hm.transAxes, fontsize=6, ha="left",
    )

    ax_stats = fig.add_subplot(grid[1, 2:])
    order = ["Embryonic", "Postnatal", "Adult"]
    palette = {"Embryonic": "#2CA02C", "Postnatal": "#FF8C00", "Adult": "#D62728"}
    sns.boxplot(
        data=fractions, x="stage", y="fraction_among_fibroblasts",
        order=order, hue="stage", palette=palette, legend=False,
        showfliers=False, width=0.55, linewidth=0.7, ax=ax_stats,
    )
    sns.stripplot(
        data=fractions, x="stage", y="fraction_among_fibroblasts",
        order=order, hue="stage", palette=palette, legend=False,
        size=3.0, jitter=0.16, edgecolor="black", linewidth=0.25, ax=ax_stats,
    )
    embryo_adult = stats.loc[stats.comparison.eq("Embryonic vs Adult")].iloc[0]
    ax_stats.set(title="Rat spatial Fib_K14-like fraction", xlabel="",
                 ylabel="Fib_K14-like / fibroblasts")
    ax_stats.yaxis.labelpad = 2
    ax_stats.text(
        0.02, 0.98,
        f"Kruskal–Wallis P={kw.pvalue:.3g}\n"
        f"Embryonic vs adult: P={embryo_adult.p_value:.3g}, "
        f"BH q={embryo_adult.q_value_BH:.3g}",
        transform=ax_stats.transAxes, ha="left", va="top", fontsize=5.8,
    )
    sns.despine(ax=ax_stats)
    fig.suptitle(
        "Cross-species and developmental comparison of the Fib_K14-like program",
        fontsize=8, fontweight="bold", y=0.992,
    )
    save(fig, figures / "FigureS5_G_cross_species_developmental_conservation")


def copy_f(source_root: str, package: Path) -> None:
    source = Path(source_root)
    destination = package / "figures" / "panels"
    destination.mkdir(parents=True, exist_ok=True)
    for ext in [".pdf", ".png", ".svg"]:
        src = source / f"FigureS4_F_all_human_spatial_joint_density{ext}"
        dst = destination / f"FigureS5_F_human_spatial_joint_density_reference{ext}"
        shutil.copy2(src, dst)


def place_pdf(page: fitz.Page, source: Path, rect: fitz.Rect) -> None:
    doc = fitz.open(source)
    page.show_pdf_page(rect, doc, 0, keep_proportion=True, overlay=True)
    doc.close()


def assemble(package: Path) -> None:
    panels = package / "figures" / "panels"
    outdir = package / "figures"
    width = 7.08 * 72
    height = 15.15 * 72
    doc = fitz.open()
    page = doc.new_page(width=width, height=height)
    margin, gap, cursor = 14, 7, 18
    layout = [
        ("AB", 167), ("C", 95), ("D", 100), ("E", 150), ("F", 235), ("G", 260)
    ]
    for key, row_h in layout:
        if key == "AB":
            a_w = 285
            place_pdf(page, panels / "FigureS5_A_rat_fibroblast_reclustered_umap.pdf",
                      fitz.Rect(margin, cursor, margin + a_w, cursor + row_h))
            place_pdf(page, panels / "FigureS5_B_rat_QC_and_doublet_reanalysis.pdf",
                      fitz.Rect(margin + a_w + gap, cursor, width - margin, cursor + row_h))
            page.insert_text((7, cursor + 12), "A", fontsize=14, fontname="helv")
            page.insert_text((margin + a_w + 1, cursor + 12), "B", fontsize=14, fontname="helv")
        else:
            filename = {
                "C": "FigureS5_C_rat_marker_featureplots.pdf",
                "D": "FigureS5_D_rat_PAGA.pdf",
                "E": "FigureS5_E_rat_monocle2_trajectory.pdf",
                "F": "FigureS5_F_human_spatial_joint_density_reference.pdf",
                "G": "FigureS5_G_cross_species_developmental_conservation.pdf",
            }[key]
            place_pdf(page, panels / filename,
                      fitz.Rect(margin, cursor, width - margin, cursor + row_h))
            page.insert_text((7, cursor + 12), key, fontsize=14, fontname="helv")
        cursor += row_h + 8
    pdf = outdir / "FigureS5_rat_omics_revised_v02.pdf"
    doc.save(pdf, garbage=4, deflate=True)
    doc.close()
    rendered = fitz.open(pdf)
    pix = rendered[0].get_pixmap(matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False)
    pix.save(outdir / "FigureS5_rat_omics_revised_v02_600dpi.png")
    rendered.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rat-scrna", required=True)
    parser.add_argument("--rat-spatial", required=True)
    parser.add_argument("--human-scrna", required=True)
    parser.add_argument("--human-f-source", required=True)
    parser.add_argument("--package", required=True)
    args = parser.parse_args()
    setup()
    package = Path(args.package)
    build_b(args.rat_scrna, package)
    copy_f(args.human_f_source, package)
    build_g(args.rat_spatial, args.human_scrna, args.rat_scrna, package)
    assemble(package)


if __name__ == "__main__":
    main()
