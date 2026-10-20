"""Build Figure S4 panel B and assemble the human-only supplement."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import fitz
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--scrna", required=True)
    p.add_argument("--package", required=True)
    return p.parse_args()


def setup() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "axes.linewidth": 0.7,
    })
    sns.set_style("white")


def build_b(scrna_path: str, package: Path) -> None:
    figures = package / "figures" / "panels"
    tables = package / "tables"
    adata = ad.read_h5ad(scrna_path, backed="r")
    obs = adata.obs.copy()
    scores = pd.read_csv(tables / "S4B_scrublet_cell_scores.csv.gz")
    scores = scores.set_index("cell_id").reindex(obs.index.astype(str))
    labels = obs["sub_labels"].astype(str)
    keep = labels.str.startswith("KC_") | labels.str.startswith("Fib_")
    group = np.where(labels.eq("Fib_K14"), "Fib_K14",
                     np.where(labels.str.startswith("Fib_"), "Other Fib",
                              "Keratinocyte"))
    df = pd.DataFrame({
        "Group": group[keep.to_numpy()],
        "Number of Gene": obs.loc[keep, "nFeature_RNA"].to_numpy(),
        "Number of UMI": obs.loc[keep, "nCount_RNA"].to_numpy(),
        "Mitochondrial (%)": obs.loc[keep, "mt_percent"].to_numpy(),
        "Doublet score": scores.loc[keep.to_numpy(), "doublet_score"].to_numpy(),
        "Predicted doublet": scores.loc[keep.to_numpy(), "predicted_doublet"].to_numpy(),
    })
    df.to_csv(tables / "S4B_QC_doublet_plot_data.csv.gz", index=False)
    order = ["Keratinocyte", "Other Fib", "Fib_K14"]
    colors = ["#BDBDBD", "#4682B4", "#E31A1C"]
    rng = np.random.default_rng(17)
    plot_idx = []
    for group_name in order:
        idx = np.flatnonzero(df["Group"].eq(group_name).to_numpy())
        plot_idx.extend(rng.choice(idx, min(3500, len(idx)), replace=False))
    plot = df.iloc[plot_idx]

    fig, axes = plt.subplots(2, 2, figsize=(3.05, 2.75))
    for ax, metric in zip(axes.ravel(), [
        "Number of Gene", "Number of UMI", "Mitochondrial (%)", "Doublet score",
    ]):
        sns.violinplot(data=plot, x="Group", y=metric, order=order, ax=ax,
                       palette=colors, cut=0, inner="quartile", linewidth=0.55)
        ax.set_title(metric, pad=2)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=28)
        sns.despine(ax=ax)
    frac = df.groupby("Group")["Predicted doublet"].mean().reindex(order)
    axes[1, 1].text(
        0.02, 0.98,
        (
            f"Pred. doublets: KC {frac['Keratinocyte']:.1%}; "
            f"other Fib {frac['Other Fib']:.1%}\n"
            f"Fib_K14 {frac['Fib_K14']:.1%}"
        ),
        transform=axes[1, 1].transAxes, ha="left", va="top", fontsize=4.3,
    )
    fig.subplots_adjust(left=0.16, right=0.98, top=0.94, bottom=0.20,
                        wspace=0.48, hspace=0.62)
    stem = figures / "FigureS4_B_QC_and_doublet_reanalysis"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def place_pdf(out_page, source: Path, rect: fitz.Rect) -> None:
    doc = fitz.open(source)
    out_page.show_pdf_page(rect, doc, 0, keep_proportion=True, overlay=True)
    doc.close()


def assemble(package: Path) -> None:
    panels = package / "figures" / "panels"
    outdir = package / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    width_pt = 7.08 * 72
    height_pt = 11.45 * 72
    doc = fitz.open()
    page = doc.new_page(width=width_pt, height=height_pt)
    margin = 14
    gap = 7
    cursor = 20

    # Row 1: A+B equal height, full A4 width.
    row_h = 167
    a_w = 285
    place_pdf(page, panels / "FigureS4_A_fibroblast_reclustered_umap.pdf",
              fitz.Rect(margin, cursor, margin + a_w, cursor + row_h))
    place_pdf(page, panels / "FigureS4_B_QC_and_doublet_reanalysis.pdf",
              fitz.Rect(margin + a_w + gap, cursor, width_pt - margin, cursor + row_h))
    page.insert_text((7, cursor + 12), "A", fontsize=14, fontname="helv")
    page.insert_text((margin + a_w + 1, cursor + 12), "B", fontsize=14, fontname="helv")
    cursor += row_h + 8

    row_h = 95
    place_pdf(page, panels / "FigureS4_C_marker_featureplots.pdf",
              fitz.Rect(margin, cursor, width_pt - margin, cursor + row_h))
    page.insert_text((7, cursor + 12), "C", fontsize=14, fontname="helv")
    cursor += row_h + 8

    row_h = 85
    place_pdf(page, panels / "FigureS4_D_proxy_statistics.pdf",
              fitz.Rect(margin, cursor, width_pt - margin, cursor + row_h))
    page.insert_text((7, cursor + 12), "D", fontsize=14, fontname="helv")
    cursor += row_h + 8

    row_h = 136
    place_pdf(page, panels / "FigureS4_E_monocle2_trajectory.pdf",
              fitz.Rect(margin, cursor, width_pt - margin, cursor + row_h))
    page.insert_text((7, cursor + 12), "E", fontsize=14, fontname="helv")
    cursor += row_h + 8

    row_h = height_pt - cursor - 14
    place_pdf(page, panels / "FigureS4_F_all_human_spatial_joint_density.pdf",
              fitz.Rect(margin, cursor, width_pt - margin, cursor + row_h))
    page.insert_text((7, cursor + 12), "F", fontsize=14, fontname="helv")

    pdf = outdir / "FigureS4_human_revised_v01.pdf"
    doc.save(pdf, garbage=4, deflate=True)
    doc.close()
    rendered = fitz.open(pdf)
    pix = rendered[0].get_pixmap(matrix=fitz.Matrix(600 / 72, 600 / 72),
                                 alpha=False)
    pix.save(outdir / "FigureS4_human_revised_v01_600dpi.png")
    rendered.close()


def main() -> None:
    args = parse_args()
    package = Path(args.package)
    setup()
    build_b(args.scrna, package)
    assemble(package)


if __name__ == "__main__":
    main()
