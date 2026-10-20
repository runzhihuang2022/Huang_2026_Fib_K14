"""Build an evidence-safe A4 review draft for Figure 4.

This script intentionally leaves unsupported panels as source-data placeholders.
It must not be used to claim that blocked analyses have been completed.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "review_draft"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def add_panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=14, fontweight="bold", va="top", ha="left")


def image_panel(fig, rect, path, crop=None, title=None):
    ax = fig.add_axes(rect)
    img = Image.open(path).convert("RGB")
    if crop is not None:
        img = img.crop(crop)
    ax.imshow(img)
    ax.set_axis_off()
    if title:
        ax.set_title(title, pad=2, fontweight="bold")
    return ax


def placeholder(fig, rect, title, body):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box = FancyBboxPatch(
        (0.02, 0.05),
        0.96,
        0.90,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        facecolor="#F5F5F5",
        edgecolor="#B33A3A",
        linewidth=1.0,
    )
    ax.add_patch(box)
    ax.text(0.05, 0.78, title, fontsize=9, fontweight="bold", color="#8B1E1E")
    ax.text(0.05, 0.58, body, fontsize=7.5, va="top", linespacing=1.25)
    return ax


def working_model(fig, rect):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(Rectangle((0.04, 0.12), 0.24, 0.72, facecolor="#E8F0F7", edgecolor="#4C78A8"))
    ax.add_patch(Rectangle((0.40, 0.33), 0.20, 0.32, facecolor="#FFF1D6", edgecolor="#E39C37"))
    ax.add_patch(Rectangle((0.72, 0.12), 0.24, 0.72, facecolor="#E6F5EC", edgecolor="#2A8C62"))
    ax.text(0.16, 0.60, "Fib_K14\nstate", ha="center", va="center", fontweight="bold")
    ax.text(0.50, 0.49, "KLF4", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(0.84, 0.61, "Altered\nmorphology", ha="center", va="center", fontweight="bold")
    ax.text(0.84, 0.33, "Keratin-associated\nfeatures", ha="center", va="center", fontsize=7.5)
    ax.annotate("", xy=(0.40, 0.49), xytext=(0.28, 0.49), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.annotate("", xy=(0.72, 0.49), xytext=(0.60, 0.49), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.annotate(
        "undefined upstream signal",
        xy=(0.50, 0.67),
        xytext=(0.50, 0.94),
        ha="center",
        fontsize=7,
        arrowprops=dict(arrowstyle="->", lw=1.0, linestyle=":"),
    )
    ax.text(0.50, 0.03, "Evidence-bounded model; fate conversion not asserted", ha="center", fontsize=6.8)


def build():
    # A4 portrait; all inter-panel gaps are approximately 3-4 mm.
    fig = plt.figure(figsize=(8.2677, 11.6929), facecolor="white")

    comp = ROOT / "assets" / "computational"
    exp = ROOT / "inputs" / "experimental"
    atac = ROOT / "inputs" / "atac"

    # Row 1
    add_panel_letter(fig, 0.025, 0.985, "A")
    image_panel(fig, [0.055, 0.755, 0.43, 0.215], comp / "Figure4A_KLF4_prioritization_spatial_gradient.png")
    add_panel_letter(fig, 0.505, 0.985, "B")
    image_panel(fig, [0.535, 0.805, 0.44, 0.155], comp / "Figure4B_OSKM_spatial_featureplots.png")
    image_panel(fig, [0.535, 0.755, 0.44, 0.065], comp / "Figure4B_OSKM_x_axis_gradients.png")

    # Row 2
    add_panel_letter(fig, 0.025, 0.735, "C")
    image_panel(
        fig,
        [0.055, 0.575, 0.43, 0.145],
        comp / "FigureS4A_B_KLF4_regulon_proxy_FibK14_vs_fibroblasts.png",
        title="KLF4 expression and target-gene proxy",
    )
    add_panel_letter(fig, 0.505, 0.735, "D")
    if_src = exp / "IF_source_400pct.png"
    if if_src.exists():
        im = Image.open(if_src)
        w, h = im.size
        image_panel(fig, [0.535, 0.575, 0.44, 0.145], if_src, crop=(0, 0, w, int(h * 0.35)))
    else:
        placeholder(fig, [0.535, 0.575, 0.44, 0.145], "SOURCE EXPORT REQUIRED", "Export assets/experimental/1.ai at high resolution.")

    # Row 3
    add_panel_letter(fig, 0.025, 0.555, "E")
    oe_src = exp / "KLF4_OE_source_400pct.png"
    if oe_src.exists():
        im = Image.open(oe_src)
        w, h = im.size
        image_panel(fig, [0.055, 0.380, 0.60, 0.16], oe_src, crop=(int(w * 0.22), 0, w, int(h * 0.32)))
    else:
        placeholder(fig, [0.055, 0.380, 0.60, 0.16], "SOURCE EXPORT REQUIRED", "Export assets/experimental/3.ai at high resolution.")
    add_panel_letter(fig, 0.675, 0.555, "F")
    if oe_src.exists():
        im = Image.open(oe_src)
        w, h = im.size
        image_panel(fig, [0.705, 0.380, 0.27, 0.16], oe_src, crop=(0, 0, int(w * 0.25), int(h * 0.30)))
    else:
        placeholder(fig, [0.705, 0.380, 0.27, 0.16], "SOURCE EXPORT REQUIRED", "KRT label check pending.")

    # Row 4
    add_panel_letter(fig, 0.025, 0.360, "G")
    placeholder(
        fig,
        [0.055, 0.205, 0.28, 0.14],
        "SOURCE DATA REQUIRED",
        "RNA-seq DEG table + ranked list + GSEA results.\nPPT crop intentionally excluded.",
    )
    add_panel_letter(fig, 0.355, 0.360, "H")
    image_panel(fig, [0.385, 0.245, 0.18, 0.10], atac / "UMAP_predictedGroup.png", title="Rat scATAC labels")
    image_panel(fig, [0.575, 0.245, 0.18, 0.10], atac / "UMAP_Klf4_gene_score.png", title="Klf4 gene score")
    image_panel(fig, [0.385, 0.205, 0.37, 0.045], atac / "BrowserTrack_Krt14.png")
    add_panel_letter(fig, 0.775, 0.360, "I")
    working_model(fig, [0.805, 0.205, 0.17, 0.14])

    fig.text(
        0.03,
        0.175,
        "REVIEW DRAFT - blocked panels are deliberately not synthesized from PowerPoint. Arial 8-14 pt; A4 portrait; panel identities follow the formal A-I legend.",
        fontsize=7,
        color="#8B1E1E",
    )
    fig.savefig(OUT / "Figure4_evidence_safe_review.png", dpi=300, facecolor="white")
    fig.savefig(OUT / "Figure4_evidence_safe_review_600dpi.png", dpi=600, facecolor="white")
    fig.savefig(OUT / "Figure4_evidence_safe_review_600dpi.tiff", dpi=600, facecolor="white")
    fig.savefig(OUT / "Figure4_evidence_safe_review.pdf", facecolor="white")
    fig.savefig(OUT / "Figure4_evidence_safe_review.svg", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    build()
