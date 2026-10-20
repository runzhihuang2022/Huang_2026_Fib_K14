"""Build a one-page A4 portrait source-availability draft for Figure S4 A-Q."""

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "review_draft"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})


def panel_letter(fig, x, y, label):
    fig.text(x, y, label, fontsize=12, fontweight="bold", va="top")


def image_panel(fig, rect, path, crop=None, title=None):
    ax = fig.add_axes(rect)
    im = Image.open(path).convert("RGB")
    if crop:
        im = im.crop(crop)
    ax.imshow(im)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=8, fontweight="bold", pad=1)


def blocked(fig, rect, title, detail):
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02, 0.05), 0.96, 0.90, boxstyle="round,pad=0.01", facecolor="#F6F6F6", edgecolor="#B33A3A", linewidth=0.8))
    ax.text(0.06, 0.75, title, color="#8B1E1E", fontweight="bold", fontsize=8)
    ax.text(0.06, 0.55, detail, fontsize=7, va="top", linespacing=1.15)


def build():
    fig = plt.figure(figsize=(8.2677, 11.6929), facecolor="white")
    comp = ROOT / "assets" / "computational"
    exp = ROOT / "inputs" / "experimental"
    atac = ROOT / "inputs" / "atac"
    gap = 0.012

    # Row 1: A-D
    xs = [0.04, 0.285, 0.53, 0.775]
    for x, lab in zip(xs, "ABCD"):
        panel_letter(fig, x - 0.022, 0.98, lab)
    blocked(fig, [xs[0], 0.80, 0.22, 0.16], "TRUE SCENIC SOURCE REQUIRED", "KLF4 regulon heat map\nAUC/loom matrix not found")
    blocked(fig, [xs[1], 0.80, 0.22, 0.16], "TRUE SCENIC SOURCE REQUIRED", "KLF4 regulon UMAP\nproxy is not substituted")
    image_panel(fig, [xs[2], 0.80, 0.22, 0.16], comp / "FigureS4C_D_KLF4_correlations_epithelialization_vs_scar.png", title="Live Spearman rho")
    image_panel(fig, [xs[3], 0.80, 0.19, 0.16], comp / "Figure4B_OSKM_x_axis_gradients.png", title="OSKM expression")

    # Row 2: E-G
    xs2 = [0.04, 0.365, 0.69]
    for x, lab in zip(xs2, "EFG"):
        panel_letter(fig, x - 0.022, 0.785, lab)
    image_panel(fig, [xs2[0], 0.63, 0.30, 0.14], comp / "Figure4B_OSKM_spatial_featureplots.png", title="Spatial OSKM features")
    blocked(fig, [xs2[1], 0.63, 0.30, 0.14], "qPCR SOURCE REQUIRED", "Raw values, biological n,\nnormalization and test")
    blocked(fig, [xs2[2], 0.63, 0.27, 0.14], "WESTERN BLOT SOURCE REQUIRED", "Original blot, molecular weights,\nreplicates and quantification")

    # Row 3: H-J from source AI
    xs3 = [0.04, 0.47, 0.76]
    widths = [0.40, 0.26, 0.20]
    for x, lab in zip(xs3, "HIJ"):
        panel_letter(fig, x - 0.022, 0.615, lab)
    oe = exp / "KLF4_OE_source_400pct.png"
    if oe.exists():
        im = Image.open(oe); w, h = im.size
        image_panel(fig, [xs3[0], 0.44, widths[0], 0.16], oe, crop=(int(w*0.22), 0, w, int(h*0.32)), title="KLF4-OE time lapse")
        image_panel(fig, [xs3[1], 0.44, widths[1], 0.16], oe, crop=(int(w*0.45), int(h*0.15), w, int(h*0.40)), title="Length and aspect ratio")
        image_panel(fig, [xs3[2], 0.44, widths[2], 0.16], oe, crop=(0, 0, int(w*0.25), int(h*0.30)), title="Keratin endpoint - label check")
    else:
        for x, wd in zip(xs3, widths):
            blocked(fig, [x, 0.44, wd, 0.16], "SOURCE EXPORT REQUIRED", "Experimental AI is present; high-resolution export missing")

    # Row 4: K-M
    xs4 = [0.04, 0.365, 0.69]
    for x, lab in zip(xs4, "KLM"):
        panel_letter(fig, x - 0.022, 0.425, lab)
    blocked(fig, [xs4[0], 0.28, 0.30, 0.13], "RNA-SEQ TABLE REQUIRED", "Volcano plot cannot be rebuilt\nfrom the PowerPoint crop")
    blocked(fig, [xs4[1], 0.28, 0.30, 0.13], "POSITIVE GSEA TABLE REQUIRED", "Ranked list, NES/ES,\nP and FDR")
    blocked(fig, [xs4[2], 0.28, 0.27, 0.13], "NEGATIVE GSEA TABLE REQUIRED", "Ranked list, NES/ES,\nP and FDR")

    # Row 5: N-Q
    xs5 = [0.04, 0.285, 0.53, 0.775]
    for x, lab in zip(xs5, "NOPQ"):
        panel_letter(fig, x - 0.022, 0.265, lab)
    image_panel(fig, [xs5[0], 0.10, 0.22, 0.15], atac / "UMAP_predictedGroup.png", title="Rat RNA-label-transferred scATAC")
    blocked(fig, [xs5[1], 0.10, 0.22, 0.15], "ATAC QC SOURCE REQUIRED", "Fragment count, TSS enrichment\nand nucleosome signal")
    image_panel(fig, [xs5[2], 0.10, 0.22, 0.15], atac / "BrowserTrack_Krt14.png", title="Rat Krt14 accessibility")
    blocked(fig, [xs5[3], 0.10, 0.19, 0.15], "MOTIF OUTPUT INVALID", "Current KLF4 motif table has zero hits;\nno footprinting output")

    fig.text(0.04, 0.06, "FIGURE S4 REVIEW DRAFT - A4 portrait, Arial 8-12 pt. Panel identities A-Q are frozen; gray/red cards are unresolved source-data gates.", fontsize=7, color="#8B1E1E")
    fig.savefig(OUT / "FigureS4_evidence_safe_review.png", dpi=300, facecolor="white")
    fig.savefig(OUT / "FigureS4_evidence_safe_review_600dpi.png", dpi=600, facecolor="white")
    fig.savefig(OUT / "FigureS4_evidence_safe_review_600dpi.tiff", dpi=600, facecolor="white")
    fig.savefig(OUT / "FigureS4_evidence_safe_review.pdf", facecolor="white")
    fig.savefig(OUT / "FigureS4_evidence_safe_review.svg", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    build()
