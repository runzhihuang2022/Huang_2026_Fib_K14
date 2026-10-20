from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import linregress

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "inputs" / "raw_panels" / "from_slide3"
EXPERIMENTAL = ROOT / "inputs" / "experimental"
DATA = ROOT / "source_data"
PANELS = ROOT / "panels"
FINAL = ROOT / "final"
QA = ROOT / "qa"
for p in (PANELS, FINAL, QA):
    p.mkdir(parents=True, exist_ok=True)

COLORS = {"KC_Spinous_Mig": "#27358E", "SAC_SG_Progenitor": "#16A6B6", "Fib_K14": "#D84A4A"}
REGION_COLORS = {"Un-epi": "#D84A4A", "Epi-Front": "#3B67A2", "Newly-epi": "#D9B92E"}
REGIONS = ["Un-epi", "Epi-Front", "Newly-epi"]

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "lines.linewidth": 1.2,
})


def save_panel(fig, name):
    fig.savefig(PANELS / f"{name}.svg", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(PANELS / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(PANELS / f"{name}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def clean_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=0.6, length=2)


def panel_c():
    path = DATA / "Figure2C_trajectory_centroids.csv"
    df = pd.read_csv(path)
    df = df[df.Cell_Type.isin(COLORS)]
    fig, axes = plt.subplots(1, 2, figsize=(4.4, 2.05), constrained_layout=True)
    for state, color in COLORS.items():
        d = df[df.Cell_Type.eq(state)].sort_values("Day")
        axes[0].plot(d.Day, d.Median_X / 1000, "o-", color=color, ms=3, label=state)
        axes[1].plot(d.Day, d.Median_Y, "o-", color=color, ms=3, label=state)
    axes[0].axhline(0, color="#777777", lw=.6, ls="--")
    axes[0].set(xlabel="Healing time (days)", ylabel="Distance to front (mm)", title="Horizontal trajectory")
    axes[1].invert_yaxis()
    axes[1].set(xlabel="Healing time (days)", ylabel="Relative depth", title="Vertical trajectory")
    for ax in axes:
        ax.set_xticks([5, 12, 19, 26, 60]); clean_axis(ax)
    axes[1].legend(frameon=False, loc="best", handlelength=1.4)
    save_panel(fig, "Figure2C_orthogonal_trajectory")


def panel_e():
    df = pd.read_csv(DATA / "Figure2E_Fib_K14_region_proportions.csv").set_index("region").loc[REGIONS].reset_index()
    p = df.proportion.to_numpy(); n = df.region_bins.to_numpy()
    se = np.sqrt(p * (1 - p) / n)
    fig, ax = plt.subplots(figsize=(2.2, 2.05), constrained_layout=True)
    ax.bar(range(3), p * 100, yerr=1.96 * se * 100, color=[REGION_COLORS[r] for r in REGIONS], width=.72, capsize=2, linewidth=.5, edgecolor="black")
    for i, row in df.iterrows():
        ax.text(i, row.proportion * 100 + 1.2, f"{row.proportion*100:.1f}%", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(range(3), REGIONS, rotation=25, ha="right")
    ax.set_ylabel("Fib_K14 bins (%)")
    ax.set_title("19 dpb_p1")
    ymax = max(p * 100)
    for x1, x2, y in [(0, 1, ymax + 3.0), (1, 2, ymax + 6.0)]:
        ax.plot([x1, x1, x2, x2], [y-.45, y, y, y-.45], color="black", lw=.65, clip_on=False)
        ax.text((x1+x2)/2, y+.15, "***", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, ymax + 9.0)
    clean_axis(ax)
    save_panel(fig, "Figure2E_Fib_K14_proportion")


def panel_s2a():
    df = pd.read_csv(DATA / "Figure2B_S2AB_19dpb_p1_region_enrichment.csv")
    order = df.groupby("cell_type").log2_enrichment.max().sort_values(ascending=False).index
    mat = df.pivot(index="cell_type", columns="region", values="log2_enrichment").loc[order, REGIONS]
    fdr = df.pivot(index="cell_type", columns="region", values="fisher_fdr").loc[order, REGIONS]
    fig, ax = plt.subplots(figsize=(5.7, 4.2), constrained_layout=True)
    vmax = np.nanmax(np.abs(mat.to_numpy()))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(3), REGIONS); ax.set_yticks(range(len(mat)), mat.index)
    for i in range(len(mat)):
        for j in range(3):
            q = fdr.iloc[i, j]
            star = "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""
            if star: ax.text(j, i, star, ha="center", va="center", fontsize=5, color="black")
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.02); cb.set_label("log2 enrichment")
    ax.set_title("19 dpb_p1 regional enrichment (Fisher FDR)")
    save_panel(fig, "FigureS2A_all_states_heatmap")


def panel_s2b():
    targets = list(COLORS)
    df = pd.read_csv(DATA / "Figure2B_S2AB_19dpb_p1_region_enrichment.csv")
    fig, axes = plt.subplots(1, 3, figsize=(4.6, 1.8), sharey=True, constrained_layout=True)
    for ax, state in zip(axes, targets):
        d = df[df.cell_type.eq(state)].set_index("region").loc[REGIONS]
        ax.bar(range(3), d.log2_enrichment, color=[REGION_COLORS[r] for r in REGIONS], width=.72)
        ax.axhline(0, color="#555", lw=.6)
        ax.set_xticks(range(3), ["Un", "Front", "New"], rotation=25)
        ax.set_title(state.replace("_", "\n"), fontsize=7)
        for j, (_, row) in enumerate(d.iterrows()):
            q = row.fisher_fdr
            star = "***" if q < .001 else "**" if q < .01 else "*" if q < .05 else ""
            if star:
                y = row.log2_enrichment + (.12 if row.log2_enrichment >= 0 else -.12)
                ax.text(j, y, star, ha="center", va="bottom" if y >= 0 else "top", fontsize=6)
        clean_axis(ax)
    axes[0].set_ylabel("log2 enrichment")
    save_panel(fig, "FigureS2B_target_enrichment")


def correlation_panel(xname, yname, method, name):
    df = pd.read_csv(DATA / "FigureS2CD_EpiFront_mapping_scores.csv")
    stat = pd.read_csv(DATA / "FigureS2CD_correlations.csv")
    row = stat[(stat.state_a.eq(xname)) & (stat.state_b.eq(yname)) & (stat.method.eq(method))].iloc[0]
    rng = np.random.default_rng(11)
    show = df.iloc[rng.choice(len(df), min(3500, len(df)), replace=False)]
    fig, ax = plt.subplots(figsize=(2.35, 2.05), constrained_layout=True)
    ax.scatter(show[xname], show[yname], s=2.5, alpha=.2, color="#4D6F8D", rasterized=True)
    slope, intercept, *_ = linregress(df[xname], df[yname])
    xx = np.linspace(df[xname].min(), df[xname].max(), 100)
    ax.plot(xx, intercept + slope * xx, color="#D84A4A", lw=1)
    ax.set_xlabel(xname.replace("_", " ")); ax.set_ylabel(yname.replace("_", " "))
    symbol = "r" if method == "Pearson" else "rho"
    ax.text(.04, .96, f"{method} {symbol}={row.statistic:.3f}\nP={row.p_value:.2g}; n={int(row.n_bins):,}", transform=ax.transAxes, va="top", fontsize=6)
    clean_axis(ax); save_panel(fig, name)


def panel_s2e():
    df = pd.read_csv(DATA / "FigureS2E_nearest_neighbor_summary.csv")
    labels = ["Fib-KC", "Fib-SAC", "KC-SAC"]
    fig, ax = plt.subplots(figsize=(3.4, 2.05), constrained_layout=True)
    x = np.arange(3)
    ax.errorbar(x - .12, df.permuted_mean_um, yerr=[df.permuted_mean_um-df.perm_q025_um, df.perm_q975_um-df.permuted_mean_um], fmt="o", color="#777", capsize=3, label="Permuted")
    ax.scatter(x + .12, df.observed_mean_nn_um, color="#D84A4A", s=22, label="Observed", zorder=3)
    for i, row in df.iterrows():
        ax.text(i, row.observed_mean_nn_um + 18, f"P(shorter)={row.permutation_p:.2f}", ha="center", fontsize=6)
    ax.set_xticks(x, labels); ax.set_ylabel("Mean nearest-neighbor distance (µm)")
    ax.legend(frameon=False, ncol=2); clean_axis(ax)
    save_panel(fig, "FigureS2E_nearest_neighbor")


def panel_s2f():
    """Baseline 0-day MIF controls extracted only from the experimental composites."""
    sources = [
        (EXPERIMENTAL / "staining_slide1_composite.png", "KRT77 / SCGB2A2 / S100A7"),
        (EXPERIMENTAL / "staining_slide2_composite.png", "SOX9 / MMP3 / KRT6B"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(6.8, 2.45), constrained_layout=True)
    for i, (path, program) in enumerate(sources):
        im = Image.open(path).convert("RGB")
        w, h = im.size
        # The source composite contains two fields per condition. Use the Merge
        # channel from the prespecified S-0d and D-0d rows without recoloring.
        if i == 0:
            crops = [(0.377, .031, .501, .184), (.879, .031, 1.0, .184),
                     (.377, .535, .501, .684), (.879, .535, 1.0, .684)]
        else:
            crops = [(0.375, .031, .502, .184), (.875, .031, 1.0, .184),
                     (.375, .514, .502, .672), (.875, .514, 1.0, .672)]
        for j, box in enumerate(crops):
            px = tuple(int(v * (w if k % 2 == 0 else h)) for k, v in enumerate(box))
            axes[i, j].imshow(im.crop(px)); axes[i, j].axis("off")
            if i == 0:
                axes[i, j].set_title(["S-0d field 1", "S-0d field 2", "D-0d field 1", "D-0d field 2"][j], fontsize=6, pad=1)
        axes[i, 0].set_ylabel(program, fontsize=7, fontweight="bold", labelpad=3)
    save_panel(fig, "FigureS2F_healthy_baseline_MIF")


def add_image(fig, rect, path, crop=None, title=None):
    ax = fig.add_axes(rect)
    im = Image.open(path).convert("RGB")
    if crop is not None: im = im.crop(crop)
    ax.imshow(im); ax.axis("off")
    if title: ax.set_title(title, fontsize=8, fontweight="bold", pad=1.5)
    return ax


def add_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=14, fontweight="bold", va="top", ha="left")


def main_figure():
    fig = plt.figure(figsize=(11.6929, 8.2677), facecolor="white")
    # A: equal-width trajectory rows
    add_image(fig, [.025, .790, .57, .165], RAW / "image15.png")
    add_image(fig, [.025, .615, .57, .165], RAW / "image16.png")
    add_image(fig, [.025, .440, .57, .165], RAW / "image17.png")
    add_letter(fig, .008, .982, "A")
    # B
    add_image(fig, [.615, .445, .36, .51], RAW / "image19.png")
    add_letter(fig, .598, .982, "B")
    # C, D, E
    add_image(fig, [.025, .215, .29, .205], PANELS / "Figure2C_orthogonal_trajectory.png")
    add_letter(fig, .008, .435, "C")
    add_image(fig, [.33, .31, .39, .105], RAW / "image29.png")
    add_image(fig, [.33, .205, .39, .105], RAW / "image28.png")
    add_letter(fig, .315, .435, "D")
    add_image(fig, [.75, .215, .225, .205], PANELS / "Figure2E_Fib_K14_proportion.png")
    add_letter(fig, .733, .435, "E")
    # F: microscopy left, enlarged quantification right; original pixels are not altered.
    add_image(fig, [.025, .015, .285, .175], RAW / "image22.png", crop=(0, 30, 1570, 558), title="SAC_SG_Progenitor: KRT77 / SCGB2A2 / S100A7")
    add_image(fig, [.315, .015, .175, .175], RAW / "image22.png", crop=(1570, 0, 2202, 558))
    add_image(fig, [.510, .015, .285, .175], RAW / "image25.png", crop=(0, 30, 1570, 570), title="KC_Spinous_Mig: SOX9 / MMP3 / KRT6B")
    add_image(fig, [.800, .015, .175, .175], RAW / "image25.png", crop=(1570, 0, 2208, 570))
    add_letter(fig, .008, .205, "F")
    fig.savefig(FINAL / "Figure2_final_editable.svg", format="svg")
    fig.savefig(FINAL / "Figure2_final.pdf", format="pdf", dpi=600)
    fig.savefig(FINAL / "Figure2_final_600dpi.png", dpi=600)
    fig.savefig(FINAL / "Figure2_final_600dpi.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def supplement_figure():
    # Cell-style supplementary page: A4 portrait, 6–14 pt Arial, and
    # approximately 3–4 mm vertical gaps. No oversized page title.
    fig = plt.figure(figsize=(8.2677, 11.6929), facecolor="white")
    add_image(fig, [.035, .625, .93, .35], PANELS / "FigureS2A_all_states_heatmap.png")
    add_letter(fig, .012, .987, "A")
    add_image(fig, [.035, .500, .93, .105], PANELS / "FigureS2B_target_enrichment.png")
    add_letter(fig, .012, .617, "B")
    add_image(fig, [.035, .345, .445, .135], PANELS / "FigureS2C_Pearson.png")
    add_letter(fig, .012, .492, "C")
    add_image(fig, [.520, .345, .445, .135], PANELS / "FigureS2D_Spearman.png")
    add_letter(fig, .497, .492, "D")
    add_image(fig, [.035, .225, .93, .10], PANELS / "FigureS2E_nearest_neighbor.png")
    add_letter(fig, .012, .337, "E")
    add_image(fig, [.035, .075, .93, .13], PANELS / "FigureS2F_healthy_baseline_MIF.png")
    add_letter(fig, .012, .217, "F")
    # Panel G is kept in its legend-defined position, but remains explicitly
    # non-submission-ready until biological-replicate counts are provided.
    rect = [.035, .012, .93, .043]
    ax = fig.add_axes(rect); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="#FAFAFA", edgecolor="#8A8A8A", linestyle="--", linewidth=.8))
    ax.text(.5, .52, "Healthy versus Epi-Front quantification requires raw positive-cell counts, denominators and biological-replicate IDs.",
            ha="center", va="center", fontsize=7, color="#444444")
    add_letter(fig, .012, .067, "G")
    fig.savefig(FINAL / "FigureS2_review_editable.svg", format="svg")
    fig.savefig(FINAL / "FigureS2_review.pdf", format="pdf", dpi=600)
    fig.savefig(FINAL / "FigureS2_review_600dpi.png", dpi=600)
    fig.savefig(FINAL / "FigureS2_review_600dpi.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def write_notes():
    (ROOT / "STATISTICAL_AUDIT.md").write_text(textwrap.dedent("""
    # Figure 2 / Figure S2 statistical audit

    - Analysis object: `BW15D_1_D1_SDSDB_19dpb_Annotated.h5ad` (19dpb_p1), 47,042 Stereo-seq bins.
    - Region bins: Un-epi 23,016; Epi-Front 9,595; Newly-epi 14,431.
    - Fib_K14: 31.64%, 34.36%, and 18.94%, respectively. Pairwise Fisher tests support a higher Epi-Front proportion than both other regions.
    - KC_Spinous_Mig is positive-enriched at the Epi-Front but has a still higher proportion in Newly-epi.
    - SAC_SG_Progenitor is not enriched at the Epi-Front in this categorical-label analysis.
    - Within Epi-Front bins, Fib_K14 vs KC_Spinous_Mig mapping scores show Pearson r=-0.268; Fib_K14 vs SAC_SG_Progenitor show Spearman rho=-0.024.
    - Categorical-label nearest-neighbor distances are longer, not shorter, than label-permuted expectations. The fixed Figure S2C-E legend therefore requires author/statistician adjudication before submission.
    - Figure S2F uses the matched 0-day baseline fields (S-0d and D-0d) found in `时空组学/染色.pptx`; the experimental pixels are preserved and only the Merge fields are cropped for layout.
    - No raw healthy-vs-front positive-cell count table with denominators and biological-replicate IDs was located. Figure S2G therefore remains explicitly non-submission-ready; no significance symbol is fabricated.
    """).strip()+"\n", encoding="utf-8")


def main():
    panel_c(); panel_e(); panel_s2a(); panel_s2b()
    correlation_panel("Fib_K14", "KC_Spinous_Mig", "Pearson", "FigureS2C_Pearson")
    correlation_panel("Fib_K14", "SAC_SG_Progenitor", "Spearman", "FigureS2D_Spearman")
    panel_s2e(); panel_s2f(); main_figure(); supplement_figure(); write_notes()


if __name__ == "__main__":
    main()
