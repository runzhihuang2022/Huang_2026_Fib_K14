# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo


BASE = Path(r"data/external/local_a82f375325f6")
EXT = Path(r"data/external/local_535c3822c456")
ROOT = BASE / "main_figure_KLF4_conservation"
FIG, TAB, SCR, REP = [ROOT / x for x in ["figures", "tables", "scripts", "reports"]]

GROUPS = ["Human normal", "Human wound", "Rat embryo", "Rat adult", "Reindeer velvet", "Reindeer back"]
DISPLAY_LABELS = ["Human\nnormal", "Human\nwound", "Rat\nembryo", "Rat\nadult", "Reindeer\nvelvet", "Reindeer\nback"]
COLORS = {"Human": "#4C78A8", "Rat": "#F58518", "Reindeer": "#72B7B2"}
GROUP_COLORS = [COLORS["Human"], COLORS["Human"], COLORS["Rat"], COLORS["Rat"], COLORS["Reindeer"], COLORS["Reindeer"]]

mpl.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
})


def ensure_dirs():
    for p in [FIG, TAB, SCR, REP]:
        p.mkdir(parents=True, exist_ok=True)


def panel_letter(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")


def draw_tree(ax):
    tree = Phylo.read(str(BASE / "resources" / "sequences" / "KLF4_strict_NJ_bootstrap.newick"), "newick")
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_func=lambda c: c.name.replace("_", " ") if c.is_terminal() else None,
        branch_labels=lambda c: int(c.confidence) if c.confidence is not None and c.confidence >= 50 else None,
    )
    ax.set_title("Strict KLF4 protein orthology", pad=7)
    ax.set_xlabel("Amino-acid distance")
    ax.set_ylabel("")
    panel_letter(ax, "A")


def draw_expression(ax, summary):
    x = np.arange(len(GROUPS))
    vals = summary.set_index("group").loc[GROUPS, "positive_pct"].to_numpy(float)
    ns = summary.set_index("group").loc[GROUPS, "n"].to_numpy(int)
    ax.bar(x, vals, color=GROUP_COLORS, edgecolor="#4D4D4D", linewidth=.6)
    for xi, value, n in zip(x, vals, ns):
        ax.text(xi, value + 2.2, f"{value:.1f}%", ha="center", va="bottom", fontsize=7)
        ax.text(xi, 2.2, f"n={n:,}", ha="center", va="bottom", fontsize=6.5,
                color="white" if value > 18 else "#333333", rotation=90)
    ax.set_xticks(x, DISPLAY_LABELS)
    ax.set_ylabel("KLF4-positive Fib_K14 cells (%)")
    ax.set_ylim(0, 105)
    ax.set_title("KLF4 expression in Fib_K14")
    panel_letter(ax, "B")


def draw_scenic(ax, auc, stats_table):
    sample = auc.copy()
    grouped = sample.groupby("group")["KLF4_regulon_AUC_mean"].agg(["mean", "sem"]).reindex(GROUPS)
    x = np.arange(len(GROUPS))
    ax.bar(x, grouped["mean"], yerr=grouped["sem"], color=GROUP_COLORS,
           edgecolor="#4D4D4D", linewidth=.6,
           error_kw={"lw": .8, "capthick": .8, "capsize": 2})
    rng = np.random.default_rng(20260713)
    for i, group in enumerate(GROUPS):
        values = sample.loc[sample.group == group, "KLF4_regulon_AUC_mean"].to_numpy()
        ax.scatter(i + rng.uniform(-.12, .12, len(values)), values, s=15,
                   facecolor="white", edgecolor="#4D4D4D", linewidth=.6, zorder=3)
    ax.set_xticks(x, DISPLAY_LABELS)
    ax.set_ylabel("KLF4 regulon AUC")
    ax.set_title("Motif-pruned KLF4 regulon activity")
    ymax = sample.KLF4_regulon_AUC_mean.max() * 1.30
    ax.set_ylim(0, ymax * 1.12)
    for k, row in stats_table.reset_index(drop=True).iterrows():
        i, j = GROUPS.index(row.group_1), GROUPS.index(row.group_2)
        y = ymax * (1 + .04 * k)
        ax.plot([i, i, j, j], [y-.003, y, y, y-.003], color="black", lw=.8)
        ax.text((i+j)/2, y+.002, row.significance, ha="center", va="bottom", fontsize=8)
    panel_letter(ax, "C")


def draw_atac(ax, peaks):
    peaks = peaks.copy()
    peaks["mlog10"] = -np.log10(peaks.FDR.clip(lower=np.finfo(float).tiny))
    significant = peaks.FDR < .05
    colors = np.where(~significant, "#BDBDBD", np.where(peaks.partial_spearman > 0, "#D55E00", "#0072B2"))
    ax.scatter(peaks.partial_spearman, peaks.mlog10, c=colors, s=np.where(significant, 18, 11), alpha=.78, linewidths=0)
    ax.axhline(-np.log10(.05), color="#555555", ls="--", lw=.8)
    ax.axvline(0, color="#555555", lw=.7)
    labels = ["Areg", "Col17a1", "Grhl3", "Tacstd2", "Krt5", "Krt17", "Col1a1", "Cd34"]
    offsets = {"Areg": (4, 3), "Col17a1": (4, 3), "Grhl3": (4, 3), "Tacstd2": (4, -9),
               "Krt5": (4, 3), "Krt17": (4, -9), "Col1a1": (4, 3), "Cd34": (4, -9)}
    for gene in labels:
        hit = peaks[peaks.target_gene.eq(gene)].sort_values("FDR").head(1)
        if len(hit):
            row = hit.iloc[0]
            ax.annotate(gene, (row.partial_spearman, row.mlog10), xytext=offsets[gene],
                        textcoords="offset points", fontsize=7)
    ax.set_xlabel("Partial Spearman rho")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("Rat Fib_K14 KLF4 motif-linked peaks")
    ax.text(.98, .98, f"60/124 peaks FDR<0.05", transform=ax.transAxes,
            ha="right", va="top", fontsize=7)
    panel_letter(ax, "D")


def write_support_tables(summary):
    inclusion = pd.DataFrame([
        ["Human normal", "include", "KLF4 detected in 60.0% of 865 Fib_K14 cells; sample-level SCENIC available"],
        ["Human wound", "include", "KLF4 detected in 49.5% of 3,869 Fib_K14 cells; sample-level SCENIC available"],
        ["Rat embryo", "include", "KLF4 detected in 84.1% of 3,850 Fib_K14 cells; strict one-to-one ortholog and SCENIC available"],
        ["Rat adult", "include", "KLF4 detected in 92.4% of 5,794 Fib_K14 cells; strict one-to-one ortholog, SCENIC and ATAC available"],
        ["Reindeer velvet", "include", "KLF4 detected in 36.5% of 1,563 Fib_K14 cells; cervid orthology and SCENIC available"],
        ["Reindeer back", "include", "KLF4 detected in 72.8% of 1,355 Fib_K14 cells; cervid orthology and SCENIC available"],
        ["Acomys wound", "exclude", "KLF4 feature unavailable in the selected 3,000-HVG object"],
        ["Acomys unwounded", "supplement only", "KLF4 detected, but only 6 Fib_K14 cells were available"],
        ["Planarian time points", "exclude", "no validated KLF4 feature/protein ortholog in the current mapping"],
    ], columns=["group", "main_figure_decision", "reason"])
    inclusion.to_csv(TAB / "KLF4_group_inclusion_audit.csv", index=False)
    summary.to_csv(TAB / "main_figure_KLF4_expression_groups.csv", index=False)

    evidence = pd.DataFrame([
        ["Strict protein orthology", "Human-rat-mouse and cervid/rodent proxies form a high-identity mammalian KLF4 tree", "yes", "strong", "direct sequence conservation"],
        ["Fib_K14 expression", "KLF4 is detected in all six retained human, rat and reindeer groups", "yes", "strong", "same target cell context across species"],
        ["Focused SCENIC", "motif-pruned KLF4 regulon activity is present in human, rat and reindeer Fib_K14", "yes", "strong", "shared regulatory program"],
        ["SCENIC group contrasts", "rat embryo/adult and reindeer back/velvet differ at FDR=0.0425; human wound/normal is n.s.", "partial", "moderate", "supports context dependence, not presence alone"],
        ["Rat chromVAR/peak links", "60 of 124 motif-linked peaks pass FDR<0.05", "indirect", "strong mechanism support", "single-species chromatin evidence"],
        ["Acomys KLF4", "wound feature unavailable; unwounded n=6", "no for main figure", "weak", "insufficient feature coverage/cell number"],
        ["Planarian KLF4", "no validated KLF4 ortholog or matched feature", "no", "unsupported", "must not be claimed as conserved"],
    ], columns=["analysis", "result", "supports_KLF4_cross_species_conservation", "evidence_strength", "interpretation"])
    evidence.to_csv(TAB / "KLF4_conservation_evidence_screening.csv", index=False)


def write_report():
    text = """# KLF4 cross-species conservation main figure

## Main conclusion

KLF4 conservation is supported at three cross-species levels: mammalian protein orthology, reproducible expression in Fib_K14 cells from human/rat/reindeer, and a motif-pruned KLF4 regulon detectable in the same cell context. Rat scATAC provides an independent chromatin-mechanism validation but is not, by itself, cross-species evidence.

## Panel interpretation

- **A, strict protein orthology:** KLF4 is highly conserved among human, rat, mouse and the included rodent/cervid references. The available Cervus KLF4 record is partial and is excluded from the full-length tree; the full-length musk-deer sequence is explicitly labeled as a cervid proxy.
- **B, Fib_K14 expression:** all six retained groups express KLF4. The bars show the proportion of KLF4-positive Fib_K14 cells. Cross-species expression magnitudes are not directly tested because normalization and assay depth differ between datasets.
- **C, SCENIC:** all three species show KLF4 regulon activity in Fib_K14. Rat embryo versus adult and reindeer back versus velvet have BH-adjusted FDR=0.0425; human normal versus wound is not significant.
- **D, ATAC/motif:** in rat Fib_K14, 60 of 124 KLF4 motif-linked target peaks pass FDR<0.05. This supports a functional chromatin mechanism and highlights epithelial/repair targets.

## Excluded groups

- Acomys wound: KLF4 was unavailable in the selected HVG feature space.
- Acomys unwounded: KLF4 was detected, but only six Fib_K14 cells were available; retained as supplementary evidence only.
- Planarian: no validated KLF4 feature or strict protein ortholog was available, so all planarian groups were removed from the KLF4 main figure.

## Recommended manuscript statement

“KLF4 showed conserved protein sequence, Fib_K14-restricted expression and motif-pruned regulon activity across human, rat and reindeer datasets. Independent rat scATAC-seq analysis further identified significant KLF4 motif accessibility and motif-linked epithelial/repair target peaks, supporting conservation of a KLF4-associated repair regulatory program.”
"""
    (REP / "KLF4_main_figure_interpretation.md").write_text(text, encoding="utf-8")


def main():
    ensure_dirs()
    group_summary = pd.read_csv(EXT / "tables" / "KLF4_FibK14_IGFL1_allcell_group_summary.csv")
    group_summary = group_summary[(group_summary.regulator == "KLF4") & group_summary.group.isin(GROUPS)].copy()
    auc = pd.read_csv(BASE / "tables" / "KLF4_SCENIC_AUC_by_sample.csv")
    scenic_stats = pd.read_csv(BASE / "tables" / "KLF4_SCENIC_AUC_group_statistics.csv")
    peaks = pd.read_csv(BASE / "tables" / "KLF4_target_peak_partial_spearman.csv")

    fig = plt.figure(figsize=(8.1, 6.7))
    gs = fig.add_gridspec(2, 2, left=.09, right=.98, bottom=.09, top=.96, wspace=.32, hspace=.38)
    draw_tree(fig.add_subplot(gs[0, 0]), )
    draw_expression(fig.add_subplot(gs[0, 1]), group_summary)
    draw_scenic(fig.add_subplot(gs[1, 0]), auc, scenic_stats)
    draw_atac(fig.add_subplot(gs[1, 1]), peaks)
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"KLF4_cross_species_conservation_main_figure.{ext}",
                    dpi=600 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    write_support_tables(group_summary)
    write_report()
    shutil.copy2(Path(__file__), SCR / Path(__file__).name)
    print({"output": str(ROOT), "groups": GROUPS, "panels": 4})


if __name__ == "__main__":
    main()
