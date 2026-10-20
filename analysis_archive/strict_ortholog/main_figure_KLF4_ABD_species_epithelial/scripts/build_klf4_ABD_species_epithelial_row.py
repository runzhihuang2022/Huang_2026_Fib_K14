from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

BASE = Path(r"data/external/local_a82f375325f6")
EXT = Path(r"data/external/local_535c3822c456")
ROOT = BASE / "main_figure_KLF4_ABD_species_epithelial"
FIG, TAB, SCR, REP = [ROOT / x for x in ["figures", "tables", "scripts", "reports"]]

WIDTH_MM, HEIGHT_MM = 140, 60
SPECIES = ["Human", "Rat", "Reindeer"]
COLORS = {"Human": "#4C78A8", "Rat": "#F58518", "Reindeer": "#72B7B2"}

mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": .7,
})


def draw_tree(ax):
    tree = Phylo.read(BASE / "resources" / "sequences" / "KLF4_strict_NJ_bootstrap.newick", "newick")
    labels = {
        "Guinea_pig_proxy": "G. pig*",
        "Naked_mole-rat_proxy": "N. mole-rat*",
        "Musk_deer_(cervid_proxy)": "Musk deer*",
    }
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_func=lambda c: labels.get(c.name, c.name.replace("_", " ")) if c.is_terminal() else None,
        branch_labels=None,
    )
    ax.set_title("KLF4 orthology", pad=4)
    ax.set_xlabel("Protein distance")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.text(.98, .98, "*proxy", transform=ax.transAxes, ha="right",
            va="top", fontsize=7, color="#555555")


def aggregate_species(summary):
    d = summary[(summary.regulator == "KLF4") & summary.group.str.startswith(tuple(SPECIES))].copy()
    d["species"] = d.group.str.split().str[0]
    d["estimated_positive_n"] = d["n"] * d["positive_pct"] / 100
    rows = []
    for species in SPECIES:
        x = d[d.species == species]
        total_n = int(x.n.sum())
        rows.append({
            "species": species,
            "n_cells": total_n,
            "estimated_positive_n": float(x.estimated_positive_n.sum()),
            "pooled_positive_pct": float(x.estimated_positive_n.sum() / total_n * 100),
            "source_groups": "; ".join(x.group.tolist()),
            "aggregation": "cell-count-weighted pooled percentage across conditions",
        })
    return pd.DataFrame(rows)


def draw_expression(ax, species_summary):
    x = np.arange(len(SPECIES))
    vals = species_summary.set_index("species").loc[SPECIES, "pooled_positive_pct"].to_numpy()
    bars = ax.bar(x, vals, width=.66, color=[COLORS[s] for s in SPECIES],
                  edgecolor="#4D4D4D", linewidth=.55)
    for xi, val in zip(x, vals):
        ax.text(xi, val + 2.2, f"{val:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x, SPECIES)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_ylabel("KLF4+ Fib_K14 (%)")
    ax.set_ylim(0, 105)
    ax.set_title("KLF4 expression", pad=4)


def best_peak(peaks, gene):
    hit = peaks[peaks.target_gene.eq(gene)].sort_values("FDR").head(1)
    return None if hit.empty else hit.iloc[0]


def draw_atac(ax, peaks):
    d = peaks.copy()
    d["mlog10"] = -np.log10(d.FDR.clip(lower=np.finfo(float).tiny))
    sig = d.FDR < .05
    point_colors = np.where(~sig, "#BDBDBD", np.where(d.partial_spearman > 0, "#D55E00", "#0072B2"))
    ax.scatter(d.partial_spearman, d.mlog10, c=point_colors,
               s=np.where(sig, 14, 8), alpha=.78, linewidths=0)
    ax.axhline(-np.log10(.05), color="#555555", ls="--", lw=.7)
    ax.axvline(0, color="#555555", lw=.65)

    annotations = {
        "Krt10": (-.160, 21.5),
        "Col17a1": (-.128, 15.0),
        "Grhl3": (-.092, 10.3),
        "Tacstd2": (-.105, 1.8),
        "Krt15": (.096, 9.3),
        "Krt5": (.101, 6.5),
        "Krt17": (.091, 4.2),
        "Krt1": (.030, 3.0),
    }
    selected_rows = []
    for gene, label_xy in annotations.items():
        row = best_peak(d, gene)
        if row is None:
            continue
        selected_rows.append(row.to_dict())
        ax.annotate(
            gene,
            (row.partial_spearman, row.mlog10),
            xytext=label_xy,
            textcoords="data",
            fontsize=7.3,
            arrowprops={"arrowstyle": "-", "lw": .45, "color": "#555555"},
        )

    ax.set_xlabel("Partial Spearman rho")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title("KLF4 epithelial targets", pad=4)
    ax.set_xlim(d.partial_spearman.min() - .012, d.partial_spearman.max() + .018)
    ax.set_ylim(-.5, d.mlog10.max() * 1.12)
    ax.text(.98, .98, f"{int(sig.sum())}/{len(d)} significant",
            transform=ax.transAxes, ha="right", va="top", fontsize=7)
    ax.text(.98, .91, "Epcam: n.d.",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#555555")
    return pd.DataFrame(selected_rows)


def main():
    for folder in [FIG, TAB, SCR, REP]:
        folder.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(EXT / "tables" / "KLF4_FibK14_IGFL1_allcell_group_summary.csv")
    species_summary = aggregate_species(summary)
    peaks = pd.read_csv(BASE / "tables" / "KLF4_target_peak_partial_spearman.csv")

    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")
    gs = fig.add_gridspec(
        1, 3, left=.055, right=.985, bottom=.16, top=.93,
        width_ratios=[1.02, .86, 1.35], wspace=.48,
    )
    draw_tree(fig.add_subplot(gs[0, 0]))
    draw_expression(fig.add_subplot(gs[0, 1]), species_summary)
    selected = draw_atac(fig.add_subplot(gs[0, 2]), peaks)

    stem = FIG / "KLF4_conservation_ABD_species_epithelial_row"
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)

    species_summary.to_csv(TAB / "KLF4_species_pooled_expression.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(TAB / "KLF4_selected_epithelial_motif_linked_peaks.csv", index=False, encoding="utf-8-sig")
    shutil.copy2(Path(__file__), SCR / Path(__file__).name)
    report = {
        "size_mm": [WIDTH_MM, HEIGHT_MM],
        "panels": ["A strict orthology", "B pooled species expression", "C epithelial motif-linked peaks"],
        "layout": "one horizontal row",
        "expression_aggregation": "cell-count-weighted pooled percentage across conditions within species",
        "epcam_status": "not detected in the current 124 tested motif-linked peaks",
        "font": "Arial, 7-9 pt",
        "editable_pdf": "TrueType text, pdf.fonttype=42",
    }
    (REP / "figure_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"output_root": str(ROOT), **report}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
