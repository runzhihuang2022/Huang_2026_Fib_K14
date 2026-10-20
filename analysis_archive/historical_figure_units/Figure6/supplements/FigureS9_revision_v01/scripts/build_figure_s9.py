from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
REP = ROOT / "reports"
QC = ROOT / "QC"
for p in [FIG, TAB, REP, QC]:
    p.mkdir(parents=True, exist_ok=True)

SPECIES_COLORS = {
    "Human": "#4C78A8",
    "Rat": "#F58518",
    "Reindeer": "#72B7B2",
    "Acomys": "#B279A2",
    "Planarian": "#9D755D",
}
DATASET_COLORS = {
    "Reindeer": SPECIES_COLORS["Reindeer"],
    "SpinyWound": SPECIES_COLORS["Acomys"],
    "SpinySteady": "#D4A6C8",
    "MouseEpi": SPECIES_COLORS["Rat"],
    "RatSpatial": "#FFB05A",
}
DATASET_LABELS = {
    "Reindeer": "Reindeer",
    "SpinyWound": "Acomys wound",
    "SpinySteady": "Acomys steady",
    "MouseEpi": "Mouse epithelial",
    "RatSpatial": "Rat spatial",
}


def bh_mw_statistics(raw: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("Human wound", "Human normal"),
        ("Rat adult", "Rat embryo"),
        ("Reindeer velvet", "Reindeer back"),
        ("Acomys wound", "Acomys unwounded"),
        ("Planarian 3dpa1", "Planarian 0hpa1"),
        ("Planarian 7dpa1", "Planarian 0hpa1"),
    ]
    rows = []
    for a, b in pairs:
        x = raw.loc[raw["condition"].eq(a), "double_positive_percent"].dropna()
        y = raw.loc[raw["condition"].eq(b), "double_positive_percent"].dropna()
        if len(x) >= 2 and len(y) >= 2:
            u, p = mannwhitneyu(x, y, alternative="two-sided")
        else:
            u, p = np.nan, np.nan
        rows.append({
            "group_a": a, "group_b": b, "n_samples_a": len(x), "n_samples_b": len(y),
            "mean_percent_a": x.mean(), "mean_percent_b": y.mean(),
            "mann_whitney_u": u, "p_value": p,
        })
    out = pd.DataFrame(rows)
    valid = out["p_value"].notna()
    out["p_adj_bh"] = np.nan
    if valid.any():
        out.loc[valid, "p_adj_bh"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    out["inference_status"] = np.where(
        out["p_adj_bh"].notna(),
        np.where(out["p_adj_bh"] < 0.05, "sample-level significant", "not significant after BH"),
        "descriptive only: insufficient sample replication",
    )
    return out


def marker_matrices(marker: pd.DataFrame):
    genes = ["PDGFRA", "KRT14", "KLF4", "IGF1", "COL1A1", "COL3A1", "POSTN", "LUM", "MMP2", "BMP2", "SOX9"]
    datasets = ["Reindeer", "SpinyWound", "SpinySteady", "MouseEpi", "RatSpatial"]
    mean = marker.pivot(index="dataset", columns="gene", values="mean_expression").reindex(index=datasets, columns=genes)
    pct = marker.pivot(index="dataset", columns="gene", values="percent_positive").reindex(index=datasets, columns=genes)
    available = marker.pivot(index="dataset", columns="gene", values="matched_symbol").reindex(index=datasets, columns=genes).notna()
    return genes, datasets, mean, pct, available


def add_panel_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=14, fontweight="bold", va="top", ha="left")


def panel_a(fig, raw, summary, stats_df):
    ax = fig.add_axes([0.075, 0.735, 0.86, 0.22])
    order = [
        "Human normal", "Human wound", "Rat embryo", "Rat adult",
        "Reindeer velvet", "Reindeer back", "Acomys wound", "Acomys unwounded",
        "Planarian 0hpa1", "Planarian 3dpa1", "Planarian 7dpa1",
    ]
    s = summary.set_index("condition").reindex(order).reset_index()
    x = np.arange(len(s))
    colors = [SPECIES_COLORS.get(v, "#777777") for v in s["condition"].map(
        lambda z: "Human" if z.startswith("Human") else "Rat" if z.startswith("Rat")
        else "Reindeer" if z.startswith("Reindeer") else "Acomys" if z.startswith("Acomys") else "Planarian"
    )]
    ax.bar(x, s["sample_mean_percent"], color=colors, edgecolor="#333333", linewidth=0.45, width=0.68)
    rng = np.random.default_rng(20260728)
    for i, condition in enumerate(order):
        vals = raw.loc[raw["condition"].eq(condition), "double_positive_percent"].dropna().to_numpy()
        if len(vals):
            ax.scatter(rng.normal(i, 0.055, len(vals)), vals, s=8, color="black", alpha=0.55, zorder=3)
    short = ["H normal", "H wound", "Rat embryo", "Rat adult", "Velvet", "Back",
             "A. wound", "A. steady", "Plan. 0 h", "Plan. 3 d", "Plan. 7 d"]
    ax.set_xticks(x, short, rotation=32, ha="right")
    ax.set_ylabel("Positive observations (%)")
    ax.set_title("Cross-context QC and Figure 6-positive fraction audit", fontweight="bold", pad=3)
    ax.set_ylim(0, max(15, np.nanmax(raw["double_positive_percent"]) * 1.08))
    h = stats_df[(stats_df.group_a == "Human wound") & (stats_df.group_b == "Human normal")].iloc[0]
    r = stats_df[(stats_df.group_a == "Rat adult") & (stats_df.group_b == "Rat embryo")].iloc[0]
    ax.text(0.01, 0.96, f"Sample-level: human BH-FDR={h.p_adj_bh:.3g}; rat BH-FDR={r.p_adj_bh:.3g}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.text(0.99, 0.96, "Bars: sample means; points: samples", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color="#444444")
    handles = [Line2D([0], [0], marker="s", linestyle="", color=c, label=k, markersize=6)
               for k, c in SPECIES_COLORS.items()]
    ax.legend(handles=handles, ncol=5, frameon=False, loc="upper right",
              bbox_to_anchor=(0.995, 0.88), handletextpad=0.35, columnspacing=0.8)


def panel_b(fig, orth):
    ax = fig.add_axes([0.075, 0.405, 0.40, 0.205])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_title("Ortholog harmonization and evidence classes", fontweight="bold", pad=2)
    boxes = [
        (0.02, "Native gene\nsymbols", "#EAF0F8"),
        (0.27, "Case-normalized\nsymbol audit", "#EAF5F3"),
        (0.52, "Strict / missing\nclassification", "#FFF3E8"),
        (0.77, "Shared marker\nmatrix", "#F4ECF4"),
    ]
    for x, txt, fc in boxes:
        patch = FancyBboxPatch((x, 0.55), 0.19, 0.20, boxstyle="round,pad=0.015,rounding_size=0.02",
                               ec="#555555", fc=fc, lw=0.7)
        ax.add_patch(patch); ax.text(x + 0.095, 0.65, txt, ha="center", va="center", fontsize=8)
    for x1, x2 in [(0.21, 0.27), (0.46, 0.52), (0.71, 0.77)]:
        ax.add_patch(FancyArrowPatch((x1, 0.65), (x2, 0.65), arrowstyle="-|>", mutation_scale=8, lw=0.8))
    count = orth.groupby(["dataset_id", "mapping_class"]).size().unstack(fill_value=0)
    strict = int((orth["mapping_class"] == "strict_symbol_match").sum())
    missing = int(orth["mapping_class"].str.contains("unmapped", na=False).sum())
    ax.text(0.02, 0.39, f"Audited records: {len(orth)}  |  strict symbol matches: {strict}  |  unresolved: {missing}",
            fontsize=8, fontweight="bold")
    ax.text(0.02, 0.25,
            "Strict Fib_K14 requires a fibroblast-lineage marker plus KRT14.\n"
            "Rat spatial and planarian data lack an auditable strict marker pair.\n"
            "Their Figure 6 signals are labelled surrogate/signature evidence.",
            fontsize=8, va="top", color="#333333")


def panel_c(fig, mean, datasets):
    ax = fig.add_axes([0.555, 0.405, 0.38, 0.205])
    x = np.log1p(mean.fillna(0).to_numpy())
    x = (x - x.mean(axis=0, keepdims=True)) / (x.std(axis=0, keepdims=True) + 1e-8)
    emb = PCA(n_components=2).fit_transform(x)
    for i, ds in enumerate(datasets):
        ax.scatter(emb[i, 0], emb[i, 1], s=85, color=DATASET_COLORS[ds], edgecolor="black", linewidth=0.5)
        dy = -0.24 if ds == "Reindeer" else 0.18
        va = "top" if ds == "Reindeer" else "bottom"
        ha = "left" if ds == "SpinyWound" else "center"
        ax.text(emb[i, 0], emb[i, 1] + dy, DATASET_LABELS[ds], ha=ha, va=va, fontsize=8)
    ax.axhline(0, color="#DDDDDD", lw=0.6); ax.axvline(0, color="#DDDDDD", lw=0.6)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title("Exploratory shared-marker embedding\nof dataset contexts", fontweight="bold", pad=2)
    ax.text(0.01, 0.02, "Dataset-level means; not a cell-level co-embedding", transform=ax.transAxes,
            color="#A33A2B", fontsize=8, va="bottom")


def panel_d(fig, genes, datasets, mean, pct, available):
    ax = fig.add_axes([0.115, 0.085, 0.345, 0.245])
    z = np.log1p(mean.copy())
    z = (z - z.mean(axis=0)) / z.std(axis=0).replace(0, np.nan)
    for yi, ds in enumerate(datasets):
        for xi, gene in enumerate(genes):
            if not available.loc[ds, gene]:
                continue
            size = 12 + float(pct.loc[ds, gene]) * 2.2
            val = float(z.loc[ds, gene]) if pd.notna(z.loc[ds, gene]) else 0
            ax.scatter(xi, yi, s=size, c=val, cmap="RdBu_r", vmin=-2, vmax=2,
                       edgecolor="#333333", linewidth=0.25)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right", fontstyle="italic")
    ax.set_yticks(range(len(datasets)), [DATASET_LABELS[x] for x in datasets])
    ax.set_xlim(-0.6, len(genes)-0.4); ax.set_ylim(len(datasets)-0.5, -0.5)
    ax.set_title("Audited marker expression across evaluable contexts", fontweight="bold", pad=2)
    ax.text(0.01, -0.34, "Dot size: % positive; color: within-gene scaled mean;\nblank: unresolved after orthology audit",
            transform=ax.transAxes, fontsize=8)


def panel_e(fig, genes, datasets, pct, available, overlap):
    ax = fig.add_axes([0.595, 0.085, 0.335, 0.245])
    programs = {
        "Lineage / epithelial": ["PDGFRA", "KRT14", "KLF4"],
        "ECM / repair": ["COL1A1", "COL3A1", "POSTN", "LUM", "MMP2"],
        "Growth / repair": ["IGF1", "BMP2", "SOX9"],
    }
    mat = pct.copy()
    mat = (mat - mat.mean(axis=0)) / mat.std(axis=0).replace(0, np.nan)
    ordered = sum(programs.values(), [])
    arr = mat[ordered].to_numpy()
    arr[~available[ordered].to_numpy()] = np.nan
    im = ax.imshow(arr, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(ordered)), ordered, rotation=45, ha="right", fontstyle="italic")
    ax.set_yticks(range(len(datasets)), [DATASET_LABELS[x] for x in datasets])
    ax.set_title("Conserved and context-dependent repair programs", fontweight="bold", pad=2)
    pos = 0
    for label, glist in programs.items():
        center = pos + (len(glist)-1)/2
        ax.text(center, -0.78, label, ha="center", va="bottom", fontsize=8, fontweight="bold")
        if pos > 0: ax.axvline(pos - 0.5, color="white", lw=1.4)
        pos += len(glist)
    cax = ax.inset_axes([1.02, 0.14, 0.025, 0.65])
    cb = plt.colorbar(im, cax=cax); cb.set_label("Scaled % positive", fontsize=8)
    best = overlap.sort_values("p_adj_bh").iloc[0]
    ax.text(0.0, -0.34,
            f"Reindeer-Acomys wound top-200 overlap: {int(best.overlap)} genes;\n"
            f"hypergeometric BH-FDR={best.p_adj_bh:.2e}.",
            transform=ax.transAxes, fontsize=8)


def save_outputs(fig):
    base = FIG / "FigureS9_cross_species_conservation_v01"
    fig.savefig(base.with_suffix(".svg"), facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=600, facecolor="white")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    return base


def main():
    raw = pd.read_csv(SRC / "cross_context_sample_positive_fraction.csv")
    summary = pd.read_csv(SRC / "cross_context_group_summary.csv")
    marker = pd.read_csv(SRC / "audited_marker_expression.csv")
    orth = pd.read_csv(SRC / "ortholog_audit.csv")
    overlap = pd.read_csv(SRC / "marker_overlap_statistics.csv")
    neg = pd.read_csv(SRC / "negative_findings.csv")
    stats_df = bh_mw_statistics(raw)
    stats_df.to_csv(TAB / "sample_level_pairwise_statistics.csv", index=False)
    genes, datasets, mean, pct, available = marker_matrices(marker)
    mean.to_csv(TAB / "marker_mean_expression_matrix.csv")
    pct.to_csv(TAB / "marker_percent_positive_matrix.csv")
    available.to_csv(TAB / "marker_availability_matrix.csv")
    orth.groupby(["dataset_id", "mapping_class"]).size().unstack(fill_value=0).to_csv(
        TAB / "ortholog_mapping_class_counts.csv")
    neg.to_csv(TAB / "negative_and_inconclusive_findings.csv", index=False)

    fig = plt.figure(figsize=(8.27, 11.69), facecolor="white")
    add_panel_letter(fig, 0.025, 0.975, "A")
    add_panel_letter(fig, 0.025, 0.625, "B")
    add_panel_letter(fig, 0.500, 0.625, "C")
    add_panel_letter(fig, 0.025, 0.350, "D")
    add_panel_letter(fig, 0.510, 0.350, "E")
    panel_a(fig, raw, summary, stats_df)
    panel_b(fig, orth)
    panel_c(fig, mean, datasets)
    panel_d(fig, genes, datasets, mean, pct, available)
    panel_e(fig, genes, datasets, pct, available, overlap)
    base = save_outputs(fig)
    plt.close(fig)

    hashes = []
    for p in sorted(SRC.glob("*")):
        hashes.append({"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
    pd.DataFrame(hashes).to_csv(QC / "source_sha256.csv", index=False)
    report = {
        "figure": "Figure S9",
        "version": "v01",
        "canvas": "A4 portrait",
        "font": "Arial 8-14 pt",
        "positive_results": [
            "Human wound positive fraction exceeds human normal at sample level (BH-FDR < 0.01).",
            "Reindeer and Acomys wound top-200 marker lists overlap by 47 genes (hypergeometric BH-FDR 4.23e-51).",
            "ECM/repair genes are auditable across all five evaluated marker matrices.",
        ],
        "boundaries": [
            "Rat adult versus embryo is a trend after sample-level BH correction (FDR about 0.068).",
            "Reindeer velvet versus back is not significant and the mean is not higher in velvet.",
            "Acomys control and planarian time points lack sample replication for sample-level inference.",
            "Panel C is a dataset-level focused marker embedding, not a cell-level co-embedding.",
            "Rat spatial and planarian signals are surrogate/signature evidence, not strict PDGFRA+KRT14 Fib_K14 calls.",
        ],
        "outputs": [str(base.with_suffix(x)) for x in [".svg", ".pdf", ".png", ".tiff"]],
    }
    (REP / "analysis_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
