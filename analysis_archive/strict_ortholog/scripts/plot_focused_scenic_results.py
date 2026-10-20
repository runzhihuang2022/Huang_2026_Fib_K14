# -*- coding: utf-8 -*-
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ctxcore.genesig import GeneSignature
from ctxcore.rnkdb import FeatherRankingDatabase
from pyscenic.utils import load_motif_annotations
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"data/external/local_a82f375325f6")
RES, TAB, FIG = ROOT / "resources" / "scenic", ROOT / "tables", ROOT / "figures"
ORDER = ["Human normal", "Human wound", "Rat embryo", "Rat adult", "Reindeer velvet", "Reindeer back"]
COLORS = {"Human": "#4C78A8", "Rat": "#F58518", "Reindeer": "#72B7B2"}
mpl.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                     "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})


def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."


def motif_rank_audit():
    db = FeatherRankingDatabase(str(RES / "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"), "hg38")
    ranks = db.load(GeneSignature("IGFL1", {"IGFL1": 1.0}))
    annotations = load_motif_annotations(str(RES / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"))
    klf4 = annotations.loc["KLF4"].reset_index().drop_duplicates("MotifID")
    out = klf4.merge(ranks.reset_index().rename(columns={"motifs": "MotifID", "IGFL1": "IGFL1_rank"}), on="MotifID")
    # Feather ranks are stored as small integers; cast before multiplication to avoid overflow.
    gene_universe = 27090
    out["rank_percentile"] = 100.0 * out["IGFL1_rank"].astype(float) / gene_universe
    out = out.sort_values("IGFL1_rank")
    out.to_csv(TAB / "IGFL1_KLF4_cisTarget_motif_rank_audit.csv", index=False)
    return out


def scenic_plot():
    auc = pd.read_csv(RES / "KLF4_regulon_AUC.csv").rename(columns={"Cell": "cell_id"})
    meta = pd.read_csv(RES / "focused_FibK14_cell_metadata.csv")
    dat = meta.merge(auc, on="cell_id", how="inner")
    sample = dat.groupby(["group", "sample_id"], as_index=False)["KLF4(+)"].mean().rename(columns={"KLF4(+)": "KLF4_regulon_AUC_mean"})
    sample.to_csv(TAB / "KLF4_SCENIC_AUC_by_sample.csv", index=False)
    pairs = [("Human normal", "Human wound"), ("Rat embryo", "Rat adult"), ("Reindeer back", "Reindeer velvet")]
    tests = []
    for a, b in pairs:
        x = sample.loc[sample.group == a, "KLF4_regulon_AUC_mean"].to_numpy()
        y = sample.loc[sample.group == b, "KLF4_regulon_AUC_mean"].to_numpy()
        p = stats.mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) and len(y) else np.nan
        tests.append({"group_1": a, "group_2": b, "n_1": len(x), "n_2": len(y), "test": "two-sided Mann-Whitney U on sample means", "p_value": p})
    tdf = pd.DataFrame(tests)
    ok = tdf.p_value.notna()
    tdf.loc[ok, "FDR"] = multipletests(tdf.loc[ok, "p_value"], method="fdr_bh")[1]
    tdf["significance"] = tdf.FDR.map(lambda p: stars(p) if pd.notna(p) else "NA")
    tdf.to_csv(TAB / "KLF4_SCENIC_AUC_group_statistics.csv", index=False)

    expr = pd.read_csv(RES / "focused_FibK14_human_harmonized_expression.csv", index_col=0)
    targets = pd.read_csv(TAB / "KLF4_SCENIC_motif_pruned_targets.csv")
    positive = targets[targets.regulon.eq("KLF4(+)")].sort_values("weight", ascending=False).head(14)["target_gene"].tolist()
    genes = [g for g in ["KLF4"] + positive if g in expr.columns]
    group_for_cell = meta.set_index("cell_id").loc[expr.index, "group"]
    means = expr[genes].groupby(group_for_cell).mean().reindex(ORDER)
    z = (means - means.mean()) / means.std(ddof=0).replace(0, 1)
    z.to_csv(TAB / "KLF4_SCENIC_target_group_zscores.csv")

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.3), gridspec_kw={"width_ratios": [1.2, 1.35]})
    summaries = sample.groupby("group")["KLF4_regulon_AUC_mean"].agg(["mean", "sem"]).reindex(ORDER)
    xpos = np.arange(len(ORDER))
    colors = [COLORS["Human"]] * 2 + [COLORS["Rat"]] * 2 + [COLORS["Reindeer"]] * 2
    axes[0].bar(xpos, summaries["mean"], yerr=summaries["sem"], color=colors, edgecolor="#555555", linewidth=.6,
                error_kw={"lw": .8, "capthick": .8, "capsize": 2})
    rng = np.random.default_rng(17)
    for i, group in enumerate(ORDER):
        vals = sample.loc[sample.group == group, "KLF4_regulon_AUC_mean"].to_numpy()
        axes[0].scatter(i + rng.uniform(-.12, .12, len(vals)), vals, s=18, facecolor="white", edgecolor="#444444", lw=.7, zorder=3)
    axes[0].set_xticks(xpos, [x.replace(" ", "\n", 1) for x in ORDER], rotation=0, fontsize=7)
    axes[0].set_ylabel("KLF4 regulon AUC")
    axes[0].set_title("Fib_K14 sample-level regulon activity")
    ymax = max(sample.KLF4_regulon_AUC_mean.max(), summaries["mean"].max()) * 1.18
    axes[0].set_ylim(0, ymax * 1.12)
    for k, row in tdf.iterrows():
        i, j = ORDER.index(row.group_1), ORDER.index(row.group_2)
        y = ymax * (1 + .045 * k)
        axes[0].plot([i, i, j, j], [y-.003, y, y, y-.003], color="black", lw=.8)
        axes[0].text((i+j)/2, y, row.significance, ha="center", va="bottom", fontsize=8)

    im = axes[1].imshow(z.T.to_numpy(), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    axes[1].set_xticks(range(len(ORDER)), [x.replace(" ", "\n", 1) for x in ORDER])
    axes[1].set_yticks(range(len(genes)), genes)
    axes[1].set_title("Motif-pruned KLF4 target program")
    cb = fig.colorbar(im, ax=axes[1], fraction=.035, pad=.03)
    cb.set_label("Group mean z-score")
    fig.tight_layout(w_pad=1.2)
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"03_KLF4_SCENIC_regulon_validation.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def igfl1_plot(motif):
    adj = pd.read_csv(RES / "KLF4_human_cervid_IGFL1_adjacencies.tsv", sep="\t").sort_values("importance", ascending=False).reset_index(drop=True)
    hit = adj[adj.target.eq("IGFL1")].copy()
    hit["GRNBoost2_rank"] = hit.index + 1
    hit.to_csv(TAB / "IGFL1_KLF4_GRNBoost2_edge_audit.csv", index=False)
    targets = pd.read_csv(TAB / "KLF4_human_cervid_IGFL1_SCENIC_motif_pruned_targets.csv")
    summary = pd.DataFrame([
        {"evidence": "IGFL1 measured in human+cervid Fib_K14", "value": 1, "interpretation": "present in joint expression matrix"},
        {"evidence": "KLF4->IGFL1 GRNBoost2 edge", "value": int(len(hit) > 0), "interpretation": f"rank {int(hit.GRNBoost2_rank.iloc[0])} of {len(adj)}" if len(hit) else "absent"},
        {"evidence": "KLF4 motif near IGFL1 in cisTarget", "value": int((motif.IGFL1_rank <= 1500).any()), "interpretation": f"best rank {int(motif.IGFL1_rank.min())}; top {motif.rank_percentile.min():.2f}%"},
        {"evidence": "IGFL1 retained in motif-pruned KLF4 regulon", "value": int(targets.target_gene.str.upper().eq("IGFL1").any()), "interpretation": "NES>=3 regulon target criterion"},
    ])
    summary.to_csv(TAB / "IGFL1_KLF4_regulatory_evidence_summary.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.3, 2.3))
    y = np.arange(len(summary))[::-1]
    colors = summary.value.map({1: "#009E73", 0: "#BDBDBD"})
    ax.scatter(summary.value, y, s=120, c=colors, edgecolor="#444444", lw=.7)
    ax.set_yticks(y, summary.evidence)
    ax.set_xticks([0, 1], ["Not retained", "Supported"])
    ax.set_xlim(-.2, 2.5)
    for yi, txt in zip(y, summary.interpretation):
        ax.text(1.25, yi, txt, ha="left", va="center", fontsize=7)
    ax.set_title("IGFL1 as a candidate KLF4 effector")
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"04_IGFL1_KLF4_regulatory_evidence.{ext}", dpi=600 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    motif = motif_rank_audit()
    scenic_plot()
    igfl1_plot(motif)


if __name__ == "__main__":
    main()
