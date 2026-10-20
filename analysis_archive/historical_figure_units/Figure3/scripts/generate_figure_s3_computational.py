"""Generate Figure S3 computational validation panels from the live h5ad.

Panels A and C-H are source-generated here. Panel B reports the three QC
metrics present in the registered h5ad and explicitly flags the absent doublet
score instead of inventing one. Panel I is assembled from additional,
source-generated spatial maps using the registered wound-section orientation.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "supplements" / "figures" / "panels"
TABLES = ROOT / "supplements" / "tables"
H5AD = Path(r"data/external/local_decd995edda5")
SOURCE_ROOT = Path(
    r"data/external/local_91ac4b6618ea"
    r"\时空组学\stereoseq\pdf_output"
)
ENRICHMENT = SOURCE_ROOT / (
    r"单细胞测序对比图\Codex_scrna_publication_analysis_20260625"
    r"\compact_white_stats\tables\Fib_K14_curated_enrichment.csv"
)


FIB_ORDER = [
    "Fib_K14", "Fib_Prolif", "Fib_Myo", "Fib_CD45", "Fib_Inflama",
    "Fib_EN1", "Fib_Fasci", "Fib_SFRP2", "Fib_Papi",
]
PALETTE = dict(zip(FIB_ORDER, sns.color_palette("tab10", len(FIB_ORDER))))
GENES = [
    "KRT14", "KRT5", "KRT17", "PDGFRA", "COL1A1", "COL3A1", "DCN",
    "LUM", "MKI67", "TOP2A", "PTPRC", "CD3D", "CD79A", "NKG7",
]


def style():
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.7,
    })
    sns.set_style("white")


def save(fig, name, w, h):
    fig.set_size_inches(w, h)
    fig.savefig(OUT / f"{name}.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def bh_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def expression_matrix(adata, indices, genes):
    present = [g for g in genes if g in adata.var_names]
    x = adata[indices, present].X
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x), present


def panel_a(adata, obs):
    labels = obs["sub_labels"].astype(str)
    mask = labels.str.startswith("Fib_").to_numpy()
    xy = np.asarray(adata.obsm["X_umap"])[mask]
    lab = labels.to_numpy()[mask]
    fig, ax = plt.subplots()
    for name in FIB_ORDER:
        m = lab == name
        if m.any():
            ax.scatter(xy[m, 0], xy[m, 1], s=0.35, c=[PALETTE[name]], label=name, rasterized=True)
    ax.set(xlabel="UMAP1", ylabel="UMAP2", title="Fibroblast subclusters")
    ax.legend(markerscale=8, frameon=False, ncol=2, loc="best")
    sns.despine(ax=ax)
    save(fig, "FigureS3_A_fibroblast_umap", 3.1, 2.5)


def panel_b(obs):
    labels = obs["sub_labels"].astype(str)
    group = np.where(labels.eq("Fib_K14"), "Fib_K14",
                     np.where(labels.str.startswith("Fib_"), "Other Fib", "Keratinocyte"))
    keep = labels.eq("Fib_K14") | labels.str.startswith("Fib_") | labels.str.startswith("KC_")
    df = pd.DataFrame({
        "Group": group[keep],
        "nFeature_RNA": obs.loc[keep, "nFeature_RNA"].to_numpy(),
        "nCount_RNA": obs.loc[keep, "nCount_RNA"].to_numpy(),
        "mt_percent": obs.loc[keep, "mt_percent"].to_numpy(),
    })
    # Reproducible plotting sample only; all cells remain in summary tables.
    plot = df.groupby("Group", group_keys=False).sample(n=min(4000, df.groupby("Group").size().min()), random_state=7)
    order = ["Keratinocyte", "Other Fib", "Fib_K14"]
    fig, axes = plt.subplots(1, 4, gridspec_kw={"width_ratios": [1, 1, 1, 0.9]})
    for ax, metric, title in zip(axes[:3], ["nFeature_RNA", "nCount_RNA", "mt_percent"],
                                 ["Features", "UMIs", "Mitochondrial %"]):
        sns.violinplot(data=plot, x="Group", y=metric, order=order, ax=ax, inner="quartile", cut=0,
                       palette=["#BDBDBD", "#5AAE61", "#D73027"], linewidth=0.6)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=35)
        sns.despine(ax=ax)
    axes[3].axis("off")
    axes[3].text(0.5, 0.62, "Doublet score\nnot present in\nregistered h5ad", ha="center", va="center",
                 fontsize=7, color="#A50F15", weight="bold",
                 bbox=dict(boxstyle="round,pad=0.5", fc="#FFF5F0", ec="#CB181D", lw=0.8))
    axes[3].text(0.5, 0.25, "Panel remains pending\nraw doublet metadata", ha="center", va="center", fontsize=6)
    df.groupby("Group")[["nFeature_RNA", "nCount_RNA", "mt_percent"]].describe().to_csv(TABLES / "S3B_QC_summary.csv")
    save(fig, "FigureS3_B_QC_doublet_audit", 6.8, 2.3)


def fib_expression(adata, obs):
    labels = obs["sub_labels"].astype(str)
    mask = labels.isin(FIB_ORDER).to_numpy()
    idx = np.flatnonzero(mask)
    x, genes = expression_matrix(adata, idx, GENES)
    labs = labels.to_numpy()[idx]
    means, pct = [], []
    for name in FIB_ORDER:
        m = labs == name
        means.append(np.asarray(x[m].mean(axis=0)).ravel())
        pct.append(np.asarray((x[m] > 0).mean(axis=0)).ravel())
    means = pd.DataFrame(means, index=FIB_ORDER, columns=genes)
    pct = pd.DataFrame(pct, index=FIB_ORDER, columns=genes)
    z = (means - means.mean(axis=0)) / means.std(axis=0).replace(0, 1)
    means.to_csv(TABLES / "S3C_G_marker_means.csv")
    pct.to_csv(TABLES / "S3C_G_marker_percent.csv")
    return means, pct, z


def panel_c(z):
    fig, ax = plt.subplots()
    sns.heatmap(z, cmap="vlag", center=0, vmin=-2, vmax=2, ax=ax,
                cbar_kws={"label": "Row z score"}, linewidths=0.2, linecolor="white")
    ax.set_title("Fib_K14 identity and contamination-marker audit")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45)
    save(fig, "FigureS3_C_marker_heatmap", 5.0, 2.5)


def panel_d(adata, obs):
    labels = obs["sub_labels"].astype(str).to_numpy()
    sample = obs["sample_id"].astype(str).to_numpy()
    harmony = np.asarray(adata.obsm["X_harmony"])
    kc = np.char.startswith(labels.astype(str), "KC_")
    kc_centroid = np.median(harmony[kc, :20], axis=0)
    rows = []
    for name in FIB_ORDER:
        for sid in np.unique(sample[labels == name]):
            m = (labels == name) & (sample == sid)
            if m.sum() < 20:
                continue
            centroid = np.median(harmony[m, :20], axis=0)
            rows.append((name, sid, np.linalg.norm(centroid - kc_centroid), int(m.sum())))
    df = pd.DataFrame(rows, columns=["subcluster", "sample_id", "distance", "n_cells"])
    scale = df["distance"].median()
    df["KC_proxy"] = np.exp(-df["distance"] / scale)
    comparisons = []
    ref = df.loc[df.subcluster.eq("Fib_K14"), "KC_proxy"]
    for name in FIB_ORDER[1:]:
        other = df.loc[df.subcluster.eq(name), "KC_proxy"]
        p = stats.mannwhitneyu(ref, other, alternative="two-sided").pvalue if len(ref) and len(other) else np.nan
        comparisons.append((name, p))
    valid = [p for _, p in comparisons if np.isfinite(p)]
    adj = multipletests(valid, method="fdr_bh")[1] if valid else []
    k = 0
    records = []
    for name, p in comparisons:
        q = float(adj[k]) if np.isfinite(p) else np.nan
        k += int(np.isfinite(p))
        records.append((name, p, q, bh_stars(q) if np.isfinite(q) else "NA"))
    pd.DataFrame(records, columns=["comparison_to_Fib_K14", "p", "FDR", "stars"]).to_csv(TABLES / "S3D_proxy_statistics.csv", index=False)
    df.to_csv(TABLES / "S3D_sample_level_proxy.csv", index=False)
    fig, ax = plt.subplots()
    sns.boxplot(data=df, x="subcluster", y="KC_proxy", order=FIB_ORDER, ax=ax, showfliers=False,
                color="white", linewidth=0.7)
    sns.stripplot(data=df, x="subcluster", y="KC_proxy", order=FIB_ORDER, ax=ax, size=2,
                  palette=[PALETTE[x] for x in FIB_ORDER], jitter=0.18)
    star_map = {r[0]: r[3] for r in records}
    ymax = df.KC_proxy.max()
    for i, name in enumerate(FIB_ORDER[1:], 1):
        ax.text(i, ymax * 1.03, star_map.get(name, ""), ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, ymax * 1.15)
    ax.set(ylabel="KC proxy score", xlabel="", title="Sample-level KC proximity")
    ax.tick_params(axis="x", rotation=45)
    sns.despine(ax=ax)
    save(fig, "FigureS3_D_proxy_statistics", 4.6, 2.5)


def panel_e_f(adata, obs):
    states = ["Fib_K14", "SAC_SG_Progenitor", "KC_Basal", "KC_Basal_Mig", "KC_Spinous_Mig", "KC_Spinous"]
    labels_all = obs["sub_labels"].astype(str).to_numpy()
    rng = np.random.default_rng(17)
    chosen = []
    for state in states:
        idx = np.flatnonzero(labels_all == state)
        chosen.extend(rng.choice(idx, size=min(2500, len(idx)), replace=False))
    chosen = np.asarray(chosen)
    coords = np.asarray(adata.obsm["X_harmony"])[chosen, :20]
    labs = labels_all[chosen]
    mini = ad.AnnData(X=coords, obs=pd.DataFrame({"state": pd.Categorical(labs, categories=states)}))
    sc.pp.neighbors(mini, n_neighbors=15, use_rep="X")
    sc.tl.paga(mini, groups="state")
    conn = mini.uns["paga"]["connectivities"].toarray()
    cats = list(mini.obs["state"].cat.categories)
    graph = mini.obsp["connectivities"].tocsr()
    rows = []
    for source in ["Fib_K14", "SAC_SG_Progenitor"]:
        sm = labs == source
        row_idx = np.flatnonzero(sm)
        for target in states[2:]:
            tm = labs == target
            obs_edges = float(graph[row_idx][:, tm].sum())
            all_edges = float(graph[row_idx].sum())
            expected_p = tm.mean()
            # Effective integer edge count for an enrichment significance audit.
            n_eff = max(1, int(round(all_edges)))
            k_eff = min(n_eff, int(round(obs_edges)))
            p = stats.binomtest(k_eff, n_eff, expected_p, alternative="greater").pvalue
            rows.append((source, target, conn[cats.index(source), cats.index(target)], p))
    df = pd.DataFrame(rows, columns=["source", "target", "PAGA_connectivity", "p"])
    df["FDR"] = multipletests(df.p, method="fdr_bh")[1]
    df["stars"] = df.FDR.map(bh_stars)
    df.to_csv(TABLES / "S3E_PAGA_connectivity_statistics.csv", index=False)
    fig, ax = plt.subplots()
    sns.barplot(data=df, y="target", x="PAGA_connectivity", hue="source", ax=ax,
                palette=["#D73027", "#00A6CA"], edgecolor="black", linewidth=0.4)
    for patch, star in zip(ax.patches, list(df[df.source.eq("Fib_K14")].stars) + list(df[df.source.eq("SAC_SG_Progenitor")].stars)):
        if np.isfinite(patch.get_width()):
            ax.text(patch.get_width() + 0.01, patch.get_y() + patch.get_height()/2, star,
                    va="center", fontsize=6)
    ax.set(title="PAGA connectivity to keratinocyte states", xlabel="PAGA connectivity", ylabel="")
    ax.legend(frameon=False, title="")
    sns.despine(ax=ax)
    save(fig, "FigureS3_E_PAGA_connectivity", 4.4, 2.4)

    # DPT places the selected states on a diffusion trajectory. Root is the
    # Fib_K14 cell furthest from the KC centroid in Harmony space.
    kc_centroid = np.median(coords[np.char.startswith(labs.astype(str), "KC_")], axis=0)
    fib_idx = np.flatnonzero(labs == "Fib_K14")
    mini.uns["iroot"] = int(fib_idx[np.argmax(np.linalg.norm(coords[fib_idx] - kc_centroid, axis=1))])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.diffmap(mini)
        sc.tl.dpt(mini)
    dpt = mini.obs["dpt_pseudotime"].to_numpy()
    dpt_df = pd.DataFrame({"state": labs, "pseudotime": dpt})
    dpt_df.to_csv(TABLES / "S3F_DPT_pseudotime.csv", index=False)
    fig, axes = plt.subplots(1, 2, gridspec_kw={"width_ratios": [1.25, 1]})
    umap = np.asarray(adata.obsm["X_umap"])[chosen]
    sca = axes[0].scatter(umap[:, 0], umap[:, 1], c=dpt, cmap="viridis", s=1, rasterized=True)
    axes[0].set(xlabel="UMAP1", ylabel="UMAP2", title="DPT pseudotime")
    fig.colorbar(sca, ax=axes[0], fraction=0.046, pad=0.02)
    sns.violinplot(data=dpt_df, x="state", y="pseudotime", order=states, ax=axes[1], cut=0,
                   color="#80CDC1", linewidth=0.5)
    axes[1].set(ylabel="Pseudotime", xlabel="", title="State position")
    axes[1].tick_params(axis="x", rotation=45)
    sns.despine(ax=axes[0]); sns.despine(ax=axes[1])
    fig.subplots_adjust(wspace=0.45)
    save(fig, "FigureS3_F_pseudotime", 6.4, 2.5)


def panel_g(means, pct, z):
    rows = []
    for state in FIB_ORDER:
        for gene in means.columns:
            rows.append((state, gene, float(z.loc[state, gene]), float(pct.loc[state, gene])))
    df = pd.DataFrame(rows, columns=["state", "gene", "z", "pct"])
    fig, ax = plt.subplots()
    xmap = {g: i for i, g in enumerate(means.columns)}
    ymap = {s: i for i, s in enumerate(FIB_ORDER)}
    sca = ax.scatter(df.gene.map(xmap), df.state.map(ymap), s=8 + df.pct * 65,
                     c=df.z, cmap="viridis", vmin=-2, vmax=2, edgecolors="none")
    ax.set_xticks(range(len(xmap)), list(xmap), rotation=45, ha="right")
    ax.set_yticks(range(len(ymap)), list(ymap))
    ax.invert_yaxis()
    ax.set_title("Epithelial, fibroblast, proliferation and contamination markers")
    fig.colorbar(sca, ax=ax, label="Average expression z score", fraction=0.03, pad=0.02)
    sns.despine(ax=ax, left=False, bottom=False)
    save(fig, "FigureS3_G_epithelialization_dotplot", 5.5, 2.5)


def panel_h():
    df = pd.read_csv(ENRICHMENT)
    # tolerate previously used column names
    term_col = next(c for c in df.columns if c.lower() in {"term", "pathway", "description"})
    p_col = next(c for c in df.columns if "adj" in c.lower() or "fdr" in c.lower())
    df = df.copy()
    df["score"] = -np.log10(df[p_col].clip(lower=np.finfo(float).tiny))
    df = df.sort_values("score")
    fig, ax = plt.subplots()
    colors = sns.color_palette("viridis", len(df))
    ax.barh(df[term_col], df.score, color=colors, edgecolor="black", linewidth=0.4)
    for y, (_, row) in enumerate(df.iterrows()):
        ax.text(row.score + 0.08, y, bh_stars(row[p_col]), va="center", fontsize=6)
    ax.set(xlabel="-log10(FDR)", ylabel="", title="Fib_K14-enriched processes")
    sns.despine(ax=ax)
    save(fig, "FigureS3_H_enrichment", 4.4, 2.5)


def panel_i():
    files = [
        SOURCE_ROOT / r"CellType_Individual_White\Fib_K14\Fib_K14_19dpb_p1.png",
        SOURCE_ROOT / r"CellType_Individual_White\Fib_K14\Fib_K14_26dpb_p1.png",
        SOURCE_ROOT / r"CellType_Individual_White\Fib_K14\Fib_K14_26dpb_p2.png",
    ]
    fig, axes = plt.subplots(1, 3)
    for ax, path, title in zip(axes, files, ["19 dpb p1", "26 dpb p1", "26 dpb p2"]):
        ax.imshow(plt.imread(path))
        ax.set_title(title)
        ax.axis("off")
    save(fig, "FigureS3_I_additional_spatial_maps", 7.0, 2.2)


def main():
    print("S3: initialize", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    style()
    adata = ad.read_h5ad(H5AD, backed="r")
    print("S3: h5ad open", flush=True)
    obs = adata.obs.copy()
    panel_a(adata, obs)
    print("S3: A", flush=True)
    panel_b(obs)
    print("S3: B", flush=True)
    means, pct, z = fib_expression(adata, obs)
    print("S3: expression", flush=True)
    panel_c(z)
    print("S3: C", flush=True)
    panel_d(adata, obs)
    print("S3: D", flush=True)
    panel_e_f(adata, obs)
    print("S3: E-F", flush=True)
    panel_g(means, pct, z)
    print("S3: G", flush=True)
    panel_h()
    print("S3: H", flush=True)
    panel_i()
    print("S3: I", flush=True)
    print(OUT)


if __name__ == "__main__":
    main()
