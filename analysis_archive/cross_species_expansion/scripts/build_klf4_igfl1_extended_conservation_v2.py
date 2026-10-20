from __future__ import annotations

import json
import shutil
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests

WORK = Path(__file__).resolve().parent
sys.path.insert(0, str(WORK))
from build_fibk14_density_reference_style import marker_call, subset_condition  # noqa: E402
from build_fibk14_umap_plus_conservation_fullwidth import specs11  # noqa: E402
from build_panelA_annotation_human_rat_flexible_others_planarian_alias import (  # noqa: E402
    fibk14_annotation_call,
)

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

ROOT = Path(r"data/external/local_535c3822c456")
FIG, TAB, SCR, REP = [ROOT / x for x in ["figures", "tables", "scripts", "reports"]]
BACKUP = ROOT / "figures_before_context_stats_update"

ORDER = [
    "Human normal",
    "Human wound",
    "Rat embryo",
    "Rat adult",
    "Reindeer velvet",
    "Reindeer back",
    "Acomys wound",
    "Acomys unwounded",
    "Planarian 0 h",
    "Planarian 36 h",
    "Planarian 3 d",
    "Planarian 5 d",
]

MODULES = {
    "KLF4 repair proxy": [
        "CDKN1A", "KRT14", "KRT5", "KRT17", "EPCAM", "TACSTD2", "CDH1",
        "CLDN1", "DSP", "JUP", "ITGA6", "ITGB4", "TP63", "NOTCH1", "IVL",
    ],
    "ECM niche": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "FN1", "POSTN", "LUM", "DCN", "MMP2"],
    "Growth-factor niche": ["FGF7", "FGF10", "HGF", "IGF1", "TGFB1", "PDGFA", "VEGFA", "WNT5A"],
    "Inflammatory niche": ["IL6", "CXCL8", "CXCL2", "CCL2", "TNF", "IL1B", "PTGS2"],
    "Epithelial repair": ["KRT14", "KRT5", "KRT17", "EPCAM", "TACSTD2", "TP63", "ITGA6", "ITGB4"],
}


def ensure_dirs() -> None:
    for p in [FIG, TAB, SCR, REP, BACKUP]:
        p.mkdir(parents=True, exist_ok=True)
    if not any(BACKUP.iterdir()):
        for p in FIG.glob("*.*"):
            shutil.copy2(p, BACKUP / p.name)


def display_label(label: str) -> str:
    return {
        "Acomys control": "Acomys unwounded",
        "Planarian 0hpa1": "Planarian 0 h",
        "Planarian 36hpa1": "Planarian 36 h",
        "Planarian 3dpa1": "Planarian 3 d",
        "Planarian 5dpa1": "Planarian 5 d",
    }.get(label, label)


def star(p: float) -> str:
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def sample_column(a: ad.AnnData) -> str:
    for c in ["sample_id", "samples", "sample", "h5ad_sample", "orig.ident", "orig.ident.x", "source_accession"]:
        if c in a.obs:
            return c
    return ""


def gene_lookup(a: ad.AnnData) -> dict[str, str]:
    lookup = {str(x).upper(): str(x) for x in a.var_names}
    for c in ["human_gene_symbol", "gene_symbol", "original_gene_name"]:
        if c not in a.var:
            continue
        for idx, value in zip(a.var_names, a.var[c].astype(str)):
            if value and value.lower() != "nan":
                lookup.setdefault(value.upper(), str(idx))
    return lookup


def read_vectors(a: ad.AnnData, requested: list[str]) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    lookup = gene_lookup(a)
    matched = {g: lookup[g.upper()] for g in requested if g.upper() in lookup}
    ordered = list(dict.fromkeys(matched.values()))
    if not ordered:
        return {}, matched
    x = a[:, ordered].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    col = {g: i for i, g in enumerate(ordered)}
    return {g: np.asarray(x[:, col[s]]).ravel() for g, s in matched.items()}, matched


def contexts(a: ad.AnnData, spec) -> tuple[np.ndarray, np.ndarray, str]:
    keep = subset_condition(a, spec) if spec.cohort != "Planarian" else np.ones(a.n_obs, dtype=bool)
    if spec.cohort == "Planarian":
        gid = "SMESG000046289.1"
        if gid in a.var_names:
            x = a[:, gid].X
            x = np.asarray(x.toarray()).ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
            fib = x > 0
            source = "KRT14-like direct ortholog surrogate"
        else:
            fib = np.zeros(a.n_obs, dtype=bool)
            source = "KRT14-like ortholog unavailable"
    elif spec.label.startswith(("Human", "Rat embryo", "Rat adult")):
        fib, _ = fibk14_annotation_call(a)
        source = "existing Fib_K14 annotation"
    elif "Fib_K14" in a.obs:
        fib = a.obs["Fib_K14"].astype(str).to_numpy() == "Fib_K14"
        source = "existing Fib_K14 annotation"
    else:
        fib = marker_call(a, allow_krt5_fallback=False)[0]
        source = "PDGFRA/KRT14 flexible double-positive call"
    return keep & fib, keep, source


def module_score(vectors: dict[str, np.ndarray], genes: list[str], idx: np.ndarray, exclude: str = "") -> tuple[np.ndarray | None, int]:
    avail = [vectors[g][idx] for g in genes if g in vectors and g != exclude]
    if not avail:
        return None, 0
    return np.mean(np.vstack(avail), axis=0), len(avail)


def audit_status(group: str, n_vars: int, matched: dict[str, str], gene: str) -> tuple[str, str]:
    if gene in matched:
        return "direct/mapped feature", f"Matched as {matched[gene]}"
    if group.startswith("Rat") and gene == "IGFL1":
        return "ortholog unresolved", "Igfl1 is absent from the 25,726-feature object; Ensembl has no rat IGFL1 ortholog in the queried release"
    if group.startswith("Acomys"):
        return "HVG-limited absence", f"Absent from the {n_vars:,}-HVG object; no raw matrix is embedded, so biological absence cannot be concluded"
    if group.startswith("Planarian"):
        return "mapping unresolved", "No IGFL1/KLF4 human-symbol match in the selected feature space; local planarian-to-human mapping has no IGFL1 candidate"
    return "feature-space absence", f"No {gene} symbol/alias match in the available feature space"


def summarize() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows, availability, correlations, cell_summaries = [], [], [], []
    requested = sorted(set(["KLF4", "IGFL1"] + sum(MODULES.values(), [])))
    rng = np.random.default_rng(20260713)
    for spec in specs11():
        group = display_label(spec.label)
        a = ad.read_h5ad(spec.path, backed="r")
        fib_idx, all_idx, annotation_source = contexts(a, spec)
        vectors, matched = read_vectors(a, requested)
        sc = sample_column(a)
        sid = a.obs[sc].astype(str).to_numpy() if sc else np.repeat("pooled", a.n_obs)
        for gene, context_name in [("KLF4", "Fib_K14"), ("IGFL1", "all cells/spots")]:
            status, reason = audit_status(group, a.n_vars, matched, gene)
            availability.append(
                {
                    "group": group,
                    "cohort": spec.cohort,
                    "gene": gene,
                    "analysis_context": context_name,
                    "available": gene in vectors,
                    "matched_feature": matched.get(gene, ""),
                    "feature_space_n": a.n_vars,
                    "audit_status": status,
                    "audit_reason": reason,
                    "annotation_source": annotation_source,
                    "source": str(spec.path),
                }
            )
        contexts_by_reg = {"KLF4": fib_idx, "IGFL1": all_idx}
        for sample in np.unique(sid[all_idx]):
            rec = {"group": group, "cohort": spec.cohort, "sample": sample}
            for regulator, base_idx in contexts_by_reg.items():
                idx = base_idx & (sid == sample)
                rec[f"{regulator}_context"] = "Fib_K14" if regulator == "KLF4" else "all cells/spots"
                rec[f"{regulator}_n"] = int(idx.sum())
                if regulator in vectors and idx.sum():
                    rec[f"{regulator}_mean"] = float(vectors[regulator][idx].mean())
                    rec[f"{regulator}_positive_pct"] = float((vectors[regulator][idx] > 0).mean() * 100)
                else:
                    rec[f"{regulator}_mean"] = np.nan
                    rec[f"{regulator}_positive_pct"] = np.nan
                for mod, gl in MODULES.items():
                    score, ng = module_score(vectors, gl, idx)
                    rec[f"{regulator} context | {mod}"] = float(score.mean()) if score is not None and len(score) else np.nan
                    rec[f"{regulator} context | {mod} n_genes"] = ng
            rows.append(rec)

        for regulator, base_idx in contexts_by_reg.items():
            if regulator not in vectors:
                continue
            idx = np.where(base_idx)[0]
            if len(idx) > 5000:
                idx = rng.choice(idx, 5000, replace=False)
            for mod, gl in MODULES.items():
                score, ng = module_score(vectors, gl, idx, exclude=regulator)
                if score is None or len(idx) < 3 or np.nanstd(vectors[regulator][idx]) == 0 or np.nanstd(score) == 0:
                    rho, p = np.nan, np.nan
                else:
                    rho, p = stats.spearmanr(vectors[regulator][idx], score)
                correlations.append(
                    {
                        "group": group,
                        "regulator": regulator,
                        "analysis_context": "Fib_K14" if regulator == "KLF4" else "all cells/spots",
                        "module": mod,
                        "rho": rho,
                        "p_value": p,
                        "n_cells_or_spots": len(idx),
                        "n_module_genes": ng,
                    }
                )
        for regulator, base_idx in contexts_by_reg.items():
            if regulator in vectors and base_idx.sum():
                cell_summaries.append(
                    {
                        "group": group,
                        "regulator": regulator,
                        "context": "Fib_K14" if regulator == "KLF4" else "all cells/spots",
                        "n": int(base_idx.sum()),
                        "mean_expression": float(vectors[regulator][base_idx].mean()),
                        "positive_pct": float((vectors[regulator][base_idx] > 0).mean() * 100),
                    }
                )
        try:
            a.file.close()
        except Exception:
            pass

    d = pd.DataFrame(rows)
    av = pd.DataFrame(availability)
    co = pd.DataFrame(correlations)
    cs = pd.DataFrame(cell_summaries)
    co["p_adj_bh"] = np.nan
    for regulator in ["KLF4", "IGFL1"]:
        m = co.regulator.eq(regulator) & co.p_value.notna()
        if m.any():
            co.loc[m, "p_adj_bh"] = multipletests(co.loc[m, "p_value"], method="fdr_bh")[1]
    co["significance"] = co.p_adj_bh.map(star)
    d.to_csv(TAB / "sample_level_expression_and_modules.csv", index=False, encoding="utf-8-sig")
    av.to_csv(TAB / "KLF4_IGFL1_feature_space_and_orthology_audit.csv", index=False, encoding="utf-8-sig")
    co.to_csv(TAB / "regulator_niche_module_correlations.csv", index=False, encoding="utf-8-sig")
    cs.to_csv(TAB / "KLF4_FibK14_IGFL1_allcell_group_summary.csv", index=False, encoding="utf-8-sig")
    return d, av, co, cs


def effects(d: pd.DataFrame) -> pd.DataFrame:
    pairs = [
        ("Human", "Human wound", "Human normal"),
        ("Rat", "Rat adult", "Rat embryo"),
        ("Reindeer", "Reindeer velvet", "Reindeer back"),
        ("Acomys", "Acomys wound", "Acomys unwounded"),
        ("Planarian", "Planarian 36 h", "Planarian 0 h"),
        ("Planarian", "Planarian 3 d", "Planarian 0 h"),
        ("Planarian", "Planarian 5 d", "Planarian 0 h"),
    ]
    metrics = ["KLF4_mean", "IGFL1_mean"] + [f"KLF4 context | {m}" for m in MODULES]
    rows = []
    for cohort, ga, gb in pairs:
        for metric in metrics:
            xa = d.loc[d.group.eq(ga), metric].dropna().to_numpy()
            xb = d.loc[d.group.eq(gb), metric].dropna().to_numpy()
            effect = np.nan
            if len(xa) and len(xb):
                pooled = np.sqrt(((len(xa)-1)*np.var(xa, ddof=1)+(len(xb)-1)*np.var(xb, ddof=1))/max(len(xa)+len(xb)-2, 1)) if len(xa) > 1 and len(xb) > 1 else np.nan
                effect = (xa.mean()-xb.mean())/pooled if np.isfinite(pooled) and pooled > 0 else np.log2((xa.mean()+1e-6)/(xb.mean()+1e-6))
            p = stats.mannwhitneyu(xa, xb, alternative="two-sided").pvalue if len(xa) > 1 and len(xb) > 1 else np.nan
            rows.append({"cohort": cohort, "contrast": f"{ga} vs {gb}", "metric": metric, "effect": effect, "p_value": p, "n_a": len(xa), "n_b": len(xb), "mean_a": xa.mean() if len(xa) else np.nan, "mean_b": xb.mean() if len(xb) else np.nan})
    e = pd.DataFrame(rows)
    e["p_adj_bh"] = np.nan
    m = e.p_value.notna()
    if m.any():
        e.loc[m, "p_adj_bh"] = multipletests(e.loc[m, "p_value"], method="fdr_bh")[1]
    e["significance"] = e.p_adj_bh.map(star)
    e.to_csv(TAB / "cross_species_effect_sizes.csv", index=False, encoding="utf-8-sig")
    return e


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIG / f"{stem}.png", dpi=600, facecolor="white")
    fig.savefig(FIG / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def plot_expression(cs: pd.DataFrame, av: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(180/25.4, 105/25.4), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], left=.17, right=.98, bottom=.12, top=.91, wspace=.42)
    vmax = {g: max(cs.loc[cs.regulator.eq(g), "mean_expression"].max(), 1e-6) for g in ["KLF4", "IGFL1"]}
    for j, (gene, title) in enumerate([("KLF4", "KLF4 in Fib_K14"), ("IGFL1", "IGFL1 in all cells/spots")]):
        ax = fig.add_subplot(gs[0, j])
        q = cs[cs.regulator.eq(gene)].set_index("group")
        for i, group in enumerate(ORDER):
            if group not in q.index:
                ax.scatter(0, i, marker="x", s=28, color="#777777", linewidth=.8)
                continue
            row = q.loc[group]
            ax.scatter(0, i, s=18 + 2.2*row.positive_pct, c=[row.mean_expression], cmap="viridis", vmin=0, vmax=vmax[gene], edgecolor="#333333", linewidth=.4)
        ax.set_xlim(-.45, .45); ax.set_xticks([0], [gene]); ax.set_ylim(len(ORDER)-.5, -.5)
        ax.set_yticks(range(len(ORDER)), ORDER if j == 0 else [])
        ax.set_title(title, loc="left", pad=7)
        ax.spines["bottom"].set_visible(False); ax.tick_params(axis="x", length=0)
    ax = fig.add_subplot(gs[0, 2])
    status_order = ["direct/mapped feature", "HVG-limited absence", "ortholog unresolved", "mapping unresolved", "feature-space absence"]
    cmap = mpl.colors.ListedColormap(["#3B8BC2", "#E4A853", "#C76D5E", "#8C79B5", "#D9D9D9"])
    mat = np.zeros((len(ORDER), 2))
    for i, group in enumerate(ORDER):
        for j, gene in enumerate(["KLF4", "IGFL1"]):
            q = av[(av.group == group) & (av.gene == gene)]
            mat[i, j] = status_order.index(q.audit_status.iloc[0]) if len(q) else 4
    ax.imshow(mat, cmap=cmap, vmin=-.5, vmax=4.5, aspect="auto")
    ax.set_xticks([0, 1], ["KLF4", "IGFL1"]); ax.set_yticks(range(len(ORDER)), [])
    ax.set_title("Feature/orthology audit", loc="left", pad=7)
    for i in range(len(ORDER)):
        for j in range(2):
            ax.text(j, i, ["yes", "HVG", "orth.", "map", "abs."][int(mat[i,j])], ha="center", va="center", fontsize=7, color="white" if mat[i,j] in [0,2,3] else "#333333")
    handles = [mpl.patches.Patch(color=cmap(i), label=l) for i, l in enumerate(status_order[:4])]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.67, .005), ncol=2, frameon=False, fontsize=7, columnspacing=.8, handlelength=.9)
    save(fig, "01_KLF4_IGFL1_expression_availability")


def plot_effects(e: pd.DataFrame) -> None:
    metrics = ["KLF4_mean", "IGFL1_mean"] + [f"KLF4 context | {m}" for m in MODULES]
    labels = ["KLF4 (Fib_K14)", "IGFL1 (all cells)"] + list(MODULES)
    contrasts = ["Human wound vs normal", "Rat adult vs embryo", "Reindeer velvet vs back", "Acomys wound vs unwounded", "Planarian 36 h vs 0 h", "Planarian 3 d vs 0 h", "Planarian 5 d vs 0 h"]
    mat = e.pivot(index="metric", columns="contrast", values="effect").reindex(index=metrics)
    colmap = dict(zip(e.contrast.unique(), contrasts))
    mat = mat.rename(columns=colmap).reindex(columns=contrasts)
    pmat = e.pivot(index="metric", columns="contrast", values="p_adj_bh").reindex(index=metrics).rename(columns=colmap).reindex(columns=contrasts)
    lim = np.nanmax(np.abs(mat.to_numpy())) if np.isfinite(mat.to_numpy()).any() else 1
    lim = max(lim, 1)
    fig, ax = plt.subplots(figsize=(180/25.4, 88/25.4), facecolor="white")
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(contrasts)), [x.replace(" vs ", "\nvs ") for x in contrasts], rotation=32, ha="right")
    ax.set_yticks(range(len(metrics)), labels)
    for i in range(len(metrics)):
        for j in range(len(contrasts)):
            v, p = mat.iloc[i,j], pmat.iloc[i,j]
            ax.text(j, i, "NA" if not np.isfinite(v) else f"{v:.2f}{star(p)}", ha="center", va="center", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=.024, pad=.025); cb.ax.set_title("Effect", fontsize=8, pad=4)
    ax.set_title("Sample-level effects; stars show BH-adjusted Mann-Whitney tests", loc="left")
    fig.subplots_adjust(left=.22, right=.94, bottom=.30, top=.89)
    save(fig, "02_cross_species_effect_size_heatmap")


def plot_correlations(co: pd.DataFrame) -> None:
    mods = list(MODULES)
    fig = plt.figure(figsize=(190/25.4, 108/25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1, .035], left=.17, right=.96, bottom=.22, top=.93, hspace=.36, wspace=.05)
    last_im = None
    for row, (reg, title) in enumerate([("KLF4", "KLF4-module association in Fib_K14"), ("IGFL1", "IGFL1-module association in all cells/spots")]):
        ax = fig.add_subplot(gs[row, 0])
        q = co[co.regulator.eq(reg)]
        mat = q.pivot(index="module", columns="group", values="rho").reindex(index=mods, columns=ORDER)
        pmat = q.pivot(index="module", columns="group", values="p_adj_bh").reindex(index=mods, columns=ORDER)
        last_im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_yticks(range(len(mods)), mods); ax.set_title(title, loc="left", pad=7)
        ax.set_xticks(range(len(ORDER)))
        if row == 1:
            ax.set_xticklabels([x.replace(" ", "\n", 1) for x in ORDER], rotation=36, ha="right")
        else:
            ax.set_xticklabels([])
        for i in range(len(mods)):
            for j in range(len(ORDER)):
                v, p = mat.iloc[i,j], pmat.iloc[i,j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.2f}{star(p)}", ha="center", va="center", fontsize=6.2)
    cax = fig.add_subplot(gs[:, 1]); cb = fig.colorbar(last_im, cax=cax); cb.set_label("Spearman rho")
    save(fig, "03_regulator_niche_module_correlations")


def plot_modules(d: pd.DataFrame) -> None:
    mods = list(MODULES)
    fig = plt.figure(figsize=(165/25.4, 105/25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1, .04], left=.18, right=.94, bottom=.17, top=.92, hspace=.32, wspace=.05)
    for row, (reg, title) in enumerate([("KLF4", "Fib_K14 context"), ("IGFL1", "All-cell/spot context")]):
        cols = [f"{reg} context | {m}" for m in mods]
        mat = d.groupby("group")[cols].mean().reindex(ORDER)
        mat.columns = mods
        z = (mat - mat.mean(axis=0)) / mat.std(axis=0).replace(0, np.nan)
        ax = fig.add_subplot(gs[row, 0]); im = ax.imshow(z, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
        ax.set_yticks(range(len(ORDER)), ORDER); ax.set_title(title, loc="left", pad=6)
        ax.set_xticks(range(len(mods)))
        ax.set_xticklabels([m.replace(" niche", "") for m in mods], rotation=28, ha="right" if row else "center")
        if row == 0: ax.set_xticklabels([])
    cax = fig.add_subplot(gs[:, 1]); cb = fig.colorbar(im, cax=cax); cb.ax.set_title("z-score", fontsize=8, pad=4)
    save(fig, "04_KLF4_target_and_niche_program_heatmap")


def update_rat_spatial() -> pd.DataFrame:
    rat_root = Path(r"data/external/local_fb7c1ee2a2f4")
    samples = [
        ("Embryo", "E16.5 6H", "E16.5_6H_C02938E1_Annotated_Spatial.h5ad"),
        ("P5", "P5 3D", "P5_3D_B02621C2_Annotated_Spatial.h5ad"),
        ("Adult", "Adult 10d", "Adult_10d_A02988C5_Annotated_Spatial.h5ad"),
    ]
    rows = []; region_values = {}
    fig = plt.figure(figsize=(180/25.4, 58/25.4), facecolor="white")
    gs = fig.add_gridspec(1, 4, width_ratios=[1,1,1,1.45], left=.025, right=.985, bottom=.20, top=.88, wspace=.34)
    palette = {"Center":"#7A9CC6", "Front":"#D95F5F", "Peri":"#B5B5B5"}
    for i, (stage, label, name) in enumerate(samples):
        a = ad.read_h5ad(rat_root/name, backed="r")
        key = "X_rot_spatial" if "X_rot_spatial" in a.obsm else "spatial"; c = np.asarray(a.obsm[key])[:,:2]
        fib = a.obs["sub_labels"].astype(str).to_numpy() == "Fib_K14"
        vectors, matched = read_vectors(a, ["KLF4"]); v = vectors.get("KLF4", np.zeros(a.n_obs))
        x, y = c[:,0], c[:,1]; bins=150
        ix=np.clip(((x-x.min())/(np.ptp(x)+1e-9)*(bins-1)).astype(int),0,bins-1); iy=np.clip(((y-y.min())/(np.ptp(y)+1e-9)*(bins-1)).astype(int),0,bins-1)
        h=np.zeros((bins,bins)); np.add.at(h,(iy,ix),np.where(fib,v,0)); h=gaussian_filter(h,1.3); score=h[iy,ix]
        hi=np.quantile(score[score>0],.99) if np.any(score>0) else 1; score=np.clip(score/hi,0,1)
        ax=fig.add_subplot(gs[0,i]); idx=np.argsort(score); idx=idx[np.linspace(0,len(idx)-1,min(35000,len(idx)),dtype=int)]
        ax.scatter(c[idx,0],c[idx,1],c=score[idx],s=.12,cmap="magma",vmin=0,vmax=1,linewidths=0); ax.set_title(label); ax.axis("off"); ax.set_aspect("equal")
        left,right=np.quantile(x,[.25,.75]); band=(right-left)*.10; front=(np.abs(x-left)<=band)|(np.abs(x-right)<=band); center=(x>left+band)&(x<right-band); region=np.where(front,"Front",np.where(center,"Center","Peri"))
        for r in ["Center","Front","Peri"]:
            k=fib&(region==r); region_values[(stage,r)] = v[k]
            rows.append({"stage":stage,"sample":label,"region":r,"n_Fib_K14":int(k.sum()),"KLF4_mean":float(v[k].mean()) if k.sum() else np.nan,"KLF4_positive_pct":float((v[k]>0).mean()*100) if k.sum() else np.nan})
        try:a.file.close()
        except Exception:pass
    d=pd.DataFrame(rows); stat_rows=[]
    for stage in ["Embryo","P5","Adult"]:
        base=region_values[(stage,"Center")]
        for r in ["Front","Peri"]:
            cur=region_values[(stage,r)]
            p=stats.mannwhitneyu(cur,base,alternative="two-sided").pvalue if len(cur)>1 and len(base)>1 else np.nan
            stat_rows.append({"stage":stage,"contrast":f"{r} vs Center","p_value":p,"n_region":len(cur),"n_center":len(base),"unit":"Fib_K14 spots within one selected slice"})
    stat_df=pd.DataFrame(stat_rows); stat_df["p_adj_bh"]=np.nan; sm=stat_df.p_value.notna()
    if sm.any():stat_df.loc[sm,"p_adj_bh"]=multipletests(stat_df.loc[sm,"p_value"],method="fdr_bh")[1]
    stat_df["significance"]=stat_df.p_adj_bh.map(star)
    ax=fig.add_subplot(gs[0,3]); width=.22
    for i,stage in enumerate(["Embryo","P5","Adult"]):
        for j,r in enumerate(["Center","Front","Peri"]):
            val=d[(d.stage==stage)&(d.region==r)].KLF4_mean.iloc[0]; ax.bar(i+(j-1)*width,val,width,color=palette[r],edgecolor="#444",linewidth=.4,label=r if i==0 else None)
            if r != "Center":
                sig=stat_df[(stat_df.stage==stage)&(stat_df.contrast==f"{r} vs Center")].significance.iloc[0]
                if sig:ax.text(i+(j-1)*width,val+.009,sig,ha="center",va="bottom",fontsize=7)
    ax.set_xticks(range(3),["Emb.","P5","Adult"]); ax.set_ylabel("KLF4 mean in Fib_K14"); ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,1.18),columnspacing=.6,handlelength=.8)
    d.to_csv(TAB/"rat_spatial_KLF4_region_summary.csv",index=False,encoding="utf-8-sig")
    ax.set_ylim(0,max(d.KLF4_mean.max()*1.25,.1))
    stat_df.to_csv(TAB/"rat_spatial_KLF4_region_statistics.csv",index=False,encoding="utf-8-sig")
    save(fig,"05_rat_spatial_KLF4_region_validation")
    return d


def update_reindeer_time() -> None:
    spec = next(s for s in specs11() if s.label == "Reindeer velvet")
    a = ad.read_h5ad(spec.path, backed="r")
    fib = a.obs["Fib_K14"].astype(str).to_numpy() == "Fib_K14"
    tissue = a.obs["tissue"].astype(str).to_numpy(); time = a.obs["timepoint"].astype(str).to_numpy(); sample = a.obs["sample"].astype(str).to_numpy()
    vectors, _ = read_vectors(a,["KLF4","IGFL1"]); rows=[]
    for ti,tlabel in [("antler_velvet","Velvet"),("back_skin","Back")]:
        for tp in ["0","3","7","14"]:
            for gene in ["KLF4","IGFL1"]:
                context = fib if gene=="KLF4" else np.ones(a.n_obs,dtype=bool)
                for sid in np.unique(sample[(tissue==ti)&(time==tp)]):
                    k=context&(tissue==ti)&(time==tp)&(sample==sid); v=vectors[gene]
                    rows.append({"tissue":tlabel,"time_day":int(tp),"sample":sid,"gene":gene,"context":"Fib_K14" if gene=="KLF4" else "all cells","n":int(k.sum()),"mean_expression":float(v[k].mean()) if k.sum() else np.nan,"positive_pct":float((v[k]>0).mean()*100) if k.sum() else np.nan})
    d=pd.DataFrame(rows); stat_rows=[]
    for gene in ["KLF4","IGFL1"]:
        for tissue_name in ["Velvet","Back"]:
            base=d[(d.gene==gene)&(d.tissue==tissue_name)&(d.time_day==0)].mean_expression.dropna().to_numpy()
            for tp in [3,7,14]:
                cur=d[(d.gene==gene)&(d.tissue==tissue_name)&(d.time_day==tp)].mean_expression.dropna().to_numpy()
                p=stats.mannwhitneyu(cur,base,alternative="two-sided").pvalue if len(cur)>1 and len(base)>1 else np.nan
                stat_rows.append({"gene":gene,"tissue":tissue_name,"contrast":f"day {tp} vs day 0","p_value":p,"n_day":len(cur),"n_day0":len(base)})
    st=pd.DataFrame(stat_rows); st["p_adj_bh"]=np.nan; m=st.p_value.notna()
    if m.any():st.loc[m,"p_adj_bh"]=multipletests(st.loc[m,"p_value"],method="fdr_bh")[1]
    st["significance"]=st.p_adj_bh.map(star)
    fig,axs=plt.subplots(1,2,figsize=(130/25.4,58/25.4),facecolor="white")
    for ax,gene in zip(axs,["KLF4","IGFL1"]):
        q=d[d.gene==gene]
        for tissue_name,col,mark in [("Velvet","#2A9D8F","o"),("Back","#6C757D","s")]:
            z=q[q.tissue==tissue_name].groupby("time_day").agg(mean=("mean_expression","mean"),sem=("mean_expression","sem"),pct=("positive_pct","mean")).reset_index()
            ax.errorbar(z.time_day,z["mean"],yerr=z["sem"],marker=mark,color=col,label=tissue_name,lw=1.1,capsize=2)
            for _,rr in z.iterrows():
                sig=st[(st.gene==gene)&(st.tissue==tissue_name)&(st.contrast==f"day {int(rr.time_day)} vs day 0")].significance
                if int(rr.time_day)>0 and len(sig) and sig.iloc[0]: ax.text(rr.time_day,rr["mean"]+max(z["mean"].max()*.08,.002),sig.iloc[0],ha="center",fontsize=7)
        ax.set_title(f"{gene} in {'Fib_K14' if gene=='KLF4' else 'all cells'}"); ax.set_xlabel("Day"); ax.set_ylabel("Sample mean expression"); ax.set_xticks([0,3,7,14]); ax.legend(frameon=False)
    fig.tight_layout(); save(fig,"06_reindeer_KLF4_IGFL1_timecourse")
    d.to_csv(TAB/"reindeer_KLF4_IGFL1_timecourse.csv",index=False,encoding="utf-8-sig"); st.to_csv(TAB/"reindeer_KLF4_IGFL1_timecourse_statistics.csv",index=False,encoding="utf-8-sig")
    try:a.file.close()
    except Exception:pass


def write_ensembl_igfl1_ortholog_audit() -> None:
    url = "https://rest.ensembl.org/homology/symbol/homo_sapiens/IGFL1?content-type=application/json;type=orthologues"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            payload = json.load(response)
        rows = []
        for item in payload.get("data", [{}])[0].get("homologies", []):
            target = item.get("target", {})
            rows.append(
                {
                    "query_symbol": "IGFL1",
                    "target_species": target.get("species", ""),
                    "target_taxon_id": target.get("taxon_id", ""),
                    "target_gene_id": target.get("id", ""),
                    "homology_type": item.get("type", ""),
                    "percent_identity": target.get("perc_id", np.nan),
                    "taxonomy_level": item.get("taxonomy_level", ""),
                    "query_url": url,
                    "retrieved_date": "2026-07-13",
                }
            )
        pd.DataFrame(rows).sort_values("target_species").to_csv(TAB/"IGFL1_Ensembl_ortholog_query.csv", index=False, encoding="utf-8-sig")
        (TAB/"IGFL1_Ensembl_ortholog_query_raw.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    except Exception as exc:
        pd.DataFrame([{"query_symbol":"IGFL1", "query_url":url, "error":str(exc), "retrieved_date":"2026-07-13"}]).to_csv(TAB/"IGFL1_Ensembl_ortholog_query.csv", index=False, encoding="utf-8-sig")


def main() -> None:
    ensure_dirs()
    d,av,co,cs=summarize(); e=effects(d)
    plot_expression(cs,av); plot_effects(e); plot_correlations(co); plot_modules(d)
    update_rat_spatial(); update_reindeer_time()
    write_ensembl_igfl1_ortholog_audit()
    shutil.copy2(Path(__file__), SCR/Path(__file__).name)
    report={
        "group_order":ORDER,
        "KLF4_context":"Fib_K14 only",
        "IGFL1_context":"all cells/spots in each selected group",
        "correlation_statistics":"Spearman correlation with BH correction within each regulator; cell/spot-level descriptive inference",
        "effect_statistics":"sample-level Mann-Whitney where both groups have at least two sample summaries; global BH correction",
        "IGFL1_audit":"Human and reindeer directly measurable. Rat Igfl1 is absent and no Ensembl rat ortholog was returned. Acomys objects are HVG-limited, so absence is unresolved. Planarian mapping did not identify IGFL1.",
        "editable_pdf":"Arial with PDF Type 42 text",
    }
    (REP/"report_context_stats_update.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"output_root":str(ROOT),"figures":6,"groups":ORDER},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
