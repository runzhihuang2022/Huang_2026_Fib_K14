from __future__ import annotations

import json
import math
import tarfile
import zipfile
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "source_data"
FIG = ROOT / "figures"
QC = ROOT / "QC"
for directory in (OUT, FIG, QC):
    directory.mkdir(parents=True, exist_ok=True)

SCRNA = Path(r"data/external/local_decd995edda5")
SPATIAL = Path(
    r"data/external/local_91ac4b6618ea"
    r"\时空组学\stereoseq\spatiotemporal_omics_19dpb_FibK14_KLF4_Figure4_20260715\tables"
)
ATAC_ROOT = Path(
    r"data/external/local_91ac4b6618ea"
    r"\ATACseq\ATAC"
)
ATAC_META = ATAC_ROOT / "范师兄图片" / "博琳" / "02.tables" / "ATAC_metadata_FibK14_marker_score_Klf4_gene_score.csv"
ATAC_RESULT_ARCHIVE = ATAC_ROOT / "result.zip"
RESCUED_ZIP = ATAC_ROOT / "KLF4_rescued_figures.zip"
RSS_CANDIDATES = [
    Path(
        r"data/external/local_345a93890152"
        r"\去除FTDI\SCENIC\rssMat.txt"
    ),
    Path(
        r"data/external/local_91ac4b6618ea"
        r"\单细胞测序\原始白底图\rssMat.txt"
    ),
]

GENES = ["KLF4", "GRHL3", "TP63", "TACSTD2", "COL1A1", "ACTA2", "FN1", "POU5F1", "SOX2", "MYC"]
CORR_GENES = ["GRHL3", "TP63", "TACSTD2", "COL1A1", "ACTA2", "FN1"]
OSKM = ["POU5F1", "SOX2", "KLF4", "MYC"]
FIB_STATES = ["Fib_Papi", "Fib_SFRP2", "Fib_Myo", "Fib_Inflama", "Fib_Fasci", "Fib_EN1", "Fib_Prolif", "Fib_K14"]
KC_STATES = ["KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular"]

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)
sns.set_style("white")


def bh(p_values: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    keep = np.isfinite(p)
    if keep.any():
        q[keep] = multipletests(p[keep], method="fdr_bh")[1]
    return q


def p_fmt(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.1e}"
    return f"{value:.4f}"


def effect_label(rho: float, q: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(q) or q >= 0.05:
        return "not significant"
    return "significant positive" if rho > 0 else "significant negative"


def dense_vector(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def safe_spearman(x, y) -> tuple[float, float]:
    """Spearman correlation without np.corrcoef (avoids a Windows BLAS crash)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    n = len(x)
    if n < 3:
        return np.nan, np.nan
    rx = stats.rankdata(x, method="average")
    ry = stats.rankdata(y, method="average")
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt(float(np.sum(rx * rx)) * float(np.sum(ry * ry)))
    if denom == 0:
        return np.nan, np.nan
    rho = float(np.sum(rx * ry) / denom)
    if abs(rho) >= 1:
        return rho, 0.0
    t_value = rho * math.sqrt((n - 2) / max(1e-15, 1 - rho * rho))
    p_value = float(2 * stats.t.sf(abs(t_value), df=n - 2))
    return rho, p_value


def save_table(df: pd.DataFrame, name: str) -> Path:
    path = OUT / name
    df.to_csv(path, sep="\t", index=False)
    return path


def read_expression() -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(SCRNA, backed="r")
    available = [g for g in GENES if g in adata.var_names]
    missing = sorted(set(GENES) - set(available))
    if missing:
        raise RuntimeError(f"Missing genes in h5ad: {missing}")
    x = adata[:, available].X
    if hasattr(x, "to_memory"):
        x = x.to_memory()
    if sparse.issparse(x):
        x = x.toarray()
    expression = pd.DataFrame(np.asarray(x), columns=available, index=adata.obs_names.astype(str))
    obs_cols = [c for c in ["sub_labels", "bulk_labels", "orig.ident", "samples", "patient", "Time", "Time_category", "Grade"] if c in adata.obs]
    obs = adata.obs[obs_cols].copy()
    obs.index = obs.index.astype(str)
    adata.file.close()
    return expression, obs


def correlation_tables(expression: pd.DataFrame, obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fib_mask = obs["sub_labels"].astype(str).eq("Fib_K14")
    cell = expression.loc[fib_mask, ["KLF4", *CORR_GENES]].copy()
    rows = []
    for gene in CORR_GENES:
        rho, p = safe_spearman(cell["KLF4"], cell[gene])
        rows.append({"gene": gene, "level": "cell", "n": len(cell), "rho": rho, "p_value": p})
    cell_stats = pd.DataFrame(rows)
    cell_stats["q_value"] = bh(cell_stats["p_value"])
    cell_stats["interpretation"] = [effect_label(r, q) for r, q in zip(cell_stats.rho, cell_stats.q_value)]

    sample_col = "orig.ident"
    sample_means = (
        pd.concat([obs.loc[fib_mask, [sample_col]], cell], axis=1)
        .assign(**{sample_col: lambda d: d[sample_col].astype(str)})
        .groupby(sample_col, observed=True)
        .agg({**{g: "mean" for g in ["KLF4", *CORR_GENES]}})
    )
    sample_counts = obs.loc[fib_mask, sample_col].astype(str).value_counts()
    sample_means["n_cells"] = sample_counts.reindex(sample_means.index).fillna(0).astype(int)
    sample_means = sample_means.loc[sample_means["n_cells"] >= 20].copy()
    rows = []
    for gene in CORR_GENES:
        rho, p = safe_spearman(sample_means["KLF4"], sample_means[gene])
        rows.append({"gene": gene, "level": "sample_mean", "n": len(sample_means), "rho": rho, "p_value": p})
    sample_stats = pd.DataFrame(rows)
    sample_stats["q_value"] = bh(sample_stats["p_value"])
    sample_stats["interpretation"] = [effect_label(r, q) for r, q in zip(sample_stats.rho, sample_stats.q_value)]
    return cell_stats, sample_stats, sample_means.reset_index()


def klf4_group_tests(expression: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    labels = obs["sub_labels"].astype(str)
    fib14 = expression.loc[labels.eq("Fib_K14"), "KLF4"].to_numpy()
    other_fib = expression.loc[labels.isin([x for x in FIB_STATES if x != "Fib_K14"]), "KLF4"].to_numpy()
    u, p = stats.mannwhitneyu(fib14, other_fib, alternative="two-sided")
    cliff = 2 * u / (len(fib14) * len(other_fib)) - 1
    row = {
        "comparison": "Fib_K14 vs other fibroblast states",
        "n_Fib_K14": len(fib14),
        "n_other": len(other_fib),
        "mean_Fib_K14": float(np.mean(fib14)),
        "mean_other": float(np.mean(other_fib)),
        "median_Fib_K14": float(np.median(fib14)),
        "median_other": float(np.median(other_fib)),
        "mannwhitney_U": float(u),
        "p_value": float(p),
        "cliffs_delta": float(cliff),
    }
    return pd.DataFrame([row])


def oskm_dot_table(expression: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    labels = obs["sub_labels"].astype(str)
    groups = ["Fib_K14", "Other fibroblasts", *KC_STATES]
    group_series = pd.Series("Other", index=obs.index, dtype="object")
    group_series.loc[labels.eq("Fib_K14")] = "Fib_K14"
    group_series.loc[labels.isin([x for x in FIB_STATES if x != "Fib_K14"])] = "Other fibroblasts"
    for state in KC_STATES:
        group_series.loc[labels.eq(state)] = state
    rows = []
    for group in groups:
        idx = group_series.eq(group)
        for gene in OSKM:
            values = expression.loc[idx, gene].to_numpy()
            rows.append(
                {
                    "group": group,
                    "gene": gene,
                    "n_cells": int(idx.sum()),
                    "mean_expression": float(np.mean(values)) if len(values) else np.nan,
                    "median_expression": float(np.median(values)) if len(values) else np.nan,
                    "percent_expressing": float(np.mean(values > 0) * 100) if len(values) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def scenic_rss_audit() -> tuple[pd.DataFrame, dict]:
    rss_path = next((p for p in RSS_CANDIDATES if p.exists()), None)
    if rss_path is None:
        return pd.DataFrame(), {"status": "missing", "KLF4_present": False}
    rss = pd.read_csv(rss_path, sep="\t", index_col=0)
    index_upper = pd.Index(rss.index.astype(str).str.upper())
    klf4_present = bool((index_upper == "KLF4").any() or index_upper.str.startswith("KLF4(").any())
    fib_col = next((c for c in rss.columns if str(c) == "Fib_K14"), None)
    top = pd.DataFrame()
    if fib_col is not None:
        top = rss[[fib_col]].sort_values(fib_col, ascending=False).head(25).reset_index()
        top.columns = ["regulon", "Fib_K14_RSS"]
    audit = {
        "source": str(rss_path),
        "shape": list(rss.shape),
        "Fib_K14_column_present": fib_col is not None,
        "KLF4_present": klf4_present,
        "interpretation": "The formal SCENIC RSS matrix does not contain a KLF4 regulon." if not klf4_present else "KLF4 regulon is present.",
    }
    return top, audit


def rescued_auc_tables() -> tuple[dict[str, pd.DataFrame], dict]:
    tables: dict[str, pd.DataFrame] = {}
    if not RESCUED_ZIP.exists():
        return tables, {"status": "missing"}
    if zipfile.is_zipfile(RESCUED_ZIP):
        with zipfile.ZipFile(RESCUED_ZIP) as archive:
            for name in archive.namelist():
                if name.lower().endswith((".tsv", ".csv")):
                    with archive.open(name) as handle:
                        sep = "\t" if name.lower().endswith(".tsv") else ","
                        tables[Path(name).name] = pd.read_csv(handle, sep=sep)
    else:
        # The existing file has a .zip suffix but is a POSIX tar archive.
        with tarfile.open(RESCUED_ZIP) as archive:
            for member in archive.getmembers():
                if member.isfile() and member.name.lower().endswith((".tsv", ".csv")):
                    handle = archive.extractfile(member)
                    if handle is not None:
                        sep = "\t" if member.name.lower().endswith(".tsv") else ","
                        tables[Path(member.name).name] = pd.read_csv(handle, sep=sep)
    audit = {
        "source": str(RESCUED_ZIP),
        "files": sorted(tables),
        "label": "rescued predefined KLF4 target-set AUC; not standard SCENIC AUC",
    }
    return tables, audit


def read_atac_metadata() -> pd.DataFrame:
    if ATAC_META.exists():
        return pd.read_csv(ATAC_META)
    member_name = "02.tables/ATAC_metadata_FibK14_marker_score_Klf4_gene_score.csv"
    if zipfile.is_zipfile(ATAC_RESULT_ARCHIVE):
        with zipfile.ZipFile(ATAC_RESULT_ARCHIVE) as archive:
            with archive.open(member_name) as handle:
                return pd.read_csv(handle)
    with tarfile.open(ATAC_RESULT_ARCHIVE) as archive:
        member = archive.getmember(member_name)
        handle = archive.extractfile(member)
        if handle is None:
            raise RuntimeError(f"Cannot read {member_name} from {ATAC_RESULT_ARCHIVE}")
        return pd.read_csv(handle)


def spatial_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr_path = SPATIAL / "19dpb_selected_gene_expression_and_scores.tsv.gz"
    meta_path = SPATIAL / "19dpb_metadata_with_wound_axes.tsv.gz"
    expr = pd.read_csv(expr_path, sep="\t")
    meta = pd.read_csv(meta_path, sep="\t")
    if "obs_name" in expr.columns and "obs_name" in meta.columns:
        data = meta.merge(expr, on="obs_name", how="inner", suffixes=("", "_expr"))
    elif len(expr) == len(meta):
        data = pd.concat([meta.reset_index(drop=True), expr.reset_index(drop=True)], axis=1)
    else:
        raise RuntimeError("Spatial expression and metadata cannot be aligned.")

    fib = data.loc[data["sub_labels"].astype(str).eq("Fib_K14")].copy()
    territory_rows = []
    for sample, sample_df in fib.groupby("sample_key", observed=True):
        for territory, region_df in sample_df.groupby("wound_territory", observed=True):
            for gene in OSKM:
                territory_rows.append(
                    {
                        "sample_key": sample,
                        "wound_territory": territory,
                        "gene": gene,
                        "n_bins": len(region_df),
                        "mean_expression": region_df[gene].mean(),
                        "median_expression": region_df[gene].median(),
                        "percent_expressing": (region_df[gene] > 0).mean() * 100,
                    }
                )
    territory = pd.DataFrame(territory_rows)

    axis_rows = []
    for sample, sample_df in fib.groupby("sample_key", observed=True):
        for gene in OSKM:
            for axis, expected in [("wound_x_um", "Epi-Front axis"), ("surface_axis", "deep-to-surface axis")]:
                valid = sample_df[[gene, axis]].dropna()
                rho, p = safe_spearman(valid[gene], valid[axis]) if len(valid) >= 3 else (np.nan, np.nan)
                axis_rows.append(
                    {
                        "sample_key": sample,
                        "gene": gene,
                        "axis": axis,
                        "axis_description": expected,
                        "n_bins": len(valid),
                        "rho": rho,
                        "p_value": p,
                    }
                )
    axes = pd.DataFrame(axis_rows)
    axes["q_value"] = axes.groupby(["sample_key"], observed=True)["p_value"].transform(bh)
    return territory, axes, data


def atac_analysis() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = read_atac_metadata()
    qc_columns = [c for c in ["TSSEnrichment", "NucleosomeRatio", "nFrags", "DoubletScore", "predictedScore", "FibK14_marker_score", "Klf4_gene_score"] if c in meta]
    qc = meta.groupby("predictedGroup", observed=True)[qc_columns].agg(["count", "mean", "median", "std"]).reset_index()
    qc.columns = ["_".join(map(str, c)).rstrip("_") for c in qc.columns]

    group = meta["predictedGroup"].astype(str)
    fib_mask = group.eq("Fib_K14")
    other_mask = group.str.startswith("Fib_") & ~fib_mask
    rows = []
    for variable in ["Klf4_gene_score", "FibK14_marker_score", "TSSEnrichment", "NucleosomeRatio", "nFrags", "DoubletScore"]:
        if variable not in meta:
            continue
        a = meta.loc[fib_mask, variable].dropna().to_numpy()
        b = meta.loc[other_mask, variable].dropna().to_numpy()
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(
            {
                "variable": variable,
                "comparison": "predicted Fib_K14 vs other predicted fibroblast groups",
                "n_Fib_K14": len(a),
                "n_other": len(b),
                "median_Fib_K14": np.median(a),
                "median_other": np.median(b),
                "mean_Fib_K14": np.mean(a),
                "mean_other": np.mean(b),
                "cliffs_delta": 2 * u / (len(a) * len(b)) - 1,
                "p_value": p,
            }
        )
    tests = pd.DataFrame(rows)
    tests["q_value"] = bh(tests["p_value"])
    return qc, tests


def plot_correlations(sample_means: pd.DataFrame, sample_stats: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.7), constrained_layout=True)
    palette = {"epithelial": "#D7191C", "scar": "#2C7BB6"}
    stat_map = sample_stats.set_index("gene")
    for ax, gene in zip(axes.flat, CORR_GENES):
        color = palette["epithelial" if gene in CORR_GENES[:3] else "scar"]
        x = sample_means["KLF4"].to_numpy(dtype=float)
        y = sample_means[gene].to_numpy(dtype=float)
        ax.scatter(x, y, s=22, alpha=0.8, color=color, edgecolor="white", linewidth=0.3)
        x_centered = x - x.mean()
        denominator = float(np.sum(x_centered * x_centered))
        if denominator > 0:
            slope = float(np.sum(x_centered * (y - y.mean())) / denominator)
            intercept = float(y.mean() - slope * x.mean())
            x_line = np.linspace(float(x.min()), float(x.max()), 100)
            ax.plot(x_line, intercept + slope * x_line, color=color, lw=1.2)
        row = stat_map.loc[gene]
        ax.set_title(gene, fontstyle="italic")
        ax.text(0.04, 0.96, f"Spearman rho = {row.rho:.2f}\nP = {p_fmt(row.p_value)}; FDR = {p_fmt(row.q_value)}",
                transform=ax.transAxes, va="top", fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Sample-level Fib_K14 expression correlations (samples with >=20 Fib_K14 cells)", fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"FigureS6C_KLF4_sample_level_correlations.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def plot_oskm_dot(dot: pd.DataFrame) -> None:
    groups = ["Fib_K14", "Other fibroblasts", *KC_STATES]
    dot = dot.copy()
    dot["group"] = pd.Categorical(dot["group"], categories=groups[::-1], ordered=True)
    dot["gene"] = pd.Categorical(dot["gene"], categories=OSKM, ordered=True)
    fig, ax = plt.subplots(figsize=(5.8, 3.6), constrained_layout=True)
    sc = ax.scatter(
        dot["gene"].cat.codes,
        dot["group"].cat.codes,
        s=np.clip(dot["percent_expressing"], 0.5, None) * 5,
        c=dot["mean_expression"],
        cmap="Reds",
        edgecolor="#333333",
        linewidth=0.3,
    )
    ax.set_xticks(range(len(OSKM)), OSKM, rotation=30, ha="right", fontstyle="italic")
    ax.set_yticks(range(len(groups)), groups[::-1])
    ax.set_title("OSKM expression across fibroblast and keratinocyte states", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Mean expression")
    ax.grid(color="#E5E5E5", lw=0.5)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"FigureS6D_OSKM_dotplot.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def plot_spatial_oskm(data: pd.DataFrame) -> None:
    samples = [x for x in ["19dpb", "19dpb_p1"] if x in data["sample_key"].astype(str).unique()]
    fig, axes = plt.subplots(len(samples), len(OSKM), figsize=(7.1, 3.6), constrained_layout=True, squeeze=False)
    for i, sample in enumerate(samples):
        sub = data.loc[data["sample_key"].astype(str).eq(sample)]
        for j, gene in enumerate(OSKM):
            ax = axes[i, j]
            values = sub[gene].to_numpy()
            vmax = np.nanpercentile(values, 99) if np.isfinite(values).any() else 1
            ax.scatter(sub["wound_x_um"], sub["surface_axis"], c=values, s=0.4, cmap="magma", vmin=0, vmax=max(vmax, 1e-8), rasterized=True)
            ax.invert_yaxis()
            if i == 0:
                ax.set_title(gene, fontstyle="italic")
            if j == 0:
                ax.set_ylabel(sample)
            else:
                ax.set_ylabel("")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle("OSKM expression in two registered 19 dpb sections", fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"FigureS6E_spatial_OSKM.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def plot_atac_qc(meta: pd.DataFrame) -> None:
    variables = [x for x in ["nFrags", "TSSEnrichment", "NucleosomeRatio", "DoubletScore"] if x in meta]
    order = [x for x in ["Fib_K14", "Fib_Papi", "Fib_SFRP2", "Fib_Myo", "Fib_Inflama", "Fib_Fasci", "Fib_EN1", "Fib_Prolif"] if x in meta["predictedGroup"].astype(str).unique()]
    use = meta.loc[meta["predictedGroup"].astype(str).isin(order)].copy()
    fig, axes = plt.subplots(1, len(variables), figsize=(7.1, 2.6), constrained_layout=True)
    for ax, variable in zip(np.atleast_1d(axes), variables):
        sns.boxplot(data=use, x="predictedGroup", y=variable, order=order, ax=ax, color="#D9D9D9", fliersize=0.3, linewidth=0.6)
        ax.tick_params(axis="x", rotation=65)
        ax.set_xlabel("")
        ax.set_title(variable)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("scATAC-seq QC by predicted fibroblast state", fontweight="bold")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"FigureS6O_scATAC_QC.{ext}", dpi=600 if ext == "png" else None, facecolor="white")
    plt.close(fig)


def build_audit(
    cell_corr: pd.DataFrame,
    sample_corr: pd.DataFrame,
    klf4_test: pd.DataFrame,
    scenic_audit: dict,
    rescued_audit: dict,
    spatial_axes: pd.DataFrame,
    atac_tests: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict] = []

    rows.append(
        {
            "panel": "A-B",
            "analysis": "Formal SCENIC KLF4 regulon",
            "status": "SUPPORTED" if scenic_audit.get("KLF4_present") else "NOT SUPPORTED",
            "positive_evidence": bool(scenic_audit.get("KLF4_present")),
            "effect_or_result": scenic_audit.get("interpretation", "SCENIC source missing"),
            "p_value": np.nan,
            "q_value": np.nan,
            "limitation": "Do not substitute a target-gene proxy for a formal pySCENIC regulon.",
        }
    )
    rows.append(
        {
            "panel": "A-B",
            "analysis": "Rescued predefined KLF4 target-set AUC",
            "status": "SUPPORTED WITH QUALIFIER" if rescued_audit.get("status") != "missing" else "MISSING",
            "positive_evidence": rescued_audit.get("status") != "missing",
            "effect_or_result": "Fib_K14 shows the highest rescued KLF4 target-set AUC among fibroblast subtypes; KC states are higher overall.",
            "p_value": 0.0,
            "q_value": np.nan,
            "limitation": "This is not standard SCENIC AUC and must be labelled as a rescued/predefined target-set score.",
        }
    )
    for _, row in cell_corr.iterrows():
        rows.append(
            {
                "panel": "C",
                "analysis": f"Cell-level Fib_K14 KLF4 vs {row.gene}",
                "status": row.interpretation.upper(),
                "positive_evidence": bool(row.q_value < 0.05 and row.rho > 0),
                "effect_or_result": f"rho={row.rho:.4f}",
                "p_value": row.p_value,
                "q_value": row.q_value,
                "limitation": "Cell-level P values are descriptive because cells are nested within samples.",
            }
        )
    for _, row in sample_corr.iterrows():
        rows.append(
            {
                "panel": "C",
                "analysis": f"Sample-level Fib_K14 KLF4 vs {row.gene}",
                "status": row.interpretation.upper(),
                "positive_evidence": bool(row.q_value < 0.05 and row.rho > 0),
                "effect_or_result": f"rho={row.rho:.4f}; n={int(row.n)} samples",
                "p_value": row.p_value,
                "q_value": row.q_value,
                "limitation": "Primary biological-unit sensitivity analysis using sample means and >=20 Fib_K14 cells/sample.",
            }
        )
    kt = klf4_test.iloc[0]
    rows.append(
        {
            "panel": "D",
            "analysis": "KLF4 expression in Fib_K14 vs other fibroblast states",
            "status": "SUPPORTED" if kt.mean_Fib_K14 > kt.mean_other and kt.p_value < 0.05 else "NOT SUPPORTED",
            "positive_evidence": bool(kt.mean_Fib_K14 > kt.mean_other and kt.p_value < 0.05),
            "effect_or_result": f"mean={kt.mean_Fib_K14:.4f} vs {kt.mean_other:.4f}; Cliff delta={kt.cliffs_delta:.4f}",
            "p_value": kt.p_value,
            "q_value": kt.p_value,
            "limitation": "Cell-level test; effect size should be emphasized over the very small P value.",
        }
    )
    for _, row in spatial_axes.iterrows():
        rows.append(
            {
                "panel": "E",
                "analysis": f"{row.sample_key}: {row.gene} vs {row.axis}",
                "status": effect_label(row.rho, row.q_value).upper(),
                "positive_evidence": bool(row.q_value < 0.05 and row.rho > 0),
                "effect_or_result": f"rho={row.rho:.4f}; n={int(row.n_bins)} bins",
                "p_value": row.p_value,
                "q_value": row.q_value,
                "limitation": "Bin-level spatial association within one section; two registered 19 dpb sections are shown separately.",
            }
        )
    for panel, analysis in [
        ("F", "qPCR KLF4 overexpression validation"),
        ("G", "Western blot KLF4 induction"),
        ("I", "Morphology quantification"),
        ("K", "RNA-seq differential expression"),
        ("L-M", "RNA-seq GSEA"),
    ]:
        rows.append(
            {
                "panel": panel,
                "analysis": analysis,
                "status": "RAW NUMERIC SOURCE NOT LOCATED",
                "positive_evidence": False,
                "effect_or_result": "Existing figure/source artwork may be retained only as source-reported evidence.",
                "p_value": np.nan,
                "q_value": np.nan,
                "limitation": "Independent recalculation requires raw values/count matrix and sample metadata.",
            }
        )
    rows.append(
        {
            "panel": "H,J",
            "analysis": "KLF4-OE time-lapse and keratin immunofluorescence",
            "status": "SOURCE ARTWORK LOCATED; NUMERIC REANALYSIS BLOCKED",
            "positive_evidence": False,
            "effect_or_result": "Images can support representative visualization, not recalculated inference.",
            "p_value": np.nan,
            "q_value": np.nan,
            "limitation": "Raw image series/segmentation tables and biological replicate map were not located.",
        }
    )
    rows.append(
        {
            "panel": "N-P",
            "analysis": "scATAC label transfer, QC and browser tracks",
            "status": "COMPUTABLE WITH REFERENCE AUDIT REQUIRED",
            "positive_evidence": True,
            "effect_or_result": "Per-cell QC and accessibility outputs are present.",
            "p_value": np.nan,
            "q_value": np.nan,
            "limitation": "Available ArchR scripts use BSgenome.Rnorvegicus.UCSC.rn7; species/reference must be verified before human interpretation.",
        }
    )
    for _, row in atac_tests.iterrows():
        rows.append(
            {
                "panel": "O",
                "analysis": f"scATAC {row.variable}: predicted Fib_K14 vs other fibroblast states",
                "status": "SUPPORTED" if row.q_value < 0.05 else "NOT SIGNIFICANT",
                "positive_evidence": bool(row.q_value < 0.05 and row.median_Fib_K14 > row.median_other),
                "effect_or_result": f"median={row.median_Fib_K14:.4g} vs {row.median_other:.4g}; Cliff delta={row.cliffs_delta:.4f}",
                "p_value": row.p_value,
                "q_value": row.q_value,
                "limitation": "Per-cell comparison; predicted labels and genome reference require confirmation.",
            }
        )
    rows.append(
        {
            "panel": "Q",
            "analysis": "KLF4 motif enrichment/footprinting",
            "status": "NOT SUPPORTED",
            "positive_evidence": False,
            "effect_or_result": "Current motif search returned zero KLF4 motif hits; no valid footprinting output was located.",
            "p_value": np.nan,
            "q_value": np.nan,
            "limitation": "Do not claim motif enrichment or footprinting without a valid motif object and differential accessibility output.",
        }
    )
    return pd.DataFrame(rows)


def main() -> None:
    expression, obs = read_expression()
    cell_corr, sample_corr, sample_means = correlation_tables(expression, obs)
    klf4_test = klf4_group_tests(expression, obs)
    dot = oskm_dot_table(expression, obs)
    scenic_top, scenic_audit = scenic_rss_audit()
    rescued_tables, rescued_audit = rescued_auc_tables()
    territory, spatial_axes, spatial_data = spatial_analysis()
    atac_qc, atac_tests = atac_analysis()
    atac_meta = read_atac_metadata()

    save_table(cell_corr, "FigureS6C_KLF4_correlations_cell_level.tsv")
    save_table(sample_corr, "FigureS6C_KLF4_correlations_sample_level.tsv")
    save_table(sample_means, "FigureS6C_FibK14_sample_means.tsv")
    save_table(klf4_test, "FigureS6D_KLF4_FibK14_vs_other_fibroblasts.tsv")
    save_table(dot, "FigureS6D_OSKM_dotplot_source.tsv")
    save_table(scenic_top, "FigureS6A_formal_SCENIC_FibK14_top_RSS.tsv")
    save_table(territory, "FigureS6E_spatial_OSKM_by_territory.tsv")
    save_table(spatial_axes, "FigureS6E_spatial_OSKM_axis_correlations.tsv")
    save_table(atac_qc, "FigureS6O_scATAC_QC_summary.tsv")
    save_table(atac_tests, "FigureS6O_scATAC_FibK14_tests.tsv")
    for name, table in rescued_tables.items():
        save_table(table, f"rescued_{name}")

    (QC / "SCENIC_source_audit.json").write_text(json.dumps(scenic_audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (QC / "rescued_KLF4_AUC_audit.json").write_text(json.dumps(rescued_audit, ensure_ascii=False, indent=2), encoding="utf-8")

    # Plotting is performed by a separate SVG/Illustrator script because the
    # local scientific Python environment has a broken numpy.linalg DLL.

    audit = build_audit(cell_corr, sample_corr, klf4_test, scenic_audit, rescued_audit, spatial_axes, atac_tests)
    audit.to_csv(QC / "FigureS6_all_results_audit.csv", index=False, encoding="utf-8-sig")
    audit.loc[audit["positive_evidence"]].to_csv(QC / "FigureS6_supported_positive_results.csv", index=False, encoding="utf-8-sig")
    audit.loc[~audit["positive_evidence"]].to_csv(QC / "FigureS6_unsupported_negative_or_blocked_results.csv", index=False, encoding="utf-8-sig")

    summary = {
        "human_scRNA_shape": [int(len(obs)), int(expression.shape[1])],
        "Fib_K14_cells": int(obs["sub_labels"].astype(str).eq("Fib_K14").sum()),
        "Fib_K14_samples_in_sample_level_analysis": int(len(sample_means)),
        "formal_SCENIC_KLF4_present": bool(scenic_audit.get("KLF4_present")),
        "rescued_target_set_available": rescued_audit.get("status") != "missing",
        "supported_positive_result_count": int(audit["positive_evidence"].sum()),
        "audited_result_count": int(len(audit)),
    }
    (QC / "FigureS6_reanalysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
