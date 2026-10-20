from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


INPUT_XLSX = Path(r"data/external/local_3712a2ef50e6")
RESULT_ROOT = Path(r"data/external/local_54ab21a5fa6d")
DATASET_NAME = "课题01+bulkRNAseq+KLF4_OE+Control_vs_OE+20260705"
PROJECT_DIR = RESULT_ROOT / DATASET_NAME
RANDOM_SEED = 20260705
CONTROL_LABEL = "Control"
CASE_LABEL = "OE"
KLF4_ENSEMBL = "ENSG00000136826"


def mm_to_in(mm: float) -> float:
    return mm / 25.4


def ensure_dirs(base: Path) -> dict[str, Path]:
    dirs = {
        "metadata": base / "00_metadata",
        "raw": base / "01_raw_data",
        "qc": base / "02_qc",
        "processed": base / "03_processed_objects",
        "per_omics": base / "04_per_omics" / "bulk_rnaseq",
        "integration": base / "05_integration",
        "figures": base / "06_figures",
        "tables": base / "07_tables",
        "reports": base / "08_reports",
        "logs": base / "09_logs",
        "scripts": base / "scripts",
    }
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return dirs


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def setup_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.dpi": 150,
            "savefig.dpi": 600,
        }
    )
    sns.set_theme(style="whitegrid", rc={"font.family": "Arial", "font.size": 7})


def savefig(fig: plt.Figure, stem: Path, width_mm: float, height_mm: float) -> None:
    fig.set_size_inches(mm_to_in(width_mm), mm_to_in(height_mm))
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=600)
    plt.close(fig)


def load_workbook(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = pd.read_excel(path, sheet_name="summary")
    counts = pd.read_excel(path, sheet_name="raw_counts")
    metadata = pd.read_excel(path, sheet_name="metadata")
    return summary, counts, metadata


def prepare_counts(counts: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_cols = metadata["specimen"].astype(str).tolist()
    missing = [c for c in sample_cols if c not in counts.columns]
    if missing:
        raise ValueError(f"Sample columns missing from raw_counts: {missing}")
    gene_ids = counts["gene_id"].astype(str)
    matrix = counts.loc[:, sample_cols].copy()
    matrix.index = gene_ids
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
    matrix = matrix.groupby(matrix.index).sum()
    metadata = metadata.set_index("specimen").loc[sample_cols].reset_index()
    return matrix, metadata


def compute_logcpm(count_matrix: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    library_sizes = count_matrix.sum(axis=0)
    cpm = count_matrix.div(library_sizes, axis=1) * 1_000_000
    return np.log2(cpm + 1), library_sizes


def differential_expression(logcpm: pd.DataFrame, count_matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    control_samples = metadata.loc[metadata["condition"] == CONTROL_LABEL, "specimen"].tolist()
    case_samples = metadata.loc[metadata["condition"] == CASE_LABEL, "specimen"].tolist()
    c = logcpm[control_samples]
    k = logcpm[case_samples]

    mean_control = c.mean(axis=1)
    mean_case = k.mean(axis=1)
    log2fc = mean_case - mean_control
    pvals = []
    for gene in logcpm.index:
        x = c.loc[gene].to_numpy(dtype=float)
        y = k.loc[gene].to_numpy(dtype=float)
        if np.allclose(x, x[0]) and np.allclose(y, y[0]) and np.isclose(x[0], y[0]):
            pvals.append(1.0)
        else:
            pvals.append(stats.ttest_ind(y, x, equal_var=False).pvalue)
    pvals = np.nan_to_num(np.array(pvals, dtype=float), nan=1.0, posinf=1.0, neginf=1.0)
    padj = multipletests(pvals, method="fdr_bh")[1]

    base_mean = count_matrix.mean(axis=1)
    out = pd.DataFrame(
        {
            "gene_id": logcpm.index,
            "gene_id_no_version": pd.Index(logcpm.index).str.replace(r"\.\d+$", "", regex=True),
            "baseMean_raw_count": base_mean.values,
            "mean_logCPM_Control": mean_control.values,
            "mean_logCPM_OE": mean_case.values,
            "log2FC_OE_vs_Control": log2fc.values,
            "pvalue_welch_logCPM": pvals,
            "padj_BH": padj,
        }
    )
    out["direction"] = np.select(
        [
            (out["padj_BH"] < 0.05) & (out["log2FC_OE_vs_Control"] >= 1),
            (out["padj_BH"] < 0.05) & (out["log2FC_OE_vs_Control"] <= -1),
        ],
        ["Up_in_OE", "Down_in_OE"],
        default="Not_significant",
    )
    return out.sort_values(["padj_BH", "pvalue_welch_logCPM", "log2FC_OE_vs_Control"], ascending=[True, True, False])


def plot_library_sizes(library_sizes: pd.Series, metadata: pd.DataFrame, figdir: Path, tabledir: Path) -> None:
    df = pd.DataFrame({"sample": library_sizes.index, "library_size": library_sizes.values})
    df = df.merge(metadata[["specimen", "condition", "source"]], left_on="sample", right_on="specimen", how="left")
    df.to_csv(tabledir / "source_data_Fig1_library_sizes.csv", index=False)
    fig, ax = plt.subplots()
    sns.barplot(data=df, x="sample", y="library_size", hue="condition", palette=["#4C78A8", "#D55E00"], ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("Raw library size")
    ax.set_title("Library size")
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    savefig(fig, figdir / "Fig1_library_sizes", 85, 65)


def plot_pca(logcpm_filtered: pd.DataFrame, metadata: pd.DataFrame, figdir: Path, tabledir: Path) -> None:
    x = logcpm_filtered.T
    x_scaled = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2, random_state=RANDOM_SEED)
    coords = pca.fit_transform(x_scaled)
    df = metadata.copy()
    df["PC1"] = coords[:, 0]
    df["PC2"] = coords[:, 1]
    df.to_csv(tabledir / "source_data_Fig2_PCA.csv", index=False)
    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x="PC1", y="PC2", hue="condition", style="source", s=55, palette=["#4C78A8", "#D55E00"], ax=ax)
    for _, row in df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["specimen"], fontsize=7, ha="left", va="bottom")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("PCA of filtered logCPM")
    savefig(fig, figdir / "Fig2_PCA_filtered_logCPM", 90, 75)


def plot_correlation(logcpm_filtered: pd.DataFrame, figdir: Path, tabledir: Path) -> None:
    corr = logcpm_filtered.corr(method="pearson")
    corr.to_csv(tabledir / "source_data_Fig3_sample_correlation.csv")
    fig, ax = plt.subplots()
    sns.heatmap(corr, vmin=0.7, vmax=1.0, cmap="viridis", annot=True, fmt=".3f", square=True, cbar_kws={"label": "Pearson r"}, ax=ax)
    ax.set_title("Sample correlation")
    savefig(fig, figdir / "Fig3_sample_correlation", 80, 70)


def plot_volcano(deg: pd.DataFrame, figdir: Path, tabledir: Path) -> None:
    df = deg.copy()
    df["neg_log10_padj"] = -np.log10(df["padj_BH"].clip(lower=1e-300))
    df.to_csv(tabledir / "source_data_Fig4_volcano.csv", index=False)
    top = df.sort_values(["padj_BH", "pvalue_welch_logCPM"]).head(12)
    palette = {"Up_in_OE": "#D55E00", "Down_in_OE": "#0072B2", "Not_significant": "#9AA0A6"}
    fig, ax = plt.subplots()
    sns.scatterplot(data=df, x="log2FC_OE_vs_Control", y="neg_log10_padj", hue="direction", palette=palette, s=7, linewidth=0, alpha=0.8, ax=ax)
    ax.axvline(1, ls="--", lw=0.7, color="black")
    ax.axvline(-1, ls="--", lw=0.7, color="black")
    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="black")
    for _, row in top.iterrows():
        ax.text(row["log2FC_OE_vs_Control"], row["neg_log10_padj"], row["gene_id_no_version"], fontsize=7)
    ax.set_xlabel("log2FC (OE vs Control)")
    ax.set_ylabel("-log10(BH-adjusted P)")
    ax.set_title("Exploratory differential expression")
    savefig(fig, figdir / "Fig4_volcano_exploratory_DEG", 105, 85)


def plot_heatmap(logcpm_filtered: pd.DataFrame, deg: pd.DataFrame, metadata: pd.DataFrame, figdir: Path, tabledir: Path) -> None:
    top_genes = deg.sort_values(["padj_BH", "pvalue_welch_logCPM"]).head(40)["gene_id"].tolist()
    mat = logcpm_filtered.loc[top_genes]
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    z.index = pd.Index(z.index).str.replace(r"\.\d+$", "", regex=True)
    z.to_csv(tabledir / "source_data_Fig5_top40_heatmap_zscore.csv")
    lut = {CONTROL_LABEL: "#4C78A8", CASE_LABEL: "#D55E00"}
    col_colors = metadata.set_index("specimen").loc[z.columns, "condition"].map(lut)
    grid = sns.clustermap(
        z,
        cmap="vlag",
        center=0,
        col_colors=col_colors,
        figsize=(mm_to_in(105), mm_to_in(120)),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Row z-score"},
    )
    grid.fig.suptitle("Top exploratory DE genes", y=1.02, fontsize=9, fontweight="bold")
    grid.savefig(figdir / "Fig5_top40_DEG_heatmap.pdf", bbox_inches="tight")
    grid.savefig(figdir / "Fig5_top40_DEG_heatmap.png", bbox_inches="tight", dpi=600)
    plt.close(grid.fig)


def plot_klf4(count_matrix: pd.DataFrame, logcpm: pd.DataFrame, metadata: pd.DataFrame, figdir: Path, tabledir: Path) -> dict[str, object]:
    stripped = pd.Index(count_matrix.index).str.replace(r"\.\d+$", "", regex=True)
    matches = count_matrix.index[stripped == KLF4_ENSEMBL].tolist()
    status: dict[str, object] = {"klf4_ensembl": KLF4_ENSEMBL, "found": bool(matches), "matched_gene_ids": matches}
    if not matches:
        return status
    gene = matches[0]
    df = metadata[["specimen", "condition", "source"]].copy()
    df["raw_count"] = count_matrix.loc[gene, df["specimen"]].to_numpy()
    df["logCPM"] = logcpm.loc[gene, df["specimen"]].to_numpy()
    df.to_csv(tabledir / "source_data_Fig6_KLF4_counts.csv", index=False)
    fig, ax = plt.subplots()
    sns.barplot(data=df, x="condition", y="logCPM", palette=["#4C78A8", "#D55E00"], errorbar="sd", ax=ax)
    sns.stripplot(data=df, x="condition", y="logCPM", color="black", size=4, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("KLF4 log2(CPM + 1)")
    ax.set_title("KLF4 expression")
    savefig(fig, figdir / "Fig6_KLF4_expression", 65, 60)
    status["table"] = str(tabledir / "source_data_Fig6_KLF4_counts.csv")
    return status


def write_report(dirs: dict[str, Path], summary: dict[str, object]) -> None:
    report = f"""# Bulk RNA-seq Analysis Report

Dataset: `{DATASET_NAME}`  
Input: `{INPUT_XLSX}`  
Generated: `{datetime.now().isoformat(timespec='seconds')}`

## Design

- Comparison: `{CASE_LABEL}` vs `{CONTROL_LABEL}`.
- Samples: Control n=2 (`C7-1`, `C7-2`), OE n=2 (`K7-1`, `K7-2`).
- Important caution: n=2 per group is underpowered for definitive bulk RNA-seq inference. Results are exploratory and should be validated with a DESeq2/edgeR/limma-voom model once the R/Bioconductor stack is available or with additional replicates.

## Methods

1. Read `summary`, `metadata`, and `raw_counts` sheets from the Excel workbook.
2. Verified metadata sample IDs match count matrix columns.
3. Removed genes with mean raw count < 10.
4. Normalized counts to CPM and transformed as `log2(CPM + 1)`.
5. Ran PCA and Pearson sample-correlation QC on filtered genes.
6. Ran exploratory gene-wise Welch t-tests on logCPM values, comparing OE against Control.
7. Adjusted P values by Benjamini-Hochberg FDR.
8. Defined exploratory DE genes as `padj < 0.05` and `abs(log2FC) >= 1`.

## CNS Figure Constraints

- Font family: Arial requested in matplotlib.
- Minimum font size: 7 pt.
- Figure widths: 65-105 mm, within the required 42-210 mm range.
- Figure heights: 60-120 mm, within the required 49.5-148.5 mm range.
- Exports: PDF and PNG at 600 dpi.
- Source data saved for every figure.

## Key Results

- Input genes: {summary['n_genes_input']:,}
- Genes retained after mean raw count >= 10: {summary['n_genes_filtered']:,}
- Up in OE, exploratory FDR < 0.05 and log2FC >= 1: {summary['n_up']:,}
- Down in OE, exploratory FDR < 0.05 and log2FC <= -1: {summary['n_down']:,}
- KLF4 Ensembl ID `{KLF4_ENSEMBL}` found: {summary['klf4_status']['found']}

## Main Outputs

- DEG table: `07_tables/DEG_exploratory_welch_logCPM_OE_vs_Control.csv`
- Normalized matrix: `03_processed_objects/logCPM_filtered_meanCount10.csv`
- Figures: `06_figures`
- Reproducibility log: `09_logs/run_summary.json`

## Not Run / Pending

- DESeq2/edgeR/limma-voom: pending because the current Windows conda/R stack did not have these Bioconductor tools fully installed at the time of this run.
- Enrichment/GSEA: not run because this workbook contains Ensembl gene IDs without a verified gene-symbol/Entrez annotation table. This avoids silently mapping wrong gene identifiers.
"""
    (dirs["reports"] / "bulk_rnaseq_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    np.random.seed(RANDOM_SEED)
    setup_plotting()
    dirs = ensure_dirs(PROJECT_DIR)

    summary_sheet, raw_counts, metadata = load_workbook(INPUT_XLSX)
    count_matrix, metadata = prepare_counts(raw_counts, metadata)
    logcpm, library_sizes = compute_logcpm(count_matrix)
    keep = count_matrix.mean(axis=1) >= 10
    count_filtered = count_matrix.loc[keep]
    logcpm_filtered = logcpm.loc[keep]
    deg = differential_expression(logcpm_filtered, count_filtered, metadata)

    summary_sheet.to_csv(dirs["metadata"] / "input_workbook_summary_sheet.csv", index=False)
    metadata.to_csv(dirs["metadata"] / "sample_metadata.csv", index=False)
    count_matrix.to_csv(dirs["raw"] / "raw_counts_matrix.csv")
    count_filtered.to_csv(dirs["processed"] / "raw_counts_filtered_meanCount10.csv")
    logcpm_filtered.to_csv(dirs["processed"] / "logCPM_filtered_meanCount10.csv")
    deg.to_csv(dirs["tables"] / "DEG_exploratory_welch_logCPM_OE_vs_Control.csv", index=False)
    deg.head(100).to_csv(dirs["tables"] / "top100_DEG_exploratory.csv", index=False)

    plot_library_sizes(library_sizes, metadata, dirs["figures"], dirs["tables"])
    plot_pca(logcpm_filtered, metadata, dirs["figures"], dirs["tables"])
    plot_correlation(logcpm_filtered, dirs["figures"], dirs["tables"])
    plot_volcano(deg, dirs["figures"], dirs["tables"])
    plot_heatmap(logcpm_filtered, deg, metadata, dirs["figures"], dirs["tables"])
    klf4_status = plot_klf4(count_matrix, logcpm, metadata, dirs["figures"], dirs["tables"])

    manifest = {
        "dataset": DATASET_NAME,
        "input_xlsx": str(INPUT_XLSX),
        "input_sha256": sha256_file(INPUT_XLSX),
        "project_dir": str(PROJECT_DIR),
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version,
        "platform": platform.platform(),
        "random_seed": RANDOM_SEED,
        "analysis": "exploratory Welch t-test on log2(CPM+1), BH FDR",
        "filter": "mean raw count >= 10",
        "figure_constraints": {
            "font": "Arial",
            "min_font_pt": 7,
            "single_width_mm_range": [42, 210],
            "single_height_mm_range": [49.5, 148.5],
            "export": "PDF + PNG 600 dpi",
        },
        "n_genes_input": int(count_matrix.shape[0]),
        "n_genes_filtered": int(count_filtered.shape[0]),
        "n_samples": int(count_matrix.shape[1]),
        "n_up": int((deg["direction"] == "Up_in_OE").sum()),
        "n_down": int((deg["direction"] == "Down_in_OE").sum()),
        "klf4_status": klf4_status,
        "packages": {
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "matplotlib": mpl.__version__,
            "seaborn": sns.__version__,
        },
    }
    (dirs["logs"] / "run_summary.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    try:
        env_export = subprocess.run(
            [r"data/external/local_3a50d8ec5ea3", "env", "export", "-n", "panomics-meta-20260705"],
            text=True,
            capture_output=True,
            check=False,
        )
        (dirs["logs"] / "conda_env_export.yml").write_text(env_export.stdout, encoding="utf-8")
    except Exception as exc:
        (dirs["logs"] / "conda_env_export_error.txt").write_text(repr(exc), encoding="utf-8")

    write_report(dirs, manifest)
    console_summary = {
        "dataset": DATASET_NAME.encode("unicode_escape").decode("ascii"),
        "project_dir": str(PROJECT_DIR).encode("unicode_escape").decode("ascii"),
        "n_genes_input": manifest["n_genes_input"],
        "n_genes_filtered": manifest["n_genes_filtered"],
        "n_up": manifest["n_up"],
        "n_down": manifest["n_down"],
        "klf4_found": manifest["klf4_status"]["found"],
    }
    print(json.dumps(console_summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
