from __future__ import annotations

import json
import math
import platform
import shutil
import sys
from datetime import datetime
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


PROJECT_DIR = Path(r"data/external/local_f2959f1905d6")
OLD_DIR = Path(r"data/external/local_136961c29c99")
OUT = PROJECT_DIR / "10_supplementary_advanced"
FIG = OUT / "figures"
TAB = OUT / "tables"
REP = OUT / "reports"
LOG = OUT / "logs"
SCRIPT_DIR = PROJECT_DIR / "scripts"
DEG_OLD = OLD_DIR / "DEG_results_all.csv"
LOGCPM = PROJECT_DIR / "03_processed_objects" / "logCPM_filtered_meanCount10.csv"
META = PROJECT_DIR / "00_metadata" / "sample_metadata.csv"
GO_CACHE = OLD_DIR / "cache" / "GO_Biological_Process_2023.gmt.json"
KEGG_CACHE = OLD_DIR / "cache" / "KEGG_2021_Human.gmt.json"
RANDOM_SEED = 20260705
DPI = 600


def mm_px(mm: float) -> int:
    return int(round(mm / 25.4 * DPI))


def font(size_pt: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "arialbd.ttf" if bold else "arial.ttf"
    p = Path(r"data/external/local_81ff0cfd85af") / name
    return ImageFont.truetype(str(p), int(round(size_pt / 72 * DPI)))


F7 = font(7)
F8 = font(8)
F9B = font(9, True)


def init() -> None:
    for p in [OUT, FIG, TAB, REP, LOG, SCRIPT_DIR]:
        p.mkdir(parents=True, exist_ok=True)


def canvas(width_mm: float, height_mm: float) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    img = Image.new("RGB", (mm_px(width_mm), mm_px(height_mm)), "white")
    return img, ImageDraw.Draw(img)


def save_img(img: Image.Image, stem: Path) -> None:
    img.save(stem.with_suffix(".png"), dpi=(DPI, DPI))
    img.save(stem.with_suffix(".pdf"), "PDF", resolution=DPI)


def text(draw: ImageDraw.ImageDraw, xy, s: str, fill=(0, 0, 0), f=F7, anchor=None) -> None:
    draw.text(xy, str(s), fill=fill, font=f, anchor=anchor)


def read_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    deg = pd.read_csv(DEG_OLD).dropna(subset=["symbol"]).copy()
    deg["symbol"] = deg["symbol"].astype(str)
    deg["padj_safe"] = deg["padj"].replace(0, 1e-300).fillna(1).clip(lower=1e-300)
    deg["neglog10_padj"] = -np.log10(deg["padj_safe"])
    expr = pd.read_csv(LOGCPM, index_col=0)
    meta = pd.read_csv(META)
    return deg, expr, meta


def copy_existing_results() -> None:
    for name in ["GO_BP_enrichment.csv", "KEGG_enrichment.csv", "GSEA_preranked_results.csv", "GSVA_ssGSEA_like_scores.csv", "analysis_report.md", "analysis_summary.json"]:
        src = OLD_DIR / name
        if src.exists():
            shutil.copy2(src, OUT / name)
    src = OLD_DIR / "GSEA_classic_top10"
    dst = OUT / "GSEA_classic_top10_integrated"
    if src.exists() and not dst.exists():
        shutil.copytree(src, dst)


def plot_volcano(deg: pd.DataFrame) -> None:
    df = deg.copy()
    df["class"] = "Other"
    df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] >= 1), "class"] = "Up"
    df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] <= -1), "class"] = "Down"
    df["yplot"] = df["neglog10_padj"].clip(upper=60)
    df.to_csv(TAB / "source_data_FigS1_enhanced_volcano.csv", index=False)
    img, d = canvas(105, 85)
    w, h = img.size
    left, right, top, bottom = mm_px(15), mm_px(6), mm_px(8), mm_px(15)
    x0, x1, y0, y1 = left, w - right, h - bottom, top
    xmin, xmax, ymax = -12, 12, 60
    def sx(x): return int(x0 + (np.clip(x, xmin, xmax) - xmin) / (xmax - xmin) * (x1 - x0))
    def sy(y): return int(y0 - np.clip(y, 0, ymax) / ymax * (y0 - y1))
    d.line((x0, y0, x1, y0), fill=(0, 0, 0), width=2)
    d.line((x0, y0, x0, y1), fill=(0, 0, 0), width=2)
    for xv in [-1, 1]:
        px = sx(xv)
        d.line((px, y0, px, y1), fill=(180, 180, 180), width=2)
    py = sy(-math.log10(0.05))
    d.line((x0, py, x1, py), fill=(200, 200, 200), width=2)
    other = df[df["class"] == "Other"].sample(min(6000, (df["class"] == "Other").sum()), random_state=RANDOM_SEED)
    for sub, col, r in [(other, (185, 185, 185), 2), (df[df["class"] == "Down"], (60, 92, 207), 3), (df[df["class"] == "Up"], (230, 75, 53), 3)]:
        for x, y in zip(sub["log2FoldChange"], sub["yplot"]):
            px, py = sx(x), sy(y)
            d.ellipse((px - r, py - r, px + r, py + r), fill=col)
    labels = pd.concat([
        df[df["class"] == "Up"].sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(8),
        df[df["class"] == "Down"].sort_values(["padj", "log2FoldChange"], ascending=[True, True]).head(8),
        df[df["symbol"].eq("KLF4")],
    ]).drop_duplicates("symbol")
    for _, row in labels.iterrows():
        px, py = sx(row["log2FoldChange"]), sy(row["yplot"])
        d.line((px, py, px + mm_px(2), py - mm_px(2)), fill=(90, 90, 90), width=1)
        text(d, (px + mm_px(2.2), py - mm_px(3)), row["symbol"], f=F7)
    text(d, ((x0 + x1) // 2, h - mm_px(6)), "log2FC (KLF4-OE / Control)", f=F8, anchor="mm")
    text(d, (mm_px(3), (y0 + y1) // 2), "-log10 adjusted P-value", f=F8)
    text(d, (x0, top // 2), "Differential expression volcano", f=F9B)
    save_img(img, FIG / "FigS1_enhanced_volcano_reference_style")


def load_sets(path: Path) -> dict[str, set[str]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {k: set(map(str.upper, v)) for k, v in raw.items()}


def ora(query: set[str], bg: set[str], sets: dict[str, set[str]], category: str) -> pd.DataFrame:
    rows = []
    q = query & bg
    for term, genes0 in sets.items():
        genes = genes0 & bg
        ov = q & genes
        if len(ov) < 3:
            continue
        p = stats.hypergeom.sf(len(ov) - 1, len(bg), len(genes), len(q))
        rows.append({"category": category, "term": term, "overlap_n": len(ov), "term_size": len(genes), "query_size": len(q), "pvalue": p, "genes": ";".join(sorted(ov)), "top8_genes": ", ".join(sorted(ov)[:8])})
    out = pd.DataFrame(rows)
    if len(out):
        out["padj"] = multipletests(out["pvalue"], method="fdr_bh")[1]
        out["neglog10_padj"] = -np.log10(out["padj"].clip(lower=1e-300))
        out = out.sort_values(["padj", "pvalue"])
    return out


def plot_bar(df: pd.DataFrame, value: str, label: str, title: str, stem: Path, width_mm=125, height_mm=110) -> None:
    img, d = canvas(width_mm, height_mm)
    w, h = img.size
    left, right, top, bottom = mm_px(55), mm_px(8), mm_px(10), mm_px(10)
    x0, x1, y0, y1 = left, w - right, h - bottom, top
    vals = df[value].astype(float).to_numpy()
    minv, maxv = float(np.nanmin(vals)), float(np.nanmax(vals))
    if minv < 0 < maxv:
        zero_x = int(x0 + (0 - minv) / (maxv - minv) * (x1 - x0))
        scale = lambda v: int(x0 + (v - minv) / (maxv - minv) * (x1 - x0))
    else:
        zero_x = x0
        vmax = max(maxv, 1)
        scale = lambda v: int(x0 + v / vmax * (x1 - x0))
    row_h = (y0 - y1) / max(len(df), 1)
    for i, (_, r) in enumerate(df.iloc[::-1].iterrows()):
        y = int(y0 - (i + 0.7) * row_h)
        xe = scale(float(r[value]))
        color = (230, 75, 53) if "Up" in str(r.get("direction", "")) or r.get("NES", 0) >= 0 else (60, 92, 207)
        d.rectangle((min(zero_x, xe), y - int(row_h * 0.3), max(zero_x, xe), y + int(row_h * 0.3)), fill=color)
        text(d, (x0 - mm_px(1), y), str(r[label])[:54], f=F7, anchor="rm")
    d.line((x0, y0, x1, y0), fill=(0, 0, 0), width=2)
    if minv < 0 < maxv:
        d.line((zero_x, y0, zero_x, y1), fill=(80, 80, 80), width=2)
    text(d, (x0, y1 - mm_px(4)), title, f=F9B)
    text(d, ((x0 + x1) // 2, h - mm_px(4)), value, f=F8, anchor="mm")
    save_img(img, stem)


def split_enrichment(deg: pd.DataFrame) -> None:
    bg = set(deg["symbol"].str.upper())
    up = set(deg.loc[(deg["padj"] < 0.05) & (deg["log2FoldChange"] >= 1), "symbol"].str.upper())
    down = set(deg.loc[(deg["padj"] < 0.05) & (deg["log2FoldChange"] <= -1), "symbol"].str.upper())
    frames = []
    for direction, genes in [("Up_in_OE", up), ("Down_in_OE", down)]:
        for cat, sets in [("GO_BP", load_sets(GO_CACHE)), ("KEGG", load_sets(KEGG_CACHE))]:
            x = ora(genes, bg, sets, cat)
            if len(x):
                x.insert(0, "direction", direction)
                frames.append(x)
    enr = pd.concat(frames, ignore_index=True)
    enr.to_csv(TAB / "split_up_down_GO_KEGG_enrichment.csv", index=False)
    plot_df = enr.sort_values("padj").groupby(["direction", "category"], group_keys=False).head(5).copy()
    plot_df["label"] = plot_df["term"].str.replace(r"\s*\(GO:\d+\)", "", regex=True)
    plot_df.to_csv(TAB / "source_data_FigS2_split_enrichment.csv", index=False)
    plot_bar(plot_df.sort_values("neglog10_padj"), "neglog10_padj", "label", "Up/down DEG functional enrichment", FIG / "FigS2_up_down_functional_enrichment", 125, 120)


def plot_gsea_gsva() -> None:
    gsea = pd.read_csv(OLD_DIR / "GSEA_preranked_results.csv")
    gsea.to_csv(TAB / "GSEA_preranked_results_integrated.csv", index=False)
    top = gsea.sort_values(["padj", "pvalue"]).head(15).copy()
    top["label"] = top["term"].str.slice(0, 52)
    top.to_csv(TAB / "source_data_FigS3_GSEA_NES.csv", index=False)
    plot_bar(top.sort_values("NES"), "NES", "label", "Pre-ranked GSEA summary", FIG / "FigS3_GSEA_NES_summary", 125, 105)
    gsva = pd.read_csv(OLD_DIR / "GSVA_ssGSEA_like_scores.csv", index_col=0)
    z = gsva.sub(gsva.mean(axis=1), axis=0).div(gsva.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    z.to_csv(TAB / "source_data_FigS4_GSVA_ssGSEA_like_heatmap_zscore.csv")
    img, d = canvas(115, 95)
    w, h = img.size
    left, topm = mm_px(55), mm_px(10)
    cell_w = (w - left - mm_px(8)) / z.shape[1]
    cell_h = (h - topm - mm_px(12)) / z.shape[0]
    lim = max(1, float(np.abs(z.to_numpy()).max()))
    for i, idx in enumerate(z.index):
        text(d, (left - mm_px(1), int(topm + (i + 0.5) * cell_h)), idx[:45], f=F7, anchor="rm")
        for j, col in enumerate(z.columns):
            v = float(z.iloc[i, j])
            t = (v + lim) / (2 * lim)
            rgb = (int(60 + 180 * t), int(80 + 80 * (1 - abs(t - 0.5) * 2)), int(210 - 160 * t))
            x = int(left + j * cell_w)
            y = int(topm + i * cell_h)
            d.rectangle((x, y, int(x + cell_w), int(y + cell_h)), fill=rgb)
    for j, col in enumerate(z.columns):
        text(d, (int(left + (j + 0.5) * cell_w), h - mm_px(7)), col, f=F7, anchor="mm")
    text(d, (left, mm_px(4)), "GSVA/ssGSEA-like pathway activity", f=F9B)
    save_img(img, FIG / "FigS4_GSVA_ssGSEA_like_heatmap")


def align_expr(expr: pd.DataFrame, deg: pd.DataFrame) -> pd.DataFrame:
    mp = deg.drop_duplicates("gene_id").set_index("gene_id")["symbol"].to_dict()
    x = expr.copy()
    x["symbol"] = [mp.get(g, g.split(".")[0]) for g in x.index]
    return x.groupby("symbol").mean()


def ml_and_network(expr_sym: pd.DataFrame, deg: pd.DataFrame, meta: pd.DataFrame) -> None:
    samples = meta["specimen"].tolist()
    y = (meta["condition"] == "OE").astype(int).to_numpy()
    features = [g for g in deg.sort_values("padj").drop_duplicates("symbol").head(30)["symbol"] if g in expr_sym.index]
    X = expr_sym.loc[features, samples].T
    Xs = StandardScaler().fit_transform(X)
    preds = []
    for train, test in LeaveOneOut().split(Xs):
        clf = LogisticRegression(solver="liblinear", random_state=RANDOM_SEED).fit(Xs[train], y[train])
        preds.append(int(clf.predict(Xs[test])[0]))
    metrics = {"n_samples": 4, "n_features": len(features), "loocv_accuracy": float(accuracy_score(y, preds)), "loocv_balanced_accuracy": float(balanced_accuracy_score(y, preds)), "note": "Exploratory only; n=4."}
    (LOG / "machine_learning_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    pd.DataFrame({"sample": samples, "condition": meta["condition"], "y_true_OE": y, "loocv_pred_OE": preds}).to_csv(TAB / "machine_learning_LOOCV_predictions.csv", index=False)
    clf = LogisticRegression(solver="liblinear", random_state=RANDOM_SEED).fit(Xs, y)
    imp = pd.DataFrame({"symbol": features, "SHAP_fallback_score": np.abs(clf.coef_[0])}).sort_values("SHAP_fallback_score", ascending=False)
    imp.to_csv(TAB / "SHAP_fallback_model_feature_importance.csv", index=False)
    plot_bar(imp.head(15).sort_values("SHAP_fallback_score"), "SHAP_fallback_score", "symbol", "Exploratory model explanation", FIG / "FigS5_SHAP_fallback_feature_importance", 95, 85)
    ctab = pd.crosstab(pd.Series(y, name="Observed"), pd.Series(preds, name="Predicted"))
    ctab.to_csv(TAB / "source_data_FigS6_machine_learning_confusion.csv")
    img, d = canvas(65, 60)
    left, topm, size = mm_px(22), mm_px(14), mm_px(12)
    for i in range(2):
        for j in range(2):
            val = int(ctab.iloc[i, j]) if i < ctab.shape[0] and j < ctab.shape[1] else 0
            shade = 245 - val * 70
            d.rectangle((left + j * size, topm + i * size, left + (j + 1) * size, topm + (i + 1) * size), fill=(shade, shade, 255), outline=(0, 0, 0))
            text(d, (left + (j + 0.5) * size, topm + (i + 0.5) * size), val, f=F8, anchor="mm")
    text(d, (left, mm_px(6)), "LOOCV classifier confusion", f=F9B)
    save_img(img, FIG / "FigS6_machine_learning_LOOCV_confusion")

    top = [g for g in deg.sort_values("padj").drop_duplicates("symbol").head(200)["symbol"] if g in expr_sym.index]
    mat = expr_sym.loc[top, samples]
    corr = mat.T.corr()
    edges = []
    for i, a in enumerate(top):
        for b in top[i + 1:]:
            r = float(corr.loc[a, b])
            if abs(r) >= 0.95:
                edges.append((a, b, r))
    edge_df = pd.DataFrame(edges, columns=["source", "target", "pearson_r"])
    edge_df.to_csv(TAB / "WGCNA_fallback_coexpression_edges_absR095.csv", index=False)
    G = nx.Graph()
    G.add_nodes_from(top)
    for a, b, r in edges:
        G.add_edge(a, b, weight=abs(r), pearson_r=r)
    hubs = pd.DataFrame({"symbol": list(G.nodes), "degree": [G.degree(n) for n in G.nodes]}).sort_values("degree", ascending=False)
    hubs.to_csv(TAB / "WGCNA_fallback_hub_genes.csv", index=False)
    nodes = hubs.head(28)["symbol"].tolist()
    H = G.subgraph(nodes).copy()
    pos = nx.circular_layout(H)
    img, d = canvas(120, 110)
    cx, cy, rad = mm_px(60), mm_px(55), mm_px(38)
    pxy = {n: (int(cx + pos[n][0] * rad), int(cy + pos[n][1] * rad)) for n in H.nodes}
    for a, b in H.edges:
        d.line((*pxy[a], *pxy[b]), fill=(180, 180, 180), width=1)
    for n, (x, yy) in pxy.items():
        rr = 5 + min(12, H.degree(n))
        d.ellipse((x - rr, yy - rr, x + rr, yy + rr), fill=(76, 120, 168), outline=(255, 255, 255))
        text(d, (x + rr + 2, yy), n, f=F7)
    text(d, (mm_px(8), mm_px(5)), "Co-expression hub network", f=F9B)
    save_img(img, FIG / "FigS7_WGCNA_fallback_coexpression_hub_network")
    edge_df.to_csv(TAB / "PPI_proxy_network_from_coexpression_NOT_STRING.csv", index=False)


def status_tables(deg: pd.DataFrame) -> None:
    tf_seeds = {"KLF4", "RELA", "IRF3", "FOXO3", "FOSL1", "FOXC1", "SOX4", "SIX1", "SIX2", "TFAP2A", "ZBTB16", "HAND2", "ARID2", "TCF21", "NFKBIB"}
    tfs = deg[deg["symbol"].isin(tf_seeds)].sort_values(["padj", "pvalue"])[["symbol", "description", "log2FoldChange", "pvalue", "padj", "direction"]]
    tfs.to_csv(TAB / "TF_regulatory_network_seed_nodes_current_data.csv", index=False)
    pd.DataFrame([{"target": r.symbol, "current_evidence": f"log2FC={r.log2FoldChange:.3g}; padj={r.padj:.3g}", "drug_target_next_database": "DGIdb/ChEMBL/DrugBank/PubChem/SwissTargetPrediction", "docking_status": "Not run: requires compound SMILES/SDF, receptor structure, docking box, and engine."} for r in tfs.itertuples()]).to_csv(TAB / "drug_targeting_and_docking_candidate_matrix.csv", index=False)
    status = pd.DataFrame([
        ("Clinical univariate regression", "not_run", "No clinical outcome/covariates; only condition/source metadata for 4 samples.", "clinical table with outcome and covariates"),
        ("Clinical multivariable regression", "not_run", "n=4 and no covariates; model would be singular.", "larger clinical cohort with predefined formula"),
        ("Canonical WGCNA", "fallback_only", "n=4 cannot support stable scale-free topology/module detection.", "prefer >=15-20 samples with traits"),
        ("SHAP", "fallback_only", "True SHAP not used because Python/R plotting/binary stack was unstable; logistic-coefficient model explanation provided.", "stable SHAP stack and larger cohort"),
        ("True TF-target regulatory network", "seed_table_only", "No TF-target/motif/ATAC/ChIP data; differential TF seeds listed.", "DoRothEA/TRRUST/ChEA/JASPAR or ATAC/ChIP data"),
        ("STRING/Cytoscape PPI", "proxy_only", "Current network is co-expression proxy, not physical PPI.", "STRING/Cytoscape export"),
        ("Molecular docking", "not_run", "No compound/receptor structures or docking configuration.", "ligands, structures, binding pocket, Vina/smina config"),
    ], columns=["module", "status", "reason", "required_input"])
    status.to_csv(TAB / "data_sufficiency_and_module_status.csv", index=False)


def write_report() -> None:
    status = pd.read_csv(TAB / "data_sufficiency_and_module_status.csv")
    md_rows = ["| module | status | reason | required_input |", "|---|---|---|---|"]
    for r in status.itertuples(index=False):
        md_rows.append(f"| {r.module} | {r.status} | {r.reason} | {r.required_input} |")
    lines = [
        "# Advanced Bulk RNA-seq Supplementary Analysis",
        "",
        f"Generated: `{datetime.now().isoformat(timespec='seconds')}`",
        "",
        "## Added Modules",
        "",
        "- FigS1 reference-style volcano plot.",
        "- FigS2 up/down GO/KEGG enrichment.",
        "- FigS3 pre-ranked GSEA NES summary.",
        "- FigS4 GSVA/ssGSEA-like heatmap.",
        "- FigS5 SHAP-fallback model explanation.",
        "- FigS6 LOOCV machine-learning confusion matrix.",
        "- FigS7 WGCNA-fallback co-expression hub network.",
        "",
        "## Data Sufficiency",
        "",
        "\n".join(md_rows),
        "",
        "All figures were generated as PNG and PDF with Arial 7 pt minimum and physical dimensions within A4 constraints.",
    ]
    (REP / "advanced_supplementary_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    np.random.seed(RANDOM_SEED)
    init()
    deg, expr, meta = read_data()
    copy_existing_results()
    plot_volcano(deg)
    split_enrichment(deg)
    plot_gsea_gsva()
    expr_sym = align_expr(expr, deg)
    ml_and_network(expr_sym, deg, meta)
    status_tables(deg)
    write_report()
    shutil.copy2(Path(__file__), SCRIPT_DIR / "run_bulk_rnaseq_advanced_supplement_pil.py")
    summary = {"output": str(OUT).encode("unicode_escape").decode("ascii"), "figures_pdf": len(list(FIG.glob("*.pdf"))), "tables_csv": len(list(TAB.glob("*.csv"))), "python": sys.version.split()[0], "platform": platform.platform()}
    (LOG / "advanced_supplement_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
