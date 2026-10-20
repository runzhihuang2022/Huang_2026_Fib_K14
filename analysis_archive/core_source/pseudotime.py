from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import pymupdf


WORK = Path(r"data/external/local_3fdd74c6b78b")
OUT = Path(r"data/external/local_68eb4cd5a8f6")
M2 = WORK / "monocle2_three_sets" / "rerooted_full_sub_labels_20260826"
M3 = WORK / "KC_SAC_Monocle3_recluster_20260903"

ANALYSES = {
    "M2: Fib_K14 + all KC + all SAC": M2 / "set1_FibK14_allKC_allSAC" / "monocle2_rerooted_coordinates.csv",
    "M2: Fib_K14 + all KC": M2 / "set2_FibK14_allKC" / "monocle2_rerooted_coordinates.csv",
    "M2: all KC + all SAC": M2 / "set3_allKC_allSAC" / "monocle2_rerooted_coordinates.csv",
    "M3: all KC + all SAC": M3 / "03_analysis" / "KC_SAC_monocle3_pseudotime.tsv.gz",
}

FIGURE_SOURCES = {
    "Figure_R11A_Monocle2_FibK14_allKC_allSAC": M2 / "strips" / "set1_FibK14_allKC_allSAC_rerooted_three_panel_strip",
    "Figure_R11B_Monocle2_FibK14_allKC": M2 / "strips" / "set2_FibK14_allKC_rerooted_three_panel_strip",
    "Figure_R11C_Monocle2_allKC_allSAC": M2 / "strips" / "set3_allKC_allSAC_rerooted_three_panel_strip",
    "Figure_R11D_Monocle3_KC_SAC_reclustered": M3 / "05_figures" / "KC_SAC_reclustered_Monocle3",
}

SUBTYPE_ORDER = [
    "Fib_K14",
    "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
    "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
    "SAC_SG_Progenitor", "SAC_SG_Clear", "SAC_SG_Dark", "SAC_SG_Ductal",
    "SAC_HF_HFSC", "SAC_HF_ORS", "SAC_HF_IRS", "SAC_HF_DP_DS", "SAC_HF_Matrix",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_analysis(name: str, path: Path) -> pd.DataFrame:
    if name.startswith("M3"):
        x = pd.read_csv(path, sep="\t")
        x = x.rename(columns={"cell": "cell_id"})
    else:
        x = pd.read_csv(path)
        x = x.rename(columns={"Pseudotime": "pseudotime", "time_bin": "observed_healing_time"})
    required = {"cell_id", "sub_labels", "pseudotime"}
    missing = required.difference(x.columns)
    if missing:
        raise ValueError(f"{name}: missing columns {sorted(missing)}")
    x = x.loc[np.isfinite(pd.to_numeric(x["pseudotime"], errors="coerce"))].copy()
    x["pseudotime"] = pd.to_numeric(x["pseudotime"])
    lo, hi = x["pseudotime"].min(), x["pseudotime"].max()
    x["pseudotime_01"] = (x["pseudotime"] - lo) / (hi - lo) if hi > lo else 0.0
    x["analysis"] = name
    return x


def copy_source_figures() -> list[Path]:
    copied: list[Path] = []
    for out_stem, src_stem in FIGURE_SOURCES.items():
        for ext in (".pdf", ".png"):
            src = src_stem.with_suffix(ext)
            dst = OUT / "individual_figures" / f"{out_stem}{ext}"
            shutil.copy2(src, dst)
            copied.append(dst)
        for ext in (".svg", ".tiff"):
            src = src_stem.with_suffix(ext)
            if src.exists():
                dst = OUT / "individual_figures" / f"{out_stem}{ext}"
                shutil.copy2(src, dst)
                copied.append(dst)
    return copied


def make_multipage_pdf() -> Path:
    out = pymupdf.open()
    toc = []
    for i, (stem, _) in enumerate(FIGURE_SOURCES.items(), start=1):
        src = OUT / "individual_figures" / f"{stem}.pdf"
        doc = pymupdf.open(src)
        start_page = out.page_count + 1
        out.insert_pdf(doc)
        toc.append([1, stem.replace("Figure_R11", "R11"), start_page])
        doc.close()
    out.set_toc(toc)
    dst = OUT / "Figure_R11_pseudotime_sensitivity_4pages_vector.pdf"
    out.save(dst, garbage=4, deflate=True)
    out.close()
    return dst


def render_individual_600dpi() -> list[Path]:
    rendered: list[Path] = []
    for stem in FIGURE_SOURCES:
        src = OUT / "individual_figures" / f"{stem}.pdf"
        doc = pymupdf.open(src)
        page = doc[0]
        pix = page.get_pixmap(dpi=600, alpha=False)
        image = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        dst = OUT / "individual_figures" / f"{stem}_600dpi.tiff"
        image.save(dst, dpi=(600, 600), compression="tiff_lzw")
        rendered.append(dst)
        doc.close()
    return rendered


def make_contact_sheet() -> tuple[Path, Path]:
    pngs = [OUT / "individual_figures" / f"{stem}.png" for stem in FIGURE_SOURCES]
    opened = [Image.open(p).convert("RGB") for p in pngs]
    cell_w, cell_h = 2400, 1200
    margin, gap = 120, 80
    sheet = Image.new("RGB", (margin * 2 + cell_w * 2 + gap, margin * 2 + cell_h * 2 + gap), "white")
    draw = ImageDraw.Draw(sheet)
    font_path = Path(r"data/external/local_28ac8473f586")
    font = ImageFont.truetype(str(font_path), 66) if font_path.exists() else ImageFont.load_default()
    labels = ["A", "B", "C", "D"]
    for i, (im, label) in enumerate(zip(opened, labels)):
        im.thumbnail((cell_w - 20, cell_h - 20), Image.Resampling.LANCZOS)
        x0 = margin + (i % 2) * (cell_w + gap)
        y0 = margin + (i // 2) * (cell_h + gap)
        x = x0 + (cell_w - im.width) // 2
        y = y0 + (cell_h - im.height) // 2
        sheet.paste(im, (x, y))
        draw.text((x0, y0 - 10), label, fill="black", font=font)
    png = OUT / "Figure_R11_overview_contact_sheet.png"
    tif = OUT / "Figure_R11_overview_contact_sheet_600dpi.tiff"
    sheet.save(png, dpi=(300, 300), optimize=True)
    sheet.save(tif, dpi=(600, 600), compression="tiff_lzw")
    return png, tif


def make_consistency_figure(medians: pd.DataFrame, corr: pd.DataFrame) -> list[Path]:
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    pivot = medians.pivot(index="sub_labels", columns="analysis", values="median_pseudotime_01")
    ordered = [x for x in SUBTYPE_ORDER if x in pivot.index]
    pivot = pivot.reindex(ordered)
    analysis_order = list(ANALYSES)
    pivot = pivot.reindex(columns=analysis_order)

    cmap1 = mpl.colormaps["viridis"].copy()
    cmap1.set_bad("#f2f2f2")
    fig = plt.figure(figsize=(7.1, 9.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.25, 1.0], hspace=0.56)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(pivot.values, cmap=cmap1, vmin=0, vmax=1, aspect="auto")
    ax.set_title("Subtype median pseudotime", fontweight="bold", pad=8)
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    short = ["M2: Fib+KC+SAC", "M2: Fib+KC", "M2: KC+SAC", "M3: KC+SAC"]
    ax.set_xticks(range(len(short)), short, rotation=55, ha="right")
    ax.set_xlabel("Within-analysis normalized pseudotime")
    for s in ax.spines.values():
        s.set_linewidth(0.8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.025)
    cbar.set_label("Median pseudotime (0–1)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax2 = fig.add_subplot(gs[1, 0])
    cmat = corr.reindex(index=analysis_order, columns=analysis_order)
    im2 = ax2.imshow(cmat.values, cmap="coolwarm", vmin=-1, vmax=1, aspect="equal")
    ax2.set_title("Shared-subtype rank agreement", fontweight="bold", pad=8)
    ax2.set_xticks(range(len(short)), short, rotation=55, ha="right")
    ax2.set_yticks(range(len(short)), short)
    for i in range(len(short)):
        for j in range(len(short)):
            v = cmat.iloc[i, j]
            if np.isfinite(v):
                ax2.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                         color="white" if abs(v) > 0.55 else "black")
    for s in ax2.spines.values():
        s.set_linewidth(0.8)
    cbar2 = fig.colorbar(im2, ax=ax2, fraction=0.028, pad=0.025)
    cbar2.set_label("Spearman rho", fontsize=8)
    cbar2.ax.tick_params(labelsize=7)
    fig.text(0.01, 0.985, "E", ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(0.5, 0.012,
             "Descriptive robustness comparison only; pseudotime is an inferred transcriptional ordering, not lineage tracing.",
             ha="center", va="bottom", fontsize=7)
    fig.subplots_adjust(left=0.31, right=0.91, top=0.95, bottom=0.16)

    base = OUT / "Figure_R11E_pseudotime_consistency"
    paths = []
    for ext, kwargs in [("pdf", {}), ("svg", {}), ("png", {"dpi": 600}), ("tiff", {"dpi": 600})]:
        p = base.with_suffix(f".{ext}")
        fig.savefig(p, facecolor="white", **kwargs)
        paths.append(p)
    plt.close(fig)
    return paths


def main() -> None:
    (OUT / "individual_figures").mkdir(parents=True, exist_ok=True)
    copied = copy_source_figures()
    copied.extend(render_individual_600dpi())

    frames = []
    audit_rows = []
    for name, path in ANALYSES.items():
        x = load_analysis(name, path)
        frames.append(x)
        audit_rows.append({
            "analysis": name,
            "source_file": str(path),
            "n_cells": len(x),
            "n_sub_labels": x["sub_labels"].nunique(),
            "pseudotime_min": x["pseudotime"].min(),
            "pseudotime_max": x["pseudotime"].max(),
        })
    all_cells = pd.concat(frames, ignore_index=True)
    medians = (all_cells.groupby(["analysis", "sub_labels"], observed=True)
               .agg(n_cells=("cell_id", "size"),
                    median_pseudotime=("pseudotime", "median"),
                    median_pseudotime_01=("pseudotime_01", "median"),
                    q25_pseudotime_01=("pseudotime_01", lambda x: x.quantile(0.25)),
                    q75_pseudotime_01=("pseudotime_01", lambda x: x.quantile(0.75)))
               .reset_index())
    medians.to_csv(OUT / "subtype_median_pseudotime.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(audit_rows).to_csv(OUT / "analysis_input_audit.csv", index=False, encoding="utf-8-sig")

    piv = medians.pivot(index="sub_labels", columns="analysis", values="median_pseudotime_01")
    corr = pd.DataFrame(np.nan, index=list(ANALYSES), columns=list(ANALYSES), dtype=float)
    metric_rows = []
    for a in ANALYSES:
        for b in ANALYSES:
            if a == b:
                pair = piv[[a]].dropna()
                rho = 1.0 if len(pair) >= 1 else np.nan
            else:
                pair = pd.concat([piv[a].rename("x"), piv[b].rename("y")], axis=1).dropna()
                rho = float(spearmanr(pair["x"], pair["y"]).statistic) if len(pair) >= 3 else np.nan
            corr.loc[a, b] = rho
            metric_rows.append({"analysis_1": a, "analysis_2": b, "n_shared_sub_labels": len(pair), "spearman_rho": rho})
    pd.DataFrame(metric_rows).to_csv(OUT / "pairwise_subtype_rank_agreement.csv", index=False, encoding="utf-8-sig")
    corr.to_csv(OUT / "pairwise_subtype_rank_agreement_matrix.csv", encoding="utf-8-sig")

    generated = []
    generated.extend(make_consistency_figure(medians, corr))
    generated.append(make_multipage_pdf())
    generated.extend(make_contact_sheet())

    root_rows = []
    for set_name in ("set1_FibK14_allKC_allSAC", "set2_FibK14_allKC", "set3_allKC_allSAC"):
        p = M2 / set_name / "root_direction_validation.csv"
        z = pd.read_csv(p)
        z.insert(0, "analysis", set_name)
        root_rows.append(z)
    roots = pd.concat(root_rows, ignore_index=True)
    roots.to_csv(OUT / "monocle2_root_direction_audit.csv", index=False, encoding="utf-8-sig")
    shutil.copy2(M3 / "03_analysis" / "KC_SAC_monocle3_root_cells.tsv", OUT / "monocle3_root_cells.tsv")

    summary = {
        "output_directory": str(OUT),
        "analyses": audit_rows,
        "pairwise_correlations": metric_rows,
        "important_boundary": "Pseudotime is an inferred transcriptional ordering and does not prove migration, lineage conversion, or transdifferentiation.",
        "rna_velocity_status": "Not generated from pbmc_final.h5ad because spliced/unspliced layers are absent.",
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest_rows = []
    for path in sorted([p for p in OUT.rglob("*") if p.is_file() and p.name != "manifest_sha256.csv"]):
        manifest_rows.append({"relative_path": str(path.relative_to(OUT)), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest_rows).to_csv(OUT / "manifest_sha256.csv", index=False, encoding="utf-8-sig")

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
