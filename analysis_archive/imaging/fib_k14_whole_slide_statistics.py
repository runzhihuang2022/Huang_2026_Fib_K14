from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from scipy import ndimage as ndi
from scipy.stats import kruskal, mannwhitneyu
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import regionprops
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed
from statsmodels.stats.multitest import multipletests


ROOT = Path(r"data/external/local_93c497bb3a06")
INPUT = ROOT / "01_ome_tiff"
OUT = ROOT / "07_whole_slide_statistics"
PIXEL_UM = 0.25284451246261597 * 4  # pyramid level 2
LEVEL = 2
TILE = 1024
MARGIN = 32
MIN_AREA_UM2 = 20.0
MAX_AREA_UM2 = 300.0
MARKER_RADIUS_UM = 5.0
BG_SIGMA_UM = 15.0

META = {
    "003": (10, "Mode 2"), "007": (7, "Mode 1"), "32-1": (6, "Mode 1"),
    "BS1": (12, "Mode 2"), "BW3": (19, "Mode 3"), "BW7": (22, "Mode 3"),
    "BW8": (22, "Mode 3"), "BW9": (26, "Mode 3"), "W14-1": (14, "Mode 2"),
    "001A": (None, "UNMAPPED"),
}


def sample_id(path: Path) -> str:
    s = path.name.upper()
    for key in ["W14-1", "32-1", "001A", "003", "007", "BS1", "BW3", "BW7", "BW8", "BW9"]:
        if key in s:
            return key
    return re.sub(r"\s+PDGF.*", "", path.stem, flags=re.I)


def open_level(path: Path):
    store = tifffile.imread(path, series=0, level=LEVEL, aszarr=True)
    return store, zarr.open(store, mode="r")


def robust_otsu(values: np.ndarray) -> float:
    x = values[np.isfinite(values)]
    if x.size == 0:
        return float("nan")
    lo, hi = np.percentile(x, [0.5, 99.5])
    x = x[(x >= lo) & (x <= hi)]
    return float(threshold_otsu(x)) if np.unique(x).size > 1 else float(x[0])


def segment_tile(dapi: np.ndarray) -> np.ndarray:
    sm = gaussian(dapi.astype(np.float32), sigma=0.8, preserve_range=True)
    positive = sm[sm > 0]
    if positive.size < 100:
        return np.zeros(sm.shape, np.int32)
    thr = threshold_otsu(positive)
    mask = remove_small_objects(sm > thr, min_size=max(3, int(MIN_AREA_UM2 / PIXEL_UM**2 / 3)))
    dist = ndi.distance_transform_edt(mask)
    coords = peak_local_max(dist, min_distance=max(2, int(2.0 / PIXEL_UM)), labels=mask)
    markers = np.zeros(mask.shape, np.int32)
    if coords.size:
        markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    else:
        markers, _ = ndi.label(mask)
    return watershed(-dist, markers, mask=mask)


def p90_disk(img: np.ndarray, cy: float, cx: float, radius: int) -> float:
    y0, y1 = max(0, int(cy) - radius), min(img.shape[0], int(cy) + radius + 1)
    x0, x1 = max(0, int(cx) - radius), min(img.shape[1], int(cx) + radius + 1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    m = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
    v = img[y0:y1, x0:x1][m]
    return float(np.percentile(v, 90)) if v.size else float("nan")


def analyze_slide(path: Path, cell_dir: Path) -> dict:
    sid = sample_id(path)
    store, arr = open_level(path)
    _, h, w = arr.shape
    min_px = MIN_AREA_UM2 / PIXEL_UM**2
    max_px = MAX_AREA_UM2 / PIXEL_UM**2
    radius = max(1, round(MARKER_RADIUS_UM / PIXEL_UM))
    bg_sigma = BG_SIGMA_UM / PIXEL_UM
    out_path = cell_dir / f"{sid}_whole_slide_cells.tsv.gz"
    n = 0
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["sample", "cell_id", "centroid_x_level2_px", "centroid_y_level2_px", "area_um2", "PDGFRA_p90_bgsub", "KRT14_p90_bgsub"])
        for y0 in range(0, h, TILE):
            for x0 in range(0, w, TILE):
                cy0, cy1 = y0, min(h, y0 + TILE)
                cx0, cx1 = x0, min(w, x0 + TILE)
                ry0, ry1 = max(0, cy0 - MARGIN), min(h, cy1 + MARGIN)
                rx0, rx1 = max(0, cx0 - MARGIN), min(w, cx1 + MARGIN)
                block = np.asarray(arr[[0, 1, 4], ry0:ry1, rx0:rx1])
                labels = segment_tile(block[0])
                pdgfra = block[1].astype(np.float32)
                krt14 = block[2].astype(np.float32)
                pdgfra = np.maximum(pdgfra - gaussian(pdgfra, sigma=bg_sigma, preserve_range=True), 0)
                krt14 = np.maximum(krt14 - gaussian(krt14, sigma=bg_sigma, preserve_range=True), 0)
                for r in regionprops(labels):
                    if not (min_px <= r.area <= max_px):
                        continue
                    yy, xx = r.centroid
                    gy, gx = yy + ry0, xx + rx0
                    if not (cy0 <= gy < cy1 and cx0 <= gx < cx1):
                        continue
                    n += 1
                    wr.writerow([sid, n, f"{gx:.3f}", f"{gy:.3f}", f"{r.area * PIXEL_UM**2:.3f}",
                                 f"{p90_disk(pdgfra, yy, xx, radius):.6g}", f"{p90_disk(krt14, yy, xx, radius):.6g}"])
    store.close()
    return {"sample": sid, "source_file": path.name, "width_level2_px": w, "height_level2_px": h,
            "pixel_size_um": PIXEL_UM, "nuclei_count": n, "cell_table": str(out_path)}


def load_marker_tables(cell_dir: Path) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep="\t", dtype={"sample": str}) for p in sorted(cell_dir.glob("*_whole_slide_cells.tsv.gz"))]
    return pd.concat(frames, ignore_index=True)


def stats_tables(per_slide: pd.DataFrame, stat_dir: Path):
    mapped = per_slide[per_slide["mode"].str.startswith("Mode")].copy()
    groups = [mapped.loc[mapped["mode"] == m, "double_positive_ratio_pct"].to_numpy() for m in ["Mode 1", "Mode 2", "Mode 3"]]
    tests = []
    if all(len(x) > 0 for x in groups):
        h, p = kruskal(*groups)
        tests.append({"test_family": "omnibus", "contrast": "Mode 1 vs Mode 2 vs Mode 3", "test": "Kruskal-Wallis", "statistic": h, "p_raw": p, "p_BH": p})
    pair_rows = []
    for i, a in enumerate(["Mode 1", "Mode 2", "Mode 3"]):
        for b in ["Mode 1", "Mode 2", "Mode 3"][i + 1:]:
            x = mapped.loc[mapped["mode"] == a, "double_positive_ratio_pct"].to_numpy()
            y = mapped.loc[mapped["mode"] == b, "double_positive_ratio_pct"].to_numpy()
            if len(x) and len(y):
                u, p = mannwhitneyu(x, y, alternative="two-sided", method="auto")
                pair_rows.append({"test_family": "pairwise", "contrast": f"{a} vs {b}", "test": "Mann-Whitney U", "statistic": u, "p_raw": p})
    if pair_rows:
        adj = multipletests([x["p_raw"] for x in pair_rows], method="fdr_bh")[1]
        for row, q in zip(pair_rows, adj):
            row["p_BH"] = q
    pd.DataFrame(tests + pair_rows).to_csv(stat_dir / "Fib_K14_statistical_tests.tsv", sep="\t", index=False)
    summary = mapped.groupby("mode", observed=True)["double_positive_ratio_pct"].agg(["count", "mean", "std", "median", "min", "max"]).reset_index()
    summary.to_csv(stat_dir / "Fib_K14_mode_summary.tsv", sep="\t", index=False)


def plot_line(per_slide: pd.DataFrame, fig_dir: Path):
    order = ["Mode 1", "Mode 2", "Mode 3"]
    df = per_slide[per_slide["mode"].isin(order)].copy()
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(42/25.4, 28/25.4))
    rng = np.random.default_rng(20260825)
    for i, mode in enumerate(order):
        y = df.loc[df["mode"] == mode, "double_positive_ratio_pct"].to_numpy()
        ax.scatter(i + rng.uniform(-.055, .055, len(y)), y, s=16, facecolors="white", edgecolors="#D95F4C", linewidths=.8, zorder=3)
    g = df.groupby("mode", observed=True)["double_positive_ratio_pct"].agg(["mean", "sem"]).reindex(order)
    ax.plot(range(3), g["mean"], color="#D95F4C", marker="o", markersize=3.5, linewidth=1.2)
    ax.fill_between(range(3), (g["mean"] - g["sem"]).to_numpy(), (g["mean"] + g["sem"]).to_numpy(), color="#D95F4C", alpha=.13, linewidth=0)
    ax.set_xticks(range(3), ["Mode\n1", "Mode\n2", "Mode\n3"])
    ax.set_ylabel("PDGFRA+ KRT14+\ncell ratio (%)", fontsize=8, labelpad=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=.8, length=3, labelsize=8)
    fig.subplots_adjust(left=.34, right=.96, bottom=.27, top=.95)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"Fib_K14_whole_slide_lineplot.{ext}", dpi=600 if ext == "png" else None, transparent=False)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*", help="Optional sample IDs for testing/resume")
    args = ap.parse_args()
    cell_dir = OUT / "01_per_cell"
    stat_dir = OUT / "02_statistics"
    fig_dir = OUT / "03_figures"
    audit_dir = OUT / "04_audit"
    for d in [cell_dir, stat_dir, fig_dir, audit_dir]: d.mkdir(parents=True, exist_ok=True)
    paths = sorted(INPUT.glob("*.ome.tiff"))
    if args.samples:
        paths = [p for p in paths if sample_id(p) in set(args.samples)]
    records = []
    for p in paths:
        sid = sample_id(p)
        table = cell_dir / f"{sid}_whole_slide_cells.tsv.gz"
        if table.exists():
            n = sum(1 for _ in gzip.open(table, "rt", encoding="utf-8")) - 1
            records.append({"sample": sid, "source_file": p.name, "pixel_size_um": PIXEL_UM, "nuclei_count": n, "cell_table": str(table)})
        else:
            print(f"ANALYZE {sid}", flush=True)
            records.append(analyze_slide(p, cell_dir))
    all_cells = load_marker_tables(cell_dir)
    pdgfra_thr = robust_otsu(all_cells["PDGFRA_p90_bgsub"].to_numpy())
    krt14_thr = robust_otsu(all_cells["KRT14_p90_bgsub"].to_numpy())
    all_cells["PDGFRA_positive"] = all_cells["PDGFRA_p90_bgsub"] > pdgfra_thr
    all_cells["KRT14_positive"] = all_cells["KRT14_p90_bgsub"] > krt14_thr
    all_cells["double_positive"] = all_cells["PDGFRA_positive"] & all_cells["KRT14_positive"]
    rows = []
    for sid, x in all_cells.groupby("sample", sort=True):
        day, mode = META.get(sid, (None, "UNMAPPED"))
        rows.append({"sample": sid, "day_dpb": day, "mode": mode, "nuclei_count": len(x),
                     "PDGFRA_positive_n": int(x["PDGFRA_positive"].sum()), "KRT14_positive_n": int(x["KRT14_positive"].sum()),
                     "double_positive_n": int(x["double_positive"].sum()),
                     "double_positive_ratio_pct": float(x["double_positive"].mean() * 100)})
    per_slide = pd.DataFrame(rows)
    per_slide.to_csv(stat_dir / "Fib_K14_whole_slide_per_slide.tsv", sep="\t", index=False)
    stats_tables(per_slide, stat_dir)
    plot_line(per_slide, fig_dir)
    config = {"analysis_scope": "whole slide; ROI not used", "pyramid_level": LEVEL, "pixel_size_um": PIXEL_UM,
              "channels": {"DAPI": 0, "PDGFRA": 1, "KRT14": 4}, "nuclear_area_um2": [MIN_AREA_UM2, MAX_AREA_UM2],
              "marker_measure": f"90th percentile in {MARKER_RADIUS_UM} um centroid disk after Gaussian background subtraction",
              "background_sigma_um": BG_SIGMA_UM, "PDGFRA_locked_threshold": pdgfra_thr, "KRT14_locked_threshold": krt14_thr,
              "threshold_scope": "one pooled automatic Otsu threshold per marker across all detected nuclei and slides",
              "inference_unit": "slide", "mode_definition": {"Mode 1": "0-7 dpb", "Mode 2": "8-14 dpb", "Mode 3": "15-28 dpb"},
              "status": "exploratory; automated segmentation/thresholds require image-review validation"}
    (audit_dir / "analysis_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT / "README_统计说明.md").write_text(
        "# Fib_K14整张切片统计\n\n统计对象为整张切片，不使用Zoom-in ROI。DAPI分割细胞核后，在核中心5 µm圆盘内测量背景扣除后的PDGFRA和KRT14信号；每个marker在全队列中只使用一个锁定Otsu阈值。统计重复单位为切片。因各Mode切片数较少，检验结果仅作探索性结果，不能把细胞数当作生物学重复。001A缺少可靠时间映射，保留整片计数但排除Mode检验。\n",
        encoding="utf-8")
    print(per_slide.to_string(index=False))
    print(json.dumps({"PDGFRA": pdgfra_thr, "KRT14": krt14_thr}, ensure_ascii=False))


if __name__ == "__main__":
    main()
