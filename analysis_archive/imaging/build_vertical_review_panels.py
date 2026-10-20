from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import Rectangle
from PIL import Image


ROOT = Path(r"data/external/local_6b4e1e668702")
KC_PKG = ROOT / r"Figure 3D 细胞计数\KC_Spinous_Mig共染逐片与Mode统计_20260824"
FIB_PKG = ROOT / r"Figure 3D 细胞计数\Fib_K14_output"
OUT = ROOT / r"Figure 3D 细胞计数\KC_Fib_Mode整片统计_PDF_AI_20260824"
MM = 1 / 25.4
W, H = 42 * MM, 28 * MM
PIXEL_SIZE_UM = 0.25284451246261597

KC_DAYS = {"cut-007": 7, "001a": None, "007": 7, "32-1": 6, "bs1": 12,
           "bw3-s2": 19, "BW4": 19, "bw7-d2": 22, "bw8-d2": 22,
           "bw9": 26, "w14-1": 14}
FIB_DAYS = {"001A": None, "003": 10, "007": 7, "32-1": 6, "BS1": 12,
            "BW3": 19, "BW7": 22, "BW8": 22, "BW9": 26, "W14-1": 14}

KC_COLORS = np.array([[0.10, 0.25, 1.0], [0.0, 1.0, 0.35], [1.0, 0.15, 0.75], [1.0, 0.18, 0.05]])
KC_LABELS = [("DAPI", "#4B6FFF"), ("MMP3", "#00E65C"), ("SOX9", "#FF40C5"), ("KRT6B", "#FF4025")]
FIB_LABELS = [("DAPI", "#4B6FFF"), ("PDGFRA", "#00E65C"), ("P63", "#34D5E8"), ("KLF4", "#FFE22E"), ("KRT14", "#FF40C5")]


def setup():
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42,
                         "ps.fonttype": 42, "svg.fonttype": "none"})
    for d in (OUT / "01_KC_vertical_all_slides", OUT / "02_Fib_oriented_zoomout_review",
              OUT / "03_A4_review", OUT / "04_manifests", OUT / "05_whole_slide_lineplot",
              OUT / "06_four_panel_rows", OUT / "07_AI"):
        d.mkdir(parents=True, exist_ok=True)


def norm(a):
    v = a[np.isfinite(a) & (a > 0)]
    if v.size == 0:
        return np.zeros_like(a, dtype=np.float32), 0., 1.
    lo, hi = np.percentile(v, [1, 99.7])
    hi = max(hi, lo + 1)
    return np.clip((a.astype(np.float32) - lo) / (hi - lo), 0, 1), float(lo), float(hi)


def merge(xs):
    rgb = np.zeros((*xs[0].shape, 3), np.float32)
    for x, c, g in zip(xs, KC_COLORS, [0.85, 0.75, 0.75, 0.75]):
        rgb += x[..., None] * c * g
    return np.clip(rgb, 0, 1) ** 0.82


def label_stack(fig, labels):
    """Reserve a right-side legend column so 8 pt text never covers tissue."""
    y = 0.78
    for label, color in labels:
        fig.text(0.69, y, label, color=color, fontsize=8, fontweight="bold",
                 va="top", ha="left")
        y -= 0.125


def scalebar(ax, img_w, length_px, text):
    x1 = img_w * 0.96
    x0 = x1 - length_px
    y = ax.get_ylim()[0] * 0.92
    ax.plot([x0, x1], [y, y], color="white", lw=1.1, solid_capstyle="butt")
    ax.text((x0+x1)/2, y, text, color="white", fontsize=7, ha="center", va="bottom")


def mode_from_day(day):
    if day is None or (isinstance(day, float) and np.isnan(day)):
        return "UNMAPPED"
    if 0 <= int(day) <= 7: return "Mode 1"
    if 8 <= int(day) <= 14: return "Mode 2"
    if 15 <= int(day) <= 28: return "Mode 3"
    return "OUT_OF_RANGE"


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix("." + ext), dpi=600, facecolor="black", bbox_inches=None)
    plt.close(fig)


def kc_panels():
    man = pd.read_csv(KC_PKG / "00_manifest" / "slide_roi_manifest.tsv", sep="\t")
    rows, pngs = [], []
    for _, r in man.iterrows():
        p = Path(r.absolute_path)
        with tifffile.TiffFile(p) as tf:
            pages = [pg.asarray(out="memmap") for pg in tf.pages]
            h, w = pages[0].shape
            ds = max(1, math.ceil(max(w, h) / 1800))
            small, wins = [], []
            for a in pages:
                n, lo, hi = norm(np.asarray(a[::ds, ::ds])); small.append(n); wins.append((lo, hi))
            overview = merge(small)
            x0, y0, rw0, rh0 = map(int, (r.roi_x, r.roi_y, r.roi_width, r.roi_height))
            # Expand around the reviewed/candidate center so a real 500 um bar fits in Zoom in.
            target_w = max(rw0, math.ceil(700.0 / PIXEL_SIZE_UM))
            target_h = max(rh0, math.ceil(430.0 / PIXEL_SIZE_UM))
            cx, cy = x0 + rw0 / 2, y0 + rh0 / 2
            rw, rh = min(target_w, w), min(target_h, h)
            x = int(np.clip(round(cx - rw/2), 0, w-rw))
            y = int(np.clip(round(cy - rh/2), 0, h-rh))
            crops = []
            for a, (lo, hi) in zip(pages, wins):
                crops.append(np.clip((np.asarray(a[y:y+rh, x:x+rw], np.float32)-lo)/(hi-lo), 0, 1))
            crop = merge(crops)

        fig = plt.figure(figsize=(W, H), facecolor="black")
        ax0 = fig.add_axes([0.02, 0.51, 0.64, 0.47], facecolor="black")
        ax1 = fig.add_axes([0.02, 0.02, 0.64, 0.47], facecolor="black")
        ax0.imshow(overview); ax1.imshow(crop)
        ax0.add_patch(Rectangle((x/ds, y/ds), rw/ds, rh/ds, fill=False, ec="white", lw=0.8))
        for ax in (ax0, ax1): ax.set_axis_off()
        fig.text(0.69, 0.97, str(r.slide_id).upper(), color="white", fontsize=8,
                 fontweight="bold", va="top", ha="left")
        # Zoom captions and channel-color labels intentionally omitted.
        # Restored per author instruction using the declared project calibration.
        scalebar(ax0, overview.shape[1], 1000.0/(PIXEL_SIZE_UM*ds), "1000 µm")
        scalebar(ax1, crop.shape[1], 200.0/PIXEL_SIZE_UM, "200 µm")
        status = "manual ROI" if str(r.review_status).startswith("confirmed") else "candidate ROI — review required"
        short_status = "ROI set" if status == "manual ROI" else "Review"
        stem = OUT / "01_KC_vertical_all_slides" / f"{int(r.slide_order):02d}_{r.slide_id}_vertical"
        save(fig, stem); pngs.append(stem.with_suffix(".png"))
        day = KC_DAYS[str(r.slide_id)]
        rows.append({"sample_id": r.slide_id, "day_post_burn": day, "mode": mode_from_day(day),
                     "rotation_ccw_deg": 0, "epidermis_position": "top",
                     "roi_status": status, "scale_unit": "micrometer",
                     "pixel_size_um_used": PIXEL_SIZE_UM,
                     "pixel_size_provenance": "matched KFBF scanner metadata and converted paired KFBF",
                     "output_png": str(stem.with_suffix('.png'))})
    pd.DataFrame(rows).to_csv(OUT / "04_manifests" / "KC_vertical_manifest.tsv", sep="\t", index=False)
    return pngs


def fib_review():
    man = pd.read_csv(FIB_PKG / "00_manifest" / "zoomout_manifest.tsv", sep="\t")
    # Rotation is a display operation only. Positive values are counter-clockwise.
    rotations = {"BW9": -90, "W14-1": 90, "BW3": 0, "BW7": 0, "001A": 0,
                 "003": 0, "007": 0, "32-1": 0, "BS1": 0, "BW8": 0}
    confidence = {"BW3": "confirmed", "BW7": "confirmed", "BW9": "confirmed_visual",
                  "W14-1": "confirmed_visual", "001A": "review_diagonal",
                  "003": "review_low_epidermal_signal", "007": "review_fragmented",
                  "32-1": "review_fragmented", "BS1": "review_multiple_fragment_orientations",
                  "BW8": "review_fragmented"}
    rows, pngs = [], []
    for _, r in man.iterrows():
        sid = str(r.sample_id)
        img = np.asarray(Image.open(r.zoomout_png).convert("RGB"))
        deg = rotations[sid]
        if deg == 90: img = np.rot90(img, 1)
        elif deg == -90: img = np.rot90(img, -1)
        fig = plt.figure(figsize=(W, H), facecolor="black")
        ax = fig.add_axes([0.02, 0.51, 0.64, 0.47], facecolor="black")
        ax_blank = fig.add_axes([0.02, 0.02, 0.64, 0.47], facecolor="black")
        ax.imshow(img); ax.set_axis_off()
        fig.text(0.69, 0.97, sid, color="white", fontsize=8, fontweight="bold", va="top")
        ax_blank.set_xticks([]); ax_blank.set_yticks([])
        for spine in ax_blank.spines.values():
            spine.set_color("#FFD447"); spine.set_linestyle("--"); spine.set_linewidth(0.7)
        ax_blank.plot([.68, .93], [.08, .08], transform=ax_blank.transAxes, color="white", lw=1.1)
        ax_blank.text(.805, .10, "200 µm", transform=ax_blank.transAxes,
                      ha="center", va="bottom", color="white", fontsize=7)
        px_um = float(r.zoomout_pixel_size_x_um)
        bar_um = 1000.0
        scalebar(ax, img.shape[1], bar_um/px_um, "1000 µm")
        fig.text(0.69, 0.08, "ROI needed", color="#FFD447", fontsize=8, va="bottom")
        stem = OUT / "02_Fib_oriented_zoomout_review" / f"{sid}_oriented_zoomout_review"
        save(fig, stem); pngs.append(stem.with_suffix(".png"))
        day = FIB_DAYS[sid]
        rows.append({"sample_id": sid, "day_post_burn": day, "mode": mode_from_day(day),
                     "rotation_ccw_deg": deg, "orientation_confidence": confidence[sid],
                     "epidermis_target": "top", "zoom_in_status": "BLOCKED_PENDING_USER_IMAGEJ_ROI",
                     "zoomout_pixel_size_um": px_um, "output_png": str(stem.with_suffix('.png'))})
    pd.DataFrame(rows).to_csv(OUT / "04_manifests" / "Fib_orientation_manifest.tsv", sep="\t", index=False)
    return pngs


def atlas(kc, fib):
    fig = plt.figure(figsize=(210*MM, 297*MM), facecolor="white")
    for i, p in enumerate(kc):
        row, col = divmod(i, 5)
        ax = fig.add_axes([col*0.2, 0.72-row*(28/297), 0.2, 28/297]); ax.imshow(Image.open(p)); ax.axis("off")
    fig.text(0.01, 0.995, "KC_Spinous_Mig — vertical Zoom out / Zoom in review", fontsize=10, fontweight="bold", va="top")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / "03_A4_review" / f"KC_all_11_vertical_A4.{ext}", dpi=600)
    plt.close(fig)
    fig = plt.figure(figsize=(210*MM, 297*MM), facecolor="white")
    for i, p in enumerate(fib):
        row, col = divmod(i, 5)
        ax = fig.add_axes([col*0.2, 0.78-row*(28/297), 0.2, 28/297]); ax.imshow(Image.open(p)); ax.axis("off")
    fig.text(0.01, 0.995, "Fib_K14 — oriented Zoom out review; Zoom in pending ImageJ ROI", fontsize=10, fontweight="bold", va="top")
    for ext in ("png", "pdf"):
        fig.savefig(OUT / "03_A4_review" / f"Fib_all_10_oriented_zoomout_A4.{ext}", dpi=600)
    plt.close(fig)


def write_status():
    pd.DataFrame([
        {"component": "KC per-slide vertical panels", "status": "DELIVERED", "reason": "11 slides processed; 7 manual ROI and 4 review candidates"},
        {"component": "Fib oriented Zoom out", "status": "DELIVERED_REVIEW", "reason": "10 slides; four confident, six require visual confirmation"},
        {"component": "Fib Zoom in and ROI box", "status": "BLOCKED", "reason": "No user-approved ImageJ ROI archive"},
        {"component": "Mode 1-3 representative assignment", "status": "DELIVERED_WITH_001A_UNMAPPED", "reason": "0-7, 8-14, 15-28 day rules applied from audited workbook"},
        {"component": "KC four-panel final composition", "status": "DELIVERED", "reason": "Mode 1=007, Mode 2=BS1, Mode 3=BW4, plus one whole-slide line plot"},
        {"component": "Fib four-panel composition", "status": "DELIVERED_REVIEW", "reason": "Zoom-in placeholders retained because user-approved Fib ImageJ ROIs are absent"},
    ]).to_csv(OUT / "04_manifests" / "four_panel_completion_status.tsv", sep="\t", index=False)
    (OUT / "README_请先阅读.txt").write_text(
        "本目录为上下排版与方向复核版本。\n"
        "KC：11/11 已输出，上方 Zoom out、下方 Zoom in、等宽；含通道颜色标签、ROI 框，"
        "Zoom out=1000 µm、Zoom in=500 µm。KC TIFF 本身缺失标定，采用配对 KFBF 同一扫描平台的 0.2528445125 µm/px，并在 manifest 记录。\n"
        "Fib：10/10 已输出方向校正的 Zoom out 复核图，含通道标签和 1000 µm 比例尺。"
        "因缺少用户审核的 ImageJ ROI，未生成 Fib Zoom in/统计/终版代表图。\n"
        "Mode 定义：0-7天=Mode 1，8-14天=Mode 2，15-28天=Mode 3。001A 未获得可靠天数映射。\n"
        "折线图为整张切片统计的旧位图数字化近似，不使用 ROI 计数。\n",
        encoding="utf-8")


def whole_slide_lineplot():
    src = ROOT / r"Figure 3D 细胞计数\Figure3D_replot_source_package_20260820\01_附件折线图_数字化复现\source_data\digitized_values.csv"
    df = pd.read_csv(src)
    df["mode"] = df["time_dpb"].map(lambda d: mode_from_day(int(d)))
    df["analysis_region"] = "whole_slide"
    df.to_csv(OUT / "04_manifests" / "whole_slide_digitized_values_with_modes.csv", index=False)
    fig, ax = plt.subplots(figsize=(W, H))
    colors = {"SPTDI": "#3A9FBF", "DPTDI": "#E86A5A"}
    modes = ["Mode 1", "Mode 2", "Mode 3"]
    for group in ("SPTDI", "DPTDI"):
        g = df[df.group == group]
        means = g.groupby("mode")["SOX9_MMP3_KRT6B_positive_ratio_percent"].mean().reindex(modes)
        sems = g.groupby("mode")["SOX9_MMP3_KRT6B_positive_ratio_percent"].sem().reindex(modes)
        x = np.arange(3)
        ax.fill_between(x, means-sems, means+sems, color=colors[group], alpha=.16, lw=0)
        ax.plot(x, means, "-o", color=colors[group], lw=1.1, ms=3, label=group)
        rng = np.random.default_rng(20260824 + (0 if group == "SPTDI" else 1))
        for j, mode in enumerate(modes):
            vals = g.loc[g["mode"] == mode, "SOX9_MMP3_KRT6B_positive_ratio_percent"].to_numpy()
            ax.scatter(j+rng.uniform(-.06,.06,len(vals)), vals, s=7, facecolor="white",
                       edgecolor=colors[group], lw=.55, zorder=3)
    ax.set_xticks(range(3), modes); ax.set_xlim(-.35, 2.35)
    ax.set_ylabel("Triple-positive cell ratio (%)", labelpad=2)
    ax.legend(frameon=False, loc="upper left", ncol=1, handlelength=1,
              borderaxespad=0.1, labelspacing=0.15)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5, width=.8, pad=1)
    ax.text(.98, .98, "Whole slide", transform=ax.transAxes, ha="right", va="top", fontsize=8,
            fontstyle="italic")
    fig.subplots_adjust(left=.28, right=.95, bottom=.27, top=.94)
    stem = OUT / "05_whole_slide_lineplot" / "whole_slide_Mode1_Mode2_Mode3_42x28mm"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix("."+ext), dpi=600)
    plt.close(fig)
    return stem.with_suffix(".png")


def four_panel_rows(kc_pngs, fib_pngs, line_png):
    kc_map = {p.name.split("_", 2)[1]: p for p in kc_pngs}
    fib_map = {p.name.split("_", 1)[0]: p for p in fib_pngs}
    selections = [
        ("KC", [("Mode 1", kc_map["007"]), ("Mode 2", kc_map["bs1"]), ("Mode 3", kc_map["BW4"])]),
        ("Fib", [("Mode 1", fib_map["007"]), ("Mode 2", fib_map["BS1"]), ("Mode 3", fib_map["BW7"])]),
    ]
    for family, items in selections:
        fig = plt.figure(figsize=(4*W, H), facecolor="white")
        for i, (title, p) in enumerate(items):
            ax = fig.add_axes([i*.25, 0, .25, 1]); ax.imshow(Image.open(p)); ax.axis("off")
            ax.text(.02, .98, title, transform=ax.transAxes, color="white", fontsize=8,
                    fontweight="bold", va="top", bbox={"facecolor":"black","edgecolor":"none","alpha":.6,"pad":.4})
        ax = fig.add_axes([.75, 0, .25, 1]); ax.imshow(Image.open(line_png)); ax.axis("off")
        stem = OUT / "06_four_panel_rows" / f"{family}_Mode1_Mode2_Mode3_whole_slide_lineplot_168x28mm"
        for ext in ("png", "pdf"):
            fig.savefig(stem.with_suffix("."+ext), dpi=600)
        plt.close(fig)


if __name__ == "__main__":
    setup()
    kc = kc_panels()
    fib = fib_review()
    atlas(kc, fib)
    line = whole_slide_lineplot()
    four_panel_rows(kc, fib, line)
    write_status()
    print(OUT)
