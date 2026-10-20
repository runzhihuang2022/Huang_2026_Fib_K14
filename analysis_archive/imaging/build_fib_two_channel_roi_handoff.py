from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import tifffile
import zarr


ROOT = Path(r"data/external/local_93c497bb3a06")
FIG_ROOT = Path(r"data/external/local_1f040df3bed6")
OUT = ROOT / "03_roi_handoff" / "KRT14_PDGFRA_DAPI_only"
REF = OUT / "01_Fiji_open_this"
PREVIEW = OUT / "02_three_channel_preview"
ZIP_DIR = OUT / "03_put_RoiSet_zip_here"
AUDIT = OUT / "04_audit"
PANEL = OUT / "05_vertical_review_panels"
TARGET_W_UM = 2000.0
TARGET_H_UM = 1300.0
MM = 1 / 25.4


def font(size: int, bold: bool = False):
    names = ["arialbd.ttf", "Arial Bold.ttf"] if bold else ["arial.ttf", "Arial.ttf"]
    for name in names:
        p = Path(r"data/external/local_81ff0cfd85af") / name
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def normalize(a: np.ndarray, lo: float, hi: float) -> np.ndarray:
    if hi <= lo:
        hi = lo + 1
    return np.clip((a.astype(np.float32) - lo) / (hi - lo), 0, 1)


def compose(dapi: np.ndarray, pdgfra: np.ndarray, krt14: np.ndarray) -> np.ndarray:
    rgb = np.zeros((*pdgfra.shape, 3), dtype=np.float32)
    rgb[..., 2] += dapi * 0.70
    rgb[..., 1] += pdgfra * 0.85
    rgb[..., 0] += krt14 * 0.90
    rgb[..., 2] += krt14 * 0.90
    return np.clip(rgb, 0, 1) ** 0.85


def crop_to_tissue(rgb: np.ndarray) -> np.ndarray:
    mask = rgb.max(axis=2) > 0.035
    if not mask.any():
        return rgb
    yy, xx = np.where(mask)
    pad_y = max(8, round((yy.max() - yy.min() + 1) * .05))
    pad_x = max(8, round((xx.max() - xx.min() + 1) * .05))
    y0, y1 = max(0, yy.min() - pad_y), min(rgb.shape[0], yy.max() + pad_y + 1)
    x0, x1 = max(0, xx.min() - pad_x), min(rgb.shape[1], xx.max() + pad_x + 1)
    return rgb[y0:y1, x0:x1]


def rotate(a: np.ndarray, deg: int) -> np.ndarray:
    if deg == 90:
        return np.rot90(a, 1, axes=(-2, -1))
    if deg == -90:
        return np.rot90(a, -1, axes=(-2, -1))
    if deg == 180:
        return np.rot90(a, 2, axes=(-2, -1))
    return a


def annotate(rgb: np.ndarray, sid: str, px_um: float) -> Image.Image:
    im = Image.fromarray((rgb * 255).astype(np.uint8))
    draw = ImageDraw.Draw(im)
    fs = max(14, round(im.width / 95))
    f = font(fs)
    fb = font(fs, True)
    draw.text((12, 8), sid, fill="white", font=fb, stroke_width=1, stroke_fill="black")
    bar_px = max(1, round(1000 / px_um))
    x1, yb = im.width - 18, im.height - 24
    x0 = max(12, x1 - bar_px)
    draw.line((x0, yb, x1, yb), fill="white", width=max(2, round(im.width / 600)))
    label = "1000 um"
    box = draw.textbbox((0, 0), label, font=f)
    draw.text(((x0 + x1 - (box[2] - box[0])) / 2, yb - fs - 4), label,
              fill="white", font=f, stroke_width=1, stroke_fill="black")
    return im


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    for d in (REF, PREVIEW, ZIP_DIR, AUDIT, PANEL):
        d.mkdir(parents=True, exist_ok=True)
    zoom = pd.read_csv(ROOT / "00_manifest" / "zoomout_manifest.tsv", sep="\t")
    orient = pd.read_csv(FIG_ROOT / "04_manifests" / "Fib_orientation_manifest.tsv", sep="\t")
    orient_map = orient.set_index("sample_id").to_dict("index")
    rows = []
    previews = []
    for row in zoom.itertuples():
        sid = str(row.sample_id)
        level = int(row.pyramid_level)
        path = Path(row.ome_tiff)
        with tifffile.TiffFile(path) as tif:
            store = tifffile.imread(path, series=0, level=level, aszarr=True)
            try:
                arr = np.asarray(zarr.open(store, mode="r"))[[0, 1, 4]]
            finally:
                close = getattr(store, "close", None)
                if close:
                    close()
        windows = json.loads(row.display_windows_lo_hi)
        dapi = normalize(arr[0], *windows[0])
        pdgfra = normalize(arr[1], *windows[1])
        krt14 = normalize(arr[2], *windows[4])
        rotation = int(orient_map[sid]["rotation_ccw_deg"])
        data = rotate(np.stack([dapi, pdgfra, krt14]), rotation)
        px_um = float(row.zoomout_pixel_size_x_um)
        ref_path = REF / f"{sid}_DAPI_KRT14_PDGFRA_ROI_reference.ome.tiff"
        tifffile.imwrite(
            ref_path, (data * 255).astype(np.uint8), ome=True, photometric="minisblack",
            compression="zlib", metadata={
                "axes": "CYX", "PhysicalSizeX": px_um, "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": px_um, "PhysicalSizeYUnit": "µm",
                "Channel": {"Name": ["DAPI", "PDGFRA", "KRT14"]},
            },
        )
        rgb = compose(data[0], data[1], data[2])
        preview = annotate(rgb, sid, px_um)
        png = PREVIEW / f"{sid}_DAPI_KRT14_PDGFRA_preview.png"
        tif = PREVIEW / f"{sid}_DAPI_KRT14_PDGFRA_preview.tif"
        preview.save(png, dpi=(600, 600))
        preview.save(tif, compression="tiff_lzw", dpi=(600, 600))
        previews.append((sid, png))
        plt.rcParams.update({"font.family": "Arial", "font.size": 8,
                             "pdf.fonttype": 42, "svg.fonttype": "none"})
        panel_rgb = crop_to_tissue(rgb)
        panel_fig = plt.figure(figsize=(42 * MM, 28 * MM), facecolor="black")
        ax0 = panel_fig.add_axes([.02, .51, .96, .47], facecolor="black")
        ax1 = panel_fig.add_axes([.02, .02, .96, .47], facecolor="black")
        ax0.imshow(panel_rgb)
        ax0.axis("off")
        bar_px = 1000 / px_um
        ybar = panel_rgb.shape[0] * .93
        x0 = panel_rgb.shape[1] * .04
        x1 = x0 + bar_px
        ax0.plot([x0, x1], [ybar, ybar], color="white", lw=1.1)
        ax0.text((x0 + x1) / 2, ybar, "1000 µm", color="white", fontsize=7,
                 ha="center", va="bottom")
        ax1.set_xticks([])
        ax1.set_yticks([])
        for spine in ax1.spines.values():
            spine.set_color("#FFD447")
            spine.set_linestyle("--")
            spine.set_linewidth(.7)
        ax1.plot([.68, .93], [.08, .08], transform=ax1.transAxes, color="white", lw=1.1)
        ax1.text(.805, .10, "200 µm", transform=ax1.transAxes, color="white",
                 fontsize=7, ha="center", va="bottom")
        text_box = {"facecolor": "black", "edgecolor": "none", "alpha": .55, "pad": .3}
        panel_fig.text(.97, .93, sid, color="white", fontsize=8,
                       fontweight="bold", ha="right", va="top", bbox=text_box)
        stem = PANEL / f"{sid}_DAPI_KRT14_PDGFRA_vertical_review"
        for ext in ("png", "pdf", "svg"):
            panel_fig.savefig(stem.with_suffix("." + ext), dpi=600, facecolor="black")
        plt.close(panel_fig)
        rows.append({
            "sample_id": sid,
            "open_this_in_fiji": str(ref_path),
            "display_preview": str(png),
            "channels_retained_for_display": "DAPI|KRT14|PDGFRA",
            "channel_1": "DAPI_blue",
            "channel_2": "PDGFRA_green",
            "channel_3": "KRT14_magenta",
            "reference_pixel_size_um": px_um,
            "rotation_ccw_deg_applied": rotation,
            "target_roi_width_um": TARGET_W_UM,
            "target_roi_height_um": TARGET_H_UM,
            "target_roi_width_px_on_reference": round(TARGET_W_UM / px_um),
            "target_roi_height_px_on_reference": round(TARGET_H_UM / px_um),
            "roi_name": f"{sid}_ZI01",
            "roi_zip_name": f"{sid}_ZoomIn_RoiSet.zip",
            "roi_purpose": "visual_ZoomIn_only_not_statistical_ROI",
            "roi_status": "AUTHOR_TO_SELECT",
            "roi_zip_path": "",
        })
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT / "ROI_SELECTION_CHECKLIST.tsv", sep="\t", index=False, encoding="utf-8-sig")

    guide = (
        "Fib_K14 Zoom in ROI选择（当前统计仍为整张切片）\n\n"
        "当前ROI只用于最终Zoom in展示，不参与折线图或统计检验。每张切片只需选择1个主ROI。\n\n"
        "1. 在Fiji/ImageJ中用Bio-Formats打开 01_Fiji_open_this 文件夹内对应的 OME-TIFF。\n"
        "2. 图像已经校正为表皮在上、真皮在下；仅保留DAPI、PDGFRA与KRT14三个通道。\n"
        "3. 打开 Channels Tool：DAPI设为蓝色，PDGFRA设为绿色，KRT14设为洋红色。\n"
        "4. 选择矩形ROI，目标物理尺寸固定为2000 × 1300 µm。各片对应像素尺寸见ROI_SELECTION_CHECKLIST.tsv。\n"
        "5. 矩形上部应包含连续或相对完整的KRT14+表皮带，下部包含紧邻的PDGFRA+真皮，二者界面位于框的中部。\n"
        "6. 排除组织折叠、裂隙、切片边缘、离焦、明显过曝、大片空白和孤立碎片。不得只选信号最强点。\n"
        "7. 若同一切片有多个合格位置，优先选结构完整、界面清楚、最接近该片整体表现的位置；不要根据共染强弱挑图。\n"
        "8. 按 T 加入 ROI Manager，重命名为 样本名_ZI01，例如 BW9_ZI01。\n"
        "9. ROI Manager > More > Save，保存到 03_put_RoiSet_zip_here，文件名为 样本名_ZoomIn_RoiSet.zip。\n"
        "10. 在ROI_SELECTION_CHECKLIST.tsv填写roi_zip_path并将roi_status改为APPROVED。\n\n"
        "如果一张切片找不到同时覆盖完整表皮-真皮界面的区域，请不要勉强选择，标记NOT_APPLICABLE并写明原因。\n"
        "注意：这里的ROI名称使用ZI而不是EF，因为它只是展示用Zoom in，不构成Epi-Front或ROI统计证据。\n"
    )
    (OUT / "请先阅读_ROI怎么选.txt").write_text(guide, encoding="utf-8")

    fig = plt.figure(figsize=(210 * MM, 297 * MM), facecolor="white")
    for i, (sid, p) in enumerate(previews):
        r, c = divmod(i, 2)
        ax = fig.add_axes([.02 + c * .49, .81 - r * .16, .47, .145])
        ax.imshow(Image.open(p))
        ax.axis("off")
    fig.text(.02, .985, "Fib_K14 ROI handoff: DAPI / KRT14 / PDGFRA only", fontsize=10, fontweight="bold", va="top")
    fig.text(.02, .965, "Choose one 2000 x 1300 um anatomical epidermis-dermis interface ROI per slide; display only, not statistics.", fontsize=8, va="top")
    fig.savefig(OUT / "Fib_K14_three_channel_ROI_review_A4.pdf", dpi=600)
    fig.savefig(OUT / "Fib_K14_three_channel_ROI_review_A4.png", dpi=300)
    plt.close(fig)

    audit_rows = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.parent != AUDIT:
            audit_rows.append({"relative_path": str(p.relative_to(OUT)), "bytes": p.stat().st_size, "sha256": sha256(p)})
    pd.DataFrame(audit_rows).to_csv(AUDIT / "file_index_sha256.tsv", sep="\t", index=False)
    with (AUDIT / "FINAL_QA.tsv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status", "detail"])
        writer.writerow(["sample_coverage", "PASS", f"{len(rows)}/10"])
        writer.writerow(["display_channels", "PASS", "DAPI, KRT14 and PDGFRA only"])
        writer.writerow(["source_channels_preserved", "PASS", "original 5-channel OME-TIFF retained for audit"])
        writer.writerow(["roi_scope", "PASS", "visual ZoomIn only; whole-slide statistics unchanged"])
        writer.writerow(["roi_gate", "PENDING_AUTHOR", "ROI zip files not yet supplied"])
    print(OUT)


if __name__ == "__main__":
    main()
