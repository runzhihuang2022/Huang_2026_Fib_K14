from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
import zarr
from PIL import Image, ImageDraw, ImageFont


OUT = Path(r"data/external/local_93c497bb3a06")
OME = OUT / "01_ome_tiff"
ZOOM = OUT / "02_registered_zoomout"
AUDIT = OUT / "08_audit"
MM = 1 / 25.4
COLORS = np.asarray([
    [0.05, 0.20, 1.00],  # DAPI
    [0.00, 1.00, 0.20],  # PDGFRa 520
    [1.00, 0.10, 0.05],  # P63 570
    [1.00, 0.90, 0.05],  # KLF4 620
    [1.00, 0.00, 0.80],  # KRT14 690
], dtype=np.float32)
MARKERS = ["DAPI", "PDGFRα", "P63", "KLF4", "KRT14"]


def sample_id(name: str) -> str:
    m = re.match(r"(32-1|编号-?003|编号-?007|编号0*01A|BS1|BW3|BW7|BW8|BW9|W14-1)", name, re.I)
    if not m:
        return name.split()[0]
    x = m.group(1).replace("编号-", "").replace("编号", "")
    return x.upper()


def robust(a: np.ndarray) -> tuple[np.ndarray, float, float]:
    v = a[a > 0]
    if not v.size:
        return np.zeros_like(a, dtype=np.float32), 0.0, 1.0
    lo, hi = np.percentile(v, [1, 99.7])
    if hi <= lo:
        hi = lo + 1
    return np.clip((a.astype(np.float32)-lo)/(hi-lo), 0, 1), float(lo), float(hi)


def parse_ome(tif: tifffile.TiffFile) -> tuple[list[str], float, float]:
    root = ET.fromstring(tif.ome_metadata)
    ns = {"o": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
    pixels = root.find(".//o:Pixels", ns)
    channels = [c.get("Name", "") for c in pixels.findall("o:Channel", ns)]
    return channels, float(pixels.get("PhysicalSizeX")), float(pixels.get("PhysicalSizeY"))


def font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for p in [Path(r"data/external/local_16a706108b0a"), Path(r"data/external/local_c87def718719")]:
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def annotate(rgb: np.ndarray, sid: str, px_um: float) -> Image.Image:
    im = Image.fromarray((np.clip(rgb, 0, 1)*255).astype("uint8"))
    draw = ImageDraw.Draw(im)
    f = font(max(12, round(im.width/110)))
    bold = font(max(14, round(im.width/95)))
    draw.text((12, 8), sid, fill="white", font=bold, stroke_width=1, stroke_fill="black")
    label = "DAPI / PDGFRα / P63 / KLF4 / KRT14"
    draw.text((12, im.height - 2*f.size - 8), label, fill="white", font=f,
              stroke_width=1, stroke_fill="black")
    # Choose a physical scale bar close to 20% of image width.
    target_um = im.width * px_um * 0.20
    choices = np.asarray([100, 200, 500, 1000, 2000, 5000], dtype=float)
    bar_um = float(choices[np.argmin(np.abs(choices-target_um))])
    bar_px = max(1, round(bar_um/px_um))
    x1, y = im.width-15, im.height-22
    x0 = x1-bar_px
    draw.line((x0, y, x1, y), fill="white", width=max(2, round(im.width/500)))
    txt = f"{bar_um/1000:g} mm" if bar_um >= 1000 else f"{bar_um:g} µm"
    box = draw.textbbox((0,0), txt, font=f)
    draw.text(((x0+x1-(box[2]-box[0]))/2, y-f.size-4), txt, fill="white", font=f,
              stroke_width=1, stroke_fill="black")
    return im


def process_one(path: Path) -> dict:
    with tifffile.TiffFile(path) as tif:
        channels, base_x, base_y = parse_ome(tif)
        levels = tif.series[0].levels
        # Highest-resolution pyramid level whose longest side is <= 2200 px.
        eligible = [i for i, lev in enumerate(levels) if max(lev.shape[-2:]) <= 2200]
        level = min(eligible) if eligible else len(levels)-1
        store = tifffile.imread(path, series=0, level=level, aszarr=True)
        try:
            arr = np.asarray(zarr.open(store, mode="r"))
        finally:
            close = getattr(store, "close", None)
            if close: close()
    if arr.shape[0] != 5:
        raise ValueError(f"Expected 5 channels, found {arr.shape}")
    norms, windows = [], []
    for ch in arr:
        x, lo, hi = robust(ch); norms.append(x); windows.append([lo, hi])
    rgb = np.zeros((*arr.shape[1:], 3), dtype=np.float32)
    for x, c, gain in zip(norms, COLORS, [0.8, 0.7, 0.7, 0.7, 0.7]):
        rgb += x[..., None]*c*gain
    rgb = np.clip(rgb, 0, 1) ** 0.85
    scale = 2**level
    level_x, level_y = base_x*scale, base_y*scale
    sid = sample_id(path.name)
    annotated = annotate(rgb, sid, level_x)
    png = ZOOM / f"{sid}_Fib_K14_zoomout_merge.png"
    tifpath = ZOOM / f"{sid}_Fib_K14_zoomout_merge.tif"
    dapi = ZOOM / f"{sid}_DAPI_analysis_reference.ome.tiff"
    annotated.save(png, dpi=(600,600))
    annotated.save(tifpath, compression="tiff_lzw", dpi=(600,600))
    tifffile.imwrite(dapi, arr[0], ome=True, bigtiff=False, photometric="minisblack",
                     compression="zlib", metadata={"axes":"YX", "PhysicalSizeX":level_x,
                     "PhysicalSizeXUnit":"µm", "PhysicalSizeY":level_y,
                     "PhysicalSizeYUnit":"µm", "Channel":{"Name":["DAPI"]}})
    return {
        "sample_id": sid, "ome_tiff": str(path), "source_channel_names": "|".join(channels),
        "marker_mapping": "DAPI=DAPI|520=PDGFRα|570=P63|620=KLF4|690=KRT14",
        "pyramid_level": level, "zoomout_shape_cyx": "x".join(map(str, arr.shape)),
        "zoomout_pixel_size_x_um": level_x, "zoomout_pixel_size_y_um": level_y,
        "display_windows_lo_hi": json.dumps(windows), "orientation_status":"REVIEW_REQUIRED_source_orientation_preserved",
        "registration_status":"NOT_APPLICABLE_single_Fib_section_only",
        "zoomout_png":str(png), "zoomout_tif":str(tifpath), "dapi_reference":str(dapi),
    }


def atlas(rows: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(210*MM, 297*MM), facecolor="white")
    for i, row in rows.iterrows():
        r, c = divmod(i, 5)
        ax = fig.add_axes([c*(42/210), 1-.02-(r+1)*(28/297), 42/210, 28/297])
        ax.imshow(Image.open(row.zoomout_png)); ax.axis("off")
    fig.text(.01,.77,"Fib_K14 ZoomOut review | source orientation preserved; anatomical review required",fontsize=8)
    fig.savefig(ZOOM/"Fib_K14_all_10_zoomout_A4.png",dpi=600,bbox_inches=None)
    fig.savefig(ZOOM/"Fib_K14_all_10_zoomout_A4.pdf",dpi=600,bbox_inches=None)
    plt.close(fig)


def main() -> None:
    ZOOM.mkdir(parents=True, exist_ok=True); AUDIT.mkdir(parents=True, exist_ok=True)
    paths = sorted([p for p in OME.glob("*.ome.tiff") if not p.name.endswith(".partial")])
    rows = pd.DataFrame([process_one(p) for p in paths]).sort_values("sample_id").reset_index(drop=True)
    rows.to_csv(OUT/"00_manifest"/"zoomout_manifest.tsv",sep="\t",index=False)
    atlas(rows)
    config = {
        "raw_dir": str(Path(r"data/external/local_1e2e2e4ac7c4")),
        "output_dir": str(OUT), "pixel_size_um": 0.25284451246261597,
        "primary_radius_um": 10, "permutations": 9999, "minimum_rois_per_slide": 2,
        "recommended_rois_per_slide":"3-5", "roi_archive_pattern":"{sample_id}_RoiSet.zip",
        "analysis_scope":"standalone_Fib_sections; no serial-section registration performed",
        "samples":[{"sample_id":r.sample_id,"fixed_kfbf":r.ome_tiff,"moving_kfbf":None,
                    "fixed_dapi_channel":"DAPI","moving_dapi_channel":None,
                    "population_A":None,"population_B":["PDGFRα","KRT14"],
                    "registered_zoomout_image":r.dapi_reference,
                    "registration_status":"NOT_APPLICABLE_single_Fib_section_only"}
                   for r in rows.itertuples()],
        "figure":{"a4_width_mm":210,"panel_width_fraction":"1/5","panel_shape":"square",
                  "arrangement":"1x3","font_family":"Arial","font_min_pt":8,"font_max_pt":12},
    }
    (OUT/"00_manifest"/"run_config.json").write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding="utf-8")
    print(f"zoomouts={len(rows)}")


if __name__ == "__main__":
    main()
