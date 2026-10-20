from __future__ import annotations

import base64
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
SOURCE_HE = ROOT / "inputs" / "raw_panels" / "ppt_slide2_media" / "image5.png"
RASTER_DIR = ROOT / "assets" / "linked_rasters" / "Figure1A"
PANEL_DIR = ROOT / "outputs" / "panels"
VECTOR_DIR = ROOT / "outputs" / "vector"
RASTER_OUT = ROOT / "outputs" / "raster"
SOURCE_DATA = ROOT / "source_data"
QC_DIR = ROOT / "QC"

WIDTH_MM = 125
HEIGHT_MM = 92
PX_600 = round(WIDTH_MM / 25.4 * 600)
VERSION = "v02"

COHORT_METADATA = ROOT.parent.parent / "Figure1_Supplement_QC_20260721" / "cohort_metadata.csv"
HUMAN_UNITS = 51
RAT_UNITS = 83
SCRNA_CELLS = 729_587
STEREO_BINS = 1_031_523
VISIUM_SECTIONS = 2
FLIPPED_MODES = {"Mode2_one_way_epithelialization", "Mode4_scar_hyperplasia"}

# Coordinates were audited against the original 2317 x 276 H&E strip.
# Only black separators are removed; no brightness, contrast, color, or sharpness
# adjustment is applied.
CROPS = {
    "Mode1_unhealed_gradient": (2, 0, 545, 276),
    "Mode2_one_way_epithelialization": (584, 0, 1117, 276),
    "Mode3_almost_healed": (1153, 0, 1702, 276),
    "Mode4_scar_hyperplasia": (1741, 0, 2295, 276),
}


def ensure_dirs() -> None:
    for path in (RASTER_DIR, PANEL_DIR, VECTOR_DIR, RASTER_OUT, SOURCE_DATA, QC_DIR):
        path.mkdir(parents=True, exist_ok=True)


def crop_histology() -> dict[str, Path]:
    source = Image.open(SOURCE_HE).convert("RGB")
    outputs: dict[str, Path] = {}
    for name, box in CROPS.items():
        out = RASTER_DIR / f"Figure1A_HE_{name}.png"
        crop = source.crop(box)
        if name in FLIPPED_MODES:
            crop = crop.transpose(Image.Transpose.FLIP_LEFT_RIGHT)
        crop.save(out, dpi=(600, 600), optimize=True)
        outputs[name] = out
    return outputs


def derive_stage_counts() -> dict[str, int]:
    counts = {"inflammation": 0, "proliferation": 0, "remodeling": 0, "normal": 0, "unencoded": 0}
    with COHORT_METADATA.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            stage = row["repair_stage"]
            if stage in {"0-7dpi", "5dpb"}:
                counts["inflammation"] += 1
            elif stage in {"8-14dpi", "15-28dpi", "12dpb DPTDI-1", "12dpb DPTDI-2", "12dpb SPTDI-1", "12dpb SPTDI-2", "19dpb", "19dpb p1", "26dpb p1", "26dpb p2"}:
                counts["proliferation"] += 1
            elif stage in {"1-2mph", "2mph"}:
                counts["remodeling"] += 1
            elif stage == "Normal":
                counts["normal"] += 1
            else:
                counts["unencoded"] += 1
    assert counts == {"inflammation": 9, "proliferation": 22, "remodeling": 7, "normal": 11, "unencoded": 2}, counts
    return counts


def data_uri(path: Path) -> str:
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{payload}"


def text(x: float, y: float, value: str, cls: str = "body", anchor: str = "middle", extra: str = "") -> str:
    return f'<text x="{x}" y="{y}" class="{cls}" text-anchor="{anchor}" {extra}>{escape(value)}</text>'


def arrow(x1: float, y1: float, x2: float, y2: float, cls: str = "axis") -> str:
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" class="{cls}" marker-end="url(#arrow)"/>'


def axes(x: float, y: float, w: float, h: float) -> str:
    return (
        arrow(x, y + h, x + w, y + h)
        + arrow(x, y + h, x, y)
        + text(x + w / 2, y + h + 3.2, "x", "axislabel", extra='font-style="italic"')
        + text(x - 2.2, y + h / 2 + 1.0, "y", "axislabel", extra='font-style="italic"')
    )


def mode1(x: float, y: float, w: float, h: float) -> str:
    return f'''
    <g id="Mode_1_schematic">
      {axes(x, y, w, h)}
      <path d="M {x} {y+4} L {x+w} {y+4} L {x+w} {y+h} L {x} {y+h} Z" class="dermis"/>
      <path d="M {x} {y+4} C {x+7} {y+4}, {x+12} {y+5}, {x+16} {y+10} C {x+20} {y+15}, {x+23} {y+13}, {x+w} {y+5} L {x+w} {y+4} Z" class="wound"/>
      <line x1="{x+20.5}" y1="{y+10.5}" x2="{x+20.5}" y2="{y+h-5.0}" class="guide"/>
      <circle cx="{x+20.5}" cy="{y+10.5}" r="0.8" class="point"/>
      {text(x+6.5, y+h-3.0, "Superficial", "small")}
      {text(x+21.2, y+h-3.0, "Deep", "small")}
    </g>'''


def mode2(x: float, y: float, w: float, h: float) -> str:
    return f'''
    <g id="Mode_2_schematic">
      {axes(x, y, w, h)}
      <path d="M {x} {y+4} L {x+w} {y+4} L {x+w} {y+h} L {x} {y+h} Z" class="dermis"/>
      <path d="M {x} {y+4} L {x+12.5} {y+4} C {x+14.5} {y+5}, {x+15.5} {y+9}, {x+17} {y+11} C {x+19} {y+13}, {x+22} {y+12}, {x+w} {y+8} L {x+w} {y+4} Z" class="wound"/>
      <path d="M {x} {y+2.8} L {x+12.3} {y+2.8} C {x+13.4} {y+3.0}, {x+14.1} {y+3.4}, {x+14.9} {y+4.4}" class="epidermis-line"/>
      <line x1="{x+14.9}" y1="{y+4.4}" x2="{x+14.9}" y2="{y+h-5.0}" class="guide"/>
      <circle cx="{x+14.9}" cy="{y+4.4}" r="0.8" class="point"/>
      {arrow(x+19.0, y+8.0, x+24.7, y+8.0, "migration")}
      {text(x+6.4, y+h-2.2, "Newly epi.", "small")}
      {text(x+15.2, y+h-2.2, "Front", "small")}
      {text(x+23.0, y+h-2.2, "Un-epi", "small")}
    </g>'''


def mode3(x: float, y: float, w: float, h: float) -> str:
    return f'''
    <g id="Mode_3_schematic">
      {axes(x, y, w, h)}
      <path d="M {x} {y+4} L {x+w} {y+4} L {x+w} {y+h} L {x} {y+h} Z" class="dermis"/>
      <path d="M {x} {y+4} L {x+20.5} {y+4} C {x+22.5} {y+4.3}, {x+23.5} {y+8.5}, {x+w} {y+10.0} L {x+w} {y+4} Z" class="wound"/>
      <path d="M {x} {y+2.8} L {x+20.5} {y+2.8} C {x+21.8} {y+2.9}, {x+22.5} {y+3.4}, {x+23.2} {y+4.5}" class="epidermis-line"/>
      <line x1="{x+22.7}" y1="{y+4.4}" x2="{x+22.7}" y2="{y+h-5.0}" class="guide"/>
      <circle cx="{x+22.7}" cy="{y+4.4}" r="0.8" class="point"/>
      {text(x+10.2, y+h-2.2, "Healed area", "small")}
      {text(x+23.7, y+h-4.2, "Focal", "small")}
      {text(x+23.7, y+h-2.2, "unhealed", "small")}
    </g>'''


def mode4(x: float, y: float, w: float, h: float) -> str:
    return f'''
    <g id="Mode_4_schematic">
      {axes(x, y, w, h)}
      <path d="M {x} {y+4} L {x+w} {y+4} L {x+w} {y+h} L {x} {y+h} Z" class="dermis"/>
      <path d="M {x} {y+4} C {x+5} {y+2}, {x+7} {y-0.5}, {x+11} {y+0.8} C {x+15} {y+2.0}, {x+17} {y+4}, {x+w} {y+4} Z" class="wound"/>
      <path d="M {x} {y+2.6} C {x+5} {y+0.5}, {x+7} {y-2.0}, {x+11} {y-0.7} C {x+15} {y+0.5}, {x+17} {y+2.8}, {x+w} {y+2.8}" class="epidermis-line"/>
      <line x1="{x+10.5}" y1="{y+0.3}" x2="{x+10.5}" y2="{y+h-5.0}" class="guide"/>
      <circle cx="{x+10.5}" cy="{y+0.3}" r="0.8" class="point"/>
      {text(x+8.6, y+h-4.2, "Scar", "small")}
      {text(x+8.6, y+h-2.2, "hyperplasia", "small")}
      {text(x+22.0, y+h-2.2, "Scarless", "small")}
    </g>'''


def build_svg(crops: dict[str, Path], stage_counts: dict[str, int]) -> str:
    columns = [4.5, 35.3, 66.1, 96.9]
    col_w = 27.6
    schematic_y = 38.5
    schematic_h = 17.2
    he_y = 66.0
    he_h = 14.2
    names = list(CROPS)
    titles = [
        ("Mode 1", "Unhealed gradient"),
        ("Mode 2", "One-way epithelialization"),
        ("Mode 3", "Almost healed"),
        ("Mode 4", "Scar / hyperplasia"),
    ]

    svg: list[str] = [f'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
      width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">
      <defs>
        <marker id="arrow" markerWidth="5" markerHeight="5" refX="4.4" refY="2.5" orient="auto" markerUnits="strokeWidth">
          <path d="M0,0 L5,2.5 L0,5 Z" fill="#111111"/>
        </marker>
        <style><![CDATA[
          text {{ font-family: Arial, Helvetica, sans-serif; fill: #111111; }}
          .panel {{ font-size: 4.23px; font-weight: 700; }}
          .study {{ font-size: 3.18px; font-weight: 700; }}
          .stage {{ font-size: 2.65px; font-weight: 700; }}
          .time {{ font-size: 2.29px; }}
          .counts {{ font-size: 2.29px; font-weight: 700; }}
          .stage-count {{ font-size: 2.12px; font-weight: 700; }}
          .mode {{ font-size: 2.47px; font-weight: 700; }}
          .mode-title {{ font-size: 2.29px; font-weight: 700; }}
          .body {{ font-size: 2.12px; }}
          .small {{ font-size: 2.12px; }}
          .axislabel {{ font-size: 2.12px; }}
          .axis {{ stroke: #111111; stroke-width: 0.42; fill: none; }}
          .migration {{ stroke: #111111; stroke-width: 0.65; fill: none; }}
          .guide {{ stroke: #A34A45; stroke-width: 0.32; stroke-dasharray: 1.1 0.8; }}
          .point {{ fill: #F28E2B; stroke: #8F3E20; stroke-width: 0.25; }}
          .dermis {{ fill: #F6DAD8; stroke: #2B2B2B; stroke-width: 0.32; }}
          .wound {{ fill: #CC6C69; stroke: #2B2B2B; stroke-width: 0.32; }}
          .epidermis-line {{ fill: none; stroke: #7EA6BF; stroke-width: 1.15; stroke-linecap: round; }}
          .separator {{ stroke: #D5D5D5; stroke-width: 0.28; }}
          .he-border {{ fill: none; stroke: #4A4A4A; stroke-width: 0.28; }}
          .footer {{ font-size: 2.12px; }}
          .footer-strong {{ font-size: 2.12px; font-weight: 700; }}
          .scale-mask {{ fill: #FFFFFF; fill-opacity: 0.88; }}
          .scale-text {{ font-size: 2.12px; font-weight: 700; }}
          .scale-bar {{ stroke: #111111; stroke-width: 0.55; }}
          .morph-label {{ font-size: 2.12px; font-weight: 700; paint-order: stroke; stroke: #FFFFFF; stroke-width: 0.8px; stroke-linejoin: round; }}
        ]]></style>
        {''.join(f'<clipPath id="clip{i}"><rect x="{columns[i]}" y="{he_y}" width="{col_w}" height="{he_h}"/></clipPath>' for i in range(4))}
      </defs>
      <rect width="125" height="88" fill="#FFFFFF"/>
      <g id="Panel_A">
        <g id="A_timeline_and_counts">
          {text(0.8, 5.2, "A", "panel", "start")}
          {text(62.5, 4.8, "Human wound epithelialization study", "study")}
          {arrow(5.0, 15.0, 121.5, 15.0, "axis")}
          {text(18.0, 9.0, "Inflammation", "stage")}
          {text(18.0, 12.1, "0–7 dpb", "time")}
          {text(18.0, 20.5, f"n = {stage_counts['inflammation']} sequencing units", "stage-count")}
          {text(61.5, 9.0, "Proliferation", "stage")}
          {text(61.5, 12.1, "8–28 dpb", "time")}
          {text(61.5, 20.5, f"n = {stage_counts['proliferation']} sequencing units", "stage-count")}
          {text(103.0, 9.0, "Remodeling", "stage")}
          {text(103.0, 12.1, "1–2 mph", "time")}
          {text(103.0, 20.5, f"n = {stage_counts['remodeling']} sequencing units", "stage-count")}
          <line x1="4.5" y1="24.0" x2="124.5" y2="24.0" class="separator"/>
        </g>
        <g id="A_mode_labels">
    ''']

    for i, (mode, title) in enumerate(titles):
        cx = columns[i] + col_w / 2
        svg.append(text(cx, 28.6, mode, "mode"))
        svg.append(text(cx, 32.1, title, "mode-title"))

    svg.append('</g><g id="A_mode_schematics">')
    svg.append(mode1(columns[0], schematic_y, col_w, schematic_h))
    svg.append(mode2(columns[1], schematic_y, col_w, schematic_h))
    svg.append(mode3(columns[2], schematic_y, col_w, schematic_h))
    svg.append(mode4(columns[3], schematic_y, col_w, schematic_h))
    svg.append('</g><g id="A_histology">')
    svg.append(text(4.5, 63.1, "Representative H&E", "body", "start", 'font-weight="700"'))

    for i, name in enumerate(names):
        with Image.open(crops[name]) as im:
            iw, ih = im.size
        image_ratio = iw / ih
        draw_w = col_w
        draw_h = col_w / image_ratio
        draw_x = columns[i]
        draw_y = he_y + (he_h - draw_h) / 2
        svg.append(
            f'<image id="HE_Mode_{i+1}" x="{draw_x}" y="{draw_y}" width="{draw_w}" height="{draw_h}" '
            f'preserveAspectRatio="xMidYMid meet" clip-path="url(#clip{i})" href="{data_uri(crops[name])}"/>'
        )
        if i in {1, 3}:
            inset_label = "2 mm" if i == 1 else "1 mm"
            # Horizontal mirroring is applied to tissue orientation. The original
            # raster scale labels therefore become mirrored; mask them and redraw
            # equivalent vector scale bars at the same displayed locations.
            svg.append(f'<rect x="{columns[i]+20.2}" y="{he_y+10.2}" width="7.25" height="3.75" rx="0.35" class="scale-mask"/>')
            svg.append(text(columns[i]+23.8, he_y+12.3, "500 μm", "scale-text"))
            svg.append(f'<line x1="{columns[i]+22.1}" y1="{he_y+13.25}" x2="{columns[i]+25.5}" y2="{he_y+13.25}" class="scale-bar"/>')
            svg.append(f'<rect x="{columns[i]+0.15}" y="{he_y+10.2}" width="6.85" height="3.75" rx="0.35" class="scale-mask"/>')
            svg.append(text(columns[i]+3.55, he_y+12.3, inset_label, "scale-text"))
            svg.append(f'<line x1="{columns[i]+2.2}" y1="{he_y+13.25}" x2="{columns[i]+4.9}" y2="{he_y+13.25}" class="scale-bar"/>')
        svg.append(f'<rect x="{draw_x}" y="{draw_y}" width="{draw_w}" height="{draw_h}" class="he-border"/>')
        morph_labels = {
            0: [(0.20, ["Superficial"]), (0.79, ["Deep"])],
            1: [(0.18, ["Newly epi."]), (0.51, ["Epi-Front"]), (0.82, ["Un-epi"])],
            2: [(0.23, ["Healed"]), (0.77, ["Focal", "unhealed"])],
            3: [(0.25, ["Scar /", "hyperplasia"]), (0.76, ["Scarless", "healed"])],
        }
        for frac, label_lines in morph_labels[i]:
            for line_index, label in enumerate(label_lines):
                svg.append(text(draw_x + draw_w * frac, draw_y + 2.15 + 1.85 * line_index, label, "morph-label"))

    svg.append(f'''</g>
        <g id="A_dataset_scale_footer">
          <line x1="4.5" y1="82.5" x2="124.5" y2="82.5" class="separator"/>
          {text(31.0, 86.0, f"Human: {HUMAN_UNITS} sequencing units", "footer-strong")}
          {text(94.0, 86.0, f"Rat validation: {RAT_UNITS} sequencing units", "footer-strong")}
          {text(22.0, 89.6, f"scRNA-seq: {SCRNA_CELLS:,} cells", "footer")}
          {text(66.0, 89.6, f"Stereo-seq: {STEREO_BINS:,} cell bins", "footer")}
          {text(108.0, 89.6, f"10x Visium: {VISIUM_SECTIONS} sections", "footer")}
        </g>
      </g></svg>''')
    return ''.join(svg)


def write_provenance(crops: dict[str, Path]) -> None:
    out = SOURCE_DATA / "Figure1A_histology_sources.csv"
    with out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["display_order", "mode", "source_file", "crop_left", "crop_top", "crop_right", "crop_bottom", "processing"])
        for index, (name, box) in enumerate(CROPS.items(), 1):
            transform = "crop plus horizontal mirror; reversed raster scale labels masked and redrawn as vector; no intensity or color adjustment" if name in FLIPPED_MODES else "crop only; no intensity or color adjustment"
            writer.writerow([index, name, str(SOURCE_HE), *box, transform])


def write_count_provenance(stage_counts: dict[str, int]) -> None:
    out = SOURCE_DATA / "Figure1A_count_sources.csv"
    with out.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["display_item", "value", "unit", "source", "verification"])
        writer.writerow(["Inflammation, 0-7 dpb", stage_counts["inflammation"], "stage-coded human wound sequencing unit", str(COHORT_METADATA), "derived 2026-07-22"])
        writer.writerow(["Proliferation, 8-28 dpb", stage_counts["proliferation"], "stage-coded human wound sequencing unit", str(COHORT_METADATA), "derived 2026-07-22"])
        writer.writerow(["Remodeling, 1-2 mph", stage_counts["remodeling"], "stage-coded human wound sequencing unit", str(COHORT_METADATA), "derived 2026-07-22"])
        writer.writerow(["Normal", stage_counts["normal"], "human sequencing unit", str(COHORT_METADATA), "derived 2026-07-22; not assigned to wound-stage timeline"])
        writer.writerow(["Visium time unencoded", stage_counts["unencoded"], "human Visium section", str(COHORT_METADATA), "derived 2026-07-22; not assigned to wound-stage timeline"])
        writer.writerow(["Human", HUMAN_UNITS, "sequencing unit", "38 scRNA samples + 11 Stereo sections + 2 Visium sections", "verified"])
        writer.writerow(["Rat validation", RAT_UNITS, "sequencing unit", "53 scRNA samples + 30 Stereo sections", "verified from live H5AD metadata/object inventory"])
        writer.writerow(["scRNA-seq", SCRNA_CELLS, "cell", "human 279305 + rat 450282", "verified from live H5AD objects"])
        writer.writerow(["Stereo-seq", STEREO_BINS, "cell bin", "human 387061 + rat 644462", "verified from registered/live H5AD objects"])
        writer.writerow(["10x Visium", VISIUM_SECTIONS, "section", "epithelialization dataset router", "verified"])


def render(svg_path: Path) -> None:
    pdf_path = PANEL_DIR / f"Figure1A_{VERSION}.pdf"
    png_path = RASTER_OUT / f"Figure1A_{VERSION}_600dpi.png"
    tiff_path = RASTER_OUT / f"Figure1A_{VERSION}_600dpi.tiff"
    chrome_candidates = [
        Path(r"data/external/local_b98a42e7e75a"),
        Path(r"data/external/local_0647cf87e7e9"),
    ]
    browser = next((path for path in chrome_candidates if path.exists()), None)
    if browser is None:
        raise FileNotFoundError("Chrome or Edge is required to render SVG proof files")

    with tempfile.TemporaryDirectory(prefix="figure1a-render-") as temp_name:
        temp_dir = Path(temp_name)
        profile = temp_dir / "browser-profile"
        raw_png = temp_dir / "raw.png"
        html_path = temp_dir / "print.html"
        svg_text = svg_path.read_text(encoding="utf-8")
        html_path.write_text(
            "<!doctype html><html><head><meta charset='utf-8'><style>"
            f"@page {{ size: {WIDTH_MM}mm {HEIGHT_MM}mm; margin: 0; }}"
            f"html,body {{ margin:0; padding:0; width:{WIDTH_MM}mm; height:{HEIGHT_MM}mm; overflow:hidden; background:white; }}"
            "svg { display:block; margin:0; padding:0; }"
            "</style></head><body>" + svg_text + "</body></html>",
            encoding="utf-8",
        )
        common = [
            str(browser), "--headless=new", "--disable-gpu", "--hide-scrollbars",
            "--no-first-run", "--disable-extensions", f"--user-data-dir={profile}",
        ]
        subprocess.run(
            common + [f"--print-to-pdf={pdf_path}", "--no-pdf-header-footer", html_path.as_uri()],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        css_width = round(WIDTH_MM / 25.4 * 96)
        css_height = round(HEIGHT_MM / 25.4 * 96)
        subprocess.run(
            common + ["--force-device-scale-factor=6.25", f"--window-size={css_width},{css_height}", f"--screenshot={raw_png}", html_path.as_uri()],
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        target_h = round(HEIGHT_MM / 25.4 * 600)
        with Image.open(raw_png).convert("RGB") as raw:
            rendered = raw.resize((PX_600, target_h), Image.Resampling.LANCZOS)
            rendered.save(png_path, dpi=(600, 600), optimize=True)
            rendered.save(tiff_path, compression="tiff_lzw", dpi=(600, 600))


def main() -> None:
    ensure_dirs()
    if not SOURCE_HE.exists():
        raise FileNotFoundError(SOURCE_HE)
    stage_counts = derive_stage_counts()
    crops = crop_histology()
    svg = build_svg(crops, stage_counts)
    panel_svg = PANEL_DIR / f"Figure1A_{VERSION}.svg"
    panel_svg.write_text(svg, encoding="utf-8")
    shutil.copy2(panel_svg, VECTOR_DIR / panel_svg.name)
    write_provenance(crops)
    write_count_provenance(stage_counts)
    render(panel_svg)
    shutil.copy2(PANEL_DIR / f"Figure1A_{VERSION}.pdf", VECTOR_DIR / f"Figure1A_{VERSION}.pdf")
    (QC_DIR / f"Figure1A_{VERSION}_content_proof.txt").write_text(
        "Figure1A v02 content proof\n"
        "- No combined N displayed.\n"
        "- Stage-coded wound units: inflammation 9; proliferation 22; remodeling 7.\n"
        "- Normal units 11 and time-unencoded Visium units 2 are not assigned to the stage timeline.\n"
        "- Bottom scale: human 51 units; rat validation 83 units; scRNA 729587 cells; Stereo 1031523 bins; Visium 2 sections.\n"
        "- Human-arm icon and gross clinical photographs removed from the main panel.\n"
        "- No fluorescence images included.\n"
        "- Mode 4 scarless surface is flat/healed rather than depressed.\n"
        "- Mode 2 and Mode 4 H&E are horizontally mirrored.\n"
        "- Mirrored raster scale labels are masked and replaced by forward-reading vector scale bars.\n"
        "- H&E processing: crop and documented orientation change only; no intensity/color adjustment.\n"
        "- Final panel size: 125 x 92 mm; Arial; minimum nominal text size 6 pt.\n",
        encoding="utf-8",
    )
    print(panel_svg)
    print(PANEL_DIR / f"Figure1A_{VERSION}.pdf")
    print(RASTER_OUT / f"Figure1A_{VERSION}_600dpi.png")
    print(RASTER_OUT / f"Figure1A_{VERSION}_600dpi.tiff")


if __name__ == "__main__":
    main()
