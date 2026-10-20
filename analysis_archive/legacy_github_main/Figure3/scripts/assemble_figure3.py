"""Assemble Figure 3 on one A4 portrait page without geometric distortion.

Bioinformatic panels are linked from their source-generated PNG files. Only
experimental panels H-J are extracted from the original experimental-image
containers in the source PowerPoint. The output SVG keeps every subpanel as a
separate movable object; PDF and 600-dpi PNG use the identical layout.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
EXPERIMENTAL = ASSETS / "linked_rasters"
OUT_VECTOR = ROOT / "outputs" / "vector"
OUT_RASTER = ROOT / "outputs" / "raster"
OUT_ASSEMBLED = ROOT / "outputs" / "assembled"

PAGE_W_MM, PAGE_H_MM = 210.0, 297.0
DPI = 600
PT_PER_MM = 72.0 / 25.4
PX_PER_MM = DPI / 25.4
ARIAL = Path(r"data/external/local_16a706108b0a")
ARIAL_BOLD = Path(r"data/external/local_28ac8473f586")


PANELS = {
    "A": ["Figure3_A1_proxy_FibK14.png", "Figure3_A2_proxy_SAC.png"],
    "B": ["Figure3_B1_PAGA_Fib.png", "Figure3_B2_PAGA_SAC.png"],
    "C": ["Figure3_C1_UMAP_Time.png", "Figure3_C2_UMAP_Grade.png"],
    "D": ["Figure3_D1_Time_stats.png", "Figure3_D2_Grade_stats.png", "Figure3_D3_Time_line.png"],
    "E": ["Figure3_E_marker_dotplot.png"],
    "F": ["Figure3_F_enrichment.png"],
    "G": ["Figure3_G_spatial_trajectory.png"],
    "H": ["linked_rasters/Figure3_H_VIM_KRT5.png", "linked_rasters/Figure3_H_HE_inset.png"],
    "I": ["linked_rasters/Figure3_I_Pdgfra_lineage.png"],
    "J": ["linked_rasters/Figure3_J_sorted_cell_IF.png"],
    "K": ["Figure3_K_flow_gating.png"],
}


# x, y, width, height in mm. Related plot classes are assigned matched boxes.
BOXES = {
    "A": (5.5, 5.5, 96.0, 42.0),
    "C": (104.7, 5.5, 99.8, 42.0),
    "B": (5.5, 50.7, 96.0, 42.0),
    "D": (104.7, 50.7, 99.8, 42.0),
    "E": (5.5, 95.9, 124.0, 38.0),
    "F": (132.7, 95.9, 71.8, 38.0),
    "G": (5.5, 137.1, 199.0, 35.0),
    "H": (5.5, 175.3, 60.0, 48.0),
    "I": (68.7, 175.3, 135.8, 48.0),
    "J": (5.5, 226.5, 82.0, 64.5),
    "K": (90.7, 226.5, 113.8, 64.5),
}


def asset(rel: str) -> Path:
    return ASSETS / rel


def trim(im: Image.Image) -> Image.Image:
    """Remove only uniform outer margins; preserve all internal annotations."""
    im = im.convert("RGB")
    corner = im.getpixel((0, 0))
    bg = Image.new("RGB", im.size, corner)
    diff = ImageChops.difference(im, bg).convert("L")
    # Suppress tiny JPEG/antialias variations around the background.
    diff = diff.point(lambda p: 0 if p < 10 else 255)
    bbox = diff.getbbox()
    if not bbox:
        return im
    left, top, right, bottom = bbox
    pad = max(2, int(min(im.size) * 0.005))
    return im.crop((max(0, left - pad), max(0, top - pad), min(im.width, right + pad), min(im.height, bottom + pad)))


def fit_rect(img_size, box, pad_left=4.8, pad_top=1.0, pad_right=1.0, pad_bottom=1.0):
    x, y, w, h = box
    x += pad_left
    y += pad_top
    w -= pad_left + pad_right
    h -= pad_top + pad_bottom
    iw, ih = img_size
    scale = min(w / iw, h / ih)
    nw, nh = iw * scale, ih * scale
    return x + (w - nw) / 2, y + (h - nh) / 2, nw, nh


def panel_slots(label: str, box, images):
    x, y, w, h = box
    inner = (x + 4.8, y + 1.0, w - 5.8, h - 2.0)
    ix, iy, iw, ih = inner
    gap = 1.2
    n = len(images)
    if label in {"A", "B", "C"}:
        sw = (iw - gap) / 2
        return [(ix, iy, sw, ih), (ix + sw + gap, iy, sw, ih)]
    if label == "D":
        sw = (iw - 2 * gap) / 3
        return [(ix + i * (sw + gap), iy, sw, ih) for i in range(3)]
    if label == "H":
        # Main IF image with a small H&E inset at lower right.
        return [(ix, iy, iw, ih), (ix + iw * 0.64, iy + ih * 0.60, iw * 0.32, ih * 0.34)]
    return [inner]


def load_images(label: str):
    return [trim(Image.open(asset(rel))) for rel in PANELS[label]]


def svg_image_element(im: Image.Image, rect, element_id: str) -> str:
    bio = io.BytesIO()
    im.save(bio, format="PNG")
    encoded = base64.b64encode(bio.getvalue()).decode("ascii")
    x, y, w, h = rect
    return (
        f'<image id="{element_id}" x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" '
        f'height="{h:.3f}" preserveAspectRatio="xMidYMid meet" '
        f'href="data:image/png;base64,{encoded}"/>'
    )


def make_svg():
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{PAGE_W_MM}mm" height="{PAGE_H_MM}mm" viewBox="0 0 {PAGE_W_MM} {PAGE_H_MM}">',
        '<rect width="210" height="297" fill="white"/>',
        '<style>.panel-label{font-family:Arial,sans-serif;font-size:4.586px;font-weight:700;}</style>',
    ]
    for label in "ABCDEFGHIJK":
        box = BOXES[label]
        images = load_images(label)
        slots = panel_slots(label, box, images)
        parts.append(f'<g id="Panel_{label}">')
        parts.append(f'<text class="panel-label" x="{box[0]:.2f}" y="{box[1] + 4.6:.2f}">{label}</text>')
        for idx, (im, slot) in enumerate(zip(images, slots), 1):
            rect = fit_rect(im.size, slot, pad_left=0, pad_top=0, pad_right=0, pad_bottom=0)
            parts.append(svg_image_element(im, rect, f"Panel_{label}_{idx}"))
        parts.append("</g>")
    parts.append("</svg>")
    OUT_VECTOR.mkdir(parents=True, exist_ok=True)
    path = OUT_VECTOR / "Figure3_final.svg"
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def make_png():
    OUT_RASTER.mkdir(parents=True, exist_ok=True)
    canvas = Image.new("RGB", (round(PAGE_W_MM * PX_PER_MM), round(PAGE_H_MM * PX_PER_MM)), "white")
    draw = ImageDraw.Draw(canvas)
    label_font = ImageFont.truetype(str(ARIAL_BOLD), round(13 * DPI / 72))
    for label in "ABCDEFGHIJK":
        box = BOXES[label]
        draw.text((round(box[0] * PX_PER_MM), round(box[1] * PX_PER_MM)), label, fill="black", font=label_font)
        images = load_images(label)
        slots = panel_slots(label, box, images)
        for im, slot in zip(images, slots):
            rect = fit_rect(im.size, slot, pad_left=0, pad_top=0, pad_right=0, pad_bottom=0)
            x, y, w, h = [round(v * PX_PER_MM) for v in rect]
            resized = im.resize((max(1, w), max(1, h)), Image.Resampling.LANCZOS)
            canvas.paste(resized, (x, y))
    path = OUT_RASTER / "Figure3_final_600dpi.png"
    canvas.save(path, dpi=(DPI, DPI), optimize=True)
    return path


def make_pdf():
    OUT_VECTOR.mkdir(parents=True, exist_ok=True)
    pdfmetrics.registerFont(TTFont("Arial", str(ARIAL)))
    pdfmetrics.registerFont(TTFont("Arial-Bold", str(ARIAL_BOLD)))
    path = OUT_VECTOR / "Figure3_final.pdf"
    c = Canvas(str(path), pagesize=A4)
    page_h = A4[1]
    for label in "ABCDEFGHIJK":
        box = BOXES[label]
        c.setFont("Arial-Bold", 13)
        c.drawString(box[0] * PT_PER_MM, page_h - (box[1] + 4.6) * PT_PER_MM, label)
        images = load_images(label)
        slots = panel_slots(label, box, images)
        for im, slot in zip(images, slots):
            rect = fit_rect(im.size, slot, pad_left=0, pad_top=0, pad_right=0, pad_bottom=0)
            x, y, w, h = rect
            bio = io.BytesIO()
            im.save(bio, format="PNG")
            from reportlab.lib.utils import ImageReader
            c.drawImage(ImageReader(bio), x * PT_PER_MM, page_h - (y + h) * PT_PER_MM,
                        width=w * PT_PER_MM, height=h * PT_PER_MM, preserveAspectRatio=True, mask="auto")
    c.showPage()
    c.save()
    return path


def main():
    OUT_ASSEMBLED.mkdir(parents=True, exist_ok=True)
    for path in (make_svg(), make_pdf(), make_png()):
        print(path)


if __name__ == "__main__":
    main()
