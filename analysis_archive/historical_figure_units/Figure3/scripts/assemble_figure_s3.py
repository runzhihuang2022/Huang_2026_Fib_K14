"""Assemble Figure S3 as one A4 portrait review sheet (panels A-O frozen).

Unavailable experimental evidence is shown as an explicit evidence placeholder
so that no panel is silently omitted or populated with invented data. The
result is a review draft, not submission-ready until the listed raw assets are
supplied.
"""

from __future__ import annotations

import base64
import io
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas


ROOT = Path(__file__).resolve().parents[1]
PANELS_DIR = ROOT / "supplements" / "figures" / "panels"
ASSETS = ROOT / "assets"
OUT = ROOT / "supplements" / "figures"
PAGE_W_MM, PAGE_H_MM = 210.0, 297.0
DPI = 600
PX_PER_MM = DPI / 25.4
PT_PER_MM = 72 / 25.4
ARIAL = Path(r"data/external/local_16a706108b0a")
ARIAL_BOLD = Path(r"data/external/local_28ac8473f586")


FILES = {
    "A": PANELS_DIR / "FigureS3_A_fibroblast_umap.png",
    "B": PANELS_DIR / "FigureS3_B_QC_doublet_audit.png",
    "C": PANELS_DIR / "FigureS3_C_marker_heatmap.png",
    "D": PANELS_DIR / "FigureS3_D_proxy_statistics.png",
    "E": PANELS_DIR / "FigureS3_E_PAGA_connectivity.png",
    "F": PANELS_DIR / "FigureS3_F_pseudotime.png",
    "G": PANELS_DIR / "FigureS3_G_epithelialization_dotplot.png",
    "H": PANELS_DIR / "FigureS3_H_enrichment.png",
    "I": PANELS_DIR / "FigureS3_I_additional_spatial_maps.png",
    "L": ASSETS / "Figure3_K_flow_gating.png",
}


PLACEHOLDERS = {
    "J": "Independent VIM/KRT5 and VIM/KRT14 wound IF\nraw images and specimen IDs not located",
    "K": "Healthy-skin tissue IF controls\n(VIM/KRT5, VIM/KRT14, PDGFRA/KRT14) not located",
    "M": "Normal-versus-wound biological-replicate flow table\nnot supplied; no valid significance test possible",
    "N": "Short-term single-cell clone images\nnot located in registered assets",
    "O": "Clone-forming efficiency and marker-retention table\nnot located; replicate-level statistics unavailable",
}


BOXES = {
    "A": (5.5, 5.5, 64.2, 42.0), "B": (72.9, 5.5, 131.6, 42.0),
    "C": (5.5, 50.7, 64.2, 40.0), "D": (72.9, 50.7, 64.2, 40.0), "E": (140.3, 50.7, 64.2, 40.0),
    "F": (5.5, 93.9, 97.9, 42.0), "G": (106.6, 93.9, 97.9, 42.0),
    "H": (5.5, 139.1, 64.2, 35.0), "I": (72.9, 139.1, 131.6, 35.0),
    "J": (5.5, 177.3, 64.2, 52.0), "K": (72.9, 177.3, 64.2, 52.0), "L": (140.3, 177.3, 64.2, 52.0),
    "M": (5.5, 232.5, 64.2, 58.5), "N": (72.9, 232.5, 64.2, 58.5), "O": (140.3, 232.5, 64.2, 58.5),
}


def trim(im):
    im = im.convert("RGB")
    bg = Image.new("RGB", im.size, im.getpixel((0, 0)))
    diff = ImageChops.difference(im, bg).convert("L").point(lambda p: 0 if p < 10 else 255)
    bbox = diff.getbbox()
    if not bbox:
        return im
    pad = max(2, int(min(im.size) * 0.004))
    l, t, r, b = bbox
    return im.crop((max(0, l-pad), max(0, t-pad), min(im.width, r+pad), min(im.height, b+pad)))


def placeholder(text, width=900, height=650):
    im = Image.new("RGB", (width, height), "#F7F7F7")
    d = ImageDraw.Draw(im)
    d.rounded_rectangle((20, 20, width-20, height-20), radius=28, outline="#CB181D", width=5, fill="#FFF5F0")
    title = ImageFont.truetype(str(ARIAL_BOLD), 34)
    body = ImageFont.truetype(str(ARIAL), 28)
    d.text((width/2, height*0.30), "SOURCE EVIDENCE REQUIRED", fill="#A50F15", font=title, anchor="mm")
    lines = text.split("\n")
    y = height*0.50
    for line in lines:
        d.text((width/2, y), line, fill="#333333", font=body, anchor="mm")
        y += 46
    d.text((width/2, height*0.82), "Panel identity retained; no synthetic data inserted",
           fill="#555555", font=body, anchor="mm")
    return im


def load(label):
    if label in FILES:
        return trim(Image.open(FILES[label]))
    return placeholder(PLACEHOLDERS[label])


def fit(im, box):
    x, y, w, h = box
    x += 4.8; y += 1.0; w -= 5.8; h -= 2.0
    scale = min(w/im.width, h/im.height)
    nw, nh = im.width*scale, im.height*scale
    return x+(w-nw)/2, y+(h-nh)/2, nw, nh


def make_svg():
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="210mm" height="297mm" viewBox="0 0 210 297">',
        '<rect width="210" height="297" fill="white"/>',
        '<style>.label{font-family:Arial,sans-serif;font-size:4.586px;font-weight:700;}</style>',
    ]
    for label in "ABCDEFGHIJKLMNO":
        im = load(label)
        x, y, w, h = fit(im, BOXES[label])
        bio = io.BytesIO(); im.save(bio, "PNG")
        b64 = base64.b64encode(bio.getvalue()).decode("ascii")
        parts.append(f'<g id="Panel_{label}"><text class="label" x="{BOXES[label][0]}" y="{BOXES[label][1]+4.6}">{label}</text>')
        parts.append(f'<image id="Panel_{label}_content" x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" preserveAspectRatio="xMidYMid meet" href="data:image/png;base64,{b64}"/></g>')
    parts.append("</svg>")
    path = OUT / "FigureS3_review.svg"
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def make_png():
    page = Image.new("RGB", (round(210*PX_PER_MM), round(297*PX_PER_MM)), "white")
    draw = ImageDraw.Draw(page)
    font = ImageFont.truetype(str(ARIAL_BOLD), round(13*DPI/72))
    for label in "ABCDEFGHIJKLMNO":
        box = BOXES[label]
        draw.text((round(box[0]*PX_PER_MM), round(box[1]*PX_PER_MM)), label, fill="black", font=font)
        im = load(label)
        x, y, w, h = [round(v*PX_PER_MM) for v in fit(im, box)]
        page.paste(im.resize((max(1,w), max(1,h)), Image.Resampling.LANCZOS), (x,y))
    path = OUT / "FigureS3_review_600dpi.png"
    page.save(path, dpi=(DPI,DPI), optimize=True)
    return path


def make_pdf():
    pdfmetrics.registerFont(TTFont("Arial", str(ARIAL)))
    pdfmetrics.registerFont(TTFont("Arial-Bold", str(ARIAL_BOLD)))
    path = OUT / "FigureS3_review.pdf"
    c = Canvas(str(path), pagesize=A4)
    ph = A4[1]
    for label in "ABCDEFGHIJKLMNO":
        box = BOXES[label]
        c.setFont("Arial-Bold", 13)
        c.drawString(box[0]*PT_PER_MM, ph-(box[1]+4.6)*PT_PER_MM, label)
        im = load(label)
        x,y,w,h = fit(im,box)
        bio=io.BytesIO(); im.save(bio,"PNG")
        c.drawImage(ImageReader(bio), x*PT_PER_MM, ph-(y+h)*PT_PER_MM, w*PT_PER_MM, h*PT_PER_MM,
                    preserveAspectRatio=True, mask="auto")
    c.showPage(); c.save()
    return path


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for p in (make_svg(), make_pdf(), make_png()):
        print(p)


if __name__ == "__main__":
    main()
