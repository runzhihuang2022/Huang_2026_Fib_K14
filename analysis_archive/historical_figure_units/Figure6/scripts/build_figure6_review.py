from __future__ import annotations

import base64
import hashlib
import io
import shutil
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "inputs" / "current"
OUT_ASSEMBLED = ROOT / "outputs" / "assembled"
OUT_RASTER = ROOT / "outputs" / "raster"
OUT_VECTOR = ROOT / "outputs" / "vector"
SUPP = ROOT / "supplements" / "figures"

MM_TO_PT = 72.0 / 25.4
DPI = 600
PX_PER_MM = DPI / 25.4
PAGE_W_MM, PAGE_H_MM = 210.0, 297.0

ARIAL = Path(r"data/external/local_16a706108b0a")
ARIAL_BOLD = Path(r"data/external/local_28ac8473f586")
ARIAL_ITALIC = Path(r"data/external/local_d871a81475d6")

COLORS = {
    "navy": "#173B5C",
    "blue": "#3F88B7",
    "teal": "#2A9D8F",
    "coral": "#E76F51",
    "gold": "#E9A82E",
    "pink": "#F7D8D5",
    "dermis": "#FBE7E2",
    "fib": "#BFE7DF",
    "gray": "#667785",
    "light": "#F4F7F8",
    "border": "#8EA0AA",
}


def ensure_dirs() -> None:
    for p in (OUT_ASSEMBLED, OUT_RASTER, OUT_VECTOR, SUPP):
        p.mkdir(parents=True, exist_ok=True)


def data_uri(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def svg_text(x, y, text, size=8, weight="normal", anchor="start", style="normal", fill="#111111"):
    safe = (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
    return (
        f'<text x="{x}" y="{y}" font-family="Arial" font-size="{size}pt" '
        f'font-weight="{weight}" font-style="{style}" text-anchor="{anchor}" fill="{fill}">{safe}</text>'
    )


def mechanism_svg_group(x: float, y: float, w: float, h: float) -> str:
    sx, sy = w / 194.0, h / 78.0
    def X(v): return x + v * sx
    def Y(v): return y + v * sy
    parts = [f'<g id="Panel_G" font-family="Arial">']
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="#FFF9F7" stroke="#8EA0AA" stroke-width="0.35"/>')
    parts.append(svg_text(X(97), Y(7), "Injury-associated Fib_K14 epithelial-like plasticity", 11, "bold", "middle", fill=COLORS["navy"]))
    # Epidermal sheet and wound gap.
    parts.append(f'<path d="M {X(3)} {Y(17)} L {X(57)} {Y(17)} Q {X(67)} {Y(18)} {X(77)} {Y(26)} Q {X(91)} {Y(34)} {X(108)} {Y(28)} Q {X(124)} {Y(19)} {X(139)} {Y(17)} L {X(191)} {Y(17)} L {X(191)} {Y(31)} Q {X(155)} {Y(27)} {X(128)} {Y(35)} Q {X(103)} {Y(43)} {X(79)} {Y(34)} Q {X(59)} {Y(27)} {X(3)} {Y(31)} Z" fill="{COLORS["pink"]}" stroke="#D66B67" stroke-width="0.55"/>')
    parts.append(svg_text(X(151), Y(25), "Advancing epithelial tongue", 8, "bold", "middle", fill="#9F3C3A"))
    parts.append(svg_text(X(98), Y(15), "Epi-Front", 9, "bold", "middle", fill=COLORS["coral"]))
    # Dermal fibroblast and transitional state.
    parts.append(f'<ellipse cx="{X(25)}" cy="{Y(56)}" rx="{10*sx}" ry="{3.2*sy}" fill="#DDE4E7" stroke="#7E8C93" stroke-width="0.45" transform="rotate(-10 {X(25)} {Y(56)})"/>')
    parts.append(svg_text(X(25), Y(66), "Dermal fibroblast", 9, "bold", "middle"))
    parts.append(svg_text(X(25), Y(71), "PDGFRA / VIM", 8, anchor="middle", style="italic", fill=COLORS["gray"]))
    parts.append(f'<rect x="{X(70)}" y="{Y(47)}" width="{39*sx}" height="{26*sy}" rx="{3*sx}" fill="{COLORS["fib"]}" stroke="{COLORS["teal"]}" stroke-width="0.6"/>')
    parts.append(f'<ellipse cx="{X(89.5)}" cy="{Y(58)}" rx="{8*sx}" ry="{5*sy}" fill="#70C7B8" stroke="{COLORS["teal"]}" stroke-width="0.45"/>')
    parts.append(svg_text(X(89.5), Y(53), "Fib_K14", 10, "bold", "middle", fill="#116F68"))
    parts.append(svg_text(X(89.5), Y(66), "PDGFRA / VIM retained", 7.5, anchor="middle", style="italic"))
    parts.append(svg_text(X(89.5), Y(70), "KRT14 / KRT5 acquired", 7.5, anchor="middle", style="italic"))
    # KLF4-supported transition.
    parts.append(f'<line x1="{X(37)}" y1="{Y(56)}" x2="{X(67)}" y2="{Y(56)}" stroke="{COLORS["coral"]}" stroke-width="1.0" marker-end="url(#arrowCoral)"/>')
    parts.append(svg_text(X(52), Y(52), "KLF4", 9, "bold", "middle", fill=COLORS["coral"]))
    # IGFL1 inferred niche cue.
    parts.append(f'<circle cx="{X(122)}" cy="{Y(37)}" r="{1.5*sx}" fill="{COLORS["gold"]}"/><circle cx="{X(119)}" cy="{Y(41)}" r="{1.2*sx}" fill="{COLORS["gold"]}"/><circle cx="{X(124)}" cy="{Y(43)}" r="{1.1*sx}" fill="{COLORS["gold"]}"/>')
    parts.append(f'<path d="M {X(121)} {Y(44)} Q {X(113)} {Y(51)} {X(108)} {Y(55)}" fill="none" stroke="{COLORS["gold"]}" stroke-width="0.9" stroke-dasharray="2.4,1.6" marker-end="url(#arrowGold)"/>')
    parts.append(svg_text(X(126), Y(38), "IGFL1-associated cue", 8, "bold", fill="#A66F00"))
    # Association with re-epithelialization, intentionally dashed.
    parts.append(f'<path d="M {X(111)} {Y(58)} Q {X(143)} {Y(56)} {X(167)} {Y(43)}" fill="none" stroke="{COLORS["teal"]}" stroke-width="1.0" stroke-dasharray="2.4,1.6" marker-end="url(#arrowTeal)"/>')
    parts.append(svg_text(X(156), Y(54), "associated with", 8, anchor="middle", fill=COLORS["gray"]))
    parts.append(svg_text(X(169), Y(39), "Re-epithelialization", 10, "bold", "middle", fill=COLORS["navy"]))
    # Compact evidence grammar.
    parts.append(f'<line x1="{X(126)}" y1="{Y(69)}" x2="{X(139)}" y2="{Y(69)}" stroke="{COLORS["coral"]}" stroke-width="0.8"/>')
    parts.append(svg_text(X(141), Y(71), "perturbation-supported", 8, fill=COLORS["gray"]))
    parts.append(f'<line x1="{X(126)}" y1="{Y(74)}" x2="{X(139)}" y2="{Y(74)}" stroke="{COLORS["teal"]}" stroke-width="0.8" stroke-dasharray="2.4,1.6"/>')
    parts.append(svg_text(X(141), Y(76), "inferred / associated", 8, fill=COLORS["gray"]))
    parts.append('</g>')
    return "".join(parts)


def build_main_svg() -> Path:
    ab = data_uri(INPUT / "Figure6AB_cross_species_current.png")
    cd = data_uri(INPUT / "Figure6CD_spatial_current.png")
    e = data_uri(INPUT / "Figure6E_KLF4_current.png")
    defs = f'''<defs>
      <marker id="arrowCoral" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["coral"]}"/></marker>
      <marker id="arrowGold" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["gold"]}"/></marker>
      <marker id="arrowTeal" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["teal"]}"/></marker>
    </defs>'''
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="210mm" height="297mm" viewBox="0 0 210 297">', defs, '<rect width="210" height="297" fill="white"/>']
    parts.append(f'<g id="Panel_AB"><image x="8" y="8" width="194" height="70.2" xlink:href="{ab}" preserveAspectRatio="xMidYMid meet"/></g>')
    parts.append(svg_text(4.0, 12.0, "A", 12, "bold"))
    parts.append(svg_text(142.0, 12.0, "B", 12, "bold"))
    parts.append(f'<g id="Panel_CD"><image x="8" y="82" width="75.7" height="107" xlink:href="{cd}" preserveAspectRatio="xMidYMid meet"/></g>')
    parts.append(svg_text(4.0, 86.0, "C", 12, "bold"))
    parts.append(svg_text(43.0, 86.0, "D", 12, "bold"))
    parts.append(f'<g id="Panel_E"><image x="87.7" y="82" width="114.3" height="60.4" xlink:href="{e}" preserveAspectRatio="xMidYMid meet"/></g>')
    parts.append(svg_text(84.0, 86.0, "E", 12, "bold"))
    parts.append('<g id="Panel_F"><rect x="87.7" y="146.4" width="114.3" height="42.6" rx="1.5" fill="#F4F7F8" stroke="#8EA0AA" stroke-width="0.45" stroke-dasharray="2.2,1.5"/>')
    parts.append(svg_text(84.0, 150.4, "F", 12, "bold"))
    parts.append(svg_text(144.85, 160.5, "SOURCE REQUIRED", 11, "bold", "middle", fill="#B54038"))
    parts.append(svg_text(144.85, 168.0, "Krt14-CreER; Pdgfra-DreER tracing", 9, anchor="middle", style="italic"))
    parts.append(svg_text(144.85, 175.0, "strategy + representative wound evidence", 8, anchor="middle", fill=COLORS["gray"]))
    parts.append(svg_text(144.85, 181.5, "Do not submit this placeholder", 8, "bold", "middle", fill="#B54038"))
    parts.append('</g>')
    parts.append(svg_text(4.0, 197.0, "G", 12, "bold"))
    parts.append(mechanism_svg_group(8, 193, 194, 78))
    parts.append(svg_text(202, 289, "Figure 6 review proof v01 - not submission ready", 8, anchor="end", fill="#667785"))
    parts.append('</svg>')
    out = OUT_VECTOR / "Figure6_review_v01_editable.svg"
    out.write_text("".join(parts), encoding="utf-8")
    return out


def build_g_panel_svg() -> Path:
    defs = f'''<defs>
      <marker id="arrowCoral" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["coral"]}"/></marker>
      <marker id="arrowGold" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["gold"]}"/></marker>
      <marker id="arrowTeal" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="{COLORS["teal"]}"/></marker>
    </defs>'''
    out = OUT_VECTOR / "Figure6G_mechanism_v01_editable.svg"
    out.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="194mm" height="78mm" viewBox="0 0 194 78">'
        + defs + mechanism_svg_group(0, 0, 194, 78) + '</svg>', encoding="utf-8"
    )
    page = Image.new("RGB", (p(194), p(78)), "white")
    draw_pil_mechanism(page, (0, 0, 194, 78))
    page.save(OUT_RASTER / "Figure6G_mechanism_v01_600dpi.png", dpi=(DPI, DPI), optimize=True)
    return out


def register_fonts() -> None:
    pdfmetrics.registerFont(TTFont("Arial", str(ARIAL)))
    pdfmetrics.registerFont(TTFont("Arial-Bold", str(ARIAL_BOLD)))
    pdfmetrics.registerFont(TTFont("Arial-Italic", str(ARIAL_ITALIC)))


def mm(v: float) -> float:
    return v * MM_TO_PT


def pdf_y(y_top_mm: float, h_mm: float = 0) -> float:
    return mm(PAGE_H_MM - y_top_mm - h_mm)


def draw_pdf_text(c, x, y_top, text, size=8, font="Arial", color=black, anchor="left"):
    c.setFont(font, size)
    c.setFillColor(color)
    width = pdfmetrics.stringWidth(text, font, size)
    xpt = mm(x)
    if anchor == "middle": xpt -= width / 2
    elif anchor == "right": xpt -= width
    c.drawString(xpt, pdf_y(y_top) - size, text)


def draw_pdf_image(c, path, x, y, w, h):
    c.drawImage(ImageReader(str(path)), mm(x), pdf_y(y, h), mm(w), mm(h), preserveAspectRatio=True, anchor='c', mask='auto')


def draw_pdf_mechanism(c, x, y, w, h):
    X = lambda v: mm(x + v * w / 194.0)
    Y = lambda v: pdf_y(y + v * h / 78.0)
    c.setStrokeColor(HexColor(COLORS["border"])); c.setFillColor(HexColor("#FFF9F7")); c.roundRect(mm(x), pdf_y(y, h), mm(w), mm(h), mm(1.5), stroke=1, fill=1)
    draw_pdf_text(c, x+w/2, y+4, "Injury-associated Fib_K14 epithelial-like plasticity", 11, "Arial-Bold", HexColor(COLORS["navy"]), "middle")
    # Epidermal band.
    c.setFillColor(HexColor(COLORS["pink"])); c.setStrokeColor(HexColor("#D66B67")); c.setLineWidth(mm(0.45))
    p = c.beginPath(); p.moveTo(X(3),Y(17)); p.lineTo(X(57),Y(17)); p.curveTo(X(67),Y(18),X(70),Y(23),X(77),Y(26)); p.curveTo(X(91),Y(34),X(108),Y(28),X(124),Y(19)); p.curveTo(X(130),Y(17),X(136),Y(17),X(139),Y(17)); p.lineTo(X(191),Y(17)); p.lineTo(X(191),Y(31)); p.curveTo(X(155),Y(27),X(128),Y(35),X(103),Y(43)); p.curveTo(X(79),Y(34),X(59),Y(27),X(3),Y(31)); p.close(); c.drawPath(p,stroke=1,fill=1)
    draw_pdf_text(c, x+w*151/194, y+h*22/78, "Advancing epithelial tongue", 8, "Arial-Bold", HexColor("#9F3C3A"), "middle")
    draw_pdf_text(c, x+w*98/194, y+h*12/78, "Epi-Front", 9, "Arial-Bold", HexColor(COLORS["coral"]), "middle")
    # Cells and arrows.
    c.setFillColor(HexColor("#DDE4E7")); c.setStrokeColor(HexColor("#7E8C93")); c.ellipse(X(15),Y(59.2),X(35),Y(52.8),stroke=1,fill=1)
    draw_pdf_text(c, x+w*25/194, y+h*63/78, "Dermal fibroblast", 9, "Arial-Bold", black, "middle")
    draw_pdf_text(c, x+w*25/194, y+h*68/78, "PDGFRA / VIM", 8, "Arial-Italic", HexColor(COLORS["gray"]), "middle")
    c.setFillColor(HexColor(COLORS["fib"])); c.setStrokeColor(HexColor(COLORS["teal"])); c.roundRect(X(70),Y(73),X(109)-X(70),Y(47)-Y(73),mm(2),stroke=1,fill=1)
    draw_pdf_text(c, x+w*89.5/194, y+h*50/78, "Fib_K14", 10, "Arial-Bold", HexColor("#116F68"), "middle")
    draw_pdf_text(c, x+w*89.5/194, y+h*63/78, "PDGFRA / VIM retained", 7.5, "Arial-Italic", black, "middle")
    draw_pdf_text(c, x+w*89.5/194, y+h*68/78, "KRT14 / KRT5 acquired", 7.5, "Arial-Italic", black, "middle")
    c.setStrokeColor(HexColor(COLORS["coral"])); c.setLineWidth(mm(0.8)); c.line(X(37),Y(56),X(67),Y(56)); c.setFillColor(HexColor(COLORS["coral"])); c.setStrokeColor(HexColor(COLORS["coral"])); c.drawPath(_pdf_arrow(c, X(67),Y(56),mm(2.2),0),fill=1,stroke=0)
    draw_pdf_text(c, x+w*52/194, y+h*49/78, "KLF4", 9, "Arial-Bold", HexColor(COLORS["coral"]), "middle")
    c.setFillColor(HexColor(COLORS["gold"])); c.circle(X(122),Y(37),mm(1.1),fill=1,stroke=0); c.circle(X(119),Y(41),mm(0.9),fill=1,stroke=0); c.circle(X(124),Y(43),mm(0.8),fill=1,stroke=0)
    draw_pdf_text(c, x+w*126/194, y+h*34/78, "IGFL1-associated cue", 8, "Arial-Bold", HexColor("#A66F00"), "left")
    c.setStrokeColor(HexColor(COLORS["gold"])); c.setDash(mm(2),mm(1.4)); c.line(X(121),Y(44),X(109),Y(55)); c.setDash()
    c.setStrokeColor(HexColor(COLORS["teal"])); c.setDash(mm(2),mm(1.4)); c.line(X(111),Y(58),X(165),Y(43)); c.setDash()
    draw_pdf_text(c, x+w*155/194, y+h*51/78, "associated with", 8, "Arial", HexColor(COLORS["gray"]), "middle")
    draw_pdf_text(c, x+w*169/194, y+h*36/78, "Re-epithelialization", 10, "Arial-Bold", HexColor(COLORS["navy"]), "middle")
    c.setStrokeColor(HexColor(COLORS["coral"])); c.line(X(126),Y(69),X(139),Y(69)); draw_pdf_text(c,x+w*141/194,y+h*67/78,"perturbation-supported",8,"Arial",HexColor(COLORS["gray"]),"left")
    c.setStrokeColor(HexColor(COLORS["teal"])); c.setDash(mm(2),mm(1.4)); c.line(X(126),Y(74),X(139),Y(74)); c.setDash(); draw_pdf_text(c,x+w*141/194,y+h*72/78,"inferred / associated",8,"Arial",HexColor(COLORS["gray"]),"left")


def _pdf_arrow(c, x, y, size, angle):
    p = c.beginPath(); p.moveTo(x,y); p.lineTo(x-size,y+size*0.55); p.lineTo(x-size,y-size*0.55); p.close(); return p


def build_main_pdf() -> Path:
    register_fonts()
    out = OUT_VECTOR / "Figure6_review_v01.pdf"
    c = canvas.Canvas(str(out), pagesize=A4)
    draw_pdf_image(c, INPUT/"Figure6AB_cross_species_current.png",8,8,194,70.2)
    draw_pdf_text(c,4,8,"A",12,"Arial-Bold"); draw_pdf_text(c,142,8,"B",12,"Arial-Bold")
    draw_pdf_image(c, INPUT/"Figure6CD_spatial_current.png",8,82,75.7,107)
    draw_pdf_text(c,4,82,"C",12,"Arial-Bold"); draw_pdf_text(c,43,82,"D",12,"Arial-Bold")
    draw_pdf_image(c, INPUT/"Figure6E_KLF4_current.png",87.7,82,114.3,60.4)
    draw_pdf_text(c,84,82,"E",12,"Arial-Bold")
    c.setFillColor(HexColor(COLORS["light"])); c.setStrokeColor(HexColor(COLORS["border"])); c.setDash(mm(2.2),mm(1.5)); c.roundRect(mm(87.7),pdf_y(146.4,42.6),mm(114.3),mm(42.6),mm(1.5),stroke=1,fill=1); c.setDash()
    draw_pdf_text(c,84,146.4,"F",12,"Arial-Bold")
    draw_pdf_text(c,144.85,155.5,"SOURCE REQUIRED",11,"Arial-Bold",HexColor("#B54038"),"middle")
    draw_pdf_text(c,144.85,164,"Krt14-CreER; Pdgfra-DreER tracing",9,"Arial-Italic",black,"middle")
    draw_pdf_text(c,144.85,171,"strategy + representative wound evidence",8,"Arial",HexColor(COLORS["gray"]),"middle")
    draw_pdf_text(c,144.85,178,"Do not submit this placeholder",8,"Arial-Bold",HexColor("#B54038"),"middle")
    draw_pdf_text(c,4,193,"G",12,"Arial-Bold")
    draw_pdf_mechanism(c,8,193,194,78)
    draw_pdf_text(c,202,286,"Figure 6 review proof v01 - not submission ready",8,"Arial",HexColor(COLORS["gray"]),"right")
    c.showPage(); c.save()
    return out


def pil_font(size_pt, bold=False, italic=False):
    path = ARIAL_BOLD if bold else ARIAL_ITALIC if italic else ARIAL
    return ImageFont.truetype(str(path), round(size_pt * DPI / 72))


def p(v): return round(v * PX_PER_MM)


def draw_pil_text(draw, xy, text, size=8, bold=False, italic=False, fill="#111111", anchor="la"):
    draw.text((p(xy[0]),p(xy[1])),text,font=pil_font(size,bold,italic),fill=fill,anchor=anchor)


def paste_fit(page, img_path, box):
    x,y,w,h=box
    with Image.open(img_path).convert("RGB") as im:
        target=(p(w),p(h)); im.thumbnail(target,Image.Resampling.LANCZOS)
        ox=p(x)+(target[0]-im.width)//2; oy=p(y)+(target[1]-im.height)//2
        page.paste(im,(ox,oy))


def draw_pil_mechanism(page, box):
    x,y,w,h=box; d=ImageDraw.Draw(page)
    X=lambda v:p(x+v*w/194); Y=lambda v:p(y+v*h/78)
    d.rounded_rectangle((p(x),p(y),p(x+w),p(y+h)),radius=p(1.5),fill="#FFF9F7",outline=COLORS["border"],width=max(1,p(.35)))
    draw_pil_text(d,(x+w/2,y+4),"Injury-associated Fib_K14 epithelial-like plasticity",11,True,fill=COLORS["navy"],anchor="ma")
    poly=[(X(3),Y(17)),(X(57),Y(17)),(X(77),Y(26)),(X(103),Y(43)),(X(128),Y(35)),(X(155),Y(27)),(X(191),Y(31)),(X(191),Y(17)),(X(139),Y(17)),(X(124),Y(19)),(X(108),Y(28)),(X(91),Y(34)),(X(77),Y(26)),(X(57),Y(17)),(X(3),Y(31))]
    d.polygon(poly,fill=COLORS["pink"],outline="#D66B67")
    draw_pil_text(d,(x+w*151/194,y+h*22/78),"Advancing epithelial tongue",8,True,fill="#9F3C3A",anchor="ma")
    draw_pil_text(d,(x+w*98/194,y+h*12/78),"Epi-Front",9,True,fill=COLORS["coral"],anchor="ma")
    d.ellipse((X(15),Y(52.8),X(35),Y(59.2)),fill="#DDE4E7",outline="#7E8C93",width=p(.35))
    draw_pil_text(d,(x+w*25/194,y+h*63/78),"Dermal fibroblast",9,True,anchor="ma")
    draw_pil_text(d,(x+w*25/194,y+h*68/78),"PDGFRA / VIM",8,italic=True,fill=COLORS["gray"],anchor="ma")
    d.rounded_rectangle((X(70),Y(47),X(109),Y(73)),radius=p(2),fill=COLORS["fib"],outline=COLORS["teal"],width=p(.5))
    draw_pil_text(d,(x+w*89.5/194,y+h*50/78),"Fib_K14",10,True,fill="#116F68",anchor="ma")
    draw_pil_text(d,(x+w*89.5/194,y+h*63/78),"PDGFRA / VIM retained",7.5,italic=True,anchor="ma")
    draw_pil_text(d,(x+w*89.5/194,y+h*68/78),"KRT14 / KRT5 acquired",7.5,italic=True,anchor="ma")
    d.line((X(37),Y(56),X(67),Y(56)),fill=COLORS["coral"],width=p(.8)); d.polygon([(X(67),Y(56)),(X(64),Y(54.5)),(X(64),Y(57.5))],fill=COLORS["coral"])
    draw_pil_text(d,(x+w*52/194,y+h*49/78),"KLF4",9,True,fill=COLORS["coral"],anchor="ma")
    for vx,vy,r in [(122,37,1.2),(119,41,.9),(124,43,.8)]: d.ellipse((X(vx-r),Y(vy-r),X(vx+r),Y(vy+r)),fill=COLORS["gold"])
    draw_pil_text(d,(x+w*126/194,y+h*34/78),"IGFL1-associated cue",8,True,fill="#A66F00",anchor="la")
    # Dashed association lines.
    for a,b in [((121,44),(109,55)),((111,58),(165,43))]:
        x1,y1=X(a[0]),Y(a[1]); x2,y2=X(b[0]),Y(b[1]); steps=12
        for i in range(0,steps,2):
            t1=i/steps; t2=min((i+1)/steps,1); color=COLORS["gold"] if a[0]==121 else COLORS["teal"]
            d.line((x1+(x2-x1)*t1,y1+(y2-y1)*t1,x1+(x2-x1)*t2,y1+(y2-y1)*t2),fill=color,width=p(.7))
    draw_pil_text(d,(x+w*155/194,y+h*51/78),"associated with",8,fill=COLORS["gray"],anchor="ma")
    draw_pil_text(d,(x+w*169/194,y+h*36/78),"Re-epithelialization",10,True,fill=COLORS["navy"],anchor="ma")
    d.line((X(126),Y(69),X(139),Y(69)),fill=COLORS["coral"],width=p(.7)); draw_pil_text(d,(x+w*141/194,y+h*67/78),"perturbation-supported",8,fill=COLORS["gray"],anchor="la")
    for i in range(3): d.line((X(126+i*5),Y(74),X(129+i*5),Y(74)),fill=COLORS["teal"],width=p(.7))
    draw_pil_text(d,(x+w*141/194,y+h*72/78),"inferred / associated",8,fill=COLORS["gray"],anchor="la")


def build_main_raster() -> tuple[Path, Path]:
    page=Image.new("RGB",(p(PAGE_W_MM),p(PAGE_H_MM)),"white"); d=ImageDraw.Draw(page)
    paste_fit(page,INPUT/"Figure6AB_cross_species_current.png",(8,8,194,70.2)); draw_pil_text(d,(4,8),"A",12,True); draw_pil_text(d,(142,8),"B",12,True)
    paste_fit(page,INPUT/"Figure6CD_spatial_current.png",(8,82,75.7,107)); draw_pil_text(d,(4,82),"C",12,True); draw_pil_text(d,(43,82),"D",12,True)
    paste_fit(page,INPUT/"Figure6E_KLF4_current.png",(87.7,82,114.3,60.4)); draw_pil_text(d,(84,82),"E",12,True)
    d.rounded_rectangle((p(87.7),p(146.4),p(202),p(189)),radius=p(1.5),fill=COLORS["light"],outline=COLORS["border"],width=p(.35)); draw_pil_text(d,(84,146.4),"F",12,True)
    draw_pil_text(d,(144.85,155.5),"SOURCE REQUIRED",11,True,fill="#B54038",anchor="ma"); draw_pil_text(d,(144.85,164),"Krt14-CreER; Pdgfra-DreER tracing",9,italic=True,anchor="ma"); draw_pil_text(d,(144.85,171),"strategy + representative wound evidence",8,fill=COLORS["gray"],anchor="ma"); draw_pil_text(d,(144.85,178),"Do not submit this placeholder",8,True,fill="#B54038",anchor="ma")
    draw_pil_text(d,(4,193),"G",12,True); draw_pil_mechanism(page,(8,193,194,78)); draw_pil_text(d,(202,286),"Figure 6 review proof v01 - not submission ready",8,fill=COLORS["gray"],anchor="ra")
    png=OUT_RASTER/"Figure6_review_v01_600dpi.png"; tif=OUT_RASTER/"Figure6_review_v01_600dpi.tiff"; page.save(png,dpi=(DPI,DPI),optimize=True); page.save(tif,dpi=(DPI,DPI),compression="tiff_lzw")
    return png,tif


SUPP_PANELS = [
    ("A","Ortholog-mapping workflow","mapping table + code"),("B","Cross-species co-embedding UMAP","integrated object + embeddings"),("C","Conserved-marker dot plot","expression matrix + marker list"),("D","Conserved-program heat map","scaled matrix + gene order"),
    ("E","Module-score comparisons","per-cell scores + exact tests"),("F","Dual-recombinase logic schematic","construct map + genotype notation"),("G","Vehicle-treated wound control","original image + scale metadata"),("H","Induced unwounded control","original image + scale metadata"),
    ("I","Induced wound, Epi-Front","original image + channels"),("J","Confocal orthogonal Z-stack","raw stack + projection settings"),("K","Reporter-positive quantification","animal-level table + exact test"),("L","Regional localization quantification","animal/section table + exact test"),
]


def build_supp_blueprint_pdf() -> Path:
    out=SUPP/"FigureS6_layout_blueprint_v01.pdf"; c=canvas.Canvas(str(out),pagesize=A4)
    margin=8; gap=4; cols=3; rows=4; bw=(PAGE_W_MM-2*margin-gap*(cols-1))/cols; bh=(PAGE_H_MM-2*margin-gap*(rows-1))/rows
    for i,(letter,title,need) in enumerate(SUPP_PANELS):
        row,col=divmod(i,cols); x=margin+col*(bw+gap); y=margin+row*(bh+gap)
        c.setFillColor(HexColor(COLORS["light"])); c.setStrokeColor(HexColor(COLORS["border"])); c.setDash(mm(2),mm(1.5)); c.roundRect(mm(x),pdf_y(y,bh),mm(bw),mm(bh),mm(1.3),stroke=1,fill=1); c.setDash()
        draw_pdf_text(c,x+2,y+2,letter,12,"Arial-Bold"); draw_pdf_text(c,x+bw/2,y+10,title,9,"Arial-Bold",HexColor(COLORS["navy"]),"middle")
        draw_pdf_text(c,x+bw/2,y+29,"SOURCE UNRESOLVED",9,"Arial-Bold",HexColor("#B54038"),"middle"); draw_pdf_text(c,x+bw/2,y+40,need,8,"Arial",HexColor(COLORS["gray"]),"middle")
        draw_pdf_text(c,x+bw/2,y+bh-9,"Panel identity locked to supplied legend",8,"Arial",HexColor(COLORS["gray"]),"middle")
    draw_pdf_text(c,PAGE_W_MM-8,PAGE_H_MM-5,"Figure S6 A4 layout blueprint v01 - not submission ready",8,"Arial",HexColor(COLORS["gray"]),"right")
    c.showPage(); c.save(); return out


def build_supp_blueprint_svg() -> Path:
    margin=8; gap=4; cols=3; rows=4; bw=(PAGE_W_MM-2*margin-gap*(cols-1))/cols; bh=(PAGE_H_MM-2*margin-gap*(rows-1))/rows
    q=['<svg xmlns="http://www.w3.org/2000/svg" width="210mm" height="297mm" viewBox="0 0 210 297"><rect width="210" height="297" fill="white"/>']
    for i,(letter,title,need) in enumerate(SUPP_PANELS):
        row,col=divmod(i,cols); x=margin+col*(bw+gap); y=margin+row*(bh+gap)
        q.append(f'<g id="Panel_{letter}"><rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="1.3" fill="{COLORS["light"]}" stroke="{COLORS["border"]}" stroke-width="0.4" stroke-dasharray="2,1.5"/>')
        q.append(svg_text(x+2,y+5,letter,12,"bold")); q.append(svg_text(x+bw/2,y+13,title,9,"bold","middle",fill=COLORS["navy"])); q.append(svg_text(x+bw/2,y+32,"SOURCE UNRESOLVED",9,"bold","middle",fill="#B54038")); q.append(svg_text(x+bw/2,y+43,need,8,anchor="middle",fill=COLORS["gray"])); q.append(svg_text(x+bw/2,y+bh-6,"Panel identity locked to supplied legend",8,anchor="middle",fill=COLORS["gray"])); q.append('</g>')
    q.append(svg_text(202,292,"Figure S6 A4 layout blueprint v01 - not submission ready",8,anchor="end",fill=COLORS["gray"])); q.append('</svg>')
    out=SUPP/"FigureS6_layout_blueprint_v01_editable.svg"; out.write_text(''.join(q),encoding='utf-8'); return out


def build_supp_blueprint_raster() -> tuple[Path, Path]:
    page=Image.new("RGB",(p(PAGE_W_MM),p(PAGE_H_MM)),"white"); d=ImageDraw.Draw(page)
    margin=8; gap=4; cols=3; rows=4; bw=(PAGE_W_MM-2*margin-gap*(cols-1))/cols; bh=(PAGE_H_MM-2*margin-gap*(rows-1))/rows
    for i,(letter,title,need) in enumerate(SUPP_PANELS):
        row,col=divmod(i,cols); x=margin+col*(bw+gap); y=margin+row*(bh+gap)
        d.rounded_rectangle((p(x),p(y),p(x+bw),p(y+bh)),radius=p(1.3),fill=COLORS["light"],outline=COLORS["border"],width=p(.35))
        draw_pil_text(d,(x+2,y+2),letter,12,True)
        draw_pil_text(d,(x+bw/2,y+10),title,9,True,fill=COLORS["navy"],anchor="ma")
        draw_pil_text(d,(x+bw/2,y+29),"SOURCE UNRESOLVED",9,True,fill="#B54038",anchor="ma")
        draw_pil_text(d,(x+bw/2,y+40),need,8,fill=COLORS["gray"],anchor="ma")
        draw_pil_text(d,(x+bw/2,y+bh-9),"Panel identity locked to supplied legend",8,fill=COLORS["gray"],anchor="ma")
    draw_pil_text(d,(PAGE_W_MM-8,PAGE_H_MM-5),"Figure S6 A4 layout blueprint v01 - not submission ready",8,fill=COLORS["gray"],anchor="ra")
    png=SUPP/"FigureS6_layout_blueprint_v01_600dpi.png"; tif=SUPP/"FigureS6_layout_blueprint_v01_600dpi.tiff"
    page.save(png,dpi=(DPI,DPI),optimize=True); page.save(tif,dpi=(DPI,DPI),compression="tiff_lzw")
    return png,tif


def main():
    ensure_dirs()
    svg = build_main_svg()
    build_g_panel_svg()
    pdf = build_main_pdf()
    png, tif = build_main_raster()
    build_supp_blueprint_svg(); build_supp_blueprint_pdf(); build_supp_blueprint_raster()
    for src in (svg, pdf, png, tif):
        shutil.copyfile(src, OUT_ASSEMBLED / src.name)
    hash_lines = []
    for path in sorted(INPUT.glob("*.png")):
        hash_lines.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}")
    (ROOT / "QC" / "source_hashes_sha256.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
    print("Figure 6 review proof and Figure S6 blueprint built.")


if __name__ == "__main__":
    main()
