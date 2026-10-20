from __future__ import annotations

from pathlib import Path

from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\ImageJ_double_positive_spatial"
)
SOURCE_DIR = ROOT / "07_vertical_three_panel_AI_source"
PDF_PATH = SOURCE_DIR / "KRT6B_SOX9_PDGFRa_KRT14_vertical_3panel_editable.pdf"

A4_WIDTH_PT = 595.275590551
FIGURE_WIDTH_PT = A4_WIDTH_PT / 4.0
SOURCE_PIXEL_WIDTH = 10865
SOURCE_PIXEL_HEIGHT = 3477
PIXEL_SIZE_UM = 1.0113780498504639
SCALE_BAR_UM = 1000.0
UNITS_PER_PT = SOURCE_PIXEL_WIDTH / FIGURE_WIDTH_PT
IMAGE_HEIGHT_PT = SOURCE_PIXEL_HEIGHT / UNITS_PER_PT
HEADER_HEIGHT_PT = 24.0
GAP_HEIGHT_PT = 4.0
FIGURE_HEIGHT_PT = 3 * (HEADER_HEIGHT_PT + IMAGE_HEIGHT_PT) + 2 * GAP_HEIGHT_PT

PANELS = [
    {
        "image": SOURCE_DIR / "panel_01_CCW90.png",
        "title": "KRT6B+SOX9+",
        "channels": [("KRT6B", "#ff3038"), ("SOX9", "#20c968"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "image": SOURCE_DIR / "panel_02_CCW90.png",
        "title": "PDGFRα+KRT14+",
        "channels": [("PDGFRα", "#20c968"), ("KRT14", "#ff3038"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "image": SOURCE_DIR / "panel_03_CCW90.png",
        "title": "Spatial correspondence",
        "channels": [
            ("KRT6B+SOX9+", "#ff3038"),
            ("PDGFRα+KRT14+", "#20c968"),
            ("DAPI", "#5475ff"),
        ],
        "legend_size": 7.0,
    },
]


def main() -> None:
    pdfmetrics.registerFont(TTFont("Arial", r"data/external/local_16a706108b0a"))
    pdfmetrics.registerFont(TTFont("Arial-Bold", r"data/external/local_28ac8473f586"))
    c = canvas.Canvas(str(PDF_PATH), pagesize=(FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT), pageCompression=1)
    c.setTitle("Figure 3G vertical three-panel immunofluorescence")
    c.setAuthor("Fiji/ImageJ registered serial-section analysis")
    c.setFillColor(black)
    c.rect(0, 0, FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT, fill=1, stroke=0)

    top_y = FIGURE_HEIGHT_PT
    scale_bar_pt = (SCALE_BAR_UM / PIXEL_SIZE_UM) / UNITS_PER_PT
    for index, panel in enumerate(PANELS, start=1):
        header_bottom = top_y - HEADER_HEIGHT_PT
        image_bottom = header_bottom - IMAGE_HEIGHT_PT
        c.setFillColor(black)
        c.rect(0, header_bottom, FIGURE_WIDTH_PT, HEADER_HEIGHT_PT, fill=1, stroke=0)

        c.setFont("Arial-Bold", 12)
        c.setFillColor(white)
        c.drawString(2.5, top_y - 10.8, panel["title"])

        legend_size = panel["legend_size"]
        c.setFont("Arial-Bold", legend_size)
        channel_x = 2.5
        for channel_name, channel_color in panel["channels"]:
            c.setFillColor(HexColor(channel_color))
            c.drawString(channel_x, top_y - 21.0, channel_name)
            channel_x += pdfmetrics.stringWidth(channel_name, "Arial-Bold", legend_size) + 5.0

        c.drawImage(
            str(panel["image"]),
            0,
            image_bottom,
            width=FIGURE_WIDTH_PT,
            height=IMAGE_HEIGHT_PT,
            preserveAspectRatio=False,
            mask=None,
        )

        bar_x = 4.0
        bar_y = image_bottom + 5.0
        label = "1000 µm"
        label_width = pdfmetrics.stringWidth(label, "Arial", 7)
        background_width = max(scale_bar_pt, label_width) + 2.4
        c.setFillColor(Color(0, 0, 0, alpha=0.58))
        c.roundRect(bar_x - 1.2, bar_y - 2.0, background_width, 13.0, 0.6, fill=1, stroke=0)
        c.setStrokeColor(white)
        c.setLineWidth(0.75)
        c.line(bar_x, bar_y, bar_x + scale_bar_pt, bar_y)
        c.setFillColor(white)
        c.setFont("Arial", 7)
        c.drawString(bar_x, bar_y + 2.0, label)

        c.setStrokeColor(HexColor("#5c5c5c"))
        c.setLineWidth(0.35)
        c.rect(0.175, image_bottom + 0.175, FIGURE_WIDTH_PT - 0.35, IMAGE_HEIGHT_PT - 0.35, fill=0, stroke=1)

        top_y = image_bottom
        if index < len(PANELS):
            c.setFillColor(black)
            c.rect(0, top_y - GAP_HEIGHT_PT, FIGURE_WIDTH_PT, GAP_HEIGHT_PT, fill=1, stroke=0)
            top_y -= GAP_HEIGHT_PT

    c.showPage()
    c.save()
    (SOURCE_DIR / "pdf_layout_spec.txt").write_text(
        "A4 width=595.275590551 pt\n"
        f"figure width=A4/4={FIGURE_WIDTH_PT:.6f} pt (52.5 mm)\n"
        f"figure height={FIGURE_HEIGHT_PT:.6f} pt\n"
        "font=Arial; titles=12 pt bold; channel legends=7-8 pt bold; scale text=7 pt regular\n"
        "protein labels=upright bold; no gene-symbol italics used\n"
        f"pixel size={PIXEL_SIZE_UM} um/pixel; scale bar={SCALE_BAR_UM} um\n"
        "microscopy panels=embedded raster evidence; all added text, rules, borders, and scale bars=PDF native objects\n",
        encoding="utf-8",
    )
    print(f"PDF_PASS {PDF_PATH}")
    print(f"SIZE_PT={FIGURE_WIDTH_PT:.6f}x{FIGURE_HEIGHT_PT:.6f}")


if __name__ == "__main__":
    main()
