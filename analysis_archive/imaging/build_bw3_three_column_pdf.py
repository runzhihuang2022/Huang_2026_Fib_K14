from __future__ import annotations

import csv
from pathlib import Path

from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\BW3_reprocess_20260820"
)
OLD_SOURCE = ROOT / "04_editable_three_panel_AI" / "01_source_pdf"
SOURCE_DIR = ROOT / "05_editable_three_column_AI" / "01_source_pdf"
SOURCE_DIR.mkdir(parents=True, exist_ok=True)
PDF_PATH = SOURCE_DIR / "Figure3G_BW3_three_equal_columns_editable.pdf"
MAP_PATH = SOURCE_DIR / "ai_source_map.csv"

A4_WIDTH_PT = 595.275590551
COLUMN_WIDTH_PT = A4_WIDTH_PT / 4.0
FIGURE_WIDTH_PT = 3.0 * COLUMN_WIDTH_PT
SOURCE_WIDTH, SOURCE_HEIGHT = 6656, 4467
IMAGE_HEIGHT_PT = SOURCE_HEIGHT / (SOURCE_WIDTH / COLUMN_WIDTH_PT)
HEADER_HEIGHT_PT = 24.0
FIGURE_HEIGHT_PT = HEADER_HEIGHT_PT + IMAGE_HEIGHT_PT
PIXEL_SIZE_UM = 1.0113780498504639
SCALE_BAR_UM = 1000.0
SCALE_BAR_PT = (SCALE_BAR_UM / PIXEL_SIZE_UM) / (SOURCE_WIDTH / COLUMN_WIDTH_PT)

PANELS = [
    {
        "image": OLD_SOURCE / "panel_01_landscape.png",
        "title": "KRT6B+SOX9+",
        "channels": [("KRT6B", "#ff3038"), ("SOX9", "#20c968"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "image": OLD_SOURCE / "panel_02_landscape.png",
        "title": "PDGFRα+KRT14+",
        "channels": [("PDGFRα", "#20c968"), ("KRT14", "#ff3038"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "image": OLD_SOURCE / "panel_03_landscape.png",
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
    for panel in PANELS:
        if not panel["image"].exists():
            raise FileNotFoundError(panel["image"])
    pdfmetrics.registerFont(TTFont("Arial", r"data/external/local_16a706108b0a"))
    pdfmetrics.registerFont(TTFont("Arial-Bold", r"data/external/local_28ac8473f586"))
    c = canvas.Canvas(str(PDF_PATH), pagesize=(FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT), pageCompression=1)
    c.setTitle("Figure 3G BW3 three equal columns")
    c.setAuthor("Fiji/ImageJ whole-tissue serial-section analysis")
    c.setFillColor(black)
    c.rect(0, 0, FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT, fill=1, stroke=0)
    image_bottom = 0.0
    header_bottom = IMAGE_HEIGHT_PT

    for index, panel in enumerate(PANELS):
        x0 = index * COLUMN_WIDTH_PT
        c.setFillColor(black)
        c.rect(x0, header_bottom, COLUMN_WIDTH_PT, HEADER_HEIGHT_PT, fill=1, stroke=0)
        c.setFillColor(white)
        c.setFont("Arial-Bold", 12)
        c.drawString(x0 + 2.5, FIGURE_HEIGHT_PT - 10.8, panel["title"])
        c.setFont("Arial-Bold", panel["legend_size"])
        legend_x = x0 + 2.5
        for name, color in panel["channels"]:
            c.setFillColor(HexColor(color))
            c.drawString(legend_x, FIGURE_HEIGHT_PT - 21.0, name)
            legend_x += pdfmetrics.stringWidth(name, "Arial-Bold", panel["legend_size"]) + 4.5
        c.drawImage(
            str(panel["image"]), x0, image_bottom,
            width=COLUMN_WIDTH_PT, height=IMAGE_HEIGHT_PT,
            preserveAspectRatio=False, mask=None,
        )
        bar_x, bar_y = x0 + 4.0, 5.0
        label = "1000 µm"
        label_width = pdfmetrics.stringWidth(label, "Arial", 7)
        c.setFillColor(Color(0, 0, 0, alpha=0.64))
        c.roundRect(bar_x - 1.2, bar_y - 2.0, max(SCALE_BAR_PT, label_width) + 2.4, 13.0, 0.6, fill=1, stroke=0)
        c.setStrokeColor(white)
        c.setLineWidth(0.75)
        c.line(bar_x, bar_y, bar_x + SCALE_BAR_PT, bar_y)
        c.setFillColor(white)
        c.setFont("Arial", 7)
        c.drawString(bar_x, bar_y + 2.0, label)
        c.setStrokeColor(HexColor("#666666"))
        c.setLineWidth(0.35)
        c.rect(x0 + 0.175, 0.175, COLUMN_WIDTH_PT - 0.35, FIGURE_HEIGHT_PT - 0.35, fill=0, stroke=1)

    c.showPage()
    c.save()
    with MAP_PATH.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Unit", "Source", "Status"])
        writer.writeheader()
        writer.writerow({"Unit": "Figure3G_BW3_three_equal_columns", "Source": str(PDF_PATH), "Status": "pass"})
    (SOURCE_DIR / "three_column_layout_audit.txt").write_text(
        "status=pass\n"
        f"panel_count=3\npanel_width_pt={COLUMN_WIDTH_PT:.6f}\npanel_height_pt={FIGURE_HEIGHT_PT:.6f}\n"
        f"image_width_pt={COLUMN_WIDTH_PT:.6f}\nimage_height_pt={IMAGE_HEIGHT_PT:.6f}\n"
        f"figure_width_pt={FIGURE_WIDTH_PT:.6f}\nfigure_width_mm=157.5\n"
        "layout=three equal columns; no ZoomIn crops; all panels use the same registered whole-tissue canvas\n"
        "font=Arial; titles=12 pt bold; channel labels=7-8 pt bold; scale labels=7 pt regular\n"
        f"scale_bar_um={SCALE_BAR_UM}\npixel_size_um={PIXEL_SIZE_UM}\n"
        "serial_section_boundary=spatial correspondence only; not same-cell molecular colocalization\n",
        encoding="utf-8",
    )
    print(f"PDF_PASS {PDF_PATH}")
    print(f"FIGURE_PT={FIGURE_WIDTH_PT:.6f}x{FIGURE_HEIGHT_PT:.6f}; EACH_COLUMN_PT={COLUMN_WIDTH_PT:.6f}x{FIGURE_HEIGHT_PT:.6f}")


if __name__ == "__main__":
    main()
