from __future__ import annotations

import csv
from pathlib import Path

import tifffile
from PIL import Image
from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\BW3_reprocess_20260820"
)
ANALYSIS = ROOT / "03_imagej_double_positive_spatial"
SOURCE_DIR = ROOT / "04_editable_three_panel_AI" / "01_source_pdf"
PDF_PATH = SOURCE_DIR / "Figure3G_BW3_vertical_3panel_editable.pdf"
MAP_PATH = SOURCE_DIR / "ai_source_map.csv"

A4_WIDTH_PT = 595.275590551
FIGURE_WIDTH_PT = A4_WIDTH_PT / 4.0
SOURCE_PIXEL_WIDTH = 6656
SOURCE_PIXEL_HEIGHT = 4467
PIXEL_SIZE_UM = 1.0113780498504639
SCALE_BAR_UM = 1000.0
UNITS_PER_PT = SOURCE_PIXEL_WIDTH / FIGURE_WIDTH_PT
IMAGE_HEIGHT_PT = SOURCE_PIXEL_HEIGHT / UNITS_PER_PT
HEADER_HEIGHT_PT = 24.0
GAP_HEIGHT_PT = 4.0
FIGURE_HEIGHT_PT = 3 * (HEADER_HEIGHT_PT + IMAGE_HEIGHT_PT) + 2 * GAP_HEIGHT_PT

PANELS = [
    {
        "source": ANALYSIS / "02_recolored_merges" / "KC_KRT6B_red_SOX9_green_DAPI_blue_RGB.tif",
        "png": SOURCE_DIR / "panel_01_landscape.png",
        "title": "KRT6B+SOX9+",
        "channels": [("KRT6B", "#ff3038"), ("SOX9", "#20c968"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "source": ANALYSIS / "02_recolored_merges" / "Fib_KRT14_red_PDGFRa_green_DAPI_blue_RGB.tif",
        "png": SOURCE_DIR / "panel_02_landscape.png",
        "title": "PDGFRα+KRT14+",
        "channels": [("PDGFRα", "#20c968"), ("KRT14", "#ff3038"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "source": ANALYSIS / "04_spatial_correspondence" / "SerialSection_double_positive_cell_centers_RGB.tif",
        "png": SOURCE_DIR / "panel_03_landscape.png",
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
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    for panel in PANELS:
        array = tifffile.imread(panel["source"])
        if tuple(array.shape[:2]) != (SOURCE_PIXEL_HEIGHT, SOURCE_PIXEL_WIDTH):
            raise ValueError(f"Unexpected panel shape {array.shape}: {panel['source']}")
        Image.fromarray(array).save(panel["png"], compress_level=4)

    pdfmetrics.registerFont(TTFont("Arial", r"data/external/local_16a706108b0a"))
    pdfmetrics.registerFont(TTFont("Arial-Bold", r"data/external/local_28ac8473f586"))
    c = canvas.Canvas(str(PDF_PATH), pagesize=(FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT), pageCompression=1)
    c.setTitle("Figure 3G BW3 registered serial-section three-panel figure")
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
        c.setFont("Arial-Bold", panel["legend_size"])
        channel_x = 2.5
        for channel_name, channel_color in panel["channels"]:
            c.setFillColor(HexColor(channel_color))
            c.drawString(channel_x, top_y - 21.0, channel_name)
            channel_x += pdfmetrics.stringWidth(channel_name, "Arial-Bold", panel["legend_size"]) + 4.5
        c.drawImage(str(panel["png"]), 0, image_bottom, FIGURE_WIDTH_PT, IMAGE_HEIGHT_PT, mask=None)
        bar_x, bar_y = 4.0, image_bottom + 5.0
        label = "1000 µm"
        label_width = pdfmetrics.stringWidth(label, "Arial", 7)
        c.setFillColor(Color(0, 0, 0, alpha=0.64))
        c.roundRect(bar_x - 1.2, bar_y - 2.0, max(scale_bar_pt, label_width) + 2.4, 13.0, 0.6, fill=1, stroke=0)
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

    with MAP_PATH.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Unit", "Source", "Status"])
        writer.writeheader()
        writer.writerow({"Unit": "Figure3G_BW3_vertical_3panel", "Source": str(PDF_PATH), "Status": "pass"})
    (SOURCE_DIR / "layout_and_orientation_audit.txt").write_text(
        "status=pass\n"
        "orientation_action=no additional rotation; registered BW3 panels are already landscape\n"
        "reason=preserves horizontal tissue presentation matching the supplied reference layout\n"
        f"source_pixel_size_um={PIXEL_SIZE_UM}\nscale_bar_um={SCALE_BAR_UM}\n"
        f"figure_width_pt={FIGURE_WIDTH_PT:.6f}\nfigure_width_mm=52.5\nfigure_height_pt={FIGURE_HEIGHT_PT:.6f}\n"
        "font=Arial; title=12 pt bold; channel labels=7-8 pt bold; scale label=7 pt regular\n"
        "protein_labels=upright bold; no gene-symbol italics used\n"
        "serial_section_boundary=spatial correspondence only; not same-cell molecular colocalization\n",
        encoding="utf-8",
    )
    print(f"PDF_PASS {PDF_PATH}")
    print(f"SIZE_PT={FIGURE_WIDTH_PT:.6f}x{FIGURE_HEIGHT_PT:.6f}; SCALE_BAR_PT={scale_bar_pt:.6f}")


if __name__ == "__main__":
    main()
