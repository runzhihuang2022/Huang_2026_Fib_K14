from __future__ import annotations

import csv
import math
from pathlib import Path

from PIL import Image
from reportlab.lib.colors import Color, HexColor, black, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


BASE = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap"
)
ROOT = BASE / "08_editable_BW3_BW9_3x2_AI"
SOURCE_DIR = ROOT / "01_source_pdf"
CONVERTED_DIR = ROOT / "02_converted_ai"
QC_DIR = CONVERTED_DIR / "QC"
PREVIEW_DIR = QC_DIR / "previews"
DELIVERY_DIR = ROOT / "03_review_delivery"
for directory in (SOURCE_DIR, CONVERTED_DIR, QC_DIR, PREVIEW_DIR, DELIVERY_DIR):
    directory.mkdir(parents=True, exist_ok=True)

STEM = "Figure3G_BW3_BW9_3x2_full_tissue_no_crop_PDGFRA"
PDF_PATH = SOURCE_DIR / f"{STEM}_editable.pdf"
MAP_PATH = SOURCE_DIR / "ai_source_map.csv"

# Preserve the previous publication layout: each sample column is A4 width / 4.
# The two-column figure therefore occupies A4 width / 2.
A4_WIDTH_PT = 595.275590551
COLUMN_WIDTH_PT = A4_WIDTH_PT / 4.0
FIGURE_WIDTH_PT = 2.0 * COLUMN_WIDTH_PT
IMAGE_FRAME_HEIGHT_PT = COLUMN_WIDTH_PT * 4467.0 / 6656.0
PANEL_HEADER_HEIGHT_PT = 24.0
PANEL_HEIGHT_PT = PANEL_HEADER_HEIGHT_PT + IMAGE_FRAME_HEIGHT_PT
COLUMN_LABEL_HEIGHT_PT = 16.0
ROW_GAP_PT = 4.0
FIGURE_HEIGHT_PT = COLUMN_LABEL_HEIGHT_PT + 3.0 * PANEL_HEIGHT_PT + 2.0 * ROW_GAP_PT
PIXEL_SIZE_UM = 1.0113780498504639
SCALE_BAR_UM = 1000.0

# Clockwise leveling rotation. It is applied as a PDF transform to the complete
# source raster (expand-equivalent bounding box), so no pixels are resampled or cropped.
BW3_ROTATION_DEG = -14.0

ROW_SPECS = [
    {
        "title": "KRT6B+SOX9+",
        "channels": [("KRT6B", "#ff3038"), ("SOX9", "#20c968"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "title": "PDGFRA+KRT14+",
        "channels": [("PDGFRA", "#20c968"), ("KRT14", "#ff3038"), ("DAPI", "#5475ff")],
        "legend_size": 8.0,
    },
    {
        "title": "Spatial correspondence",
        "channels": [
            ("KRT6B+SOX9+", "#ff3038"),
            ("PDGFRA+KRT14+", "#20c968"),
            ("DAPI", "#5475ff"),
        ],
        "legend_size": 7.0,
    },
]

SLIDES = [
    {
        "name": "BW3",
        "images": [
            BASE / "BW3_reprocess_20260820" / "04_editable_three_panel_AI" / "01_source_pdf" / f"panel_0{i}_landscape.png"
            for i in (1, 2, 3)
        ],
        "rotation_deg": BW3_ROTATION_DEG,
        "orientation": "clockwise 14 degrees; complete raster retained by transform; no cropping",
    },
    {
        "name": "BW9",
        "images": [
            BASE / "ImageJ_double_positive_spatial" / "07_vertical_three_panel_AI_source" / f"panel_0{i}_CCW90.png"
            for i in (1, 2, 3)
        ],
        "rotation_deg": 0.0,
        "orientation": "audited BW9 CCW90 source retained; no additional rotation",
    },
]


def rotated_fit(source_width: int, source_height: int, angle_deg: float) -> tuple[float, float, float, float, float]:
    angle = math.radians(angle_deg)
    c_abs, s_abs = abs(math.cos(angle)), abs(math.sin(angle))
    bbox_width = source_width * c_abs + source_height * s_abs
    bbox_height = source_width * s_abs + source_height * c_abs
    scale = min(COLUMN_WIDTH_PT / bbox_width, IMAGE_FRAME_HEIGHT_PT / bbox_height)
    return source_width * scale, source_height * scale, bbox_width * scale, bbox_height * scale, scale


def draw_complete_rotated_image(
    c: canvas.Canvas,
    path: Path,
    source_width: int,
    source_height: int,
    frame_x: float,
    frame_y: float,
    angle_deg: float,
) -> tuple[float, float, float, float, float]:
    draw_width, draw_height, bbox_width, bbox_height, scale = rotated_fit(
        source_width, source_height, angle_deg
    )
    center_x = frame_x + COLUMN_WIDTH_PT / 2.0
    center_y = frame_y + IMAGE_FRAME_HEIGHT_PT / 2.0
    c.saveState()
    c.translate(center_x, center_y)
    c.rotate(angle_deg)
    c.drawImage(
        str(path),
        -draw_width / 2.0,
        -draw_height / 2.0,
        width=draw_width,
        height=draw_height,
        preserveAspectRatio=False,
        mask=None,
    )
    c.restoreState()
    bbox_x = center_x - bbox_width / 2.0
    bbox_y = center_y - bbox_height / 2.0
    return bbox_x, bbox_y, bbox_width, bbox_height, scale


def main() -> None:
    dimensions: dict[str, tuple[int, int]] = {}
    for slide in SLIDES:
        for path in slide["images"]:
            if not path.exists():
                raise FileNotFoundError(path)
            with Image.open(path) as image:
                dimensions[str(path)] = image.size

    pdfmetrics.registerFont(TTFont("Arial", r"data/external/local_16a706108b0a"))
    pdfmetrics.registerFont(TTFont("Arial-Bold", r"data/external/local_28ac8473f586"))
    c = canvas.Canvas(str(PDF_PATH), pagesize=(FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT), pageCompression=1)
    c.setTitle("Figure 3G BW3 and BW9 3 by 2 whole-tissue matrix")
    c.setAuthor("Fiji/ImageJ registered serial-section analysis")
    c.setFillColor(black)
    c.rect(0, 0, FIGURE_WIDTH_PT, FIGURE_HEIGHT_PT, fill=1, stroke=0)

    fit_records: list[dict[str, object]] = []
    for column_index, slide in enumerate(SLIDES):
        x0 = column_index * COLUMN_WIDTH_PT
        c.setFillColor(white)
        c.setFont("Arial-Bold", 12)
        sample_width = pdfmetrics.stringWidth(slide["name"], "Arial-Bold", 12)
        c.drawString(x0 + (COLUMN_WIDTH_PT - sample_width) / 2.0, FIGURE_HEIGHT_PT - 11.5, slide["name"])

        for row_index, (row_spec, image_path) in enumerate(zip(ROW_SPECS, slide["images"])):
            panel_top = FIGURE_HEIGHT_PT - COLUMN_LABEL_HEIGHT_PT - row_index * (PANEL_HEIGHT_PT + ROW_GAP_PT)
            header_bottom = panel_top - PANEL_HEADER_HEIGHT_PT
            frame_bottom = header_bottom - IMAGE_FRAME_HEIGHT_PT
            c.setFillColor(black)
            c.rect(x0, frame_bottom, COLUMN_WIDTH_PT, PANEL_HEIGHT_PT, fill=1, stroke=0)

            c.setFillColor(white)
            c.setFont("Arial-Bold", 12)
            c.drawString(x0 + 2.5, panel_top - 10.8, row_spec["title"])
            c.setFont("Arial-Bold", row_spec["legend_size"])
            legend_x = x0 + 2.5
            for channel_name, channel_color in row_spec["channels"]:
                c.setFillColor(HexColor(channel_color))
                c.drawString(legend_x, panel_top - 21.0, channel_name)
                legend_x += pdfmetrics.stringWidth(channel_name, "Arial-Bold", row_spec["legend_size"]) + 4.5

            source_width, source_height = dimensions[str(image_path)]
            bbox_x, bbox_y, bbox_width, bbox_height, display_scale = draw_complete_rotated_image(
                c,
                image_path,
                source_width,
                source_height,
                x0,
                frame_bottom,
                float(slide["rotation_deg"]),
            )
            fit_records.append(
                {
                    "sample": slide["name"],
                    "row": row_index + 1,
                    "source": str(image_path),
                    "source_px": f"{source_width}x{source_height}",
                    "rotation_deg": slide["rotation_deg"],
                    "bbox_pt": f"{bbox_width:.6f}x{bbox_height:.6f}",
                    "frame_pt": f"{COLUMN_WIDTH_PT:.6f}x{IMAGE_FRAME_HEIGHT_PT:.6f}",
                    "cropped": "false",
                }
            )

            scale_bar_pt = (SCALE_BAR_UM / PIXEL_SIZE_UM) * display_scale
            bar_x, bar_y = x0 + 4.0, frame_bottom + 5.0
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

            c.setStrokeColor(HexColor("#666666"))
            c.setLineWidth(0.35)
            c.rect(x0 + 0.175, frame_bottom + 0.175, COLUMN_WIDTH_PT - 0.35, PANEL_HEIGHT_PT - 0.35, fill=0, stroke=1)

    c.showPage()
    c.save()

    with MAP_PATH.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Unit", "Source", "Status"])
        writer.writeheader()
        writer.writerow({"Unit": "Figure3G_BW3_BW9_3x2", "Source": str(PDF_PATH), "Status": "review_ready"})

    with (SOURCE_DIR / "image_fit_audit.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fit_records[0].keys()))
        writer.writeheader()
        writer.writerows(fit_records)

    audit_lines = [
        "status=review_ready",
        "matrix=2 columns x 3 rows; 6 whole-tissue images",
        "columns=BW3|BW9",
        "rows=KRT6B+SOX9+|PDGFRA+KRT14+|Spatial correspondence",
        "zoom_in=none",
        f"column_width_pt={COLUMN_WIDTH_PT:.6f}",
        f"image_frame_width_pt={COLUMN_WIDTH_PT:.6f}",
        f"image_frame_height_pt={IMAGE_FRAME_HEIGHT_PT:.6f}",
        f"figure_width_pt={FIGURE_WIDTH_PT:.6f}",
        f"figure_height_pt={FIGURE_HEIGHT_PT:.6f}",
        "fit_method=complete source raster transformed and centered inside identical frame; no stretching; no clipping; no cropping",
        f"BW3_rotation_deg={BW3_ROTATION_DEG}",
        "marker_name=PDGFRA (gene-symbol replacement for PDGFRalpha; receptor identity retained)",
        "font=Arial; sample labels and panel titles=12 pt bold; channel labels=7-8 pt bold; scale labels=7 pt regular",
        f"scale_bar_um={SCALE_BAR_UM}",
        f"pixel_size_um={PIXEL_SIZE_UM}",
        "serial_section_boundary=spatial correspondence only; not same-cell molecular colocalization",
    ]
    for slide in SLIDES:
        audit_lines.append(f"orientation_{slide['name']}={slide['orientation']}")
        for row_index, path in enumerate(slide["images"], start=1):
            audit_lines.append(f"source_{slide['name']}_row{row_index}={path}")
    (SOURCE_DIR / "3x2_layout_and_source_audit.txt").write_text("\n".join(audit_lines) + "\n", encoding="utf-8")
    print(f"PDF_PASS {PDF_PATH}")
    print(f"FIGURE_PT={FIGURE_WIDTH_PT:.6f}x{FIGURE_HEIGHT_PT:.6f}; PANELS=6")


if __name__ == "__main__":
    main()
