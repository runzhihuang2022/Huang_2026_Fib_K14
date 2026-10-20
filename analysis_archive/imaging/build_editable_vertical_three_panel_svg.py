from __future__ import annotations

import html
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\ImageJ_double_positive_spatial"
)
SOURCE_DIR = ROOT / "07_vertical_three_panel_AI_source"
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

PIXEL_SIZE_UM = 1.0113780498504639
SCALE_BAR_UM = 1000.0
A4_WIDTH_PT = 595.275590551
FIGURE_WIDTH_PT = A4_WIDTH_PT / 4.0

PANELS = [
    {
        "path": ROOT / "02_recolored_merges" / "KC_KRT6B_red_SOX9_green_DAPI_blue_RGB.tif",
        "title": "KRT6B+SOX9+",
        "channels": [("KRT6B", "#ff3038"), ("SOX9", "#20c968"), ("DAPI", "#5475ff")],
        "legend_size_pt": 8.0,
    },
    {
        "path": ROOT / "02_recolored_merges" / "Fib_KRT14_red_PDGFRa_green_DAPI_blue_RGB.tif",
        "title": "PDGFRα+KRT14+",
        "channels": [("PDGFRα", "#20c968"), ("KRT14", "#ff3038"), ("DAPI", "#5475ff")],
        "legend_size_pt": 8.0,
    },
    {
        "path": ROOT / "04_spatial_correspondence" / "SerialSection_double_positive_cell_centers_RGB.tif",
        "title": "Spatial correspondence",
        "channels": [
            ("KRT6B+SOX9+", "#ff3038"),
            ("PDGFRα+KRT14+", "#20c968"),
            ("DAPI", "#5475ff"),
        ],
        "legend_size_pt": 7.0,
    },
]


def rotated_png_file(path: Path, index: int) -> tuple[Path, int, int]:
    array = tifffile.imread(path)
    if array.ndim != 3 or array.shape[-1] != 3:
        raise ValueError(f"Expected RGB TIFF: {path}, got {array.shape}")
    rotated = np.ascontiguousarray(np.rot90(array, k=1))
    image = Image.fromarray(rotated, mode="RGB")
    output_path = SOURCE_DIR / f"panel_{index:02d}_CCW90.png"
    image.save(output_path, format="PNG", compress_level=6)
    return output_path, image.width, image.height


def text_element(
    x: float,
    y: float,
    text: str,
    color: str,
    size: float,
    weight: str = "700",
    style: str = "normal",
) -> str:
    return (
        f'<text x="{x}" y="{y}" fill="{color}" font-family="Arial" '
        f'font-size="{size:.3f}" font-weight="{weight}" font-style="{style}">{html.escape(text)}</text>'
    )


def main() -> None:
    encoded = []
    for index, panel in enumerate(PANELS, start=1):
        png_path, width, height = rotated_png_file(panel["path"], index)
        encoded.append((panel, png_path, width, height))
        print(f"EMBED_READY {panel['path'].name}: {width}x{height}", flush=True)

    sizes = {(width, height) for _, _, width, height in encoded}
    if len(sizes) != 1:
        raise ValueError(f"Rotated panel sizes differ: {sizes}")
    panel_width, panel_height = next(iter(sizes))
    units_per_pt = panel_width / FIGURE_WIDTH_PT
    header_height = 24.0 * units_per_pt
    gap_height = 4.0 * units_per_pt
    title_size = 12.0 * units_per_pt
    scale_text_size = 7.0 * units_per_pt
    total_height = len(encoded) * (header_height + panel_height) + (len(encoded) - 1) * gap_height
    figure_height_pt = total_height / units_per_pt
    scale_bar_px = SCALE_BAR_UM / PIXEL_SIZE_UM

    svg = [
        '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{FIGURE_WIDTH_PT:.6f}pt" height="{figure_height_pt:.6f}pt" viewBox="0 0 {panel_width} {total_height:.3f}">',
        '<rect id="background" x="0" y="0" width="100%" height="100%" fill="#000000"/>',
        '<g id="editable_titles_legends_scale_bars">',
    ]

    y = 0
    manifest = [
        "order\ttitle\tsource\trotation\tpixel_size_um\tscale_bar_um\timage_width\timage_height"
    ]
    for index, (panel, png_path, width, height) in enumerate(encoded, start=1):
        image_y = y + header_height
        svg.append(f'<g id="panel_{index}">')
        svg.append(f'<rect x="0" y="{y:.3f}" width="{panel_width}" height="{header_height:.3f}" fill="#000000"/>')
        svg.append(text_element(2.5 * units_per_pt, y + 10.8 * units_per_pt, panel["title"], "#ffffff", title_size))

        channel_x = 2.5 * units_per_pt
        legend_size_pt = float(panel["legend_size_pt"])
        channel_size = legend_size_pt * units_per_pt
        for channel_name, channel_color in panel["channels"]:
            # Protein labels in immunofluorescence figures are upright bold; italics are reserved for gene symbols.
            svg.append(text_element(channel_x, y + 21.0 * units_per_pt, channel_name, channel_color, channel_size))
            channel_x += max(19.0 * units_per_pt, (0.50 * legend_size_pt * len(channel_name) + 5.0) * units_per_pt)

        svg.append(
            f'<image id="panel_{index}_microscopy_raster" x="0" y="{image_y}" width="{width}" height="{height}" '
            f'preserveAspectRatio="none" href="{html.escape(png_path.as_uri())}" '
            f'xlink:href="{html.escape(png_path.as_uri())}"/>'
        )

        bar_x = 4.0 * units_per_pt
        bar_y = image_y + height - 5.0 * units_per_pt
        svg.append(
            f'<rect id="panel_{index}_scale_background" x="{bar_x - 1.2 * units_per_pt:.3f}" '
            f'y="{bar_y - 11.0 * units_per_pt:.3f}" width="{scale_bar_px + 2.4 * units_per_pt:.3f}" '
            f'height="{13.0 * units_per_pt:.3f}" rx="{0.6 * units_per_pt:.3f}" fill="#000000" fill-opacity="0.58"/>'
        )
        svg.append(
            f'<line id="panel_{index}_scale_bar" x1="{bar_x}" y1="{bar_y}" '
            f'x2="{bar_x + scale_bar_px:.3f}" y2="{bar_y}" stroke="#ffffff" '
            f'stroke-width="{0.75 * units_per_pt:.3f}"/>'
        )
        svg.append(text_element(bar_x, bar_y - 2.0 * units_per_pt, "1000 µm", "#ffffff", scale_text_size, weight="400"))
        svg.append(
            f'<rect id="panel_{index}_border" x="4" y="{image_y + 4}" width="{width - 8}" height="{height - 8}" '
            f'fill="none" stroke="#5c5c5c" stroke-width="{0.35 * units_per_pt:.3f}"/>'
        )
        svg.append('</g>')
        manifest.append(
            f"{index}\t{panel['title']}\t{panel['path']}\t90 degrees counterclockwise\t"
            f"{PIXEL_SIZE_UM}\t{SCALE_BAR_UM}\t{width}\t{height}"
        )
        y = image_y + height
        if index < len(encoded):
            svg.append(f'<rect x="0" y="{y:.3f}" width="{panel_width}" height="{gap_height:.3f}" fill="#000000"/>')
            y += gap_height

    svg.extend(['</g>', '</svg>'])
    svg_path = SOURCE_DIR / "KRT6B_SOX9_PDGFRa_KRT14_vertical_3panel_editable.svg"
    svg_path.write_text("\n".join(svg), encoding="utf-8")
    (SOURCE_DIR / "source_manifest.tsv").write_text("\ufeff" + "\n".join(manifest) + "\n", encoding="utf-8")
    (SOURCE_DIR / "editability_note.txt").write_text(
        "Hybrid editable source: microscopy panels are embedded raster evidence images; titles, channel labels, "
        "scale bars, borders, and layout are vector/text objects. Scale bars use 1.0113780498504639 um/pixel.\n",
        encoding="utf-8",
    )
    print(f"SVG_PASS {svg_path}")
    print(f"ARTBOARD_PT={FIGURE_WIDTH_PT:.6f}x{figure_height_pt:.6f}")
    print(f"SCALE_BAR_PX={scale_bar_px:.6f}")


if __name__ == "__main__":
    main()
