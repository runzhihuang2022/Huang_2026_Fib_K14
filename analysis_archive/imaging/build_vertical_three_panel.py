from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\ImageJ_double_positive_spatial"
)
OUTPUT = ROOT / "07_vertical_three_panel"
OUTPUT.mkdir(parents=True, exist_ok=True)

PANELS = [
    (
        ROOT / "02_recolored_merges" / "KC_KRT6B_red_SOX9_green_DAPI_blue_RGB.tif",
        "KRT6B⁺SOX9⁺",
    ),
    (
        ROOT / "02_recolored_merges" / "Fib_KRT14_red_PDGFRa_green_DAPI_blue_RGB.tif",
        "PDGFRα⁺KRT14⁺",
    ),
    (
        ROOT / "04_spatial_correspondence" / "SerialSection_double_positive_cell_centers_RGB.tif",
        "Spatial correspondence",
    ),
]

HEADER_HEIGHT = 180
GAP_HEIGHT = 80
TITLE_COLOR = (255, 220, 40)
BACKGROUND = (0, 0, 0)


def load_rotate(path: Path) -> np.ndarray:
    image = tifffile.imread(path)
    if image.ndim != 3 or image.shape[-1] != 3:
        raise ValueError(f"Expected RGB TIFF: {path}, got {image.shape}")
    # Exact 90-degree counterclockwise rotation; no interpolation or resampling.
    return np.ascontiguousarray(np.rot90(image, k=1))


def main() -> None:
    rotated = [(load_rotate(path), title, path) for path, title in PANELS]
    shapes = {array.shape for array, _, _ in rotated}
    if len(shapes) != 1:
        raise ValueError(f"Rotated panels do not share one shape: {sorted(shapes)}")
    panel_height, panel_width, _ = rotated[0][0].shape
    canvas_height = len(rotated) * (HEADER_HEIGHT + panel_height) + (len(rotated) - 1) * GAP_HEIGHT
    canvas = np.zeros((canvas_height, panel_width, 3), dtype=np.uint8)

    font_path = Path(r"data/external/local_28ac8473f586")
    font = ImageFont.truetype(str(font_path), 92)
    canvas_image = Image.fromarray(canvas, mode="RGB")
    draw = ImageDraw.Draw(canvas_image)

    y = 0
    manifest_rows = ["order\ttitle\tsource\trotation\tpanel_width\tpanel_height"]
    for index, (array, title, source) in enumerate(rotated, start=1):
        draw.rectangle((0, y, panel_width, y + HEADER_HEIGHT), fill=BACKGROUND)
        draw.text((110, y + 38), title, font=font, fill=TITLE_COLOR)
        image_y = y + HEADER_HEIGHT
        canvas_image.paste(Image.fromarray(array, mode="RGB"), (0, image_y))
        manifest_rows.append(
            f"{index}\t{title}\t{source}\t90 degrees counterclockwise\t{panel_width}\t{panel_height}"
        )
        y = image_y + panel_height
        if index < len(rotated):
            draw.rectangle((0, y, panel_width, y + GAP_HEIGHT), fill=BACKGROUND)
            y += GAP_HEIGHT

    final = np.asarray(canvas_image)
    tiff_path = OUTPUT / "KRT6B_SOX9_PDGFRa_KRT14_spatial_vertical_3panel_CCW90.tif"
    png_path = OUTPUT / "KRT6B_SOX9_PDGFRa_KRT14_spatial_vertical_3panel_CCW90.png"
    preview_path = OUTPUT / "KRT6B_SOX9_PDGFRa_KRT14_spatial_vertical_3panel_CCW90_preview.png"
    tifffile.imwrite(
        tiff_path,
        final,
        bigtiff=True,
        photometric="rgb",
        compression="deflate",
        predictor=True,
        tile=(256, 256),
        metadata={"axes": "YXS"},
    )
    Image.fromarray(final, mode="RGB").save(png_path, compress_level=6)
    preview = Image.fromarray(final, mode="RGB")
    preview.thumbnail((2400, 2400), Image.Resampling.LANCZOS)
    preview.save(preview_path, compress_level=6)
    (OUTPUT / "panel_manifest.tsv").write_text("\ufeff" + "\n".join(manifest_rows) + "\n", encoding="utf-8")
    print(f"PASS TIFF={tiff_path}")
    print(f"PASS PNG={png_path}")
    print(f"PASS PREVIEW={preview_path}")
    print(f"FINAL_SIZE={final.shape[1]}x{final.shape[0]}")


if __name__ == "__main__":
    main()
