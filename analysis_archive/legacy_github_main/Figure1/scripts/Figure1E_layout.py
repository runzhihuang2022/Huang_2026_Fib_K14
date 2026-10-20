from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image

from figure_svg_utils import ROOT, ensure_dirs, image_data_uri, render_svg


CACHE = ROOT / "source_data/Figure1EF_registered11_scores"
META = ROOT / "source_data/Figure1EF_score_provenance.json"
ASSET = ROOT / "assets/linked_rasters/Figure1E"
OUT = ROOT / "outputs/panels"
STEM = "Figure1E_v01"
WIDTH_MM, HEIGHT_MM = 190, 38.6
ORDER = ["5dpb", "19dpb", "12dpb_SPTDI2", "2mph"]
TITLES = ["5 dpb", "19 dpb", "12 dpb SPTDI", "2 mph"]
PALETTE = np.array([
    [238, 238, 238], [198, 219, 239], [107, 174, 214],
    [33, 113, 181], [241, 105, 19], [253, 217, 118],
], dtype=float)


def colors(values: np.ndarray) -> np.ndarray:
    scaled = np.clip(values, 0, 1) * (len(PALETTE) - 1)
    lo = np.floor(scaled).astype(int)
    hi = np.minimum(lo + 1, len(PALETTE) - 1)
    frac = (scaled - lo)[..., None]
    return ((1 - frac) * PALETTE[lo] + frac * PALETTE[hi]).astype(np.uint8)


def rasterize_section(key: str, vmin: float, vmax: float) -> Path:
    data = np.load(CACHE / f"{key}.npz")
    x, y, score = data["x"], data["y"], data["score"]
    width, height = 1400, 560
    pad = 14
    xspan = max(float(x.max() - x.min()), 1.0)
    yspan = max(float(y.max() - y.min()), 1.0)
    usable_w, usable_h = width - 2 * pad, height - 2 * pad
    scale = min(usable_w / xspan, usable_h / yspan)
    draw_w, draw_h = xspan * scale, yspan * scale
    x0, y0 = (width - draw_w) / 2, (height - draw_h) / 2
    px = np.rint(x0 + (x - x.min()) * scale).astype(int)
    py = np.rint(y0 + (y.max() - y) * scale).astype(int)
    norm = np.clip((score - vmin) / max(vmax - vmin, 1e-12), 0, 1)
    grid = np.full((height, width), -1.0, dtype=np.float32)
    # A 3 x 3 footprint improves print visibility without interpolation or smoothing.
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            xx = np.clip(px + dx, 0, width - 1)
            yy = np.clip(py + dy, 0, height - 1)
            np.maximum.at(grid, (yy, xx), norm)
    rgb = np.full((height, width, 3), 255, dtype=np.uint8)
    mask = grid >= 0
    rgb[mask] = colors(grid[mask])
    out = ASSET / f"Figure1E_{key}.png"
    Image.fromarray(rgb, "RGB").save(out, dpi=(600, 600), optimize=True)
    return out


def build_svg(paths: list[Path]) -> str:
    x_positions = [5, 50.2, 95.4, 140.6]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH_MM}mm" height="{HEIGHT_MM}mm" viewBox="0 0 {WIDTH_MM} {HEIGHT_MM}">',
        '<rect width="190" height="38.6" fill="white"/>',
        '<defs><linearGradient id="score" x1="0" y1="1" x2="0" y2="0">'
        '<stop offset="0" stop-color="#EEEEEE"/><stop offset="0.25" stop-color="#C6DBEF"/>'
        '<stop offset="0.5" stop-color="#6BAED6"/><stop offset="0.67" stop-color="#2171B5"/>'
        '<stop offset="0.84" stop-color="#F16913"/><stop offset="1" stop-color="#FDD976"/></linearGradient></defs>',
        '<style>text{font-family:Arial,sans-serif;fill:#111}.panel{font-size:5px;font-weight:700}.head{font-size:3.2px;font-weight:700}.title{font-size:2.5px;font-weight:700}.small{font-size:2.12px}</style>',
        '<text class="panel" x="0.8" y="5">E</text>',
        '<text class="head" x="7" y="4.8">Epithelialization module score</text>',
    ]
    for x, title, path in zip(x_positions, TITLES, paths):
        parts.append(f'<text class="title" x="{x+21}" y="9" text-anchor="middle">{title}</text>')
        parts.append(f'<image href="{image_data_uri(path)}" x="{x}" y="10.5" width="42" height="24" preserveAspectRatio="xMidYMid meet"/>')
    parts.extend([
        '<rect x="186.2" y="12" width="1.7" height="19" fill="url(#score)" stroke="#777" stroke-width="0.2"/>',
        '<text class="small" x="185.4" y="34" text-anchor="middle">Low</text>',
        '<text class="small" x="185.4" y="10.5" text-anchor="middle">High</text>',
        '<text class="small" x="5" y="37">Shared 1st–99th percentile scale across the registered 11-section subset</text>',
        '</svg>',
    ])
    return "\n".join(parts)


def main() -> None:
    ensure_dirs(); ASSET.mkdir(parents=True, exist_ok=True)
    meta = json.loads(META.read_text(encoding="utf-8"))
    vmin, vmax = meta["epithelialization_display_limits_registered11_percentile_1_99"]
    paths = [rasterize_section(key, vmin, vmax) for key in ORDER]
    out = OUT / f"{STEM}.svg"
    out.write_text(build_svg(paths), encoding="utf-8")
    render_svg(out, WIDTH_MM, HEIGHT_MM, STEM)
    (ROOT / f"QC/{STEM}_content_proof.txt").write_text(
        "Figure 1E v01\n- Recomputed on the exact registered 11-section human Stereo-seq subset.\n"
        "- Shared 1st-99th percentile scale derived only from those 11 sections.\n"
        "- Pixel aggregation uses the maximum observed score per location; no interpolation or smoothing.\n"
        "- Spatial marks are rasterized; all labels and the color scale remain editable vector objects.\n",
        encoding="utf-8",
    )
    print(out)


if __name__ == "__main__":
    main()
