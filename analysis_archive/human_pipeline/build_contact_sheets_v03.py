#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def sheet(paths: list[Path], output: Path, columns: int = 10, cell: int = 240) -> None:
    rows = (len(paths) + columns - 1) // columns
    canvas = Image.new("RGB", (columns * cell, rows * cell), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", 14) if Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf").exists() else ImageFont.load_default()
    for i, path in enumerate(paths):
        image = Image.open(path).convert("RGB"); image.thumbnail((cell - 8, cell - 35), Image.Resampling.LANCZOS)
        x = (i % columns) * cell + (cell - image.width) // 2; y = (i // columns) * cell + 2
        canvas.paste(image, (x, y)); label = path.stem[:32]
        draw.text(((i % columns) * cell + 4, (i // columns) * cell + cell - 28), label, fill="black", font=font)
    output.parent.mkdir(parents=True, exist_ok=True); canvas.save(output, dpi=(150, 150))


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--project", type=Path, required=True); args = parser.parse_args()
    root = args.project.resolve(); out = root / "10_qa/contact_sheets"
    groups = {
        "TF_network_all.png": sorted((root / "04_figures/tf_network_enrichment_reference").rglob("*_signed_network_reference_v26.png")),
        "TF_GO_all.png": sorted((root / "04_figures/tf_network_enrichment_reference").rglob("*_GO_BP_barplot_reference_v26.png")),
        "CellChat_default_chord_all.png": sorted([x for x in (root / "04_figures/cellchat_default_chord").glob("*_CellChat_default_LR_chord.png") if "raw_default" not in x.name]),
        "CellChat_targeted_all.png": sorted((root / "04_figures/cellchat_targeted_locked_v27").glob("*.png")),
    }
    for name, paths in groups.items(): sheet(paths, out / name)
    report = {name: len(paths) for name, paths in groups.items()}
    (out / "CONTACT_SHEET_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__": main()
