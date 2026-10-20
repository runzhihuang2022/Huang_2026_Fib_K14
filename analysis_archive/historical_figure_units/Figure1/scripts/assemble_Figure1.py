from __future__ import annotations

import re
import shutil
from pathlib import Path

from figure_svg_utils import ROOT, ensure_dirs, render_svg


PANELS = {
    "A": ROOT / "outputs/panels/Figure1A_v02.svg",
    "B": ROOT / "outputs/panels/Figure1B_v02.svg",
    "C": ROOT / "outputs/panels/Figure1C_v02.svg",
    "D": ROOT / "outputs/panels/Figure1D_v02.svg",
    "E": ROOT / "outputs/panels/Figure1E_v01.svg",
    "F": ROOT / "outputs/panels/Figure1F_v01.svg",
}
STEM = "Figure1_A4_v02"


def inner_svg(path: Path) -> tuple[str, str]:
    text = path.read_text(encoding="utf-8")
    root = re.search(r"<svg\b([^>]*)>", text, flags=re.S)
    match = re.search(r"<svg\b[^>]*>(.*)</svg>\s*$", text, flags=re.S)
    if not match or not root:
        raise ValueError(f"Cannot parse SVG: {path}")
    viewbox_match = re.search(r'viewBox="([^"]+)"', root.group(1))
    if not viewbox_match:
        raise ValueError(f"SVG lacks viewBox: {path}")
    return viewbox_match.group(1), match.group(1)


def main() -> None:
    ensure_dirs()
    missing = [str(path) for path in PANELS.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing panels: " + "; ".join(missing))
    placements = [
        ("A", 10.0, 10.0, 125.0, 92.0),
        ("B", 139.0, 10.0, 61.0, 92.0),
        ("C", 10.0, 105.7, 190.0, 39.0),
        ("D", 10.0, 148.4, 190.0, 53.0),
        ("E", 10.0, 205.1, 190.0, 38.6),
        ("F", 10.0, 247.4, 190.0, 38.6),
    ]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="210mm" height="297mm" viewBox="0 0 210 297">',
        '<rect width="210" height="297" fill="white"/>',
    ]
    for name, x, y, width, height in placements:
        viewbox, content = inner_svg(PANELS[name])
        parts.append(f'<svg id="Figure1{name}" x="{x}" y="{y}" width="{width}" height="{height}" viewBox="{viewbox}" preserveAspectRatio="xMidYMid meet">')
        parts.append(content)
        parts.append('</svg>')
    parts.append('</svg>')
    out = ROOT / f"outputs/panels/{STEM}.svg"
    out.write_text("\n".join(parts), encoding="utf-8")
    render_svg(out, 210, 297, STEM)
    assembled = ROOT / "outputs/assembled"
    assembled.mkdir(parents=True, exist_ok=True)
    for source in (
        ROOT / f"outputs/panels/{STEM}.svg",
        ROOT / f"outputs/panels/{STEM}.pdf",
        ROOT / f"outputs/raster/{STEM}_600dpi.png",
        ROOT / f"outputs/raster/{STEM}_600dpi.tiff",
    ):
        shutil.copy2(source, assembled / source.name)
    (ROOT / f"QC/{STEM}_layout_proof.txt").write_text(
        "Figure 1 A4 v01\n- Physical page size: 210 x 297 mm.\n"
        "- Content margins: 10 mm left/right; at least 10 mm top/bottom.\n"
        "- A and B share the first row; C-F use full-width rows.\n"
        "- Source SVG panels remain individually editable and are grouped by panel ID in the assembled SVG.\n",
        encoding="utf-8",
    )
    print(out)


if __name__ == "__main__":
    main()
