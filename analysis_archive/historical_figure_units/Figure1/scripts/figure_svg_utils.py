from __future__ import annotations

import base64
import shutil
import subprocess
import tempfile
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]


def ensure_dirs() -> None:
    for rel in (
        "assets/linked_rasters",
        "outputs/panels",
        "outputs/vector",
        "outputs/raster",
        "source_data",
        "QC",
    ):
        (ROOT / rel).mkdir(parents=True, exist_ok=True)


def image_data_uri(path: Path) -> str:
    mime = "image/png" if path.suffix.lower() == ".png" else "image/jpeg"
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def crop_png(source: Path, target: Path, box: tuple[int, int, int, int]) -> Path:
    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source).convert("RGB") as image:
        image.crop(box).save(target, dpi=(600, 600), optimize=True)
    return target


def render_svg(svg_path: Path, width_mm: float, height_mm: float, stem: str) -> None:
    panel_dir = ROOT / "outputs/panels"
    vector_dir = ROOT / "outputs/vector"
    raster_dir = ROOT / "outputs/raster"
    pdf_path = panel_dir / f"{stem}.pdf"
    png_path = raster_dir / f"{stem}_600dpi.png"
    tiff_path = raster_dir / f"{stem}_600dpi.tiff"
    browser = next(
        (
            p
            for p in (
                Path(r"data/external/local_b98a42e7e75a"),
                Path(r"data/external/local_0647cf87e7e9"),
            )
            if p.exists()
        ),
        None,
    )
    if browser is None:
        raise FileNotFoundError("Chrome or Edge is required for proof rendering")
    svg_text = svg_path.read_text(encoding="utf-8")
    with tempfile.TemporaryDirectory(prefix=f"{stem}-") as temp_name:
        temp = Path(temp_name)
        html = temp / "print.html"
        raw_png = temp / "raw.png"
        profile = temp / "profile"
        html.write_text(
            "<!doctype html><html><head><meta charset='utf-8'><style>"
            f"@page{{size:{width_mm}mm {height_mm}mm;margin:0}}"
            f"html,body{{margin:0;width:{width_mm}mm;height:{height_mm}mm;overflow:hidden;background:#fff}}"
            "svg{display:block}</style></head><body>" + svg_text + "</body></html>",
            encoding="utf-8",
        )
        common = [
            str(browser),
            "--headless=new",
            "--disable-gpu",
            "--hide-scrollbars",
            "--no-first-run",
            "--disable-extensions",
            f"--user-data-dir={profile}",
        ]
        subprocess.run(
            common + [f"--print-to-pdf={pdf_path}", "--no-pdf-header-footer", html.as_uri()],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        css_w = round(width_mm / 25.4 * 96)
        css_h = round(height_mm / 25.4 * 96)
        subprocess.run(
            common
            + [
                "--force-device-scale-factor=6.25",
                f"--window-size={css_w},{css_h}",
                f"--screenshot={raw_png}",
                html.as_uri(),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        px_w = round(width_mm / 25.4 * 600)
        px_h = round(height_mm / 25.4 * 600)
        with Image.open(raw_png).convert("RGB") as image:
            proof = image.resize((px_w, px_h), Image.Resampling.LANCZOS)
            proof.save(png_path, dpi=(600, 600), optimize=True)
            proof.save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
    shutil.copy2(svg_path, vector_dir / svg_path.name)
    shutil.copy2(pdf_path, vector_dir / pdf_path.name)

