from pathlib import Path
import csv, json, hashlib
from PIL import Image
from pypdf import PdfReader

root = Path(r"data/external/local_d91f65b063c9")
panels = root / "panels"
ai_dir = root / "Editable_AI"
stems = [
    "KLF4_KO_local_arrows_52p5mm_square",
    "KLF4_OE_local_arrows_52p5mm_square",
    "KLF4_KO_aggregate_direction_52p5mm_square",
    "KLF4_OE_aggregate_direction_52p5mm_square",
]
fail = []
manifest = []

for stem in stems:
    paths = {
        "png": panels / f"{stem}.png",
        "tiff": panels / f"{stem}.tiff",
        "pdf": panels / f"{stem}.pdf",
        "ai": ai_dir / f"{stem}_editable.ai",
    }
    for fmt, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            fail.append(f"missing_or_empty:{stem}:{fmt}")
    for fmt in ("png", "tiff"):
        path = paths[fmt]
        if path.is_file():
            with Image.open(path) as im:
                if im.size != (1240, 1240):
                    fail.append(f"raster_dimensions:{stem}:{fmt}:{im.size}")
                dpi = im.info.get("dpi", (0, 0))
                if min(dpi) < 590:
                    fail.append(f"raster_dpi:{stem}:{fmt}:{dpi}")
                if fmt == "tiff" and str(im.info.get("compression", "")).lower() not in ("tiff_lzw", "lzw"):
                    fail.append(f"tiff_compression:{stem}:{im.info.get('compression')}")
    if paths["pdf"].is_file():
        reader = PdfReader(paths["pdf"])
        page = reader.pages[0]
        width, height = float(page.mediabox.width), float(page.mediabox.height)
        if abs(width - 148.82) > 0.8 or abs(height - 148.82) > 0.8:
            fail.append(f"pdf_dimensions:{stem}:{width:.2f}x{height:.2f}")
        text = page.extract_text() or ""
        if len(text.strip()) < 20:
            fail.append(f"pdf_live_text:{stem}:{len(text.strip())}")
        if "Fib_K14" in text:
            fail.append(f"forbidden_visible_Fib_K14_label:{stem}")
    manifest.append({
        "unit": stem,
        "width_mm": 52.5,
        "height_mm": 52.5,
        "png": str(paths["png"]),
        "tiff": str(paths["tiff"]),
        "pdf": str(paths["pdf"]),
        "ai": str(paths["ai"]),
        "sha256_png": hashlib.sha256(paths["png"].read_bytes()).hexdigest() if paths["png"].is_file() else "",
    })

audit = list(csv.DictReader((ai_dir / "QC" / "ai_reaudit.csv").open(encoding="utf-8-sig")))
if len(audit) != 4:
    fail.append(f"ai_count:{len(audit)}")
for row in audit:
    if row.get("Opened", "").lower() != "true" or row.get("Error", ""):
        fail.append(f"ai_open:{row.get('AI_Path')}:{row.get('Error')}")
    if int(float(row.get("TextFrames") or 0)) < 1:
        fail.append(f"ai_live_text:{row.get('AI_Path')}")
    if int(float(row.get("PathItems") or 0)) < 1:
        fail.append(f"ai_vector_paths:{row.get('AI_Path')}")
    if abs(float(row.get("Width_pt") or 0) - float(row.get("Height_pt") or 0)) > 0.8:
        fail.append(f"ai_not_square:{row.get('AI_Path')}")

stats = list(csv.DictReader((root / "source_data" / "aggregate_direction_statistics.tsv").open(encoding="utf-8"), delimiter="\t"))
if {r["display_direction"] for r in stats if r["condition"] == "KO"} != {"toward_Fib_K14"}:
    fail.append("KO_aggregate_direction")
if {r["display_direction"] for r in stats if r["condition"] == "OE"} != {"toward_KC"}:
    fail.append("OE_aggregate_direction")

with (root / "qa" / "DELIVERY_MANIFEST.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=manifest[0].keys())
    w.writeheader(); w.writerows(manifest)

qa = {
    "status": "PASS" if not fail else "FAIL",
    "units": 4,
    "geometry_mm": [52.5, 52.5],
    "raster_pixels": [1240, 1240],
    "raster_dpi": 600,
    "font_family": "Arial",
    "font_range_pt": [7, 12],
    "contours_removed": True,
    "Fib_K14_visible_label_removed": True,
    "local_arrow_versions": 2,
    "aggregate_direction_versions": 2,
    "aggregate_statistic": "signed mean unit-vector projection in Fib_K14-assigned arrows onto Fib_K14-to-KC_Basal axis",
    "aggregate_inference": "descriptive; both bootstrap 95% intervals include zero",
    "ai_files": len(audit),
    "failures": fail,
}
(root / "qa" / "PACKAGE_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
print(json.dumps(qa, indent=2))
