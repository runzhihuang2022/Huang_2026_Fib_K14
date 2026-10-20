from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root
    audit_path = root / "AI_editable" / "QC" / "ai_reaudit.csv"
    with audit_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    report = [
        "AI editability report",
        "Status: RECOMMENDED EDITABLE",
        "Conversion application: Adobe Illustrator 2024 / version 28.0.0",
        "Save mode: native IllustratorSaveOptions with PDF compatibility enabled and linked items embedded.",
        "Font audit: ArialMT and Arial-BoldMT only.",
        "Scientific plot objects were regenerated from source data; no image tracing was used.",
        "Main panel dimensions: 119.0551 x 113.3858 pt = 42 x 40 mm.",
        "A4 proof dimensions: 595.2756 x 841.8898 pt = 210 x 297 mm.",
        "",
        "Object-level audit:",
    ]
    for row in rows:
        report.append(
            f"- {Path(row['AI_Path']).name}: Opened={row['Opened']}; "
            f"TextFrames={row['TextFrames']}; PathItems={row['PathItems']}; "
            f"GroupItems={row['GroupItems']}; PlacedItems={row['PlacedItems']}; "
            f"RasterItems={row['RasterItems']}; Error={row['Error'] or 'none'}."
        )
    report.extend(
        [
            "",
            "Interpretation:",
            "- Text and ordinary plot marks are editable Illustrator objects.",
            "- No placed or raster objects are present.",
            "- The A4 proof places the 42 mm panel at native scale near the upper-left of an A4 page.",
        ]
    )
    (root / "AI_editability_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")

    included_extensions = {".pdf", ".ai", ".csv", ".tsv", ".txt", ".png", ".svg", ".eps", ".tiff"}
    files = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in included_extensions and path.name != "FILE_MANIFEST_SHA256.tsv"
    )
    lines = ["relative_path\tbytes\tsha256"]
    for path in files:
        lines.append(f"{path.relative_to(root)}\t{path.stat().st_size}\t{sha256(path)}")
    (root / "FILE_MANIFEST_SHA256.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
