#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pymupdf
from PIL import Image


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    folders = [
        root / "04_figures/tf_network_enrichment_reference",
        root / "04_figures/cellchat_targeted_locked_v27",
        root / "04_figures/cellchat_default_chord",
        root / "04_figures/cellchat_lr_enrichment_locked_v27",
    ]
    records = []
    for folder in folders:
        for pdf in sorted(folder.rglob("*.pdf")):
            doc = pymupdf.open(pdf)
            if doc.page_count != 1:
                records.append({"pdf": str(pdf.relative_to(root)), "status": "MULTIPAGE_SKIPPED", "pages": doc.page_count})
                doc.close()
                continue
            png = pdf.with_suffix(".png")
            if not png.exists():
                pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(600 / 72, 600 / 72), alpha=False)
                pix.save(png)
            doc.close()
            with Image.open(png) as image:
                rgb = image.convert("RGB")
                rgb.save(pdf.with_suffix(".tiff"), compression="tiff_lzw", dpi=(600, 600))
            records.append({"pdf": str(pdf.relative_to(root)), "status": "PASS", "png": str(png.relative_to(root)), "tiff": str(pdf.with_suffix('.tiff').relative_to(root))})
    report = {"status": "PASS" if records and all(x["status"] in {"PASS", "MULTIPAGE_SKIPPED"} for x in records) else "FAIL", "records": len(records), "single_page": sum(x["status"] == "PASS" for x in records), "multipage_skipped": sum(x["status"] == "MULTIPAGE_SKIPPED" for x in records)}
    (root / "10_qa/LOCKED_RENDER_PDF_COMPANION_EXPORT_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
