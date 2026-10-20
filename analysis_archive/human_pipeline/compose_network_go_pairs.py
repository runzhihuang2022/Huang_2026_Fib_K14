#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import fitz
from PIL import Image


MM_TO_PT = 72 / 25.4


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    figure_root = root / "04_figures/tf_network_enrichment_reference"
    records = []
    for network_pdf in sorted(figure_root.rglob("*_signed_network_reference_v26.pdf")):
        prefix = network_pdf.name.removesuffix("_signed_network_reference_v26.pdf")
        go_pdf = network_pdf.with_name(prefix + "_GO_BP_barplot_reference_v26.pdf")
        network_png = network_pdf.with_suffix(".png")
        go_png = go_pdf.with_suffix(".png")
        if not go_pdf.exists() or not network_png.exists() or not go_png.exists():
            records.append({"instance": prefix, "status": "FAIL", "reason": "missing locked source panel"})
            continue
        out_pdf = network_pdf.with_name(prefix + "_network_GO_paired_150x70mm.pdf")
        doc = fitz.open()
        page = doc.new_page(width=150 * MM_TO_PT, height=70 * MM_TO_PT)
        left = fitz.Rect(5 * MM_TO_PT, 8.75 * MM_TO_PT, 57.5 * MM_TO_PT, 61.25 * MM_TO_PT)
        right = fitz.Rect(92.5 * MM_TO_PT, 8.75 * MM_TO_PT, 145 * MM_TO_PT, 61.25 * MM_TO_PT)
        nsrc = fitz.open(network_pdf); gsrc = fitz.open(go_pdf)
        page.show_pdf_page(left, nsrc, 0); page.show_pdf_page(right, gsrc, 0)
        page.draw_line(fitz.Point(75 * MM_TO_PT, 4 * MM_TO_PT), fitz.Point(75 * MM_TO_PT, 66 * MM_TO_PT), color=(0, 0, 0), width=1.2)
        doc.save(out_pdf); doc.close(); nsrc.close(); gsrc.close()

        nimg = Image.open(network_png).convert("RGB").resize((1240, 1240), Image.Resampling.LANCZOS)
        gimg = Image.open(go_png).convert("RGB").resize((1240, 1240), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (3543, 1654), "white")
        canvas.paste(nimg, (118, 207)); canvas.paste(gimg, (2185, 207))
        for x in range(1767, 1776):
            for y in range(95, 1559):
                canvas.putpixel((x, y), (0, 0, 0))
        out_png = out_pdf.with_suffix(".png")
        canvas.save(out_png, dpi=(600, 600))
        canvas.save(out_pdf.with_suffix(".tiff"), compression="tiff_lzw", dpi=(600, 600))
        records.append({"instance": prefix, "status": "PASS", "pdf": str(out_pdf.relative_to(root)), "png": str(out_png.relative_to(root))})
    report = {"status": "PASS" if records and all(x["status"] == "PASS" for x in records) else "FAIL", "instances": len(records), "records": records}
    (root / "10_qa/NETWORK_GO_PAIRED_COMPOSITION_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "instances": len(records)}, indent=2))


if __name__ == "__main__":
    main()
