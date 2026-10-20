"""Assemble the already-rendered Figure S4/S5 panel PDFs on one A4 page."""

from __future__ import annotations

import argparse
from pathlib import Path

import fitz


A4_WIDTH_PT = 210 / 25.4 * 72
A4_HEIGHT_PT = 297 / 25.4 * 72


def place_pdf(page: fitz.Page, source: Path, rectangle: fitz.Rect) -> None:
    document = fitz.open(source)
    page.show_pdf_page(rectangle, document, 0, keep_proportion=True, overlay=True)
    document.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", choices=["human", "rat"], required=True)
    parser.add_argument("--package", required=True)
    args = parser.parse_args()
    package = Path(args.package)
    prefix = "FigureS4" if args.species == "human" else "FigureS5"
    panels = package / "figures" / "panels"
    names = {
        "A": f"{prefix}_A_fibroblast_umap_v02.pdf",
        "B": f"{prefix}_B_qc_doublets_v02.pdf",
        "C": f"{prefix}_C_markers_merged_v02.pdf",
        "D": f"{prefix}_D_compact_v02.pdf",
        "E": f"{prefix}_E_monocle_triplet_v02.pdf",
        "F": f"{prefix}_F_spatial_time_v03.pdf",
    }
    if args.species == "human":
        names["G"] = "FigureS4_G_VIM_KRT14_staining_v01.pdf"
    document = fitz.open()
    page = document.new_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
    arial_font = Path(r"data/external/local_16a706108b0a")
    panel_font = "Arial"
    page.insert_font(fontname=panel_font, fontfile=str(arial_font))
    margin, gap = 18, 10
    usable_width = A4_WIDTH_PT - 2 * margin
    cursor = 20
    row_heights = {"AB": 142, "C": 88, "DE": 138}

    half_width = (usable_width - gap) / 2
    place_pdf(
        page,
        panels / names["A"],
        fitz.Rect(margin, cursor, margin + half_width, cursor + row_heights["AB"]),
    )
    place_pdf(
        page,
        panels / names["B"],
        fitz.Rect(
            margin + half_width + gap,
            cursor,
            A4_WIDTH_PT - margin,
            cursor + row_heights["AB"],
        ),
    )
    page.insert_text((7, cursor + 12), "A", fontsize=14, fontname=panel_font)
    page.insert_text(
        (margin + half_width + 1, cursor + 12), "B", fontsize=14, fontname=panel_font
    )
    cursor += row_heights["AB"] + gap

    place_pdf(
        page,
        panels / names["C"],
        fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, cursor + row_heights["C"]),
    )
    page.insert_text((7, cursor + 12), "C", fontsize=14, fontname=panel_font)
    cursor += row_heights["C"] + gap

    third_width = (usable_width - gap) / 3
    place_pdf(
        page,
        panels / names["D"],
        fitz.Rect(margin, cursor, margin + third_width, cursor + row_heights["DE"]),
    )
    place_pdf(
        page,
        panels / names["E"],
        fitz.Rect(
            margin + third_width + gap,
            cursor,
            A4_WIDTH_PT - margin,
            cursor + row_heights["DE"],
        ),
    )
    page.insert_text((7, cursor + 12), "D", fontsize=14, fontname=panel_font)
    page.insert_text(
        (margin + third_width + 1, cursor + 12), "E", fontsize=14, fontname=panel_font
    )
    cursor += row_heights["DE"] + gap

    if args.species == "human":
        f_height = 225
        place_pdf(
            page,
            panels / names["F"],
            fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, cursor + f_height),
        )
        page.insert_text((7, cursor + 12), "F", fontsize=14, fontname=panel_font)
        cursor += f_height + gap
        place_pdf(
            page,
            panels / names["G"],
            fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, A4_HEIGHT_PT - margin),
        )
        page.insert_text((7, cursor + 12), "G", fontsize=14, fontname=panel_font)
    else:
        place_pdf(
            page,
            panels / names["F"],
            fitz.Rect(margin, cursor, A4_WIDTH_PT - margin, A4_HEIGHT_PT - margin),
        )
        page.insert_text((7, cursor + 12), "F", fontsize=14, fontname=panel_font)

    output = package / "figures" / f"{prefix}_revised_layout_v03.pdf"
    document.save(output, garbage=4, deflate=True)
    document.close()
    rendered = fitz.open(output)
    pixmap = rendered[0].get_pixmap(
        matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False
    )
    pixmap.save(
        package / "figures" / f"{prefix}_revised_layout_v03_600dpi.png"
    )
    rendered.close()
    print(output)


if __name__ == "__main__":
    main()
