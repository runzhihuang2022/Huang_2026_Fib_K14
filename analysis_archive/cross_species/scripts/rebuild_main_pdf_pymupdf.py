from pathlib import Path

import fitz


MM = 72.0 / 25.4
ROOT = Path(r"data/external/local_90817591e523")
OUT = ROOT / "figures" / "Figure6_computational_no_dual_recombinase_compact.pdf"

SOURCES = [
    (
        Path(r"data/external/local_81e76774f482"),
        (5.526315789451257, 112.0, 198.9473684210975, 72.0),
    ),
    (
        Path(r"data/external/local_a7eefaecf6bc"),
        (4.0, 21.72857142858722, 61.0, 86.27142857141277),
    ),
    (
        Path(r"data/external/local_455e5f5e19d6"),
        (67.66666666667244, 49.0, 137.66666666665512, 59.0),
    ),
    (
        Path(r"data/external/local_6af6a7b177da"),
        (74.82764327517721, 4.0, 123.34471344964558, 43.0),
    ),
]

LABELS = [
    ("A", 2.0, 5.0),
    ("B", 142.0, 5.0),
    ("C", 2.0, 82.0),
    ("D", 33.0, 82.0),
    ("E", 66.0, 82.0),
    ("F", 66.0, 143.0),
]


def mm(value: float) -> float:
    return value * MM


def main() -> None:
    page_w_mm, page_h_mm = 210.0, 190.0
    output = fitz.open()
    page = output.new_page(width=mm(page_w_mm), height=mm(page_h_mm))

    opened = []
    try:
        for source_path, (x, y_bottom, width, height) in SOURCES:
            if not source_path.exists():
                raise FileNotFoundError(source_path)
            source = fitz.open(source_path)
            opened.append(source)
            top = page_h_mm - (y_bottom + height)
            rect = fitz.Rect(mm(x), mm(top), mm(x + width), mm(top + height))
            page.show_pdf_page(rect, source, 0, keep_proportion=False, overlay=True)

        font = "data/external/local_83f91b55b69d"
        for label, x, top_y in LABELS:
            page.insert_text(
                (mm(x), mm(top_y + 4.2)),
                label,
                fontsize=12,
                fontfile=font,
                color=(0, 0, 0),
                overlay=True,
            )

        OUT.parent.mkdir(parents=True, exist_ok=True)
        output.save(OUT, garbage=4, deflate=True, clean=True)
    finally:
        for source in opened:
            source.close()
        output.close()

    print(OUT)
    print(OUT.stat().st_size)


if __name__ == "__main__":
    main()
