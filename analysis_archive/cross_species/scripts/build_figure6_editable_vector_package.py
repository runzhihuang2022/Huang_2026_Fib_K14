from __future__ import annotations

import csv
import hashlib
import json
import shutil
from io import BytesIO
from pathlib import Path

import pypdfium2 as pdfium
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


MM = 72.0 / 25.4
DATE = "20260714"
ROOT = Path(r"data/external/local_54ab21a5fa6d") / f"课题01+Fib_K14_KLF4+跨物种主附图AI无双重重组酶+{DATE}"
FIG = ROOT / "figures"
SCR = ROOT / "scripts"
REP = ROOT / "reports"

MAIN_AB = Path(r"data/external/local_81e76774f482")
MAIN_CD = Path(r"data/external/local_a7eefaecf6bc")
MAIN_E = Path(r"data/external/local_455e5f5e19d6")
MAIN_F = Path(r"data/external/local_6af6a7b177da")

SUPP = [
    Path(r"data/external/local_7aa1ed4a02d7"),
    Path(r"data/external/local_2f1f522d5c25"),
    Path(r"data/external/local_83622dec31b3"),
    Path(r"data/external/local_01cd1bf55e0c"),
    Path(r"data/external/local_6af6a7b177da"),
    Path(r"data/external/local_d05e9550534e"),
]


def mm(value: float) -> float:
    return value * MM


def ensure_inputs() -> None:
    for path in [MAIN_AB, MAIN_CD, MAIN_E, MAIN_F, *SUPP]:
        if not path.exists():
            raise FileNotFoundError(path)
    for folder in [FIG, SCR, REP]:
        folder.mkdir(parents=True, exist_ok=True)


def add_source_page(target, source_path: Path, x_mm: float, top_mm: float,
                    box_w_mm: float, box_h_mm: float, page_h_mm: float) -> dict:
    reader = PdfReader(str(source_path))
    source = reader.pages[0]
    if source.rotation:
        source.transfer_rotation_to_content()
    sw, sh = float(source.mediabox.width), float(source.mediabox.height)
    scale = min(mm(box_w_mm) / sw, mm(box_h_mm) / sh)
    placed_w, placed_h = sw * scale, sh * scale
    x = mm(x_mm) + (mm(box_w_mm) - placed_w) / 2
    y = mm(page_h_mm - top_mm) - placed_h
    target.merge_transformed_page(
        source,
        Transformation().scale(scale).translate(x, y),
        over=True,
        expand=False,
    )
    return {
        "source": str(source_path),
        "source_pt": [sw, sh],
        "scale": scale,
        "placed_mm": [x / MM, y / MM, placed_w / MM, placed_h / MM],
    }


def label_overlay(page_w_mm: float, page_h_mm: float, labels: list[tuple[str, float, float]]) -> object:
    packet = BytesIO()
    c = canvas.Canvas(packet, pagesize=(mm(page_w_mm), mm(page_h_mm)))
    arial = Path(r"data/external/local_16a706108b0a")
    arial_bold = Path(r"data/external/local_28ac8473f586")
    if arial.exists() and "Arial" not in pdfmetrics.getRegisteredFontNames():
        pdfmetrics.registerFont(TTFont("Arial", str(arial)))
    if arial_bold.exists() and "Arial-Bold" not in pdfmetrics.getRegisteredFontNames():
        pdfmetrics.registerFont(TTFont("Arial-Bold", str(arial_bold)))
    c.setFont("Arial-Bold" if arial_bold.exists() else "Helvetica-Bold", 12)
    for text, x_mm, top_mm in labels:
        c.drawString(mm(x_mm), mm(page_h_mm - top_mm), text)
    c.save()
    packet.seek(0)
    return PdfReader(packet).pages[0]


def write_page(path: Path, page) -> None:
    writer = PdfWriter()
    writer.add_page(page)
    with path.open("wb") as handle:
        writer.write(handle)


def build_main() -> tuple[Path, dict]:
    page_w, page_h = 210.0, 190.0
    writer = PdfWriter()
    page = writer.add_blank_page(width=mm(page_w), height=mm(page_h))
    placements = []
    placements.append(add_source_page(page, MAIN_AB, 4, 6, 202, 72, page_h))
    placements.append(add_source_page(page, MAIN_CD, 4, 82, 61, 104, page_h))
    placements.append(add_source_page(page, MAIN_E, 67, 82, 139, 59, page_h))
    placements.append(add_source_page(page, MAIN_F, 67, 143, 139, 43, page_h))
    labels = [
        ("A", 2.0, 5.0), ("B", 142.0, 5.0),
        ("C", 2.0, 82.0), ("D", 33.0, 82.0),
        ("E", 66.0, 82.0), ("F", 66.0, 143.0),
    ]
    page.merge_page(label_overlay(page_w, page_h, labels), over=True)
    page.compress_content_streams()
    out = FIG / "Figure6_computational_no_dual_recombinase_editable.pdf"
    with out.open("wb") as handle:
        writer.write(handle)
    return out, {"size_mm": [page_w, page_h], "placements": placements, "labels": labels}


def build_supplement() -> tuple[Path, dict]:
    page_w, page_h = 210.0, 297.0
    writer = PdfWriter()
    page = writer.add_blank_page(width=mm(page_w), height=mm(page_h))
    boxes = [
        (4, 7, 99, 64), (107, 7, 99, 64),
        (4, 75, 99, 64), (107, 75, 99, 64),
        (4, 143, 202, 57), (4, 204, 202, 86),
    ]
    placements = [
        add_source_page(page, path, *box, page_h)
        for path, box in zip(SUPP, boxes)
    ]
    labels = [
        ("A", 2.0, 6.0), ("B", 105.0, 6.0),
        ("C", 2.0, 74.0), ("D", 105.0, 74.0),
        ("E", 2.0, 142.0), ("F", 2.0, 203.0),
    ]
    page.merge_page(label_overlay(page_w, page_h, labels), over=True)
    page.compress_content_streams()
    out = FIG / "FigureS6_computational_support_editable.pdf"
    with out.open("wb") as handle:
        writer.write(handle)
    return out, {"size_mm": [page_w, page_h], "placements": placements, "labels": labels}


def render_pdf(pdf_path: Path, png_path: Path, dpi: int = 220) -> None:
    doc = pdfium.PdfDocument(str(pdf_path))
    page = doc[0]
    bitmap = page.render(scale=dpi / 72.0, rotation=0)
    bitmap.to_pil().save(png_path, dpi=(dpi, dpi))
    page.close()
    doc.close()


def write_legends() -> None:
    text = """# Updated legends for the computational-only vector package

## Figure 6 | Cross-species conservation and spatial localization of Fib_K14-like cells reveal a KLF4-centered epithelialization program

**(A)** Fib_K14-like cell density projected onto the native UMAP of each human, rat, reindeer, Acomys and planarian dataset. Density is scaled within each dataset.

**(B)** Pooled proportions of Fib_K14-like cells or spots. Brackets show prespecified two-sided Fisher's exact tests on positive and negative cell/spot counts. These comparisons are descriptive at the cell/spot level and are not animal-level inference. *P* < 0.05 and ***P*** < 0.001.

**(C)** Representative spatial density maps of Fib_K14-like spots in rat embryonic, P5 and adult wound sections. Percentages indicate positive spots in the displayed section; the grouped bar chart summarizes wound-center, epithelialization-front and peri-wound regions. Scale bars, 1 mm.

**(D)** Spatial density maps of Fib_K14-like spots in planarians at 0 h, 36 h, 3 d and 5 d after amputation. Percentages indicate positive spots in the displayed section. Scale bars, 500 um.

**(E)** Left, full-length KLF4 protein neighbor-joining tree; asterisks denote explicitly labelled proxy species. Middle, cell-count-weighted proportion of KLF4-positive Fib_K14 cells pooled within human, rat and reindeer. Right, KLF4 motif-linked peaks in rat Fib_K14 cells plotted by partial Spearman correlation and FDR. Orange and blue denote significant positive and negative associations; grey denotes non-significant associations. Sixty of 124 peaks reached FDR < 0.05; Epcam was not detected among tested peaks.

**(F)** Sample-level KLF4 regulon activity and motif-pruned KLF4 target program across human, rat and reindeer Fib_K14 cells. Points represent sample-level summaries. Prespecified contrasts were tested using two-sided Mann-Whitney U tests with Benjamini-Hochberg correction. ns, not significant; *FDR* < 0.05.

## Figure S6 | Computational controls and multi-layer validation supporting Figure 6

**(A)** Availability and orthology audit of KLF4 and IGFL1 in the Fib_K14 and all-cell/spot feature spaces.

**(B)** Cross-species effect-size summary for KLF4, IGFL1 and associated repair modules across prespecified injury or regenerative contrasts.

**(C)** Spearman associations between KLF4 or IGFL1 and extracellular-matrix, growth-factor, inflammatory and epithelial-repair modules. Coefficients and significance are displayed in the cells.

**(D)** Group-level KLF4 repair, niche and epithelial-repair programs in Fib_K14 and all-cell/spot contexts. Blank cells indicate unresolved or unavailable features rather than confirmed biological absence.

**(E)** Sample-level KLF4 regulon activity and the motif-pruned KLF4 target program across human, rat and reindeer Fib_K14 cells.

**(F)** Detailed rat Fib_K14 KLF4 chromatin evidence, including motif activity and motif-linked target-peak associations. FDR correction was applied as indicated in the panels.
"""
    (REP / "Figure6_and_FigureS6_computational_legends.md").write_text(text, encoding="utf-8")


def manifest() -> None:
    rows = []
    manifest_path = REP / "file_manifest_sha256.csv"
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == manifest_path:
            continue
        rows.append({
            "relative_path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    with manifest_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ensure_inputs()
    main_pdf, main_report = build_main()
    supp_pdf, supp_report = build_supplement()
    render_pdf(main_pdf, FIG / "Figure6_computational_no_dual_recombinase_preview.png")
    render_pdf(supp_pdf, FIG / "FigureS6_computational_support_preview.png")
    write_legends()
    shutil.copy2(Path(__file__), SCR / Path(__file__).name)
    report = {
        "backend": "Python",
        "archetype": "asymmetric mixed-modality main figure plus quantitative supplementary grid",
        "core_conclusion": "Fib_K14-like states recur across injury and regenerative contexts and are supported by conserved KLF4 evolutionary, expression, chromatin and regulon evidence.",
        "dual_recombinase_results": "excluded",
        "main": main_report,
        "supplement": supp_report,
        "editable_master": "vector PDF; converted separately to native Illustrator AI",
    }
    (REP / "vector_package_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest()
    print(json.dumps({"output_root": str(ROOT), **report}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
