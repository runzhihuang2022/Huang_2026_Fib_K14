#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path

import fitz
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MM_PT = 72 / 25.4
MM_IN = 1 / 25.4

MODES = [
    ("Mode1", "D0_7", "Mode 1", "0–7 dpb", "Unhealed gradient"),
    ("Mode2", "D15_28", "Mode 2", "15–28 dpb", "One-way epithelialization"),
    ("Mode3", "D8_14", "Mode 3", "8–14 dpb", "Almost healed"),
    ("Mode4", "M1_2", "Mode 4", "1–2 mph", "Scar / hyperplasia"),
]
INTERFACES = ["SAC", "KC", "Endo", "Immune"]
DISPLAY = {"SAC": "Fib–SAC", "KC": "Fib–KC", "Endo": "Fib–Endo", "Immune": "Fib–Immune"}
REQUIRED = {
    "SAC": ["Fib_Papi", "Fib_K14", "SAC_SG_Progenitor", "SAC_SG_Ductal", "SAC_SG_Dark"],
    "KC": ["Fib_Papi", "Fib_K14", "KC_Basal_Mig", "KC_Basal", "KC_Spinous_Mig"],
    "Endo": ["Fib_Prolif", "Fib_Myo", "Endo_Capillary", "Endo_Venous", "Endo_Arterial"],
    "Immune": ["Fib_Inflama", "Fib_CD45", "M1_Macrophage", "M2_Macrophage", "Mast"],
}


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 7,
        "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def add_fonts(page: fitz.Page) -> None:
    page.insert_font(fontname="Arial", fontfile=r"data/external/local_16a706108b0a")
    page.insert_font(fontname="Arial-Bold", fontfile=r"data/external/local_28ac8473f586")


def add_text(page: fitz.Page, x_mm: float, y_mm: float, text: str, size: float = 7,
             bold: bool = False, color=(0, 0, 0), align=0, width_mm: float | None = None) -> None:
    font = "Arial-Bold" if bold else "Arial"
    if width_mm is None:
        page.insert_text((x_mm * MM_PT, y_mm * MM_PT), text, fontsize=size, fontname=font, color=color)
    else:
        rect = fitz.Rect(x_mm * MM_PT, (y_mm - size * 0.38) * MM_PT,
                         (x_mm + width_mm) * MM_PT, (y_mm + 5) * MM_PT)
        page.insert_textbox(rect, text, fontsize=size, fontname=font, color=color, align=align)


def place_pdf(page: fitz.Page, source: Path, x_mm: float, y_mm: float, w_mm: float, h_mm: float) -> None:
    src = fitz.open(source)
    page.show_pdf_page(
        fitz.Rect(x_mm * MM_PT, y_mm * MM_PT, (x_mm + w_mm) * MM_PT, (y_mm + h_mm) * MM_PT),
        src, 0, keep_proportion=True,
    )
    src.close()


def rasterize(pdf: Path, png: Path, dpi: int = 300) -> None:
    doc = fitz.open(pdf)
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(dpi / 72, dpi / 72), alpha=False)
    png.parent.mkdir(parents=True, exist_ok=True)
    pix.save(png)
    doc.close()


def make_chord_panels(root: Path, palette: dict[str, str], audit: pd.DataFrame) -> None:
    out_dir = root / "04_figures" / "A_chords" / "final_vector"
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = root / "04_figures" / "A_chords" / "raw_extended5"
    for mode, group, mode_label, time_label, pattern in MODES:
        for interface in INTERFACES:
            row = audit[(audit.group == group) & (audit.interface == f"Fib_{interface}")].iloc[0]
            out = out_dir / f"{mode}_{interface}_required_LR_chord_52p5x74mm.pdf"
            doc = fitz.open()
            page = doc.new_page(width=52.5 * MM_PT, height=74 * MM_PT)
            add_fonts(page)
            add_text(page, 1, 4.2, f"{mode_label} | {time_label} | {DISPLAY[interface]}", 7.3, True, width_mm=50.5, align=1)
            if row.status == "PASS":
                raw = raw_dir / f"{group}_Fib_{interface}_CellChat_default_LR_chord_raw_default.pdf"
                place_pdf(page, raw, 0, 6.0, 52.5, 52.5)
                y0 = 62.0
                nodes = REQUIRED[interface]
                cols = 2
                cell_w = 25.5
                for i, node in enumerate(nodes):
                    rr, cc = divmod(i, cols)
                    x = 1.2 + cc * cell_w
                    y = y0 + rr * 3.2
                    color = palette.get(node, "#BDBDBD")[:7]
                    rgb = tuple(int(color[j:j + 2], 16) / 255 for j in (1, 3, 5))
                    page.draw_circle((x * MM_PT, (y - 0.55) * MM_PT), 0.85 * MM_PT, color=rgb, fill=rgb)
                    add_text(page, x + 1.8, y, node, 7.0)
            else:
                page.draw_rect(fitz.Rect(2 * MM_PT, 10 * MM_PT, 50.5 * MM_PT, 57 * MM_PT),
                               color=(0.72, 0.74, 0.76), fill=(0.97, 0.97, 0.97), width=0.6)
                add_text(page, 3, 28, "NOT_ESTIMABLE", 10, True, color=(0.65, 0.12, 0.12), width_mm=46.5, align=1)
                counts = str(row.reason).split(":", 1)[-1].strip()
                add_text(page, 4, 35, "Required subtype count <10", 7.0, color=(0.25, 0.27, 0.29), width_mm=44.5, align=1)
                add_text(page, 4, 39, counts, 7.0, color=(0.25, 0.27, 0.29), width_mm=44.5, align=1)
                add_text(page, 4, 48, "Panel retained; no nonsignificant LR pairs were added.", 6.2,
                         color=(0.25, 0.27, 0.29), width_mm=44.5, align=1)
            doc.save(out, garbage=4, deflate=True)
            doc.close()
            rasterize(out, out.with_suffix(".png"), 300)


def make_centrality_panels(root: Path, palette: dict[str, str], audit: pd.DataFrame) -> None:
    table_dir = root / "02_analysis" / "required_chord" / "tables_extended5"
    out_dir = root / "04_figures" / "C_centrality" / "required_exact"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for mode, group, mode_label, time_label, pattern in MODES:
        for interface in INTERFACES:
            row = audit[(audit.group == group) & (audit.interface == f"Fib_{interface}")].iloc[0]
            out = out_dir / f"{mode}_{interface}_required_in_out_centrality_52p5mm"
            fig, ax = plt.subplots(figsize=(52.5 * MM_IN, 52.5 * MM_IN))
            if row.status != "PASS":
                ax.text(0.5, 0.56, "NOT_ESTIMABLE", transform=ax.transAxes, ha="center", va="center",
                        fontsize=9, fontweight="bold", color="#A61E1E")
                ax.text(0.5, 0.42, str(row.reason).replace("; ", "\n"), transform=ax.transAxes,
                        ha="center", va="center", fontsize=6.2, color="#4A4A4A")
                ax.axis("off")
            else:
                lr = pd.read_csv(table_dir / f"{group}_Fib_{interface}_all_significant_LR.tsv", sep="\t")
                nodes = REQUIRED[interface]
                outgoing = lr.groupby("source").prob.sum().reindex(nodes, fill_value=0.0)
                incoming = lr.groupby("target").prob.sum().reindex(nodes, fill_value=0.0)
                scale = max(float(outgoing.max()), float(incoming.max()), 1e-12)
                outgoing_rel = outgoing / scale
                incoming_rel = incoming / scale
                for i, node in enumerate(nodes):
                    x, y = float(outgoing_rel[node]), float(incoming_rel[node])
                    summary_rows.append({
                        "mode": mode_label, "group": group, "time_group": time_label,
                        "interface": f"Fib_{interface}", "cell_type": node,
                        "summed_outgoing_probability": float(outgoing[node]),
                        "summed_incoming_probability": float(incoming[node]),
                        "within_panel_scale": scale,
                        "relative_outgoing_strength": x,
                        "relative_incoming_strength": y,
                    })
                    ax.scatter(x, y, s=34, facecolor=palette.get(node, "#BDBDBD")[:7],
                               edgecolor="#222222", linewidth=0.45, zorder=3)
                    dx = 4 if x < 0.72 else -4
                    ha = "left" if dx > 0 else "right"
                    if y < 0.08:
                        dy = 8
                    elif interface == "Endo":
                        dy = (-10, 10, 0, -7, 7)[i]
                    elif interface == "Immune":
                        dy = (10, -10, 10, -10, 0)[i]
                    else:
                        dy = (7, -7, 7, -7, 0)[i]
                    ax.annotate(node, (x, y), xytext=(dx, dy), textcoords="offset points",
                                ha=ha, va="center", fontsize=5.8,
                                arrowprops=dict(arrowstyle="-", lw=0.35, color="#777777"))
                ax.set_xlim(-0.08, 1.38)
                ax.set_ylim(-0.08, 1.22)
                ax.set_xticks([0, 0.5, 1.0])
                ax.set_yticks([0, 0.5, 1.0])
                ax.grid(True, color="#E5E7EA", linewidth=0.45)
                ax.spines[["top", "right"]].set_visible(False)
                ax.set_xlabel("Relative outgoing strength")
                ax.set_ylabel("Relative incoming strength")
            ax.set_title(f"{mode_label} | {time_label} | {DISPLAY[interface]}", fontsize=7.5, fontweight="bold", pad=4)
            fig.subplots_adjust(left=0.25, right=0.95, bottom=0.22, top=0.86)
            fig.savefig(out.with_suffix(".pdf"))
            fig.savefig(out.with_suffix(".png"), dpi=600)
            plt.close(fig)
    pd.DataFrame(summary_rows).to_csv(
        root / "02_analysis" / "required_chord" / "tables" / "REQUIRED_CELL_CENTRALITY_SUMMARY.tsv",
        sep="\t", index=False,
    )


def replace_go_placeholders(root: Path, audit: pd.DataFrame) -> None:
    # The enrichment renderer writes both statistically empty and low-count
    # placeholders. Preserve those files so their evidence-boundary messages
    # are not replaced during sheet composition.
    return


def source_for(root: Path, section: str, mode: str, group: str, interface: str) -> Path:
    if section == "A":
        return root / "04_figures" / "A_chords" / "final_vector" / f"{mode}_{interface}_required_LR_chord_52p5x74mm.pdf"
    if section == "B":
        return root / "04_figures" / "B_networks" / f"{mode}_{interface}_required_cell_network_42mm.pdf"
    if section == "C":
        return root / "04_figures" / "C_centrality" / "required_exact" / f"{mode}_{interface}_required_in_out_centrality_52p5mm.pdf"
    if section == "D":
        return root / "05_enrichment" / "D_GO_extended5" / f"{mode}_{interface}_required_cell_GO_42x84mm.pdf"
    raise ValueError(section)


def build_sheet(root: Path, section: str, title: str) -> Path:
    width = 240.0
    if section == "A":
        panel_w, panel_h, cell_w, cell_h = 52.5, 74.0, 57.0, 74.0
    elif section in {"B", "C"}:
        panel_w, panel_h, cell_w, cell_h = 52.5, 52.5, 57.0, 55.0
    else:
        panel_w, panel_h, cell_w, cell_h = 42.0, 84.0, 57.0, 84.0
    top, bottom, gap_y = 14.0, 6.0, 3.0
    height = top + 4 * cell_h + 3 * gap_y + bottom
    out = root / "07_figure_s3" / "sheets" / f"FigureS3_{section}_Mode1-4_vector.pdf"
    doc = fitz.open()
    page = doc.new_page(width=width * MM_PT, height=height * MM_PT)
    add_fonts(page)
    add_text(page, 4, 8, section, 12, True)
    add_text(page, 12, 8, title, 10, True)
    x0 = (width - 4 * cell_w) / 2
    for r, (mode, group, mode_label, time_label, pattern) in enumerate(MODES):
        for c, interface in enumerate(INTERFACES):
            src = source_for(root, section, mode, group, interface)
            x = x0 + c * cell_w + (cell_w - panel_w) / 2
            y = top + r * (cell_h + gap_y) + (cell_h - panel_h) / 2
            place_pdf(page, src, x, y, panel_w, panel_h)
    note = {
        "A": "CellChat::netVisual_chord_gene; 8 unique significant LR pairs for each estimable interface.",
        "B": "Line width represents summed CellChat communication probability among the prespecified required subtypes.",
        "C": "Incoming/outgoing values are sums across all significant required-cell LR pairs within each pooled mode object.",
        "D": "GO Biological Process: one-sided hypergeometric test with within-interface BH correction; stars denote q values.",
    }[section]
    add_text(page, 10, height - 2.2, note, 7.0, color=(0.30, 0.32, 0.34), width_mm=220, align=1)
    doc.save(out, garbage=4, deflate=True)
    doc.close()
    rasterize(out, out.with_suffix(".png"), 200)
    return out


def build_combined(root: Path, sheets: list[Path]) -> Path:
    widths, heights = [], []
    for sheet in sheets:
        d = fitz.open(sheet)
        widths.append(d[0].rect.width / MM_PT)
        heights.append(d[0].rect.height / MM_PT)
        d.close()
    gap = 4.0
    out = root / "07_figure_s3" / "combined" / "FigureS3_CellChat_Mode1-4_ABCD_complete_vector.pdf"
    doc = fitz.open()
    page = doc.new_page(width=max(widths) * MM_PT, height=(sum(heights) + gap * (len(sheets) - 1)) * MM_PT)
    y = 0.0
    for sheet, h in zip(sheets, heights):
        place_pdf(page, sheet, 0, y, max(widths), h)
        y += h + gap
    doc.save(out, garbage=4, deflate=True)
    doc.close()
    rasterize(out, out.with_suffix(".png"), 120)
    return out


def build_html(root: Path, audit: pd.DataFrame, sheets: list[Path], combined: Path) -> Path:
    out = root / "06_html" / "OPEN_ME_FIRST_FigureS3_CellChat_Mode1-4_chords_v03.html"
    cards = []
    for section, sheet in zip("ABCD", sheets):
        cards.append(
            f'<figure><a href="../{sheet.relative_to(root).as_posix()}"><img src="../{sheet.with_suffix(".png").relative_to(root).as_posix()}"></a>'
            f'<figcaption>Panel {section}: 16 positions retained.</figcaption></figure>'
        )
    rows = audit.to_html(index=False, classes="audit", border=0)
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>Figure S3 CellChat Mode 1–4 v03</title><style>
body{{font-family:Arial,sans-serif;margin:0;background:#f4f6f8;color:#20242a}}header,main{{max-width:1400px;margin:auto;padding:18px}}header{{background:#183153;color:white;max-width:none}}.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}}figure{{background:white;border:1px solid #d9dfe6;padding:10px;margin:0}}img{{width:100%;height:auto}}.notice{{background:#fff8dc;border-left:4px solid #c68a00;padding:10px;margin:12px 0}}table{{border-collapse:collapse;background:white;font-size:12px}}td,th{{border:1px solid #d9dfe6;padding:4px}}</style></head><body>
<header><h1>Figure S3 · Complete Mode 1–4 fibroblast-centered CellChat</h1></header><main>
<div class="notice"><b>Completeness:</b> A–D each retain all 4 modes × 4 interfaces = 16 positions. Mode 1 and Mode 4 Fib–Immune are explicitly NOT_ESTIMABLE because the required M1/M2 identities fail the locked minimum of 10 cells; no nonsignificant LR pairs were inserted.</div>
<p><a href="../{combined.relative_to(root).as_posix()}">Download combined vector PDF</a></p><div class="grid">{''.join(cards)}</div><h2>Required-cell and LR audit</h2>{rows}</main></body></html>'''
    out.write_text(page, encoding="utf-8")
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True, type=Path)
    a = p.parse_args()
    root = a.root
    configure()
    palette_t = pd.read_csv(root / "01_input_adapter" / "subtype_palette.tsv", sep="\t")
    palette = dict(zip(palette_t.subtype.astype(str), palette_t.color.astype(str)))
    audit = pd.read_csv(root / "02_analysis" / "required_chord" / "tables_extended5" / "REQUIRED_CHORD_INSTANCE_AUDIT.tsv", sep="\t")
    make_chord_panels(root, palette, audit)
    make_centrality_panels(root, palette, audit)
    replace_go_placeholders(root, audit)
    titles = {
        "A": "Required-cell ligand–receptor chord diagrams",
        "B": "Required-cell aggregate communication networks",
        "C": "Required-cell incoming and outgoing communication strength",
        "D": "Ligand–receptor gene functional enrichment",
    }
    sheets = [build_sheet(root, s, titles[s]) for s in "ABCD"]
    combined = build_combined(root, sheets)
    build_html(root, audit, sheets, combined)
    manifest = []
    for f in sorted(root.rglob("*")):
        if f.is_file() and f.name != "FILE_MANIFEST_SHA256.tsv":
            manifest.append({"relative_path": f.relative_to(root).as_posix(), "bytes": f.stat().st_size, "sha256": sha256(f)})
    pd.DataFrame(manifest).to_csv(root / "00_manifest" / "FILE_MANIFEST_SHA256.tsv", sep="\t", index=False)
    print(combined)


if __name__ == "__main__":
    main()
