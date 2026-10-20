#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pymupdf
import numpy as np
import pandas as pd
from PIL import Image


LOCKED = {
    "render_top5_tf_network_go_reference_v26.py": "7f7410d1815c1de832d482208cc0bbfdc8ef5e4a7169b34b4007b1075b0b3133",
    "render_cellchat_default_chord_v26.R": "65488c8b88e9963ba2a183b2cdb02038c647daefe2610bd504fc9af4cc4d307a",
    "compose_cellchat_default_chord_v26.py": "54943bca8b75760b39f3966fe803949b87689ce5717c02eaa2dfe991d6c06f40",
    "render_targeted_cellchat_v25.R": "cf4fd4436fc3af65d7d4ea0393172d91d2c82918762d0f43dd0f29121adf9f31",
    "render_cellchat_lr_enrichment_v25.py": "439726d3b11cad11ab470aa56ca3f5f4e1ccf9f734098e87c406262c2a52df5d",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def raster_metrics(path: Path) -> dict[str, object]:
    with Image.open(path) as image:
        rgb = np.asarray(image.convert("RGB"))
    ink = np.any(rgb < 245, axis=2)
    h, w = ink.shape
    center = ink[int(.18*h):int(.82*h), int(.18*w):int(.82*w)]
    return {
        "path": str(path), "width": w, "height": h,
        "ink_fraction": float(ink.mean()), "center_ink_fraction": float(center.mean()),
        "near_blank": bool(ink.mean() < 0.008),
    }


def text_overlap_metrics(pdf: Path) -> dict[str, object]:
    doc = pymupdf.open(pdf); page = doc[0]
    words = page.get_text("words")
    boxes = [(float(x0), float(y0), float(x1), float(y1), str(word)) for x0, y0, x1, y1, word, *_ in words if str(word).strip()]
    severe = []
    for i, a in enumerate(boxes):
        ax0, ay0, ax1, ay1, at = a; aa = max((ax1-ax0)*(ay1-ay0), 1e-9)
        for b in boxes[i+1:]:
            bx0, by0, bx1, by1, bt = b
            ix = max(0.0, min(ax1,bx1)-max(ax0,bx0)); iy = max(0.0, min(ay1,by1)-max(ay0,by0))
            if not ix or not iy: continue
            ba = max((bx1-bx0)*(by1-by0), 1e-9); ratio = ix*iy/min(aa,ba)
            if ratio > 0.55: severe.append([at, bt, ratio])
    rect = page.rect; clipped = [w[4] for w in boxes if w[0] < -0.2 or w[1] < -0.2 or w[2] > rect.width+0.2 or w[3] > rect.height+0.2]
    fonts = sorted({font[3] for font in page.get_fonts(full=True)})
    doc.close()
    return {"path": str(pdf), "severe_overlap_pairs": len(severe), "clipped_words": clipped, "fonts": fonts}


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--project", type=Path, required=True); args = parser.parse_args()
    root = args.project.resolve(); qa_dir = root / "10_qa"
    netroot = root / "04_figures/tf_network_enrichment_reference"
    chordroot = root / "04_figures/cellchat_default_chord"
    targeted = root / "04_figures/cellchat_targeted_locked_v27"
    html_path = root / "06_html/OPEN_ME_FIRST_human_epithelialization_complete_v03.html"
    counts = {
        "network_pdf": len(list(netroot.rglob("*_signed_network_reference_v26.pdf"))),
        "go_pdf": len(list(netroot.rglob("*_GO_BP_barplot_reference_v26.pdf"))),
        "paired_pdf": len(list(netroot.rglob("*_network_GO_paired_150x70mm.pdf"))),
        "cellchat_default_chord_pdf": len([x for x in chordroot.glob("*_CellChat_default_LR_chord.pdf") if "raw_default" not in x.name]),
        "cellchat_network_pdf": len(list(targeted.glob("*_network.pdf"))),
        "cellchat_centrality_pdf": len(list(targeted.glob("*_in_out_hub.pdf"))),
        "cellchat_key_lr_chord_pdf": len(list(targeted.glob("*_key_LR_chord.pdf"))),
        "cellchat_violin_pdf": len(list(targeted.glob("*_violin.pdf"))),
    }
    network_count_pass = counts["network_pdf"] == 225 and counts["go_pdf"] == 225 and counts["paired_pdf"] == 225
    lr = pd.read_csv(root / "05_tables/cellchat_default_chord/CELLCHAT_DEFAULT_CHORD_INSTANCE_AUDIT.tsv", sep="\t")
    lr_pass = bool((lr.loc[lr.status.eq("PASS"), "n_display_lr"].astype(int) >= 8).all())
    chord_count_pass = counts["cellchat_default_chord_pdf"] == int(lr.status.eq("PASS").sum())

    hash_rows = []
    for name, expected in LOCKED.items():
        path = root / "00_manifest/locked_renderer" / name; observed = sha256(path)
        hash_rows.append({"file": name, "expected": expected, "observed": observed, "status": "PASS" if expected == observed else "FAIL"})
    hash_pass = all(row["status"] == "PASS" for row in hash_rows)

    visual_paths = list(netroot.rglob("*_signed_network_reference_v26.png")) + list(netroot.rglob("*_GO_BP_barplot_reference_v26.png")) + list(netroot.rglob("*_network_GO_paired_150x70mm.png")) + [x for x in chordroot.glob("*_CellChat_default_LR_chord.png") if "raw_default" not in x.name] + list(targeted.glob("*.png"))
    visuals = [raster_metrics(path) for path in visual_paths]
    blank = [row for row in visuals if row["near_blank"]]
    chord_center_fail = [row for row in visuals if "CellChat_default_LR_chord" in row["path"] and row["center_ink_fraction"] < 0.01]

    pdf_samples = list(netroot.rglob("*_signed_network_reference_v26.pdf")) + list(netroot.rglob("*_GO_BP_barplot_reference_v26.pdf")) + [x for x in chordroot.glob("*_CellChat_default_LR_chord.pdf") if "raw_default" not in x.name]
    text = [text_overlap_metrics(path) for path in pdf_samples]
    clipped = [row for row in text if row["clipped_words"]]
    severe = [row for row in text if row["severe_overlap_pairs"] > 2]

    source = html_path.read_text(encoding="utf-8") if html_path.exists() else ""
    refs = re.findall(r'(?:href|src)="([^"]+)"', source)
    missing = []
    for ref in refs:
        if ref.startswith(("http:", "https:", "#")): continue
        target = (html_path.parent / ref).resolve()
        if not target.exists(): missing.append(ref)
    old_reuse = "CellChat_native_exact_v01" in source or (root / "04_figures/CellChat_native_exact_v01").exists()

    status = "PASS" if all([network_count_pass, lr_pass, chord_count_pass, hash_pass, not blank, not chord_center_fail, not clipped, not severe, not missing, not old_reuse]) else "FAIL"
    report = {
        "status": status, "counts": counts, "network_count_pass": network_count_pass,
        "cellchat_estimable_interfaces": int(lr.status.eq("PASS").sum()), "cellchat_not_estimable_interfaces": int(lr.status.ne("PASS").sum()),
        "cellchat_at_least_8_significant_lr_pass": lr_pass, "cellchat_chord_count_pass": chord_count_pass,
        "locked_renderer_hash_pass": hash_pass, "locked_renderer_hashes": hash_rows,
        "visual_panels_checked": len(visuals), "near_blank_failures": blank, "chord_center_failures": chord_center_fail,
        "pdf_text_panels_checked": len(text), "clipped_text_failures": clipped, "severe_text_overlap_failures": severe,
        "html_references": len(refs), "html_missing_references": missing, "old_v01_cellchat_figure_reuse": old_reuse,
        "causal_claim": False,
    }
    (qa_dir / "FINAL_V03_LOCKED_RENDER_QA.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(visuals).to_csv(qa_dir / "VISUAL_OCCUPANCY_AUDIT.tsv", sep="\t", index=False)
    pd.DataFrame(text).to_json(qa_dir / "PDF_TEXT_OVERLAP_AUDIT.json", orient="records", indent=2)
    print(json.dumps({k: report[k] for k in ["status", "counts", "cellchat_estimable_interfaces", "cellchat_not_estimable_interfaces", "visual_panels_checked", "pdf_text_panels_checked", "html_references"]}, ensure_ascii=False, indent=2))
    if status != "PASS": raise SystemExit(2)


if __name__ == "__main__":
    main()
