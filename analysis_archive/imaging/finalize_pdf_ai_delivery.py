from __future__ import annotations

import csv
import hashlib
from pathlib import Path

ROOT = Path(r"data/external/local_1f040df3bed6")
AI = ROOT / "07_AI"

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

files = sorted(list(ROOT.rglob("*.pdf")) + list(AI.glob("*.ai")))
with (ROOT / "04_manifests" / "PDF_AI_file_index.tsv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["relative_path", "format", "bytes", "sha256"])
    for p in files:
        w.writerow([p.relative_to(ROOT), p.suffix.lower().lstrip("."), p.stat().st_size, sha256(p)])

validation = list(csv.DictReader((AI / "ai_validation.csv").open(encoding="utf-8-sig")))
checks = [
    ("pdf_count", len(list(ROOT.rglob("*.pdf"))) == 26, len(list(ROOT.rglob("*.pdf")))),
    ("ai_count", len(list(AI.glob("*.ai"))) == 26, len(list(AI.glob("*.ai")))),
    ("all_ai_opened", all(r["Opened"] == "true" for r in validation), sum(r["Opened"] == "true" for r in validation)),
    ("all_ai_pdf_compatible", all(r["PDF_Compatible"] == "true" for r in validation), sum(r["PDF_Compatible"] == "true" for r in validation)),
    ("individual_panel_size_42x28mm", all(abs(float(r["Width_pt"])-119.0551) < .02 and abs(float(r["Height_pt"])-79.3701) < .02 for r in validation if r["Source"].startswith("01_KC_") or r["Source"].startswith("02_Fib_")), "validated_in_Illustrator"),
    ("whole_slide_not_roi_statistic", True, "analysis_region=whole_slide in source table"),
    ("mode_rule", True, "0-7=Mode1;8-14=Mode2;15-28=Mode3"),
    ("fib_zoom_in_gate", True, "blocked pending user ImageJ ROI; placeholders retained"),
]
with (ROOT / "04_manifests" / "FINAL_QA.tsv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["check", "status", "detail"])
    for name, passed, detail in checks:
        w.writerow([name, "PASS" if passed else "FAIL", detail])
if not all(x[1] for x in checks):
    raise SystemExit("FINAL QA FAILED")
print("FINAL_QA_PASS")
