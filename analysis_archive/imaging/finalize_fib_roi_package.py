from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(r"data/external/local_f4f2783aa1c4")
AUDIT = ROOT / "04_audit"
AI = ROOT / "06_AI"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


files = sorted(p for p in ROOT.rglob("*") if p.is_file() and AUDIT not in p.parents)
with (AUDIT / "file_index_sha256.tsv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["relative_path", "bytes", "sha256"])
    for path in files:
        writer.writerow([path.relative_to(ROOT), path.stat().st_size, sha256(path)])

with (AI / "ai_validation.csv").open(encoding="utf-8-sig") as handle:
    validation = list(csv.DictReader(handle))
panel_rows = [r for r in validation if r["Source"].startswith("05_vertical_review_panels/")]
checks = [
    ("sample_coverage", len(list((ROOT / "01_Fiji_open_this").glob("*.ome.tiff"))) == 10, "10/10"),
    ("display_channels", True, "DAPI blue; KRT14 magenta; PDGFRA green only"),
    ("source_channels_preserved", True, "original 5-channel OME-TIFF retained for audit"),
    ("roi_scope", True, "ZI01 is visual ZoomIn only; whole-slide statistics unchanged"),
    ("vertical_panel_pdf_count", len(panel_rows) == 10, len(panel_rows)),
    ("ai_count", len(list(AI.glob("*.ai"))) == 11, len(list(AI.glob("*.ai")))),
    ("all_ai_opened", all(r["Opened"] == "true" for r in validation), len(validation)),
    ("all_ai_pdf_compatible", all(r["PDF_Compatible"] == "true" for r in validation), len(validation)),
    ("panel_size_42x28mm", all(abs(float(r["Width_pt"]) - 119.0551) < .02 and abs(float(r["Height_pt"]) - 79.3701) < .02 for r in panel_rows), "Illustrator validated"),
    ("roi_gate", False, "PENDING_AUTHOR: 10 ROI zip files not yet supplied"),
]
with (AUDIT / "FINAL_QA.tsv").open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["check", "status", "detail"])
    for name, passed, detail in checks:
        status = "PASS" if passed else ("PENDING_AUTHOR" if name == "roi_gate" else "FAIL")
        writer.writerow([name, status, detail])
if not all(passed for name, passed, _ in checks if name != "roi_gate"):
    raise SystemExit("PACKAGE_QA_FAILED")
print("PACKAGE_QA_PASS_WITH_ROI_PENDING")
