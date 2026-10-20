from pathlib import Path
import gzip
import hashlib
import math
import pandas as pd

OUT = Path(r"data/external/local_565f551d74b7")

per = pd.read_csv(OUT / "04_per_slide/per_slide.tsv", sep="\t", dtype={"sample_id": str, "subject_id": str})
checks = []

def add(check, observed, expected, status, note=""):
    checks.append({"check": check, "observed": observed, "expected": expected, "status": status, "note": note})

add("per-slide rows", len(per), 21, "PASS" if len(per) == 21 else "FAIL")
add("KC rows", int((per.family == "KC").sum()), 11, "PASS" if (per.family == "KC").sum() == 11 else "FAIL")
add("Fib_K14 rows", int((per.family == "Fib_K14").sum()), 10, "PASS" if (per.family == "Fib_K14").sum() == 10 else "FAIL")
add("whole-slide inference rows", int(per.inference_included.sum()), 20, "PASS" if per.inference_included.sum() == 20 else "FAIL")
cut = per.loc[per.sample_id == "CUT_007"]
add("CUT_007 exclusion", f"region={cut.analysis_region.iloc[0]}; included={cut.inference_included.iloc[0]}", "cropped_field; False", "PASS" if len(cut) == 1 and cut.analysis_region.iloc[0] == "cropped_field" and not bool(cut.inference_included.iloc[0]) else "FAIL")

reconcile_failures = []
active_files = []
for family in ["KC", "Fib_K14"]:
    for fp in sorted((OUT / "03_per_cell" / family).glob("*_per_cell.tsv.gz")):
        active_files.append(fp)
        sid = fp.name.removesuffix("_per_cell.tsv.gz")
        n = pos = 0
        for chunk in pd.read_csv(fp, sep="\t", compression="gzip", chunksize=200000):
            n += len(chunk)
            pos += int(chunk["co_positive"].astype(bool).sum())
        row = per[(per.family == family) & (per.sample_id == sid)]
        if len(row) != 1 or int(row.denominator_n.iloc[0]) != n or int(row.numerator_n.iloc[0]) != pos:
            reconcile_failures.append(f"{family}:{sid} table={pos}/{n}")
add("active per-cell tables", len(active_files), 21, "PASS" if len(active_files) == 21 else "FAIL")
add("per-cell/per-slide reconciliation", len(reconcile_failures), 0, "PASS" if not reconcile_failures else "FAIL", "; ".join(reconcile_failures))

ratio_calc = 100.0 * per.numerator_n / per.denominator_n
ratio_err = float((ratio_calc - per.ratio_pct).abs().max())
add("ratio arithmetic max abs error", ratio_err, "<=1e-10", "PASS" if ratio_err <= 1e-10 else "FAIL")

sens = pd.read_csv(OUT / "05_sensitivity/segmentation_threshold_sensitivity.tsv", sep="\t", dtype={"sample_id": str})
add("sensitivity rows", len(sens), 189, "PASS" if len(sens) == 189 else "FAIL", "21 slides x 3 cell expansions x 3 threshold multipliers")
add("sensitivity design coverage", f"expansions={sorted(sens.expansion.unique())}; multipliers={sorted(sens.threshold_multiplier.unique())}", "[2,3,5] and [0.9,1.0,1.1]", "PASS" if set(sens.expansion) == {2,3,5} and set(round(float(x),1) for x in sens.threshold_multiplier) == {0.9,1.0,1.1} else "FAIL")

overlays = list((OUT / "02_qc/overlays").rglob("*_strict_overlay.png"))
add("overlay QC images", len(overlays), 21, "PASS" if len(overlays) == 21 else "FAIL")

thresholds = pd.read_csv(OUT / "02_thresholds/locked_thresholds.tsv", sep="\t")
exploratory = thresholds.astype(str).apply(lambda c: c.str.contains("EXPLORATORY", case=False, na=False)).any(axis=1).all()
add("threshold status", "all exploratory" if exploratory else "mixed", "all exploratory because controls absent", "PASS" if exploratory else "FAIL")
add("group mapping", "007 and 32-1 flagged REVIEW_REQUIRED", "manual confirmation required", "REVIEW_REQUIRED", "does not invalidate arithmetic; affects biological grouping")
add("segmentation biological QC", "all slides flagged REVIEW_REQUIRED", "manual overlay review", "REVIEW_REQUIRED", "automated segmentation is not a substitute for blinded biological review")
add("overall release status", "EXPLORATORY_REVIEW_REQUIRED", "no failed computational checks; controls and manual review still required", "REVIEW_REQUIRED")

qa = pd.DataFrame(checks)
(OUT / "08_audit").mkdir(exist_ok=True)
qa.to_csv(OUT / "08_audit/final_qa.tsv", sep="\t", index=False)

key_patterns = [
    "README_请先阅读.md", "00_manifest/*.tsv", "00_manifest/*.json",
    "02_thresholds/*.tsv", "04_per_slide/*.tsv", "05_statistics/*.tsv",
    "05_statistics/*.json", "05_sensitivity/*.tsv", "06_figures/*.*",
    "07_excel/*.xlsx", "08_audit/final_qa.tsv"
]
files = []
for pattern in key_patterns:
    files.extend(OUT.glob(pattern))
rows = []
for fp in sorted(set(f for f in files if f.is_file())):
    h = hashlib.sha256()
    with fp.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    rows.append({"relative_path": fp.relative_to(OUT).as_posix(), "size_bytes": fp.stat().st_size, "sha256": h.hexdigest()})
pd.DataFrame(rows).to_csv(OUT / "08_audit/key_file_sha256.tsv", sep="\t", index=False)
print(qa.to_string(index=False))
