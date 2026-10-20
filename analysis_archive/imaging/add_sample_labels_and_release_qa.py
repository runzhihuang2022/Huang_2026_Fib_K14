from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(r"data/external/local_5cc3be176607")


def add_label(path: Path, id_col: str) -> None:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "sample_label" in df.columns:
        df = df.drop(columns=["sample_label"])
    pos = df.columns.get_loc(id_col) + 1
    df.insert(pos, "sample_label", df["family"] + ":" + df[id_col])
    df.to_csv(path, sep="\t", index=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


for rel, col in [
    ("00_manifest/input_manifest.tsv", "sample_id"),
    ("00_manifest/sample_metadata.tsv", "sample"),
    ("04_per_slide/per_slide.tsv", "sample_id"),
    ("04_per_slide/KC_per_slide.tsv", "sample_id"),
    ("04_per_slide/KC_inference_input.tsv", "sample_id"),
    ("04_per_slide/Fib_K14_per_slide.tsv", "sample_id"),
    ("04_per_slide/Fib_K14_inference_input.tsv", "sample_id"),
    ("02_qc/qc_summary.tsv", "sample_id"),
]:
    add_label(ROOT / rel, col)

per_slide = pd.read_csv(ROOT / "04_per_slide/per_slide.tsv", sep="\t")
coverage = pd.DataFrame([
    ["raw_inputs", 21, len(pd.read_csv(ROOT / "00_manifest/input_manifest.tsv", sep="\t")), "PASS"],
    ["KC_images", 11, int((per_slide["family"] == "KC").sum()), "PASS"],
    ["Fib_K14_images", 10, int((per_slide["family"] == "Fib_K14").sum()), "PASS"],
    ["whole_slide_results", 20, int((per_slide["analysis_region"] == "whole_slide").sum()), "PASS"],
    ["cropped_derivatives", 1, int((per_slide["analysis_region"] == "cropped_field").sum()), "PASS"],
    ["per_cell_tables", 21, len(list((ROOT / "03_per_cell").rglob("*.tsv.gz"))), "PASS"],
])
coverage.columns = ["item", "expected_n", "observed_n", "status"]
coverage["status"] = coverage.apply(lambda r: "PASS" if r.expected_n == r.observed_n else "FAIL", axis=1)
coverage.to_csv(ROOT / "08_audit/coverage.tsv", sep="\t", index=False)

reconciliation = []
for row in per_slide.itertuples(index=False):
    cell_path = ROOT / "03_per_cell" / row.family / f"{row.sample_id}_per_cell.tsv.gz"
    cells = pd.read_csv(cell_path, sep="\t", compression="gzip", usecols=["DAPI_positive", "co_positive"])
    observed_denominator = int(cells["DAPI_positive"].astype(bool).sum())
    observed_numerator = int(cells["co_positive"].astype(bool).sum())
    expected_ratio = 100.0 * observed_numerator / observed_denominator
    passed = (
        observed_denominator == int(row.denominator_n)
        and observed_numerator == int(row.numerator_n)
        and abs(expected_ratio - float(row.ratio_pct)) < 1e-10
    )
    reconciliation.append({
        "family": row.family,
        "sample_id": row.sample_id,
        "per_cell_denominator_n": observed_denominator,
        "per_slide_denominator_n": int(row.denominator_n),
        "per_cell_numerator_n": observed_numerator,
        "per_slide_numerator_n": int(row.numerator_n),
        "status": "PASS" if passed else "FAIL",
    })
pd.DataFrame(reconciliation).to_csv(ROOT / "08_audit/count_reconciliation.tsv", sep="\t", index=False)

rows = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or path.name == "file_index_sha256.tsv":
        continue
    rows.append({
        "relative_path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    })
pd.DataFrame(rows).to_csv(ROOT / "08_audit/file_index_sha256.tsv", sep="\t", index=False)
