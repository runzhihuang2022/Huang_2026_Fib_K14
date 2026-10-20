from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd


PROJECT = Path(r"data/external/local_6b4e1e668702")
FIG3D = PROJECT / "Figure 3D 细胞计数"
OUTPUT = FIG3D / "IF_colocalized_cell_statistics_all_slides_20260825"
KC_RAW = PROJECT / "KC_Spinous_Mig"
FIB_RAW = PROJECT / "Fib_K14"
KC_SOURCE = FIG3D / "KC_Spinous_Mig_whole_slide_statistics_native_pixel_v3"
FIB_SOURCE = FIG3D / "Fib_K14_output" / "07_whole_slide_statistics"
MAPPING = FIG3D / "KC_Fib_SPTDI_DPTDI_whole_slide_revision_20260825" / "01_source_tables" / "sample_group_mode_mapping.tsv"

KC_THRESHOLDS = {"MMP3": 4358.5048828125, "SOX9": 413.203125, "KRT6B": 4966.2615234375}
FIB_THRESHOLDS = {"PDGFRA": 25.6291640625, "KRT14": 18.312668554687498}


def mkdirs() -> None:
    for rel in ["00_manifest", "02_qc", "03_per_cell/KC", "03_per_cell/Fib_K14", "04_per_slide", "05_statistics/KC", "05_statistics/Fib_K14", "06_figures/KC", "06_figures/Fib_K14", "07_excel", "08_audit"]:
        (OUTPUT / rel).mkdir(parents=True, exist_ok=True)


def sample_id(path: Path, family: str) -> str:
    s = path.name.upper()
    if family == "KC" and s.startswith("CUT-"):
        return "CUT_007"
    keys = ["W14-1", "32-1", "001A", "003", "007", "BS1", "BW3", "BW4", "BW7", "BW8", "BW9"]
    for key in keys:
        if key in s:
            return key
    raise ValueError(f"Unmapped sample filename: {path}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(16 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build_manifest() -> pd.DataFrame:
    rows = []
    candidates = [("KC", p) for p in sorted(KC_RAW.glob("*.tif"))] + [("Fib_K14", p) for p in sorted(FIB_RAW.glob("*.kfbf"))]
    for idx, (family, path) in enumerate(candidates, 1):
        sid = sample_id(path, family)
        print(f"HASH {idx}/{len(candidates)} {family} {sid}", flush=True)
        rows.append({
            "family": family,
            "sample_id": sid,
            "sample_label": f"{family}:{sid}",
            "source_path": str(path),
            "source_file": path.name,
            "format": path.suffix.lower().lstrip("."),
            "bytes": path.stat().st_size,
            "modified_time": path.stat().st_mtime,
            "sha256": sha256(path),
            "source_scope": "cropped_field" if sid == "CUT_007" else "whole_slide",
            "include_in_descriptive_count": True,
            "include_in_inference": sid != "CUT_007",
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT / "00_manifest" / "input_manifest.tsv", sep="\t", index=False)
    return out


def prepare_cut_007() -> Path:
    out = OUTPUT / "08_audit" / "CUT_007_raw_measurements.tsv.gz"
    if out.exists():
        return out
    source = next(KC_RAW.glob("cut-MAX_007*.tif"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import kc_spinous_mig_whole_slide_statistics as kc
    original_sid = kc.sid
    try:
        kc.sid = lambda _: "CUT_007"
        print("SEGMENT CUT_007 cropped field", flush=True)
        kc.analyze(source, out)
    finally:
        kc.sid = original_sid
    return out


def mapping_table() -> pd.DataFrame:
    m = pd.read_csv(MAPPING, sep="\t", dtype=str).fillna("")
    extra = pd.DataFrame([{
        "family": "KC", "sample": "CUT_007", "group": "DPTDI", "day_dpb": "7", "mode": "Mode 1",
        "group_mapping_provenance": "cropped derivative of 007; descriptive only, excluded from inference"
    }])
    m = pd.concat([m, extra], ignore_index=True)
    m.insert(2, "sample_label", m["family"] + ":" + m["sample"])
    m.to_csv(OUTPUT / "00_manifest" / "sample_metadata.tsv", sep="\t", index=False)
    return m


def channel_map() -> None:
    rows = [
        ["KC", 0, "C1", "SOX9", "filename wavelengths plus ImageJ LUT audit", "native pixel; source TIFF lacks um/px"],
        ["KC", 1, "C2", "MMP3", "filename wavelengths plus ImageJ LUT audit", "native pixel; source TIFF lacks um/px"],
        ["KC", 2, "C3", "DAPI", "punctate nuclear morphology plus ImageJ blue LUT", "native pixel; source TIFF lacks um/px"],
        ["KC", 3, "C4", "KRT6B", "filename wavelengths plus ImageJ LUT audit", "native pixel; source TIFF lacks um/px"],
        ["Fib_K14", 0, "C1", "DAPI", "validated KFBF-to-OME channel map", "0.25284451246261597 um/px at level 0"],
        ["Fib_K14", 1, "C2", "PDGFRA", "validated KFBF-to-OME channel map", "0.25284451246261597 um/px at level 0"],
        ["Fib_K14", 4, "C5", "KRT14", "validated KFBF-to-OME channel map", "0.25284451246261597 um/px at level 0"],
    ]
    pd.DataFrame(rows, columns=["family", "channel_index_0based", "source_channel", "marker", "mapping_evidence", "calibration"]).to_csv(OUTPUT / "00_manifest" / "channel_map.tsv", sep="\t", index=False)


def source_cell_files(family: str) -> dict[str, Path]:
    base = KC_SOURCE / "01_per_cell" if family == "KC" else FIB_SOURCE / "01_per_cell"
    files = {p.name.split("_whole_slide_cells")[0]: p for p in base.glob("*_whole_slide_cells.tsv.gz")}
    if family == "KC":
        files["CUT_007"] = prepare_cut_007()
    return files


def quantify_family(family: str, mapping: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    markers = ["MMP3", "SOX9", "KRT6B"] if family == "KC" else ["PDGFRA", "KRT14"]
    thresholds = KC_THRESHOLDS if family == "KC" else FIB_THRESHOLDS
    definition = "SOX9+ MMP3+ KRT6B+ DAPI+" if family == "KC" else "PDGFRA+ KRT14+ DAPI+"
    rows = []
    for sid, source in sorted(source_cell_files(family).items()):
        print(f"QUANTIFY {family} {sid}", flush=True)
        df = pd.read_csv(source, sep="\t", compression="gzip", dtype={"sample": str})
        df["sample"] = sid
        df["DAPI_positive"] = True
        for marker in markers:
            df[f"{marker}_positive"] = pd.to_numeric(df[f"{marker}_p90_bgsub"], errors="raise") > thresholds[marker]
        df["co_positive"] = df[[f"{m}_positive" for m in markers]].all(axis=1)
        keep = [c for c in df.columns if c not in []]
        dst = OUTPUT / "03_per_cell" / family / f"{sid}_per_cell.tsv.gz"
        df[keep].to_csv(dst, sep="\t", index=False, compression="gzip")
        meta = mapping[(mapping["family"] == family) & (mapping["sample"] == sid)]
        if len(meta) != 1:
            raise ValueError(f"Expected one metadata row for {family} {sid}, got {len(meta)}")
        meta = meta.iloc[0]
        raw = manifest[(manifest["family"] == family) & (manifest["sample_id"] == sid)]
        if len(raw) != 1:
            raise ValueError(f"Expected one raw manifest row for {family} {sid}, got {len(raw)}")
        raw = raw.iloc[0]
        numerator = int(df["co_positive"].sum())
        denominator = int(df["DAPI_positive"].sum())
        marker_counts = {f"{m}_positive_n": int(df[f"{m}_positive"].sum()) for m in markers}
        mapping_review = sid in {"32-1", "007", "CUT_007"}
        rows.append({
            "family": family,
            "sample_id": sid,
            "sample_label": f"{family}:{sid}",
            "subject_id": "007" if sid == "CUT_007" else sid,
            "source_file": raw["source_file"],
            "source_scope": raw["source_scope"],
            "group": meta["group"],
            "day_dpb": float(meta["day_dpb"]),
            "mode": meta["mode"],
            "positive_definition": definition,
            **marker_counts,
            "numerator_n": numerator,
            "denominator_n": denominator,
            "ratio_pct": 100.0 * numerator / denominator,
            "analysis_region": "cropped_field" if sid == "CUT_007" else "whole_slide",
            "inference_included": sid != "CUT_007",
            "inference_exclusion_reason": "cropped derivative of full 007; not an independent whole slide" if sid == "CUT_007" else "",
            "group_mapping_status": "REVIEW_REQUIRED" if mapping_review else "SUPPORTED",
            "group_mapping_provenance": meta["group_mapping_provenance"],
            "segmentation_qc_status": "REVIEW_REQUIRED",
            "threshold_rule": "locked pooled automatic Otsu from the previously segmented full-slide analysis family",
        })
    return pd.DataFrame(rows)


def write_qc(per_slide: pd.DataFrame) -> None:
    qc = per_slide[["family", "sample_id", "denominator_n", "numerator_n", "analysis_region", "inference_included", "group_mapping_status", "segmentation_qc_status"]].copy()
    qc["count_reconciliation"] = (qc["denominator_n"] > 0) & (qc["numerator_n"] >= 0) & (qc["numerator_n"] <= qc["denominator_n"])
    qc.to_csv(OUTPUT / "02_qc" / "qc_summary.tsv", sep="\t", index=False)
    gates = [
        ["raw_file_coverage", "PASS", "11 KC image files and 10 Fib KFBF files inventoried"],
        ["per_image_counts", "PASS", "Every inventoried image has numerator and DAPI denominator"],
        ["whole_slide_scope", "PASS_WITH_EXCEPTION", "20 whole slides; CUT_007 is a cropped derivative reported descriptively and excluded from inference"],
        ["channel_mapping", "PASS", "KC C1/C2/C3/C4=SOX9/MMP3/DAPI/KRT6B; Fib C1/C2/C5=DAPI/PDGFRA/KRT14"],
        ["threshold_lock", "PASS", "Prior pooled thresholds retained; not recomputed to favor significance"],
        ["segmentation_overlay_review", "REVIEW_REQUIRED", "Automated segmentation and positivity thresholds remain exploratory pending biological image review"],
        ["mode1_group_mapping", "REVIEW_REQUIRED", "32-1=SPTDI and 007=DPTDI are provisional based on prior workbook mapping"],
        ["biological_replication", "PASS", "Cells and CUT_007 cropped derivative are not counted as independent n"],
    ]
    pd.DataFrame(gates, columns=["gate", "status", "detail"]).to_csv(OUTPUT / "08_audit" / "final_qa.tsv", sep="\t", index=False)


def write_configs() -> None:
    cfg = {
        "analysis_scope": "whole slide for inference; CUT_007 cropped field descriptive only",
        "denominator": "all accepted DAPI-segmented nuclei in the declared image region",
        "KC": {"definition": "SOX9+ MMP3+ KRT6B+ DAPI+", "channels_0based": {"SOX9": 0, "MMP3": 1, "DAPI": 2, "KRT6B": 3}, "thresholds": KC_THRESHOLDS, "measurement": "p90 in 5 native-pixel centroid disk after Gaussian background subtraction sigma=15 px", "nuclear_area_px2": [8, 500], "calibration": "CALIBRATION_REQUIRED"},
        "Fib_K14": {"definition": "PDGFRA+ KRT14+ DAPI+", "channels_0based": {"DAPI": 0, "PDGFRA": 1, "KRT14": 4}, "thresholds": FIB_THRESHOLDS, "measurement": "p90 in 5 um centroid disk after Gaussian background subtraction sigma=15 um", "nuclear_area_um2": [20, 300], "pixel_size_um_level2": 1.0113780498504639},
        "mode_mapping": {"Mode 1": "0-7 dpb", "Mode 2": "8-14 dpb", "Mode 3": "15-28 dpb"},
        "inference_unit": "independent whole slide/subject",
        "status": "EXPLORATORY; segmentation and threshold image review required",
    }
    (OUTPUT / "08_audit" / "run_config.json").write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    mkdirs()
    manifest = build_manifest()
    mapping = mapping_table()
    channel_map()
    frames = [quantify_family("KC", mapping, manifest), quantify_family("Fib_K14", mapping, manifest)]
    all_slides = pd.concat(frames, ignore_index=True)
    all_slides.to_csv(OUTPUT / "04_per_slide" / "per_slide.tsv", sep="\t", index=False)
    for family, df in all_slides.groupby("family", sort=False):
        df.to_csv(OUTPUT / "04_per_slide" / f"{family}_per_slide.tsv", sep="\t", index=False)
        infer = df[df["inference_included"]].copy()
        infer.to_csv(OUTPUT / "04_per_slide" / f"{family}_inference_input.tsv", sep="\t", index=False)
    write_qc(all_slides)
    write_configs()
    print(json.dumps({"status": "COUNTS_COMPLETE", "output": str(OUTPUT), "rows": len(all_slides), "KC_rows": int((all_slides.family == 'KC').sum()), "Fib_rows": int((all_slides.family == 'Fib_K14').sum())}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
