from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path


BASE = Path(
    r"data/external/local_3413b1de7846"
    r"\20260813找曾老师指导\免疫荧光结果统计"
)
RAW_DIR = BASE / "Fib_K14"
LEGACY = BASE / "Figure 3D 细胞计数" / "Fib_K14_output"
OUT = BASE / "Figure 3D 细胞计数" / "Fib_K14_EpiFront_Control_ROI_gated_20260904"

DIRS = [
    "00_manifest",
    "01_converted_reference",
    "02_zoomout_review",
    "03_roi_handoff",
    "04_user_selected_zoomin",
    "05_cell_tables",
    "06_statistics",
    "07_publication_figures",
    "08_audit",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    for rel in DIRS:
        (OUT / rel).mkdir(parents=True, exist_ok=True)

    source_rows = read_tsv(LEGACY / "00_manifest" / "source_and_conversion_manifest.tsv")
    zoom_rows = read_tsv(LEGACY / "00_manifest" / "zoomout_manifest.tsv")
    zoom_by_sample = {row["sample_id"]: row for row in zoom_rows}

    sample_lookup = {
        "编号001A": "001A",
        "编号-003": "003",
        "编号-007": "007",
        "32-1": "32-1",
        "BS1": "BS1",
        "BW3": "BW3",
        "BW7": "BW7",
        "BW8": "BW8",
        "BW9": "BW9",
        "W14-1": "W14-1",
    }

    input_manifest: list[dict[str, object]] = []
    conversion_reference: list[dict[str, object]] = []
    sample_metadata: list[dict[str, object]] = []
    roi_manifest: list[dict[str, object]] = []
    zoomout_review: list[dict[str, object]] = []
    copied_assets: list[dict[str, object]] = []

    for row in source_rows:
        raw = Path(row["source_path"])
        if not raw.exists():
            raise FileNotFoundError(raw)
        sample_id = next((value for token, value in sample_lookup.items() if raw.name.startswith(token)), None)
        if sample_id is None:
            raise ValueError(f"Cannot map sample ID: {raw.name}")
        stat = raw.stat()
        expected_size = int(row["source_bytes"])
        size_match = stat.st_size == expected_size
        input_manifest.append(
            {
                "sample_id": sample_id,
                "source_path": str(raw),
                "format": ".kfbf",
                "size_bytes": stat.st_size,
                "mtime_local": stat.st_mtime,
                "sha256": row["source_sha256"],
                "hash_provenance": "reused_from_prior_conversion_audit_after_size_match",
                "size_match": str(size_match).lower(),
                "inclusion_state": "included_pre_roi",
            }
        )
        conversion_reference.append(
            {
                "sample_id": sample_id,
                "ome_tiff": row["output_path"],
                "ome_sha256": row["output_sha256"],
                "conversion_status": row["conversion_status"],
                "pixel_equality_status": row["pixel_equality_status"],
                "pixel_blocks_compared": row["pixel_blocks_compared"],
                "axes": row["axes"],
                "shape": row["shape"],
                "dtype": row["dtype"],
                "channel_names": row["channel_names"],
                "physical_size_x_um": row["physical_size_x_um"],
                "physical_size_y_um": row["physical_size_y_um"],
                "pyramid_levels": row["pyramid_levels"],
                "fiji_bioformats_status": "USER_OPEN_VERIFICATION_REQUIRED_THIS_RUN",
            }
        )
        sample_metadata.append(
            {
                "sample_id": sample_id,
                "input_file": raw.name,
                "subject_id": "AUTHOR_INPUT_NEEDED",
                "group": "AUTHOR_INPUT_NEEDED",
                "day_dpb": "AUTHOR_INPUT_NEEDED",
                "section_id": sample_id,
                "exclude": "false",
                "exclusion_reason": "",
            }
        )

        zrow = zoom_by_sample[sample_id]
        src_dapi = Path(zrow["dapi_reference"])
        src_png = Path(zrow["zoomout_png"])
        src_tif = Path(zrow["zoomout_tif"])
        for src in (src_dapi, src_png, src_tif):
            dst = OUT / "02_zoomout_review" / src.name
            shutil.copy2(src, dst)
            copied_assets.append(
                {
                    "sample_id": sample_id,
                    "source": str(src),
                    "copied_to": str(dst),
                    "bytes": dst.stat().st_size,
                    "sha256": sha256(dst),
                }
            )

        px_um = float(zrow["zoomout_pixel_size_x_um"])
        roi_px = round(500.0 / px_um)
        roi_manifest.append(
            {
                "sample_id": sample_id,
                "reference_image": str(OUT / "02_zoomout_review" / src_dapi.name),
                "reference_pixel_size_um": f"{px_um:.12g}",
                "suggested_500um_square_px": roi_px,
                "expected_roi_zip": f"{sample_id}_RoiSet.zip",
                "roi_zip_path": "",
                "coordinate_space": "DAPI_analysis_reference_source_orientation",
                "required_roi_names": f"{sample_id}_EF01... and {sample_id}_CTL01...",
                "minimum_rois": "2 EF + 2 CTL",
                "recommended_rois": "3-5 EF + 3-5 CTL",
                "orientation_review": "AUTHOR_INPUT_NEEDED",
                "reviewer": "",
                "status": "AUTHOR_INPUT_NEEDED",
                "note": "Choose by DAPI/anatomy, not PDGFRA or KRT14 intensity; match EF and CTL physical size within slide.",
            }
        )
        zoomout_review.append(
            {
                "sample_id": sample_id,
                "zoomout_merge": str(OUT / "02_zoomout_review" / src_png.name),
                "dapi_roi_reference": str(OUT / "02_zoomout_review" / src_dapi.name),
                "sample_mapping_review": "AUTHOR_INPUT_NEEDED",
                "channel_mapping_review": "PASS_DAPI_C1_PDGFRA_C2_P63_C3_KLF4_C4_KRT14_C5",
                "tissue_quality_review": "AUTHOR_INPUT_NEEDED",
                "orientation_review": "AUTHOR_INPUT_NEEDED",
                "approved_for_roi_drawing": "AUTHOR_INPUT_NEEDED",
                "reviewer_note": "",
            }
        )

    for name in ("Fib_K14_all_10_zoomout_A4.pdf", "Fib_K14_all_10_zoomout_A4.png"):
        src = LEGACY / "02_registered_zoomout" / name
        dst = OUT / "02_zoomout_review" / name
        shutil.copy2(src, dst)
        copied_assets.append(
            {"sample_id": "ALL", "source": str(src), "copied_to": str(dst), "bytes": dst.stat().st_size, "sha256": sha256(dst)}
        )

    write_tsv(
        OUT / "00_manifest" / "input_manifest.tsv",
        sorted(input_manifest, key=lambda x: str(x["sample_id"])),
        ["sample_id", "source_path", "format", "size_bytes", "mtime_local", "sha256", "hash_provenance", "size_match", "inclusion_state"],
    )
    write_tsv(
        OUT / "01_converted_reference" / "existing_ome_tiff_reference.tsv",
        sorted(conversion_reference, key=lambda x: str(x["sample_id"])),
        ["sample_id", "ome_tiff", "ome_sha256", "conversion_status", "pixel_equality_status", "pixel_blocks_compared", "axes", "shape", "dtype", "channel_names", "physical_size_x_um", "physical_size_y_um", "pyramid_levels", "fiji_bioformats_status"],
    )
    write_tsv(
        OUT / "00_manifest" / "sample_metadata_PENDING.tsv",
        sorted(sample_metadata, key=lambda x: str(x["sample_id"])),
        ["sample_id", "input_file", "subject_id", "group", "day_dpb", "section_id", "exclude", "exclusion_reason"],
    )
    write_tsv(
        OUT / "03_roi_handoff" / "epi_front_control_roi_manifest.tsv",
        sorted(roi_manifest, key=lambda x: str(x["sample_id"])),
        ["sample_id", "reference_image", "reference_pixel_size_um", "suggested_500um_square_px", "expected_roi_zip", "roi_zip_path", "coordinate_space", "required_roi_names", "minimum_rois", "recommended_rois", "orientation_review", "reviewer", "status", "note"],
    )
    write_tsv(
        OUT / "02_zoomout_review" / "zoomout_review_status.tsv",
        sorted(zoomout_review, key=lambda x: str(x["sample_id"])),
        ["sample_id", "zoomout_merge", "dapi_roi_reference", "sample_mapping_review", "channel_mapping_review", "tissue_quality_review", "orientation_review", "approved_for_roi_drawing", "reviewer_note"],
    )
    write_tsv(
        OUT / "08_audit" / "copied_asset_manifest.tsv",
        copied_assets,
        ["sample_id", "source", "copied_to", "bytes", "sha256"],
    )

    channel_rows = [
        {"channel_index_zero_based": 0, "source_name": "DAPI", "marker": "DAPI", "analysis_role": "nuclear segmentation and denominator", "compartment": "nucleus", "status": "locked"},
        {"channel_index_zero_based": 1, "source_name": "520", "marker": "PDGFRA", "analysis_role": "Fib_K14 positive marker", "compartment": "membrane/outer-cell ring", "status": "exploratory_without_independent_membrane_channel"},
        {"channel_index_zero_based": 2, "source_name": "570", "marker": "P63", "analysis_role": "not used in primary definition", "compartment": "nucleus", "status": "retained"},
        {"channel_index_zero_based": 3, "source_name": "620", "marker": "KLF4", "analysis_role": "not used in primary definition", "compartment": "nucleus", "status": "retained"},
        {"channel_index_zero_based": 4, "source_name": "690", "marker": "KRT14", "analysis_role": "Fib_K14 positive marker", "compartment": "cytoplasm", "status": "locked"},
    ]
    write_tsv(
        OUT / "00_manifest" / "channel_map.tsv",
        channel_rows,
        ["channel_index_zero_based", "source_name", "marker", "analysis_role", "compartment", "status"],
    )

    roi_manifest_path = OUT / "03_roi_handoff" / "epi_front_control_roi_manifest.tsv"
    metadata_path = OUT / "00_manifest" / "sample_metadata_PENDING.tsv"
    config = {
        "raw_dir": str(RAW_DIR),
        "input_dir": str(RAW_DIR),
        "output_dir": str(OUT),
        "codebook": str(metadata_path),
        "sample_metadata_tsv": str(metadata_path),
        "independent_unit": "slide",
        "subject_mapping_status": "AUTHOR_INPUT_NEEDED_before_inference",
        "analysis_region": "approved_roi",
        "approved_roi_manifest": str(roi_manifest_path),
        "pixel_size_um": 0.25284451246261597,
        "channel_map": {"DAPI": 0, "PDGFRA": 1, "P63": 2, "KLF4": 3, "KRT14": 4},
        "positive_definition": {
            "name": "PDGFRA+KRT14+ Fib_K14",
            "markers": ["PDGFRA", "KRT14"],
            "logic": "all",
            "denominator": "valid_DAPI_cells",
            "measurement": "PDGFRA_outer_cell_ring_and_KRT14_cytoplasm_assigned_to_same_cell_id",
            "threshold_rule": "pooled_outcome_blinded_locked_across_analysis_family; exploratory_without_external_controls",
        },
        "roi": {
            "size_um": 500,
            "regions": ["Epi-Front", "Control"],
            "minimum_per_region_per_slide": 2,
            "recommended_per_region_per_slide": "3-5",
            "selection_channel": "DAPI/anatomy",
            "selection_blinded_to": ["PDGFRA", "KRT14"],
            "archive_pattern": "{sample_id}_RoiSet.zip",
        },
        "orientation": {
            "anatomical_target": "epidermis_up_dermis_down",
            "one_transform_per_slide": True,
            "audit_required": True,
            "current_status": "AUTHOR_INPUT_NEEDED_source_orientation_review",
        },
        "panels": {
            "Fib_K14": {
                "markers": ["PDGFRA", "KRT14"],
                "primary_rule": "PDGFRA_and_KRT14_in_compartments_assigned_to_same_cell_id",
                "sensitivity_rule": "cell_expansion_and_threshold_sensitivity",
            }
        },
        "statistics_gate": "BLOCKED_UNTIL_ALL_APPROVED_ROI_ARCHIVES_AND_SUBJECT_MAPPING",
    }
    (OUT / "00_manifest" / "run_config_PRE_ROI.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")

    instructions = """Fib_K14 Epi-Front / Control ROI handoff（必须由作者在 ImageJ/Fiji 完成）

1. 打开 02_zoomout_review 中相应的 样本名_DAPI_analysis_reference.ome.tiff。
2. 用 Plugins > Bio-Formats > Bio-Formats Importer 打开；确认像素单位为 µm。
3. 先核对切片方向并在 manifest 的 orientation_review 填写 approved 或需要旋转的角度。
4. ROI 必须根据 DAPI 和组织解剖结构选择；选择时不要根据 PDGFRA 或 KRT14 信号强弱移动边界。
5. 每张切片同时圈选互不重叠的 Epi-Front 和 Control：每类至少2个，建议每类3–5个。
6. 建议使用相同的 500 × 500 µm 方形 ROI。每张切片对应的像素边长见 epi_front_control_roi_manifest.tsv。
7. Epi-Front 命名为 样本名_EF01、样本名_EF02……；Control 命名为 样本名_CTL01、样本名_CTL02……。
8. 同一张切片的两类 ROI 一并保存为 样本名_RoiSet.zip，例如 BW9_RoiSet.zip。
9. 把 zip 放入 03_roi_handoff/ROI_ZIPS_HERE，并在 manifest 填写 reviewer、roi_zip_path，将 status 改为 approved。
10. 请同时补全 00_manifest/sample_metadata_PENDING.tsv 中的 subject_id、group 和 day_dpb。若同一动物有多张切片，subject_id 必须相同。

Control 必须按预先规定的解剖标准选择，并在 reviewer_note 中写清楚定义；不能把“PDGFRA/KRT14信号低”本身作为 Control 的选择标准。

在 ROI 和元数据审核完成前，04–07 目录必须保持为空，不会生成细胞计数或统计结果。
"""
    (OUT / "03_roi_handoff" / "IMAGEJ_FIB_K14_EF_CONTROL_ROI_INSTRUCTIONS_CN.txt").write_text(instructions, encoding="utf-8")
    (OUT / "03_roi_handoff" / "ROI_ZIPS_HERE").mkdir(exist_ok=True)

    claim_rows = [
        {"claim": "whole-slide ZoomOut availability", "status": "SUPPORTED_FOR_REVIEW", "basis": "10 complete Fib_K14 source sections"},
        {"claim": "anatomical orientation", "status": "AUTHOR_INPUT_NEEDED", "basis": "source orientation retained; each slide requires review"},
        {"claim": "Epi-Front versus Control localization", "status": "NOT_TESTED", "basis": "approved manual ROI archives absent"},
        {"claim": "same-cell PDGFRA+KRT14+ Fib_K14", "status": "NOT_TESTED", "basis": "post-ROI compartment-aware segmentation required"},
        {"claim": "statistical significance", "status": "BLOCKED", "basis": "ROI gate and subject mapping not complete"},
        {"claim": "publication figures and editable AI", "status": "BLOCKED", "basis": "must be generated only after validated statistics"},
    ]
    write_tsv(OUT / "08_audit" / "claim_status.tsv", claim_rows, ["claim", "status", "basis"])

    forbidden_counts = {}
    for rel in ("04_user_selected_zoomin", "05_cell_tables", "06_statistics", "07_publication_figures"):
        forbidden_counts[rel] = sum(1 for p in (OUT / rel).rglob("*") if p.is_file())
    gate = {
        "status": "AUTHOR_INPUT_NEEDED",
        "sample_count": len(input_manifest),
        "all_source_sizes_match_prior_hash_manifest": all(row["size_match"] == "true" for row in input_manifest),
        "conversion_reuse": "10 prior lossless OME-TIFF conversions referenced; no raw files modified",
        "fiji_bioformats": "USER_OPEN_VERIFICATION_REQUIRED_THIS_RUN",
        "orientation": "AUTHOR_INPUT_NEEDED",
        "roi_archives_found": 0,
        "required_archive_pattern": "sample_id_RoiSet.zip",
        "required_regions": ["EF", "CTL"],
        "sample_metadata": "AUTHOR_INPUT_NEEDED subject_id/group/day_dpb",
        "forbidden_pre_roi_file_counts": forbidden_counts,
        "forbidden_outputs_empty": all(value == 0 for value in forbidden_counts.values()),
    }
    (OUT / "08_audit" / "PRE_ROI_GATE_AUDIT.json").write_text(json.dumps(gate, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"status": "AUTHOR_INPUT_NEEDED", "output": str(OUT), "samples": len(input_manifest), "forbidden_outputs_empty": gate["forbidden_outputs_empty"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
