from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\001A_reprocess_20260820"
)
LANE_ROOT = ROOT / "02_registration_and_lanes"
MANIFEST = LANE_ROOT / "00_audit" / "registration_and_lane_manifest.json"
FIXED_OME = ROOT / "01_ome_tiff" / "KC_Mig_编号001A.ome.tiff"
OUTPUT = ROOT / "03_imagej_double_positive_spatial" / "01_registered_analysis_channels"


def write_gray(path: Path, array: np.ndarray, pixel_size_um: float) -> None:
    resolution = (10000.0 / pixel_size_um, 10000.0 / pixel_size_um)
    tifffile.imwrite(
        path,
        np.asarray(array, dtype=np.uint8),
        bigtiff=True,
        photometric="minisblack",
        compression="deflate",
        predictor=True,
        tile=(256, 256),
        resolution=resolution,
        resolutionunit="CENTIMETER",
        metadata={
            "PhysicalSizeX": pixel_size_um,
            "PhysicalSizeY": pixel_size_um,
            "PhysicalSizeUnit": "um",
        },
    )


def main() -> None:
    report = json.loads(MANIFEST.read_text(encoding="utf-8"))
    level = int(report["output_level"])
    x0, y0, x1, y1 = map(int, report["crop_xyxy_level2"])
    pixel_size_um = float(report["pixel_size_um"])
    OUTPUT.mkdir(parents=True, exist_ok=True)

    with tifffile.TiffFile(FIXED_OME) as tif:
        series = tif.series[0].levels[level]
        fixed = {
            "KC_DAPI_fixed.tif": np.asarray(series.asarray(key=0), dtype=np.uint8)[y0:y1, x0:x1],
            "KC_SOX9_fixed.tif": np.asarray(series.asarray(key=2), dtype=np.uint8)[y0:y1, x0:x1],
            "KC_KRT6B_fixed.tif": np.asarray(series.asarray(key=3), dtype=np.uint8)[y0:y1, x0:x1],
        }

    for name, array in fixed.items():
        write_gray(OUTPUT / name, array, pixel_size_um)
        print(f"WROTE {name}: {array.shape[1]}x{array.shape[0]}", flush=True)

    aligned_source = LANE_ROOT / "01_aligned_channels"
    source_map = {
        "Fib_DAPI_aligned.tif": aligned_source / "DAPI_FibK14_aligned_to_KCMig.tif",
        "Fib_PDGFRa_aligned.tif": aligned_source / "PDGFRa_aligned_to_KCMig.tif",
        "Fib_KRT14_aligned.tif": aligned_source / "KRT14_aligned_to_KCMig.tif",
    }
    for output_name, source in source_map.items():
        with tifffile.TiffFile(source) as tif:
            array = np.asarray(tif.asarray(), dtype=np.uint8)
        write_gray(OUTPUT / output_name, array, pixel_size_um)
        print(f"WROTE {output_name}: {array.shape[1]}x{array.shape[0]}", flush=True)

    audit = {
        "status": "pass",
        "reference_section": "KC_Mig_编号001A",
        "moving_section": "Fib_K14_编号001A aligned to KC_Mig_编号001A",
        "shape_yx": report["output_shape_yx"],
        "pixel_size_um": pixel_size_um,
        "crop_xyxy_level2": report["crop_xyxy_level2"],
        "channel_mapping": {
            "KC": ["DAPI", "SOX9(570)", "KRT6B(690)"],
            "Fib": ["DAPI", "PDGFRa(520)", "KRT14(690)"],
        },
        "note": "No source files modified. KC has four true channels in the parent OME-TIFF, not five.",
    }
    (OUTPUT / "channel_preparation_manifest.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
