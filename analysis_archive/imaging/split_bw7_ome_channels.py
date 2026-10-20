from __future__ import annotations

import gc
import json
from pathlib import Path

import numpy as np
import tifffile


ROOT = Path(
    r"data/external/local_0fc1527f456f"
    r"\免疫荧光结果统计\Figure 3G overlap\BW7_reprocess_20260820\01_ome_tiff"
)
OUTPUT = ROOT / "02_split_channels"
PIXEL_SIZE_UM = 0.25284451246261597
PANELS = [
    (
        ROOT / "Fib_K14_BW7 PDGFRα 520+P63 570+KLF4 620+KRT14 690.ome.tiff",
        "Fib_K14_BW7",
        ["DAPI", "PDGFRa", "P63", "KLF4", "KRT14"],
    ),
    (
        ROOT / "KC_Mig_涤ⅡBW7.ome.tiff",
        "KC_Mig_BW7",
        ["DAPI", "MMP3", "SOX9", "KRT6B"],
    ),
]


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for ome_path, panel_name, markers in PANELS:
        with tifffile.TiffFile(ome_path) as tif:
            level = tif.series[0].levels[0]
            for index, marker in enumerate(markers):
                array = np.asarray(level.asarray(key=index), dtype=np.uint8)
                target = OUTPUT / f"{panel_name}_C{index + 1}_{marker}.tif"
                tifffile.imwrite(
                    target,
                    array,
                    bigtiff=True,
                    photometric="minisblack",
                    compression="deflate",
                    compressionargs={"level": 6},
                    predictor=True,
                    tile=(256, 256),
                    resolution=(10000.0 / PIXEL_SIZE_UM, 10000.0 / PIXEL_SIZE_UM),
                    resolutionunit="CENTIMETER",
                    metadata={
                        "PhysicalSizeX": PIXEL_SIZE_UM,
                        "PhysicalSizeY": PIXEL_SIZE_UM,
                        "PhysicalSizeUnit": "um",
                        "Marker": marker,
                    },
                )
                manifest.append(
                    {
                        "panel": panel_name,
                        "channel_index_1based": index + 1,
                        "marker": marker,
                        "path": str(target),
                        "shape_yx": list(array.shape),
                        "dtype": str(array.dtype),
                        "pixel_size_um": PIXEL_SIZE_UM,
                        "bytes": target.stat().st_size,
                    }
                )
                print(f"WROTE {target.name} {array.shape[1]}x{array.shape[0]}", flush=True)
                del array
                gc.collect()
    (OUTPUT / "split_channel_manifest.json").write_text(
        json.dumps({"status": "pass", "count": len(manifest), "files": manifest}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
