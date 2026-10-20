from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "scripts" / "plot_19dpb_fibk14_regional_features.py"


def main() -> None:
    output = ROOT / "outputs" / "panels" / "FibK14_19dpb_regional_expression"
    h5ad = os.environ.get("FIGURE5_STEREO_H5AD")
    metadata = os.environ.get(
        "FIGURE5_METADATA",
        str(ROOT / "source_data" / "19dpb_metadata_wound_axes.tsv.gz"),
    )
    mask_dir = os.environ.get("FIGURE5_MASK_DIR")
    if not h5ad or not mask_dir:
        raise SystemExit(
            "Set FIGURE5_STEREO_H5AD and FIGURE5_MASK_DIR before extraction. "
            "Large H5AD and mask assets are intentionally excluded from Git."
        )
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--mode",
            "extract",
            "--h5ad",
            h5ad,
            "--metadata",
            metadata,
            "--mask-dir",
            mask_dir,
            "--output",
            str(output),
        ],
        check=True,
    )
    subprocess.run(
        [sys.executable, str(SCRIPT), "--mode", "plot", "--output", str(output)],
        check=True,
    )


if __name__ == "__main__":
    main()
