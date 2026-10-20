from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main() -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_figure4_evidence_safe.py")],
        check=True,
    )


if __name__ == "__main__":
    main()
