from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"


def main() -> None:
    for name in ("compute_figure2_source_data.py", "build_figure2_package.py"):
        subprocess.run([sys.executable, str(SCRIPTS / name)], check=True)


if __name__ == "__main__":
    main()
