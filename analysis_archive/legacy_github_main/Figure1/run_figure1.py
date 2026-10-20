from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"


def run(name: str) -> None:
    subprocess.run([sys.executable, str(SCRIPTS / name)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Rebuild the Figure 1 A4 evidence unit.")
    parser.add_argument(
        "--skip-scores",
        action="store_true",
        help="Reuse source_data/Figure1EF_registered11_scores and the section-level score CSV.",
    )
    args = parser.parse_args()
    for script in ("Figure1A_layout.py", "Figure1B_layout.py", "Figure1C_layout.py", "Figure1D_layout.py"):
        run(script)
    if not args.skip_scores:
        run("Figure1EF_compute_scores.py")
    run("Figure1E_layout.py")
    run("Figure1F_temporal.py")
    run("assemble_Figure1.py")


if __name__ == "__main__":
    main()
