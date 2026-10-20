"""One-command entry point for rebuilding Figure S1 v11 source panels."""

from pathlib import Path
import runpy
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "scripts"))
runpy.run_path(str(HERE / "scripts" / "build_figure_s1_ef_v11.py"), run_name="__main__")
