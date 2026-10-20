from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parent
runpy.run_path(str(ROOT / "scripts" / "build_figure_s1_efg_v07.py"), run_name="__main__")

print(
    "Lower panels generated. Open Adobe Illustrator and run "
    "`scripts/merge_current_ai_with_efg_v07.jsx` to merge E-G into the "
    "author-supplied v06 AI while retaining A-D."
)
