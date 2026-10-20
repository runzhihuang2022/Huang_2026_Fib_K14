from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
raise SystemExit(subprocess.call([sys.executable,str(ROOT/'script/run_figure.py'),'--figure',Path(__file__).resolve().parent.name]+sys.argv[1:],cwd=ROOT))
