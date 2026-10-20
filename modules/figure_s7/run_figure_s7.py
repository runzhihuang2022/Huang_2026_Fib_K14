from pathlib import Path
import subprocess, sys, shutil, json, hashlib
root=Path(__file__).resolve().parent
sources=json.loads((root/'provenance/source_manifest.json').read_text(encoding='utf8'))
for row in sources:
    p=root/row['copy'].replace('\\','/')
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'], str(p)
subprocess.run([sys.executable,str(root/'scripts/build_figure_s7.py'),str(root)],check=True)
for p in (root/'outputs').glob('Figure_S7*'):
    if not p.is_file():continue
    dst='panels' if p.stem.startswith(tuple('Figure_S7'+x+'_' for x in 'ABCDEF')) else 'assembled'
    (root/'outputs'/dst).mkdir(exist_ok=True)
    shutil.copy2(p,root/'outputs'/dst/p.name)
    if dst=='assembled' and p.suffix in ['.pdf','.svg']:
        (root/'outputs/vector').mkdir(exist_ok=True)
        shutil.copy2(p,root/'outputs/vector'/p.name)
shutil.copy2(root/'outputs/figure_QA.json',root/'QC/figure_QA.json')
shutil.copy2(root/'outputs/Figure_S7_V12.pdf',root/'supplements/figures/Figure_S7_V12.pdf')
print('Figure S7 rebuilt from verified processed data. Differential-expression and GSEA models were not refitted.')
