"""Run real analyses in a separate output tree; source inputs remain unchanged."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.metadata as md,json,os,platform,shutil,subprocess,sys,time

ROOT=Path(__file__).resolve().parent
CORE=['scripts/audit_table_s1.py','scripts/rebuild_current_histology.py','scripts/reanalyze_bulk.py','scripts/rebuild_atlas_audits.py','scripts/spatial_protein/rebuild.py','scripts/rebuild_current_supplements.py']

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--modules',default='core,figure_s7,mechanism,derived',help='Comma-separated: core,figure_s7,mechanism,mechanism_r,gsea,derived')
 ap.add_argument('--output-dir',type=Path,default=ROOT/'outputs/rerun')
 ap.add_argument('--rscript',default=shutil.which('Rscript'))
 ap.add_argument('--external-gmt',type=Path)
 ap.add_argument('--list',action='store_true')
 args=ap.parse_args()
 if args.list:
  print('core: six processed-data/count-model tasks\nfigure_s7: rebuild V12 S7\nmechanism: patient aggregation, correlations, six mechanism plots\nmechanism_r: DESeq2 and edgeR from counts\ngsea: fgsea from locked full gene-set library\nderived: CellOracle source-table checks and cross-species/ATAC summaries');return
 out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
 if (out/'execution_report.json').exists():raise SystemExit('Choose a new output directory to preserve the previous run.')
 env=dict(os.environ,PYTHONIOENCODING='utf-8',PYTHONUTF8='1',MPLBACKEND='Agg')
 env['MPLCONFIGDIR']=str(out/'matplotlib_cache')
 records=[];modules=args.modules.split(',');start_all=time.time()
 for name in modules:
  if name not in ['core','figure_s7','mechanism','mechanism_r','gsea','derived']:raise SystemExit('Unknown module '+name)
  work=out/name;work.mkdir(parents=True,exist_ok=True)
  if name=='derived':commands=[[sys.executable,str(ROOT/'scripts/rebuild_derived.py'),str(work)]]
  else:
   src=ROOT/'modules'/('mechanism' if name in ['mechanism_r','gsea'] else name)
   for item in src.iterdir():
    if item.name in ['outputs','QC','validation','__pycache__']:continue
    target=work/item.name
    if item.is_dir():shutil.copytree(item,target,dirs_exist_ok=True)
    else:shutil.copy2(item,target)
   for d in ['outputs','validation','figures','reports','logs','analysis','tables']:(work/d).mkdir(exist_ok=True)
   if name=='core':commands=[[sys.executable,str(work/s)]for s in CORE]
   elif name=='figure_s7':commands=[[sys.executable,str(work/'scripts/build_figure_s7.py'),str(work)]]
   elif name=='mechanism':commands=[[sys.executable,str(work/'scripts/analyze_human.py')],[sys.executable,str(work/'scripts/make_figures.py')]]
   elif name=='mechanism_r':
    if not args.rscript:raise SystemExit('Supply --rscript or place Rscript on PATH.')
    commands=[[args.rscript,'scripts/run_bulk.R','.']]
   else:
    if not args.rscript or not args.external_gmt or not args.external_gmt.is_file():raise SystemExit('gsea requires --rscript and --external-gmt for the locked library.')
    # This external input is used for this local run; no redistribution permission is implied.
    shutil.copy2(args.external_gmt,work/'data/GO_KEGG_version_locked.gmt')
    commands=[[args.rscript,'scripts/run_gsea.R','.']]
  for i,command in enumerate(commands,1):
   print('RUN',name,Path(command[1]).name,flush=True);t=time.time()
   p=subprocess.run(command,cwd=work,env=env,capture_output=True,text=True,encoding='utf8',errors='replace')
   log=work/f'logs/step_{i:02}.log';log.parent.mkdir(exist_ok=True);log.write_text(p.stdout+'\n'+p.stderr,encoding='utf8')
   portable_command=[Path(command[0]).name]+[str(Path(x).relative_to(out)) if isinstance(x,str) and x.startswith(str(out)) else str(Path(x).relative_to(ROOT)) if isinstance(x,str) and x.startswith(str(ROOT)) else x for x in command[1:]]
   record=dict(module=name,step=Path(command[1]).name,command=portable_command,exit_code=p.returncode,seconds=round(time.time()-t,3),log=log.relative_to(out).as_posix(),log_sha256=digest(log),status='executed_verified' if p.returncode==0 else 'executed_failed')
   records.append(record);print('DONE',name,record['step'],'exit',p.returncode,'seconds',record['seconds'],flush=True)
   if p.returncode:
    print(p.stderr[-2000:],flush=True);break
  if records[-1]['exit_code']:break
 versions={}
 for pkg in ['numpy','pandas','scipy','matplotlib','Pillow','openpyxl','pydeseq2','anndata','h5py','PyMuPDF','statsmodels','opencv-python','formulaic','formulaic-contrasts']:
  try:versions[pkg]=md.version(pkg)
  except md.PackageNotFoundError:versions[pkg]='not installed'
 outputs=[]
 for d in [out/m for m in modules]:
  if d.exists():
   for p in d.rglob('*'):
    if p.is_file() and p.suffix in ['.csv','.tsv','.pdf','.svg','.png','.json'] and any(x in p.relative_to(d).parts for x in ['outputs','tables','figures','reports','validation']):
     outputs.append(dict(path=p.relative_to(out).as_posix(),bytes=p.stat().st_size,sha256=digest(p)))
 report=dict(date='2026-10-20',python=platform.python_version(),platform=platform.platform(),versions=versions,requested_modules=modules,steps=records,elapsed_seconds=round(time.time()-start_all,2),all_commands_passed=all(r['exit_code']==0 for r in records) and set(r['module']for r in records)==set(modules),outputs=outputs,raw_alignment_reexecuted=False,full_atlas_reclustered=False,scope='Named processed-data/count-model analyses only; archived source code is not certified by this run.')
 (out/'execution_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
 print(json.dumps({k:v for k,v in report.items()if k not in ['outputs','steps']},indent=2),flush=True)
 raise SystemExit(0 if report['all_commands_passed'] else 1)

if __name__=='__main__':main()
