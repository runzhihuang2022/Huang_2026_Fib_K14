"""Execute the available supporting analyses for one manuscript figure."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--figure',required=True);ap.add_argument('--output-dir',type=Path);a=ap.parse_args()
cfg=ROOT/a.figure/'analysis_config.json'
if not cfg.is_file():raise SystemExit('No validated current-composite task is registered for '+a.figure+'. Read its panel manifest and source archive.')
spec=json.loads(cfg.read_text('utf8'));out=(a.output_dir or ROOT/'outputs'/a.figure).resolve()
if out.exists():raise SystemExit('Choose a new --output-dir to preserve previous results.')
out.mkdir(parents=True);env=dict(os.environ,PYTHONUTF8='1',PYTHONIOENCODING='utf-8',MPLBACKEND='Agg');records=[]
if spec['core_scripts']:
 work=out/'core';shutil.copytree(ROOT/'modules/core',work);(work/'validation').mkdir(exist_ok=True)
 for step in spec['core_scripts']:
  t=time.time();p=subprocess.run([sys.executable,str(work/'scripts'/step)],cwd=work,env=env,capture_output=True,text=True,encoding='utf8',errors='replace')
  (out/(Path(step).stem+'.log')).write_text(p.stdout+'\n'+p.stderr,'utf8');records.append(dict(script=step,exit_code=p.returncode,seconds=round(time.time()-t,2)));print(step,p.returncode,flush=True)
  if p.returncode:raise SystemExit(p.returncode)
if spec['extra_groups']:
 p=subprocess.run([sys.executable,str(ROOT/'run_reproduce.py'),'--modules',','.join(spec['extra_groups']),'--output-dir',str(out/'additional')],cwd=ROOT,env=env)
 if p.returncode:raise SystemExit(p.returncode)
(out/'figure_run.json').write_text(json.dumps(dict(figure=a.figure,tasks=records,scope=spec['scope']),indent=2),'utf8')
print('Finished the registered supporting analyses. Exact composite and experimental-source status remain in the figure panel manifest.')
