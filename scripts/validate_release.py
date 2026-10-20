"""Check the packaged release without opening author workstation data."""
from pathlib import Path
import ast,csv,gzip,hashlib,json,re,sys
R=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for c in iter(lambda:f.read(1048576),b''):h.update(c)
 return h.hexdigest()
issues=[];checked=0;large=[];syntax=[]
manifest=R/'provenance/file_checksums.csv'
if manifest.exists():
 for row in csv.DictReader(manifest.open(encoding='utf-8-sig')):
  p=R/row['path'];checked+=1
  if not p.is_file() or sha(p)!=row['sha256']:issues.append({'type':'checksum_mismatch','file':row['path']})
else:issues.append({'type':'missing_checksum_manifest'})
path_re=re.compile(r'(?<![A-Za-z])\b[A-Za-z]:[\\/]|/(?:home|media|mnt)/')
credential_re=re.compile(r'''(?i)(?:api_key|password|passwd|access_token)\s*[:=]\s*["'][^"'\s]{8,}["']|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----|gh[pousr]_[A-Za-z0-9]{30,}''')
for p in R.rglob('*'):
 if not p.is_file() or p.relative_to(R).parts[0]=='outputs' or any(x in p.relative_to(R).parts for x in ['.git','__pycache__']):continue
 rel=p.relative_to(R).as_posix()
 if p.stat().st_size>=25*1024**2:large.append({'file':rel,'bytes':p.stat().st_size})
 if p.suffix in ['.py','.R','.r','.md','.txt','.json','.csv','.tsv','.ipynb','.yaml','.yml','.cff'] or p.name.endswith(('.csv.gz','.tsv.gz')):
  try:t=gzip.open(p,'rt',encoding='utf-8-sig').read() if p.suffix=='.gz' else p.read_text('utf-8-sig')
  except (UnicodeError,OSError):continue
  if rel!=Path(__file__).relative_to(R).as_posix():
   if path_re.search(t):issues.append({'type':'workstation_path','file':rel})
   if credential_re.search(t):issues.append({'type':'possible_credential','file':rel})
  if p.suffix=='.ipynb':
   try:json.loads(t)
   except Exception:issues.append({'type':'invalid_notebook_json','file':rel})
  if p.suffix=='.py' and rel.startswith(('modules/','scripts/','script/','Figure','run_')):
   try:ast.parse(t);syntax.append(rel)
   except SyntaxError as e:issues.append({'type':'active_script_syntax','file':rel,'line':e.lineno})
issues.extend({'type':'browser_upload_limit',**x}for x in large)
report={'status':'passed' if not issues else 'failed','checksummed_files':checked,'active_python_scripts_checked':len(syntax),'files_over_25_MiB':large,'issues':issues,'scope':'Automated packaging checks only; not a certification of consent, licensing, scientific validity or complete raw-data reproduction.'}
(R/'validation').mkdir(exist_ok=True)
(R/'validation/release_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2));sys.exit(0 if not issues else 1)
