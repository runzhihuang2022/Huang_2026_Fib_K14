"""Validate and redraw retained CellOracle, orthology and ATAC summary evidence."""
from pathlib import Path
import json,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=True)
T=OUT/'tables';F=OUT/'figures';T.mkdir(exist_ok=True);F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
base=ROOT/'analysis_archive/human_celloracle/source_data';checks=[]
allrows=[]
for condition in ['KO','OE']:
 path=next(base.glob(f'KLF4_{condition}_cell_vectors*.gz'))
 d=pd.read_csv(path,sep='\t');s=d.groupby('sub_labels').KC_proxy_perturbation_score.agg(n='size',median_score='median',mean_score='mean')
 s['positive_fraction']=d.assign(positive=d.KC_proxy_perturbation_score>0).groupby('sub_labels').positive.mean();s['condition']=condition
 allrows.append(s.reset_index())
actual=pd.concat(allrows,ignore_index=True);actual.to_csv(T/'CellOracle_recomputed_subtype_summary.csv',index=False)
expected=pd.read_csv(base/'KLF4_KO_OE_by_subtype_summary.tsv',sep='\t')
e=expected.merge(actual,on=['condition','sub_labels'],suffixes=('_source','_rerun'),validate='one_to_one')
assert len(e)==len(expected)==len(actual)>0
for c in ['n','median_score','mean_score','positive_fraction']:
 assert np.allclose(e[c+'_source'],e[c+'_rerun'],rtol=1e-8,atol=1e-10),c
checks.append({'analysis':'CellOracle per-cell score aggregation','source_match':True,'rows':len(e),'model_refitted':False})
fig,ax=plt.subplots(figsize=(7,3.2));p=actual.pivot(index='sub_labels',columns='condition',values='mean_score');p.plot.bar(ax=ax,color=['#368caa','#e17b63']);ax.set_ylabel('Retained CellOracle perturbation score');fig.tight_layout()
for ext in ['pdf','svg','png']:fig.savefig(F/f'CellOracle_score_source_rebuild.{ext}',dpi=300)
plt.close(fig)
for name,pattern in [('orthology','*strict_orthology_audit.csv'),('chromVAR','KLF4_chromVAR_FibK14_by_sample.csv'),('cross_species','cross_species_effect_sizes.csv')]:
 paths=list((ROOT/'analysis_archive').rglob(pattern))
 if paths:
  d=pd.read_csv(paths[0]);d.to_csv(T/f'{name}_retained_evidence.csv',index=False);checks.append({'analysis':name,'retained_table_rows':len(d),'new_inference':False})
(OUT/'reports').mkdir(exist_ok=True)
(OUT/'reports/derived_validation.json').write_text(json.dumps(checks,indent=2),encoding='utf8')
print(json.dumps(checks),flush=True)
