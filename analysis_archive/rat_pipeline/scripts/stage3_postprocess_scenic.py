#!/usr/bin/env python3
import json,re
from pathlib import Path
import numpy as np,pandas as pd,loompy
from scipy.spatial.distance import jensenshannon

ROOT=Path('data/external/local_6f6697f200d0')
OUT=ROOT/'working'/'stage3_scenic'
# pySCENIC stores structured Regulons/RegulonsAUC attributes that are readable
# by loompy but rejected by strict Loom 3 validation.  Disable validation for
# this documented pySCENIC output convention.
with loompy.connect(str(OUT/'scenic_output.loom'),mode='r',validate=False) as lf:
    cells=np.asarray(lf.ca.CellID).astype(str)
    auc=np.asarray(lf.ca.RegulonsAUC)
    names=np.asarray(lf.ca.RegulonsAUC.dtype.names if lf.ca.RegulonsAUC.dtype.names else [])
    if auc.dtype.names:
        names=np.asarray(auc.dtype.names); auc=np.vstack([auc[n] for n in names]).T
    elif 'Regulons' in lf.ca:
        names=np.asarray(lf.ca.Regulons).astype(str)
if auc.shape[0]!=len(cells) and auc.shape[1]==len(cells): auc=auc.T
if len(names)!=auc.shape[1]: names=np.array([f'Regulon_{i}' for i in range(auc.shape[1])])
meta=pd.read_csv(OUT/'cell_metadata.csv',index_col=0).reindex(cells)
aucdf=pd.DataFrame(auc,index=cells,columns=names)
aucdf.to_csv(OUT/'regulon_AUC.csv.gz',compression='gzip')

# RSS = 1 - Jensen-Shannon distance between regulon activity distribution and a cell-state indicator.
groups=meta['sub_labels_reannotated'].fillna('Unresolved').astype(str)
rss=[]
for reg in aucdf.columns:
    p=np.maximum(aucdf[reg].to_numpy(float),0)+1e-12; p=p/p.sum()
    for g in sorted(groups.unique()):
        q=(groups.to_numpy()==g).astype(float)+1e-12; q=q/q.sum()
        rss.append((reg,g,1-float(jensenshannon(p,q,base=2))))
rss=pd.DataFrame(rss,columns=['regulon','cell_state','RSS'])
rss.to_csv(OUT/'regulon_specificity_scores.csv',index=False)
rss.sort_values(['cell_state','RSS'],ascending=[True,False]).groupby('cell_state').head(20).to_csv(OUT/'top20_regulons_by_state.csv',index=False)

adj=pd.read_csv(OUT/'adjacencies.tsv',sep='\t')
adj.to_csv(OUT/'TF_target_edges.csv',index=False)
klf=adj[adj['TF'].astype(str).str.replace(r'\(.*','',regex=True).str.lower()=='klf4'].sort_values('importance',ascending=False)
klf.to_csv(OUT/'Klf4_TF_target_edges.csv',index=False)

strict=meta['Fib_K14_strict_gate'].astype(str).str.lower().isin(['true','1'])
stats=[]
for reg in aucdf.columns:
    x=aucdf[reg].to_numpy(float)
    stats.append((reg,float(np.mean(x[strict])),float(np.mean(x[~strict])),float(np.mean(x[strict])-np.mean(x[~strict]))))
pd.DataFrame(stats,columns=['regulon','mean_strict_FibK14','mean_other_fibroblasts','delta']).sort_values('delta',ascending=False).to_csv(OUT/'strict_FibK14_regulon_contrast.csv',index=False)
summary={'n_cells':len(cells),'n_regulons':aucdf.shape[1],'n_edges':len(adj),'n_Klf4_edges':len(klf),
         'formal_workflow':'pySCENIC GRNBoost2 -> cisTarget motif pruning -> AUCell; RSS calculated from AUCell matrix'}
(OUT/'stage3_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(OUT/'STAGE3_DONE').write_text('PASS\n')
print(json.dumps(summary,indent=2))
