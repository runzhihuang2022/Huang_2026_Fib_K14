from pathlib import Path
import json,pandas as pd,numpy as np,hashlib
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1];cfg=json.loads((O/'analysis/defined_modules.json').read_text('utf8'));genes=cfg['genes']
md=pd.read_csv(O/'data/rat_scrna_obs.csv',usecols=['Unnamed: 0','sample','tissue','Time','sub_labels','bulk_labels.y.y'])
md=md.rename(columns={'Unnamed: 0':'obs_id','bulk_labels.y.y':'bulk_labels'})
v=pd.read_csv(O/'data/rat_scrna_panel_logCP10K.csv',index_col=0)
assert len(md)==450282 and len(md['sample'].unique())==53 and np.array_equal(md.obs_id,v.index)
av=pd.read_csv(O/'tables/rat_scrna_gene_availability.csv').set_index('gene');gok=av[av.present==True].index
for m,gs in cfg['modules'].items():
 use=[g for g in gs if g in gok];v[m+'_score']=v[use].mean(1)if use else np.nan
rows=[]
for (sample,sub),idx in md.groupby(['sample','sub_labels']).groups.items():
 if len(idx)<20:continue
 z=v.iloc[idx];r=dict(sample=sample,sub_labels=sub,n_cells=len(idx),tissue=md.loc[idx,'tissue'].iloc[0],Time=md.loc[idx,'Time'].iloc[0])
 for g in ['KLF4','IGFL1','IGFLR1','KRT14','KRT17','TGFA','HBEGF','EGFR','IGF1','IGF1R']:
  r[g+'_mean']=z[g].mean()if g in gok else np.nan;r[g+'_pct']=100*(z[g]>0).mean()if g in gok else np.nan
 for m in cfg['modules']:r[m+'_score']=z[m+'_score'].mean()
 rows.append(r)
d=pd.DataFrame(rows);d.to_csv(O/'tables/rat_scrna_sample_subtype_expression.csv',index=False)
cor=[]
for sub in ['Fib_K14','KC_Spinous_Mig','SAC_SG_Progenitor']:
 q=d[d.sub_labels==sub]
 for gene in ['KLF4','IGFLR1']:
  for m in ['basal_epithelial','epithelial_differentiation','wound_keratin','proliferation','MAPK_response','NFKB_response','fibroblast_matrix']:
   z=q[[gene+'_mean',m+'_score']].dropna()
   if len(z)<8 or z.iloc[:,0].nunique()<3:continue
   r,p=spearmanr(z.iloc[:,0],z.iloc[:,1]);cor.append(dict(subtype=sub,gene=gene,module=m,n_libraries=len(z),rho=r,pvalue=p,inference_status='Exploratory library-level association; animal independence undocumented; no cross-species causal test'))
c=pd.DataFrame(cor);c['FDR_BH']=multipletests(c.pvalue,method='fdr_bh')[1];c.to_csv(O/'tables/rat_library_axis_module_correlations.csv',index=False)
print('RAT availability',av.loc[['KLF4','IGFL1','IGFLR1','KRT14']].to_string(),flush=True)
print('RAT associations',c[c.FDR_BH<.05].to_string(index=False),flush=True)
# Human spatial ligand/receptor geometry. Within-section null preserves tissue footprint,
# but not cell-type composition; ratio is descriptive and is not a clinical P value.
st=pd.read_csv(O/'data/human_stereo_obs.csv.gz');z=np.load(O/'data/human_stereo_panel.npz');p=z['positive'];val=z['values'];gi={g:i for i,g in enumerate(genes)};rng=np.random.default_rng(20261017)
sp=[]
for sec,idx in st.groupby('analysis_sample').groups.items():
 idx=np.array(idx);xy=st.loc[idx,['x','y']].to_numpy(float)
 for ligand,receptor in [('IGFL1','IGFLR1'),('TGFA','EGFR'),('HBEGF','EGFR')]:
  a=p[idx,gi[ligand]];b=p[idx,gi[receptor]];r=dict(section=sec,ligand=ligand,receptor=receptor,n_bins=len(idx),ligand_positive=int(a.sum()),receptor_positive=int(b.sum()),co_positive=int((a&b).sum()),analysis_status='descriptive_within_section_no_biological_replication')
  if a.sum()>=3 and b.sum()>=3:
   # Exclude co-positive bins from the ligand set to avoid self-distance zero.
   send=a&~b
   if send.sum()>=3:
    dst=np.median(cKDTree(xy[b]).query(xy[send])[0]);null=[]
    for it in range(200):
     rp=rng.choice(len(idx),int(b.sum()),replace=False);null.append(np.median(cKDTree(xy[rp]).query(xy[send])[0]))
    r.update(observed_median_distance_native_units=float(dst),null_median_distance=float(np.median(null)),distance_ratio=float(dst/np.median(null))if np.median(null)>0 else np.nan,n_permutations=200)
  sp.append(r)
pd.DataFrame(sp).to_csv(O/'tables/human_stereo_LR_spatial_distance_DESCRIPTIVE.csv',index=False)
print('SPATIAL DISTANCE DONE',flush=True)
