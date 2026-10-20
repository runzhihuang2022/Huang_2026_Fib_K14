from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parents[1]
cfg=json.loads((O/'analysis/defined_modules.json').read_text('utf8'))
d=pd.read_csv(O/'tables/human_scrna_sub_labels_sample_expression.csv')
b=pd.read_csv(O/'tables/human_scrna_bulk_labels_sample_expression.csv')
print('SUBTYPES',list(d.sub_labels.unique()),flush=True)
# A patient contributes at most one value to each subtype; multiple specimens are averaged.
vars=['KLF4','IGFL1','IGFL2','IGFL3','IGFL4','IGFLR1','EGFR','FGF7','HBEGF','TGFA','AREG','IGF1','IGF1R']
cols=[g+s for g in vars for s in ['_mean','_pct']]+[g+'_score'for g in cfg['modules']]
recs=[]
for subtype,t in d.groupby('sub_labels'):
 q=t[t.n_obs>=20].copy()
 for pid,p in q.groupby('patient'):
  r=dict(sub_labels=subtype,patient=pid,n_samples=len(p),n_cells=p.n_obs.sum(),group='|'.join(sorted(p.group.dropna().unique())))
  r.update(p[cols].mean().to_dict());recs.append(r)
p=pd.DataFrame(recs);p.to_csv(O/'tables/human_scrna_patient_subtype_expression.csv',index=False)
summ=[]
for subtype,q in p.groupby('sub_labels'):
 r=dict(sub_labels=subtype,n_patients=len(q),n_samples=q.n_samples.sum(),n_cells=q.n_cells.sum());r.update(q[cols].mean().to_dict());summ.append(r)
s=pd.DataFrame(summ).sort_values('IGFL1_mean',ascending=False);s.to_csv(O/'tables/human_scrna_donor_balanced_sources_receivers.csv',index=False)
print(s[['sub_labels','n_patients','n_cells','IGFL1_mean','IGFL1_pct','IGFLR1_mean','IGFLR1_pct','KLF4_mean']].to_string(index=False),flush=True)
rows=[];rng=np.random.default_rng(20261017)
focus=['Fib_K14','KC_Spinous_Mig','SAC_SG_Progenitor']
pairs=[('KLF4_mean','IGFL1_mean'),('KLF4_mean','IGFLR1_mean'),('IGFL1_mean','IGFLR1_mean')]+[(g+'_mean',m+'_score')for g in ['KLF4','IGFL1','IGFLR1']for m in ['basal_epithelial','epithelial_differentiation','wound_keratin','proliferation','MAPK_response','NFKB_response','fibroblast_matrix']]
for subtype in focus:
 q=p[p.sub_labels==subtype]
 for x,y in pairs:
  z=q[[x,y,'group']].dropna()
  if len(z)<8 or z[x].nunique()<3 or z[y].nunique()<3:continue
  r,pv=spearmanr(z[x],z[y]);boot=[]
  for it in range(1000):
   ix=rng.integers(0,len(z),len(z));u=z.iloc[ix]
   if u[x].nunique()>1 and u[y].nunique()>1:boot.append(spearmanr(u[x],u[y]).statistic)
  ranked=z[[x,y]].rank().to_numpy();design=pd.get_dummies(z.group,dtype=float);design=np.c_[np.ones(len(z)),design.to_numpy()[:,1:]]
  res=ranked-design@np.linalg.lstsq(design,ranked,rcond=None)[0]
  partial=float(np.corrcoef(res.T)[0,1])if np.std(res[:,0])>0 and np.std(res[:,1])>0 else np.nan
  rows.append(dict(subtype=subtype,x=x,y=y,n_patients=len(z),spearman_rho=r,pvalue=pv,bootstrap_lo=np.quantile(boot,.025),bootstrap_hi=np.quantile(boot,.975),group_adjusted_rank_correlation=partial,note='Exploratory patient-aggregated association; group-adjusted coefficient descriptive; no causal or temporal direction'))
r=pd.DataFrame(rows);r['FDR_BH']=multipletests(r.pvalue,method='fdr_bh')[1];r.to_csv(O/'tables/human_patient_axis_module_correlations.csv',index=False)
print('CORRELATIONS',r[(r.y=='IGFL1_mean')|(r.FDR_BH<.05)].to_string(index=False),flush=True)
# Observed pair expression availability, not a binding/communication probability.
lr=[]
for receptor in ['IGFLR1','EGFR','IGF1R']:
 for _,row in s.iterrows():lr.append(dict(receptor=receptor,receiver=row.sub_labels,n_patients=row.n_patients,mean_logCP10K=row[receptor+'_mean'],mean_positive_pct=row[receptor+'_pct']))
pd.DataFrame(lr).to_csv(O/'tables/human_candidate_receptor_receiver_availability.csv',index=False)
# Unfiltered bulk counts make true zero / insufficient-count status visible.
ann=pd.read_csv(O/'tables/gene_annotation_source.tsv',sep='\t');axis=[]
for nm in ['KLF4','IGFL1']:
 c=pd.read_csv(O/f'data/{nm}_counts.csv');c=c.merge(ann,on='gene_id',validate='one_to_one')
 a=c[c.gene_name.isin(cfg['axis'])].copy();a['total_raw_count']=a.drop(columns=['gene_id','gene_name']).sum(1);a['experiment']=nm
 axis.append(a)
pd.concat(axis,ignore_index=True).to_csv(O/'tables/bulk_axis_raw_counts_including_filtered_genes.csv',index=False)
