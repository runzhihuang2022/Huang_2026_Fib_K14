"""Reuse source-preserving normalized SAC expression, never reviewed replacement labels."""
from pathlib import Path
import h5py,json,numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
import shutil,hashlib
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent.parent
cfg=json.loads((ROOT/'source_data/figure_config.json').read_text('utf8'));subtypes=cfg['LINEAGES']['SAC']
DERIVED=BASE/'Figure2_SAC_AB_reclustered_20260903/derived/SAC_reclustered_marker_reviewed.h5ad'
LEGACY=BASE.parent/'20260728找曾老师指导/图片修改记录/20260809CSC代码重新整理'
obs=pd.read_csv(ROOT/'source_data/scRNA_metadata.tsv.gz',sep='\t',index_col=0)
ol=obs[obs.sub_labels.isin(subtypes)].copy()
def dec(a):return np.array([x.decode() if isinstance(x,bytes) else str(x)for x in a])
with h5py.File(DERIVED,'r') as f:
 ids=dec(f['obs'][f['obs'].attrs.get('_index','_index')][:]);genes=dec(f['var'][f['var'].attrs.get('_index','_index')][:])
 assert set(ids)==set(ol.index) and len(ids)==len(ol)
 label=f['obs']['sub_labels_source'];source_labels=dec(label['categories'][:])[label['codes'][:]]
 assert np.array_equal(source_labels,ol.loc[ids].sub_labels.to_numpy())
 q=f['X'];X=sparse.csr_matrix((q['data'][:],q['indices'][:],q['indptr'][:]),shape=tuple(q.attrs['shape']))
 order=pd.Index(ids).get_indexer(ol.index);X=X[order]
 # The derived X is already log1p(CP10k) after filtering genes in >=3 cells.
 assert np.min(X.data)>=0 and np.quantile(X.data,.99)<15
 candidates=list(dict.fromkeys(g for g in cfg['marker_programs']['SAC']+[g for st in subtypes for g in cfg['MARKERS'][st]] if g in genes))
 jj=[int(np.flatnonzero(genes==g)[0])for g in candidates];Xm=X[:,jj]
 avg=[];pct=[]
 for st in subtypes:
  z=Xm[ol.sub_labels.to_numpy()==st];avg.append(np.asarray(z.mean(0)).ravel());pct.append(np.asarray((z>0).mean(0)).ravel()*100)
 avg=np.array(avg);pct=np.array(pct);selected=[];owners=[]
 for i,st in enumerate(subtypes):
  options=[]
  if st==cfg['targets']['SAC']:
   for g in cfg['marker_programs']['SAC']:
    if g in candidates and g not in selected:selected.append(g);owners.append(st)
  for g in cfg['MARKERS'][st]:
   if g not in candidates or g in selected:continue
   j=candidates.index(g);score=np.log2((avg[i,j]+1e-4)/(np.delete(avg[:,j],i).mean()+1e-4))*pct[i,j];options.append((score,g))
  limit=5-sum(s==st for s in owners)
  for _,g in sorted(options,reverse=True)[:limit]:selected.append(g);owners.append(st)
 j=[candidates.index(g)for g in selected];av=avg[:,j];pp=pct[:,j];zz=np.clip((av-av.mean(0))/np.maximum(av.std(0),1e-8),-1,2.5)
 for name,array in [('percent_expressed',pp),('average_log_expression',av),('average_expression_zscore',zz)]:pd.DataFrame(array,index=subtypes,columns=selected).to_csv(ROOT/f'source_data/SAC_dotplot_{name}.tsv',sep='\t')
 pd.DataFrame({'gene':selected,'subtype':owners}).to_csv(ROOT/'source_data/SAC_marker_blocks.tsv',sep='\t',index=False)
sm=obs[obs.Time_category.isin(cfg['MODES'])][['sample_id','patient','Time_category']].drop_duplicates('sample_id')
counts=ol.groupby(['sample_id','sub_labels']).size().unstack(fill_value=0).reindex(columns=subtypes,fill_value=0)
counts=sm.set_index('sample_id').join(counts).fillna(0);counts['lineage_total']=counts[subtypes].sum(1);counts.to_csv(ROOT/'source_data/SAC_sample_counts.tsv',sep='\t')
long=counts.reset_index().melt(id_vars=['sample_id','patient','Time_category','lineage_total'],value_vars=subtypes,var_name='sub_labels',value_name='count');long['proportion']=long['count']/long.lineage_total.replace(0,np.nan)
patient=long.groupby(['patient','Time_category','sub_labels'],observed=True).agg(proportion=('proportion','mean'),mean_recovered_count=('count','mean')).reset_index();patient.to_csv(ROOT/'source_data/SAC_patient_proportions.tsv',sep='\t',index=False)
da=[]
for st in subtypes:
 b=counts.loc[counts.Time_category==cfg['MODES'][0],st].to_numpy(float);t=counts.loc[counts.Time_category.isin(cfg['MODES'][1:]),st].to_numpy(float);p=stats.ttest_ind(b,t,equal_var=False).pvalue
 da.append(dict(sub_labels=st,log2FC_count=float(np.log2((t.mean()+.5)/(b.mean()+.5))),P_sample_nominal=float(p),n_mode1=len(b),n_modes234=len(t),observational_unit='sequencing sample; patient non-independence remains',analysis='Welch on recovered counts; EXPLORATORY'))
pd.DataFrame(da).assign(BH_FDR_within_lineage=multipletests([r['P_sample_nominal']for r in da],method='fdr_bh')[1]).to_csv(ROOT/'source_data/SAC_recovered_count_comparison.tsv',sep='\t',index=False)
e=pd.read_csv(LEGACY/'02_embeddings_cache/SAC_harmony_umap.tsv.gz',sep='\t').set_index('cell_id');xy=e.loc[ol.index,['UMAP1','UMAP2']].to_numpy();knn=cKDTree(xy).query(xy,k=50,workers=4)[1];clinical=ol.Time_category.to_numpy();late=np.isin(clinical,cfg['MODES'][1:]);early=clinical==cfg['MODES'][0];proxy=np.log2(((late[knn].sum(1)+1)/(late.sum()+1))/((early[knn].sum(1)+1)/(early.sum()+1)))
np.savez_compressed(ROOT/'source_data/SAC_embedding_inputs.npz',xy=xy,labels=ol.sub_labels.to_numpy(str),cell_ids=ol.index.to_numpy(str),time=clinical.astype(str),neighbour_density_log2FC=proxy);ol.to_csv(ROOT/'source_data/SAC_cell_metadata.tsv.gz',sep='\t',compression='gzip')
src=BASE/'SAC_SG_Progenitor_enrichment_PanelE_20260903/source_data/SAC_SG_Progenitor_ORA_displayed_terms.tsv';shutil.copy2(src,ROOT/'source_data/SAC_enrichment_displayed.tsv');shutil.copy2(src.parent/'analysis_audit.json',ROOT/'QC/SAC_prior_marker_ORA_method.json')
(ROOT/'QC/SAC_expression_cache_audit.json').write_text(json.dumps(dict(path=str(DERIVED),source_cell_ids_and_sub_labels_source_exactly_match_primary=True,normalized_X='log1p(CP10k) after original >=3-cell gene filter; NOT normalized again',annotation_used='primary sub_labels; reviewed replacement annotations NOT used',plot_embedding='original lineage-local Harmony embedding matched by cell_id',sha256=hashlib.sha256(DERIVED.read_bytes()).hexdigest()),ensure_ascii=False,indent=2),'utf8')
print('SAC prepared from verified source-preserving expression cache:',len(ol),'cells',flush=True)
