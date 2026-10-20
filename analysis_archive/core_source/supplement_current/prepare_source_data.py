"""Source-preserving preparation. No source h5ad/annotation is written."""
from pathlib import Path
import ast,json,math,hashlib,gc,csv,shutil
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent.parent
LEGACY=BASE.parent/'20260728找曾老师指导/图片修改记录/20260809CSC代码重新整理'
PRIMARY=Path('data/external/local_c75f45c8d1a1')
for folder in ['source_data','QC','outputs/panels','outputs/assembled','assets','inputs','legends']:(ROOT/folder).mkdir(parents=True,exist_ok=True)
tree=ast.parse((LEGACY/'Figure 2A-E.py').read_text('utf8'))
cfg={}
for n in tree.body:
 if isinstance(n,ast.Assign):
  for t in n.targets:
   if isinstance(t,ast.Name) and t.id in ['LINEAGES','RAW_COLORS','MARKERS'] and t.id not in cfg:cfg[t.id]=ast.literal_eval(n.value)
MODES=['0-7dpb','8-14dpb','15-28dpb','1-2mph']
cfg.update(MODES=MODES,targets={'KC':'KC_Spinous_Mig','SAC':'SAC_SG_Progenitor'},marker_programs={'KC':['SOX9','MMP3','KRT6B'],'SAC':['SCGB2A2','S100A7','KRT77']})
(ROOT/'source_data/figure_config.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),'utf8')

def decode(v):return [x.decode() if isinstance(x,bytes) else str(x) for x in v]
def column(g,k):
 q=g[k]
 if isinstance(q,h5py.Group):
  codes=q['codes'][:];cats=np.array(decode(q['categories'][:]),dtype=object);return np.array([cats[c] if c>=0 else '' for c in codes],dtype=object)
 return np.array(decode(q[:]),dtype=object)
def frame(g,keys):return pd.DataFrame({k:column(g,k) for k in keys},index=column(g,g.attrs.get('_index','_index')))

with h5py.File(PRIMARY,'r') as f:
 obs=frame(f['obs'],['sample_id','patient','Time_category','sub_labels','bulk_labels'])
 obs['Time_category']=obs['Time_category'].str.replace('dpi','dpb',regex=False)
 obs.index.name='cell_id';obs.to_csv(ROOT/'source_data/scRNA_metadata.tsv.gz',sep='\t',compression='gzip')
 genes=np.array(decode(f['var'][f['var'].attrs.get('_index','_index')][:]))
 indptr=f['X']['indptr'][:]
 all_ids=obs.index.to_numpy();all_labels=obs.sub_labels.to_numpy();n_genes=len(genes)
 # SAC uses the separately verified source-label/logCP10k cache preparation.
 for lineage in ['KC']:
  subtypes=cfg['LINEAGES'][lineage];target=cfg['targets'][lineage]
  idx=np.flatnonzero(np.isin(all_labels,subtypes));ol=obs.iloc[idx].copy()
  emb=pd.read_csv(LEGACY/f'02_embeddings_cache/{lineage}_harmony_umap.tsv.gz',sep='\t').set_index('cell_id')
  assert emb.index.is_unique and set(ol.index).issubset(emb.index)
  xy=emb.loc[ol.index,['UMAP1','UMAP2']].to_numpy()
  print(lineage,'reading',len(idx),'cells',flush=True)
  candidates=list(dict.fromkeys(g for st in subtypes for g in cfg['MARKERS'][st] if g in genes))
  gi=np.array([np.flatnonzero(genes==g)[0] for g in candidates])
  chunks=[]
  for start in range(0,len(obs),5000):
   stop=min(len(obs),start+5000);mask=np.isin(all_labels[start:stop],subtypes)
   if not mask.any():continue
   p0,p1=int(indptr[start]),int(indptr[stop])
   x=sparse.csr_matrix((f['X']['data'][p0:p1],f['X']['indices'][p0:p1],indptr[start:stop+1]-p0),shape=(stop-start,n_genes))
   chunks.append(x[mask])
  X=sparse.vstack(chunks,format='csr');del chunks;gc.collect()
  assert X.shape[0]==len(ol)
  totals=np.asarray(X.sum(1)).ravel();pct=[];avg=[]
  Xm=X[:,gi].astype('float64');Xm=Xm.multiply((1e4/np.maximum(totals,1))[:,None]).tocsr();Xm.data=np.log1p(Xm.data)
  for st in subtypes:
   z=Xm[ol.sub_labels.to_numpy()==st]
   pct.append(np.asarray((z>0).mean(0)).ravel()*100);avg.append(np.asarray(z.mean(0)).ravel())
  avg=np.array(avg);pct=np.array(pct);owners=[];selected=[]
  for i,st in enumerate(subtypes):
   options=[]
   for g in cfg['MARKERS'][st]:
    if g not in candidates or g in selected:continue
    j=candidates.index(g);other=np.delete(avg[:,j],i).mean();score=np.log2((avg[i,j]+1e-4)/(other+1e-4))*pct[i,j]
    options.append((score,g))
   for score,g in sorted(options,reverse=True)[:5]:selected.append(g);owners.append(st)
  jj=[candidates.index(g) for g in selected];av=avg[:,jj];pp=pct[:,jj];zz=np.clip((av-av.mean(0))/np.maximum(av.std(0),1e-8),-1,2.5)
  for name,array in [('percent_expressed',pp),('average_log_expression',av),('average_expression_zscore',zz)]:
   pd.DataFrame(array,index=subtypes,columns=selected).to_csv(ROOT/f'source_data/{lineage}_dotplot_{name}.tsv',sep='\t')
  pd.DataFrame({'gene':selected,'subtype':owners}).to_csv(ROOT/f'source_data/{lineage}_marker_blocks.tsv',sep='\t',index=False)
  sm=obs[obs.Time_category.isin(MODES)][['sample_id','patient','Time_category']].drop_duplicates('sample_id')
  counts=ol.groupby(['sample_id','sub_labels']).size().unstack(fill_value=0).reindex(columns=subtypes,fill_value=0)
  counts=sm.set_index('sample_id').join(counts).fillna(0);counts['lineage_total']=counts[subtypes].sum(1)
  counts.to_csv(ROOT/f'source_data/{lineage}_sample_counts.tsv',sep='\t')
  long=counts.reset_index().melt(id_vars=['sample_id','patient','Time_category','lineage_total'],value_vars=subtypes,var_name='sub_labels',value_name='count')
  long['proportion']=long['count']/long.lineage_total.replace(0,np.nan)
  patient=long.groupby(['patient','Time_category','sub_labels'],observed=True).agg(proportion=('proportion','mean'),mean_recovered_count=('count','mean')).reset_index()
  patient.to_csv(ROOT/f'source_data/{lineage}_patient_proportions.tsv',sep='\t',index=False)
  da=[]
  for st in subtypes:
   b=counts.loc[counts.Time_category==MODES[0],st].to_numpy(float);t=counts.loc[counts.Time_category.isin(MODES[1:]),st].to_numpy(float)
   p=stats.ttest_ind(b,t,equal_var=False).pvalue
   da.append(dict(sub_labels=st,log2FC_count=float(np.log2((t.mean()+.5)/(b.mean()+.5))),P_sample_nominal=float(p),n_mode1=len(b),n_modes234=len(t),observational_unit='sequencing sample; patient non-independence remains',analysis='Welch on recovered counts; EXPLORATORY'))
  pd.DataFrame(da).assign(BH_FDR_within_lineage=multipletests([r['P_sample_nominal'] for r in da],method='fdr_bh')[1]).to_csv(ROOT/f'source_data/{lineage}_recovered_count_comparison.tsv',sep='\t',index=False)
  # Original plotting recipe used a 2D 50-neighbour density proxy, NOT Milo.
  clinical=ol.Time_category.to_numpy();knn=cKDTree(xy).query(xy,k=min(50,len(xy)),workers=4)[1]
  late=np.isin(clinical,MODES[1:]);early=clinical==MODES[0]
  proxy=np.log2(((late[knn].sum(1)+1)/(late.sum()+1))/((early[knn].sum(1)+1)/(early.sum()+1)))
  np.savez_compressed(ROOT/f'source_data/{lineage}_embedding_inputs.npz',xy=xy,labels=ol.sub_labels.to_numpy(str),cell_ids=ol.index.to_numpy(str),time=clinical.astype(str),neighbour_density_log2FC=proxy)
  ol.to_csv(ROOT/f'source_data/{lineage}_cell_metadata.tsv.gz',sep='\t',compression='gzip')
  if lineage=='KC':
   # Cell-level exploratory marker screen, not donor-level DE or efficacy test.
   X=X.astype('float64');X=X.multiply((1e4/np.maximum(totals,1))[:,None]).tocsr();X.data=np.log1p(X.data)
   mask=ol.sub_labels.to_numpy()==target;a=X[mask];b=X[~mask];n1,n0=a.shape[0],b.shape[0]
   ma=np.asarray(a.mean(0)).ravel();mb=np.asarray(b.mean(0)).ravel()
   va=(np.asarray(a.power(2).mean(0)).ravel()-ma*ma)*n1/(n1-1);vb=(np.asarray(b.power(2).mean(0)).ravel()-mb*mb)*n0/(n0-1)
   tstat,p=stats.ttest_ind_from_stats(ma,np.sqrt(np.maximum(va,0)),n1,mb,np.sqrt(np.maximum(vb,0)),n0,equal_var=False)
   p=np.nan_to_num(p,nan=1.);q=multipletests(p,method='fdr_bh')[1]
   pa=np.asarray((a>0).mean(0)).ravel();pb=np.asarray((b>0).mean(0)).ravel();lfc=np.log2((np.expm1(ma)+1e-8)/(np.expm1(mb)+1e-8))
   markers=pd.DataFrame(dict(gene=genes,mean_log_target=ma,mean_log_other=mb,pct_target=pa,pct_other=pb,log2FC=lfc,P_cell_level=p,BH_FDR=q))
   markers['selected']=(q<.05)&(lfc>=.5)&(pa>=.1)&((pa-pb)>=.1)
   markers.to_csv(ROOT/'source_data/KC_Spinous_Mig_vs_other_KC_markers.tsv.gz',sep='\t',index=False,compression='gzip')
   selectedgenes=set(markers[markers.selected].sort_values(['BH_FDR','log2FC'],ascending=[True,False]).head(500).gene)
   background=set(genes[np.asarray((X>0).mean(0)).ravel()>=.01]);selectedgenes&=background
   gmt=Path('data/external/local_58cc93d4922c')
   result=[]
   for line in gmt.read_text('utf8').splitlines():
    words=line.split('\t');term=words[0];gs=set(words[2:])&background;overlap=gs&selectedgenes
    if not gs:continue
    pval=stats.hypergeom.sf(len(overlap)-1,len(background),len(gs),len(selectedgenes))
    result.append(dict(Term=term,P_value=float(pval),overlap_count=len(overlap),gene_set_background=len(gs),selected_genes=len(selectedgenes),background_genes=len(background),Genes=';'.join(sorted(overlap))))
   ora=pd.DataFrame(result);ora['Adjusted P-value']=multipletests(ora.P_value,method='fdr_bh')[1]
   ora.to_csv(ROOT/'source_data/KC_Spinous_Mig_GO_BP_all.tsv',sep='\t',index=False)
   chosen=[];sets=[]
   for _,r in ora[(ora['Adjusted P-value']<.05)&(ora.overlap_count>=3)].sort_values(['Adjusted P-value','overlap_count'],ascending=[True,False]).iterrows():
    gs=set(r.Genes.split(';'))
    if any(len(gs&s)/len(gs|s)>=.65 for s in sets):continue
    chosen.append(r);sets.append(gs)
    if len(chosen)==6:break
   pd.DataFrame(chosen).to_csv(ROOT/'source_data/KC_enrichment_displayed.tsv',sep='\t',index=False)
   print('KC ORA',len(selectedgenes),'markers',len(chosen),'terms',flush=True)
   (ROOT/'QC/KC_marker_ORA_method.json').write_text(json.dumps(dict(source=str(PRIMARY),contrast='KC_Spinous_Mig vs other KC',method='cell-level Welch on log1p(CP10k), BH across all genes; exploratory marker screen, NOT patient-level DE',filter='BH<0.05, log2FC>=0.5, pct_target>=0.1, pct difference>=0.1; at most 500; intersect detectable background',ORA='local hypergeometric; GO BP 2025; BH across library',gmt_source=str(gmt),gmt_sha256=hashlib.sha256(gmt.read_bytes()).hexdigest()),ensure_ascii=False,indent=2),'utf8')
  else:
   src=BASE/'SAC_SG_Progenitor_enrichment_PanelE_20260903/source_data/SAC_SG_Progenitor_ORA_displayed_terms.tsv'
   shutil.copy2(src,ROOT/'source_data/SAC_enrichment_displayed.tsv')
   shutil.copy2(src.parent/'analysis_audit.json',ROOT/'QC/SAC_prior_marker_ORA_method.json')
  print(lineage,'prepared',len(idx),'cells',flush=True)
  del X,Xm;gc.collect()
print('Source preparation complete',flush=True)
