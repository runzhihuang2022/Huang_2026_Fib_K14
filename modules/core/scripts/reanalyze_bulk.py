"""Count-model sensitivity reanalysis. Does not claim to reproduce the legacy z-test P values."""
from pathlib import Path
import json,sys,importlib.metadata as md,time
import numpy as np,pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];D=R/'data/processed/bulk';O=R/'outputs/tables/bulk';F=R/'outputs/figures/bulk'
def main():
 O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True);np.random.seed(20261016)
 plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
 reports=[]
 names=pd.read_csv(D/'IGFL1_counts.csv',usecols=['gene_id','gene_name']).set_index('gene_id').gene_name
 for exp in ['KLF4','IGFL1']:
  start=time.time();raw=pd.read_csv(D/f'{exp}_counts.csv');samples=[c for c in raw.columns if c not in ['gene_id','gene_name']]
  x=raw.set_index('gene_id')[samples].T
  assert x.index.is_unique and x.columns.is_unique and np.isfinite(x.values).all() and (x.values>=0).all() and np.allclose(x.values,np.round(x.values))
  x=x.astype('int64');keep=(x>=10).sum(axis=0)>=2;x=x.loc[:,keep]
  condition=['Control' if s.startswith(('C','PBS')) else 'Treated' for s in samples]
  meta=pd.DataFrame({'condition':condition},index=samples)
  print('Fit',exp,x.shape,flush=True)
  dds=DeseqDataSet(counts=x,metadata=meta,design='~condition',refit_cooks=True,n_cpus=1,quiet=True)
  dds.deseq2();st=DeseqStats(dds,contrast=['condition','Treated','Control'],alpha=.05,cooks_filter=True,independent_filter=True,n_cpus=1,quiet=True);st.summary()
  res=st.results_df.copy();res.insert(0,'gene_name',names.reindex(res.index));res.to_csv(O/f'{exp}_pydeseq2_all_genes.csv')
  norm=pd.DataFrame(dds.layers['normed_counts'],index=x.index,columns=x.columns).T;norm.to_csv(O/f'{exp}_normalized_counts.csv')
  sig=(res.padj<.05)&(res.log2FoldChange.abs()>1)
  summary={'experiment':exp,'n_libraries':len(samples),'group_sizes':meta.condition.value_counts().to_dict(),'input_genes':len(raw),'filtered_genes':x.shape[1],'filter':'at least 10 counts in at least 2 libraries','model':'PyDESeq2 negative-binomial ~condition; Treated/Control','alpha':.05,'effect_threshold_abs_log2FC':1,'up':int((sig&(res.log2FoldChange>0)).sum()),'down':int((sig&(res.log2FoldChange<0)).sum()),'runtime_seconds':time.time()-start,'versions':{p:md.version(p) for p in ['pydeseq2','anndata','numpy','pandas','scipy']},'replicate_identity':'library-level; independent donor/culture and batch metadata not supplied','purpose':'count-model sensitivity; replaces neither experiment identity nor source-figure version audit'}
  reports.append(summary)
  fig,ax=plt.subplots(figsize=(4.8,3.5));y=-np.log10(res.padj.fillna(1).clip(lower=1e-300));colors=np.where(sig,np.where(res.log2FoldChange>0,'#c54444','#397eb5'),'#cccccc')
  ax.scatter(res.log2FoldChange,y,c=colors,s=4,alpha=.55,rasterized=True);ax.axhline(-np.log10(.05),color='grey',ls='--',lw=.5)
  ax.set(xlabel=f'log2 fold change ({exp} / control)',ylabel='−log10 BH-adjusted P',title=f'{exp}: count-model reanalysis')
  label_positions={'KLF4':(3.5,38),'KRT14':(8,18),'COL1A1':(3,10),'COL3A1':(-6,45)} if exp=='KLF4' else {'KLF4':(4,58),'TACSTD2':(-9,26),'COL1A1':(0,280),'COL3A1':(-8,285)}
  for gene,pos in label_positions.items():
   z=res[res.gene_name==gene]
   if len(z):
    row=z.iloc[0];ax.annotate(gene,(row.log2FoldChange,float(-np.log10(max(row.padj if pd.notna(row.padj) else 1,1e-300)))),xytext=pos,textcoords='data',fontsize=6,arrowprops={'arrowstyle':'-','lw':.4,'color':'#555555'})
  if exp=='IGFL1':ax.text(.98,.03,'Display floor for adjusted P: 1e−300',transform=ax.transAxes,ha='right',fontsize=5)
  fig.tight_layout()
  for ext in ['png','pdf','svg']:fig.savefig(F/f'{exp}_count_model_volcano.{ext}',dpi=600)
  plt.close(fig)
  res[res.gene_name.isin(['KLF4','KRT14','KRT5','TACSTD2','COL1A1','COL3A1','VIM','PDGFRA'])].to_csv(O/f'{exp}_key_markers.csv')
  print(json.dumps(summary),flush=True)
 (O/'analysis_summary.json').write_text(json.dumps(reports,indent=2),encoding='utf8')
if __name__=='__main__':main()
