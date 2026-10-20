"""Recompute tabular checks and redraw current KC/SAC evidence from portable processed data."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];D=R/'data/processed/supplement_current';T=R/'outputs/tables/supplement_current';F=R/'outputs/figures/supplement_current'
T.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
modes=['0-7dpb','8-14dpb','15-28dpb','1-2mph'];cfg=json.loads((D/'figure_config.json').read_text('utf8'));reports=[]
def save(fig,name):
 fig.tight_layout()
 for ext in ['png','pdf','svg']:fig.savefig(F/f'{name}.{ext}',dpi=600,bbox_inches='tight')
 plt.close(fig)
for lin in ['KC','SAC']:
 c=pd.read_csv(D/f'{lin}_sample_counts.tsv',sep='\t');cols=list(cfg['LINEAGES'][lin]);target=cfg['targets'][lin];assert (c[cols].sum(axis=1)==c.lineage_total).all()
 r=[]
 for st in cols:
  b=c.loc[c.Time_category==modes[0],st].to_numpy(float);a=c.loc[c.Time_category.isin(modes[1:]),st].to_numpy(float)
  r.append(dict(sub_labels=st,log2FC_count=float(np.log2((a.mean()+.5)/(b.mean()+.5))),P_sample_nominal=float(stats.ttest_ind(a,b,equal_var=False).pvalue)))
 rr=pd.DataFrame(r);p=rr.P_sample_nominal.to_numpy();order=np.argsort(p);qs=np.minimum.accumulate((p[order]*len(p)/(np.arange(len(p))+1))[::-1])[::-1];adj=np.empty_like(p);adj[order]=np.minimum(qs,1);rr['BH_FDR_within_lineage']=adj
 old=pd.read_csv(D/f'{lin}_recovered_count_comparison.tsv',sep='\t');j=rr.merge(old,on='sub_labels',suffixes=('_new','_source'))
 for v in ['log2FC_count','P_sample_nominal','BH_FDR_within_lineage']:assert np.allclose(j[v+'_new'],j[v+'_source'],atol=1e-10)
 rr.to_csv(T/f'{lin}_recovered_count_comparison.csv',index=False)
 summary=c.groupby('Time_category').agg(n_libraries=('sample_id','nunique'),n_patients=('patient','nunique'),recovered_lineage_cells=('lineage_total','sum'),recovered_target_cells=(target,'sum')).reindex(modes)
 summary.to_csv(T/f'{lin}_phase_units.csv')
 long=c.melt(id_vars=['sample_id','patient','Time_category','lineage_total'],value_vars=cols,var_name='state',value_name='count');long['fraction']=long['count']/long.lineage_total.replace(0,np.nan)
 pat=long.groupby(['patient','Time_category','state']).fraction.mean().reset_index();comp=pat.groupby(['Time_category','state']).fraction.mean().unstack().reindex(modes)[cols];comp.to_csv(T/f'{lin}_patient_averaged_composition.csv')
 fig,axs=plt.subplots(1,2,figsize=(7.5,3.2));bottom=np.zeros(4)
 for i,st in enumerate(cols):
  color=plt.get_cmap('tab10')(i%10);y=comp[st].fillna(0).to_numpy()*100;axs[0].bar(range(4),y,bottom=bottom,color=color,label=st,width=.7);bottom+=y
 axs[0].set(xticks=range(4),xticklabels=['0–7 d','8–14 d','15–28 d','1–2 mo'],ylabel=f'Within-{lin} proportion (%)',ylim=(0,100));axs[0].legend(fontsize=6,loc='upper left',bbox_to_anchor=(1,1),frameon=False)
 for i,m in enumerate(modes):
  vals=c.loc[c.Time_category==m,target];axs[1].scatter(np.full(len(vals),i),vals,s=15,facecolors='none',edgecolors='#25596e');axs[1].errorbar(i,vals.mean(),yerr=vals.sem(),fmt='o',c='#dd594e',capsize=3)
 axs[1].set(xticks=range(4),xticklabels=['0–7 d','8–14 d','15–28 d','1–2 mo'],ylabel=f'{target}: recovered cells / library');fig.subplots_adjust(wspace=.8);save(fig,f'{lin}_composition_and_counts')
 z=np.load(D/f'{lin}_embedding_inputs.npz',allow_pickle=False);xy=z['xy'];labels=z['labels'];time=z['time'];fig,axs=plt.subplots(2,2,figsize=(6,5))
 for ax,m in zip(axs.flat,modes):
  mask=time==m
  for i,st in enumerate(cols):
   v=mask&(labels==st);ax.scatter(xy[v,0],xy[v,1],s=.1,color=plt.get_cmap('tab10')(i%10),rasterized=True)
  ax.set_title(f'{m}: {mask.sum():,} cells; {summary.loc[m,"n_libraries"]} libraries');ax.set(xticks=[],yticks=[])
 from matplotlib.lines import Line2D
 handles=[Line2D([],[],marker='o',linestyle='None',color=plt.get_cmap('tab10')(i%10),markersize=3,label=st) for i,st in enumerate(cols)]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.01),ncol=3,fontsize=6,frameon=False)
 save(fig,f'{lin}_phase_embeddings')
 reports.append(dict(lineage=lin,contrast_table_exact=True,phase_counts=summary.reset_index().to_dict('records'),scope='retained processed embeddings and source labels; no new clustering'))
go=pd.read_csv(D/'KC_selected_enriched_functions.tsv',sep='\t');fig,ax=plt.subplots(figsize=(4.2,2.8));lab=go.display_label.str.replace('\\n','\n',regex=False);ax.barh(range(len(go)),-np.log10(go['Adjusted P-value']),color=['#bdd51c','#39b986','#249b8c','#44549b']);ax.set(yticks=range(len(go)),yticklabels=lab,xlabel='−log10 (BH-adjusted P)');ax.invert_yaxis();save(fig,'KC_selected_GO')
go.to_csv(T/'KC_selected_GO.csv',index=False)
(T/'audit.json').write_text(json.dumps({'reports':reports,'go_selection':'post hoc selected view; library-wide BH unchanged','not_reproduced':['current S2 microscopy acquisition','current S3 older marker plot exact version','raw matrices, integration or clustering'],'S2A_artwork_N_conflict':'actual libraries 8,5,9,6 and patients 4,5,5,2; original artwork labels 9,9,14,7'},indent=2),'utf8');print(json.dumps(reports))
