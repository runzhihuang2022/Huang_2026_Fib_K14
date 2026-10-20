"""Re-render retained observations and audit existing trajectories; no new integration."""
from pathlib import Path
import json,sys,itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1];D=R/'data/processed';F=R/'outputs/figures/atlas_audits';T=R/'outputs/tables/atlas_audits'
F.mkdir(parents=True,exist_ok=True);T.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(F/(name+'.'+ext),dpi=200,bbox_inches='tight')
 plt.close(fig)
d=pd.read_csv(D/'omics/human_scrna_observations.csv.gz');assert len(d)==279305 and d.sample_code.nunique()==38
colors=pd.read_csv(D/'human_tree/bulk_labels_color_key.tsv',sep='\t').set_index('bulk_labels').color.to_dict()
fig,ax=plt.subplots(figsize=(7.4,4.8))
for label,g in d.groupby('bulk_labels',sort=True):ax.scatter(g.coord_x,g.coord_y,s=.15,c=colors.get(label,'gray'),label=label,rasterized=True,linewidths=0)
ax.set(xlabel='Retained UMAP 1',ylabel='Retained UMAP 2',title='Human atlas: 279,305 cells; 38 assay samples');ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',markerscale=10,frameon=False,fontsize=7)
save(fig,'human_broad_lineage_umap')
fig,ax=plt.subplots(figsize=(6.7,4.8));fib=d[d.bulk_labels=='Fibroblast'];palette=plt.get_cmap('tab10')
for i,(label,g) in enumerate(fib.groupby('sub_labels',sort=True)):ax.scatter(g.coord_x,g.coord_y,s=.6,color=palette(i),label=label,rasterized=True,linewidths=0)
ax.set(xlabel='Retained whole-atlas UMAP 1',ylabel='Retained whole-atlas UMAP 2',title='Fibroblast states in the retained atlas embedding');ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',markerscale=6,frameon=False)
save(fig,'fibroblast_states_whole_atlas_embedding')
summ=[]
for modality,counts,genes,n in [('human_scrna','nCount_RNA','nFeature_RNA',38),('human_stereo','total_counts','n_genes_by_counts',11),('human_visium','total_counts','n_genes_by_counts',2)]:
 x=pd.read_csv(D/'omics'/f'{modality}_observations.csv.gz');assert x.sample_code.nunique()==n
 fig,axes=plt.subplots(2,1,figsize=(max(6,n*.20),5.3),sharex=True)
 samples=sorted(x.sample_code.unique())
 for ax,col,title in zip(axes,[counts,genes],['Transcripts/counts','Detected genes']):
  vals=[np.log10(1+x.loc[x.sample_code==s,col].to_numpy()) for s in samples]
  ax.boxplot(vals,showfliers=False,patch_artist=True,boxprops={'facecolor':'#76a5af'},medianprops={'color':'black'},widths=.65)
  ax.set(ylabel='log10(1 + count)',title=modality+': '+title)
 for s,g in x.groupby('sample_code'):
  summ.append(dict(modality=modality,sample_code=s,n_observations=len(g),median_counts=float(g[counts].median()),median_genes=float(g[genes].median()),observational_unit='cell' if modality=='human_scrna' else 'bin' if modality=='human_stereo' else 'spot'))
 axes[-1].set_xticks(np.arange(1,len(samples)+1),[s.rsplit('_',1)[1] for s in samples],rotation=90);axes[-1].set_xlabel('Pseudonymous assay sample');fig.tight_layout();save(fig,modality+'_retained_qc')
pd.DataFrame(summ).to_csv(T/'retained_sample_qc.csv',index=False)
# Recovered SAC marker table: deliberately not called an exact S3B reproduction.
a=pd.read_csv(D/'SAC/SAC_dotplot_average_log1p_cpm.tsv',sep='\t',index_col=0);p=pd.read_csv(D/'SAC/SAC_dotplot_percent_expressed.tsv',sep='\t',index_col=0).reindex(index=a.index,columns=a.columns)
fig,ax=plt.subplots(figsize=(11,4.5));xx,yy=np.meshgrid(np.arange(a.shape[1]),np.arange(a.shape[0]));z=ax.scatter(xx.ravel(),yy.ravel(),s=p.to_numpy().ravel()*.6,c=a.to_numpy().ravel(),cmap='Reds',vmin=0)
ax.set_xticks(np.arange(a.shape[1]),[c.split('__')[0] for c in a.columns],rotation=90);ax.set_yticks(np.arange(a.shape[0]),a.index);ax.invert_yaxis();ax.set_title('Recovered September 17 SAC source table (marker set differs from current S3B)');fig.colorbar(z,ax=ax,label='Mean log1p CPM');save(fig,'SAC_source_version_marker_plot')
# No cell/patient/sample identity is needed for this sensitivity audit.
medians={};trajcounts={}
names=['M2_Fib_KC_SAC','M2_Fib_KC','M2_KC_SAC','M3_KC_SAC']
for i,name in enumerate(names,1):
 path=D/'pseudotime'/('analysis_4.tsv.gz' if i==4 else f'analysis_{i}.csv')
 x=pd.read_csv(path,sep='\t' if i==4 else ',');pc='pseudotime' if i==4 else 'Pseudotime';x=x[np.isfinite(x[pc])]
 medians[name]=x.groupby('sub_labels')[pc].median();trajcounts[name]=len(x)
 rows=[]
for a,b in itertools.combinations(names,2):
 common=medians[a].index.intersection(medians[b].index);rho,pv=spearmanr(medians[a][common],medians[b][common]);rows.append(dict(analysis_a=a,analysis_b=b,n_shared_states=len(common),spearman_rho=rho,p_value=pv,interpretation='state_rank_sensitivity_not_patient_inference'))
pd.DataFrame(rows).to_csv(T/'pseudotime_state_rank_sensitivity.csv',index=False);pd.DataFrame(medians).to_csv(T/'pseudotime_state_medians.csv')
heat=np.eye(4)
for q in rows:ii=names.index(q['analysis_a']);jj=names.index(q['analysis_b']);heat[ii,jj]=heat[jj,ii]=q['spearman_rho']
fig,ax=plt.subplots(figsize=(5,4));im=ax.imshow(heat,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(4),names,rotation=40,ha='right');ax.set_yticks(range(4),names)
for i in range(4):
 for j in range(4):ax.text(j,i,f'{heat[i,j]:.2f}',ha='center',va='center',color='black' if abs(heat[i,j])<.7 else 'white')
fig.colorbar(im,ax=ax,label='Spearman rho of state-median pseudotime');ax.set_title('Trajectory ranking depends on input composition');save(fig,'trajectory_rank_sensitivity')
e=pd.read_csv(D/'ATAC_network/FibK14_KLF4_network_evidence_audit.csv').set_index('metric').value
fig,ax=plt.subplots(figsize=(6,3));ax.barh(['Positive sample-level FDR < 0.05','Negative cell-level FDR < 0.05','Positive cell-level FDR < 0.05'],[float(e['Positive sample-level FDR<0.05']),float(e['Negative cell-level FDR<0.05']),float(e['Positive cell-level FDR<0.05'])],color=['#888888','#4d90b6','#d9785b']);ax.set(xlabel='Number of candidate associations',title='KLF4 network: cell-level and sample-level evidence');save(fig,'KLF4_network_evidence_levels')
report=dict(human_cells=len(d),human_samples=d.sample_code.nunique(),human_states=d.sub_labels.nunique(),qc_samples=len(summ),trajectory_finite_cells=trajcounts,trajectory_rows=rows,scope='Re-rendering processed observations; no raw-read alignment, integration, clustering, trajectory fitting, or new experimental data')
(R/'validation/atlas_audits').mkdir(parents=True,exist_ok=True);(R/'validation/atlas_audits/report.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2))
