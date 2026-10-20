from pathlib import Path
import json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
O=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
files=[]
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(O/f'figures/{name}.{ext}',dpi=400,facecolor='white')
 files.append(name);plt.close(fig)
def stars(q):return '***'if q<.001 else '**'if q<.01 else '*'if q<.05 else ''
ks={n:pd.read_csv(O/f'tables/{n}_DESeq2_all.csv').set_index('gene_name')for n in ['KLF4','IGFL1']}
gs=['KLF4','KRT14','KRT17','PDGFRA','VIM','EPCAM','TACSTD2','TGFA','HBEGF','AREG','EGFR','DUSP6','NFKBIA','BIRC3','MKI67','PCNA','FGF7','STAT3']
mat=np.array([[ks[n].loc[g,'log2FoldChange']if g in ks[n].index else np.nan for n in ks]for g in gs])
fig,ax=plt.subplots(figsize=(140/25.4,145/25.4));fig.subplots_adjust(left=.25,right=.90,bottom=.13,top=.92)
im=ax.imshow(mat,cmap='RdBu_r',vmin=-5,vmax=5,aspect='auto');ax.set_yticks(range(len(gs)),gs,fontstyle='italic');ax.set_xticks([0,1],['KLF4-OE','IGFL1 treatment']);ax.set_title('Primary human dermal fibroblasts: transcriptional responses')
for j,n in enumerate(ks):
 for r,g in enumerate(gs):
  if g in ks[n].index:
   x=ks[n].loc[g];ax.text(j,r,f'{x.log2FoldChange:+.2f}'+stars(x.padj),ha='center',va='center',color='white'if abs(x.log2FoldChange)>2.8 else 'black',fontsize=8)
  else:ax.text(j,r,'<10 counts',ha='center',va='center')
fig.colorbar(im,ax=ax,fraction=.045,pad=.04,label='DESeq2 log2 fold change')
fig.text(.03,.035,'* BH FDR <0.05; ** <0.01; *** <0.001. Colors clipped at +/-5.\nKLF4: 2 vs 2; IGFL1: 2 vs 3. Replicate independence unconfirmed.',fontsize=8)
save(fig,'01_bulk_shared_response_heatmap')

d=pd.read_csv(O/'tables/human_scrna_donor_balanced_sources_receivers.csv').set_index('sub_labels')
sel=['KC_Spinous_Mig','KC_Granular','KC_Basal_Mig','KC_Basal_Prolif','Fib_K14','SAC_SG_Progenitor','Fib_CD45','CD4_T','CD8_T','Mast','cDC2','pDC']
fig,axes=plt.subplots(1,2,figsize=(190/25.4,100/25.4),sharey=True);fig.subplots_adjust(left=.27,right=.98,wspace=.12,bottom=.20,top=.88)
for ax,g,color in zip(axes,['IGFL1','IGFLR1'],['#329fbb','#ee7960']):
 ax.barh(range(len(sel)),d.loc[sel,g+'_pct'],color=color);ax.set_yticks(range(len(sel)),sel);ax.set_title(g);ax.set_xlabel('Mean expressing cells per donor (%)');ax.grid(axis='x',alpha=.15)
axes[0].invert_yaxis();fig.suptitle('Human wound atlas: ligand sources and candidate receptor receivers',fontsize=9)
fig.text(.02,.05,'Donor-balanced; >=20 cells per sample-subtype. A positive transcript is not surface receptor protein.\nRare populations have limited donor coverage; exact n is provided in the source table.',fontsize=8)
save(fig,'02_human_ligand_receptor_distribution')

sets=['Cell cycle','NF-kappa B signaling pathway','TNF signaling pathway','PI3K-Akt signaling pathway','MAPK signaling pathway','ErbB signaling pathway','JAK-STAT signaling pathway','TGF-beta signaling pathway','Notch signaling pathway','Wnt signaling pathway','Epithelium Development','Epithelial Cell Differentiation','Epithelial Cell Migration','Keratinocyte Differentiation']
fig,ax=plt.subplots(figsize=(170/25.4,115/25.4));fig.subplots_adjust(left=.49,right=.88,bottom=.16,top=.90)
norm=Normalize(-2.1,2.1);cmap=plt.get_cmap('RdBu_r')
for j,n in enumerate(['KLF4','IGFL1']):
 dd=pd.read_csv(O/f'tables/{n}_GSEA_mechanism_focus.csv')
 for r,term in enumerate(sets):
  q=dd[dd.pathway.str.split('::').str[-1].str.replace(r' \(GO:.*','',regex=True)==term]
  if q.empty:continue
  q=q.iloc[0];ax.scatter(j,r,s=100,facecolor=cmap(norm(q.NES))if q.padj<.05 else 'white',edgecolor=cmap(norm(q.NES)),linewidth=1)
  ax.text(j+.16,r,f'{q.NES:+.2f}',va='center',fontsize=8)
ax.set_ylim(len(sets)-.4,-.8);ax.set_xlim(-.4,1.8);ax.set_yticks(range(len(sets)),sets);ax.set_xticks([0,1],['KLF4-OE','IGFL1']);ax.set_title('Preranked GSEA: full-library BH correction')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.06,pad=.08);cb.set_label('NES (direction of enrichment)')
fig.text(.03,.04,'Filled: BH FDR <0.05 across all tested GO+KEGG sets; open: not significant.\nTranscriptional enrichment does not establish kinase phosphorylation or IGFLR1 dependence.',fontsize=8)
save(fig,'03_GSEA_pathway_evidence')

p=pd.read_csv(O/'tables/human_scrna_patient_subtype_expression.csv');p=p[p.sub_labels=='Fib_K14']
cs=pd.read_csv(O/'tables/human_patient_axis_module_correlations.csv')
fig,axes=plt.subplots(1,2,figsize=(180/25.4,85/25.4));fig.subplots_adjust(left=.09,right=.98,wspace=.38,bottom=.25,top=.78)
for ax,(x,y,xl,yl)in zip(axes,[('KLF4_mean','IGFL1_mean','KLF4: mean log1p(CP10K)','IGFL1: mean log1p(CP10K)'),('IGFL1_mean','wound_keratin_score','IGFL1: mean log1p(CP10K)','Wound-keratin module score')]):
 row=cs[(cs.subtype=='Fib_K14')&(cs.x==x)&(cs.y==y)].iloc[0]
 ax.scatter(p[x],p[y],s=20,edgecolor='#287ca3',facecolor='white');ax.set_xlabel(xl);ax.set_ylabel(yl)
 ax.set_title(f'n={int(row.n_patients)} patients; rho={row.spearman_rho:.2f}\nBH FDR={row.FDR_BH:.3g}')
fig.text(.03,.06,'One point per patient; samples averaged within patient. Association is not regulatory direction.\nFrozen Fib_K14 labels retained; doublet / ambient-RNA uncertainty remains.',fontsize=8)
save(fig,'04_human_FibK14_patient_associations')

sp=pd.read_csv(O/'tables/human_stereo_LR_spatial_distance_DESCRIPTIVE.csv');sp=sp[sp.ligand=='IGFL1']
fig,ax=plt.subplots(figsize=(180/25.4,75/25.4));fig.subplots_adjust(left=.38,right=.95,bottom=.24,top=.90)
names=sp.section.tolist();rat=sp.distance_ratio.to_numpy();ax.scatter(rat,np.arange(len(rat)),s=25,facecolor='white',edgecolor='#287ca3');ax.axvline(1,color='gray',ls='--',lw=.8);ax.set_yticks(range(len(names)),names);ax.invert_yaxis();ax.set_xlabel('Observed / permuted median nearest-receptor distance');ax.set_title('IGFL1-positive to IGFLR1-positive bins: descriptive geometry')
fig.text(.03,.04,'No consistent proximity enrichment across sections. Missing points: insufficient ligand detection.\n200 within-section randomizations; cell-type composition not preserved; no biological P values.',fontsize=8)
save(fig,'05_spatial_IGFL1_IGFLR1_geometry')

a=pd.read_csv(O/'tables/MOUSE_REVIEW_healing_summary_mean_SD_SEM.csv')
fig,axes=plt.subplots(1,2,figsize=(180/25.4,85/25.4),sharey=True);fig.subplots_adjust(left=.09,right=.98,bottom=.24,top=.87,wspace=.22)
for ax,n,groups in zip(axes,['KLF4','IGFL1'],[['OEC','KLF4-OE'],['NC','IGFL1']]):
 for g,c in zip(groups,['#329fbb','#ef6557']):
  q=a[(a.experiment==n)&(a.group==g)].sort_values('time_dpb');ax.plot(q.time_dpb,q.mean_percent,'o-',ms=3,lw=1.2,color=c,label=g)
 ax.set_xticks([3,5,7,10,12]);ax.set_xlabel('Days post burn (dpb)');ax.set_title(n+' H&E (review data)');ax.legend(frameon=False)
axes[0].set_ylabel('Mean re-epithelialization (%)');axes[0].set_ylim(0,105)
fig.text(.03,.05,'Slice-level descriptive means, not independent-animal efficacy tests. KLF4 D12 includes unconfirmed\nmorphology calls; IGFL1 D12 contains a source discrepancy. No mechanistic dependency is established.',fontsize=8)
save(fig,'06_mouse_HE_DESCRIPTIVE_NOT_FINAL')
(O/'reports/figure_index.json').write_text(json.dumps(files,indent=2),'utf8')
print('FIGURES',files,flush=True)
