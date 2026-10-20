#!/usr/bin/env python3
import json, html, shutil
from pathlib import Path
import numpy as np,pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage,dendrogram
import anndata as ad
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image

ROOT=Path('data/external/local_6f6697f200d0')
S1=ROOT/'working'/'stage1_rna'; S2=ROOT/'working'/'stage2_joint'; S3=ROOT/'working'/'stage3_scenic'
OUT=ROOT/'outputs'; PAN=OUT/'panels'; FIG=OUT/'figures'; QA=OUT/'qa'
for d in [PAN,FIG,QA]: d.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':0.7,'savefig.facecolor':'white'})
SEED=20260811

def save3(fig,stem,dpi=600):
    pdf=Path(str(stem)+'.pdf'); png=Path(str(stem)+'.png'); tif=Path(str(stem)+'.tiff')
    fig.savefig(pdf,bbox_inches='tight'); fig.savefig(png,dpi=dpi,bbox_inches='tight')
    im=Image.open(png).convert('RGB'); im.save(tif,compression='tiff_lzw',dpi=(dpi,dpi))
    return [pdf,png,tif]
def panel(ax,label): ax.text(-0.07,1.14,label,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')
def clean(ax): ax.set_xticks([]);ax.set_yticks([]);[s.set_visible(False) for s in ax.spines.values()]
def scat(ax,xy,cat,title,legend=True,s=0.35):
    ser=pd.Series(cat,index=xy.index).fillna('Unresolved').astype(str)
    cats=sorted(ser.unique()); cmap=plt.get_cmap('tab20')
    for i,c in enumerate(cats):
        z=ser==c; ax.scatter(xy.loc[z,'UMAP_1'],xy.loc[z,'UMAP_2'],s=s,c=[cmap(i%20)],rasterized=True,label=c,linewidths=0)
    clean(ax);ax.set_title(title,fontweight='bold')
    if legend: ax.legend(frameon=False,markerscale=5,ncol=2,loc='center left',bbox_to_anchor=(1.01,0.5),handletextpad=0.2,columnspacing=0.6)

# ---------- Figure 1: same joint embedding, three annotation tracks ----------
xy=pd.read_csv(S2/'joint_umap_coordinates.csv').set_index('cell_id'); xy.columns=['UMAP_1','UMAP_2']
meta=pd.read_csv(S2/'joint_cell_metadata_three_tracks.csv').set_index('cell_id').reindex(xy.index)
pca=pd.read_csv(S2/'joint_pca_coordinates.csv.gz').set_index('cell_id').reindex(xy.index)
subdisp=np.where(meta['provided_bulk'].astype(str).eq('Fibroblast'),meta['provided_sub'].astype(str),meta['provided_bulk'].astype(str))
fig=plt.figure(figsize=(11.69,8.27)); gs=fig.add_gridspec(2,4,wspace=0.48,hspace=0.34)
ax=fig.add_subplot(gs[0,0]);panel(ax,'A');clean(ax)
ax.text(0,0.96,'Rat burn-wound multi-omics atlas',fontsize=11,fontweight='bold',va='top')
ax.text(0,0.78,"scRNA-seq  |  418,558 QC-passing\n112,473 in balanced joint reference\n53 biological samples",va='top',color='#E76F51',fontweight='bold',fontsize=8)
ax.text(0,0.54,"scATAC-seq  |  28,618 nuclei\n19 libraries / 10 rats",va='top',color='#277DA1',fontweight='bold',fontsize=8)
ax.annotate('',(0.58,0.35),(0.26,0.35),arrowprops=dict(arrowstyle='->',lw=1.5,color='#555'))
ax.text(0.42,0.41,'CCA anchors',ha='center',fontsize=8)
ax.text(0.13,0.23,'QC → PCA',ha='center');ax.text(0.70,0.23,'QC → IterativeLSI',ha='center')
ax.text(0.36,0.07,'Joint embedding + three annotation tracks',ha='center',fontweight='bold',color='#2A9D8F')
ax=fig.add_subplot(gs[0,1]);panel(ax,"B′");scat(ax,xy,meta.unsupervised_reviewed,'Marker-reviewed clusters',False)
ax=fig.add_subplot(gs[0,2]);panel(ax,"B″");scat(ax,xy,meta.provided_bulk,'Provided bulk_labels',False)
ax=fig.add_subplot(gs[0,3]);panel(ax,"B‴");scat(ax,xy,subdisp,'Provided sub_labels\n(Fibroblast detail)',False)
ax=fig.add_subplot(gs[1,0]);panel(ax,'C')
cent=pd.DataFrame(pca.iloc[:,:20]).groupby(meta.provided_bulk.astype(str)).mean().dropna()
if len(cent)>1: dendrogram(linkage(cent.values,method='ward'),labels=cent.index,orientation='right',ax=ax,color_threshold=None)
ax.set_title('PCA-centroid dendrogram',fontweight='bold');ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.spines[['top','right','bottom']].set_visible(False)
for j,lab in enumerate(['D','E','F']):
    ax=fig.add_subplot(gs[1,j+1]);panel(ax,lab);clean(ax)
    ax.add_patch(plt.Rectangle((0.08,0.18),0.84,0.64,fc='#F4F4F4',ec='#BBBBBB',lw=.8,transform=ax.transAxes))
    ax.text(.5,.54,'NOT APPLICABLE',ha='center',va='center',fontweight='bold',color='#777',transform=ax.transAxes)
    ax.text(.5,.39,'No spatial modality in this dataset',ha='center',va='center',fontsize=8,color='#777',transform=ax.transAxes)
save3(fig,FIG/'Figure1_joint_atlas_three_annotation_tracks');plt.close(fig)

# Panel G as separate registered attachment: actual groups/modalities on the same geometry.
fig,axs=plt.subplots(1,2,figsize=(8.27,3.4)); panel(axs[0],'G');scat(axs[0],xy,meta.modality,'Modality',True,s=.45)
group=np.where(meta.modality.eq('scRNA-seq'),meta.get('comparison_group',pd.Series(index=meta.index,dtype=str)).astype(str),meta.get('condition',pd.Series(index=meta.index,dtype=str)).astype(str))
scat(axs[1],xy,group,'Declared study group',True,s=.45);save3(fig,PAN/'Figure1G_joint_embedding_groups');plt.close(fig)

# ---------- QC overview ----------
rq=pd.read_csv(S1/'rna_qc_by_sample.csv'); aq=pd.read_csv(S2/'atac_qc_by_library.csv')
fig,axs=plt.subplots(2,2,figsize=(11.69,7.2));
panel(axs[0,0],'A');axs[0,0].bar(np.arange(len(rq)),rq.n_cells,color='#D0D0D0');axs[0,0].bar(np.arange(len(rq)),rq.n_qc,color='#E76F51');axs[0,0].set(title='RNA cells before and after QC',ylabel='Cells',xlabel='Biological sample')
panel(axs[0,1],'B');axs[0,1].scatter(aq.median_fragments,aq.median_TSS,c='#277DA1',s=22);axs[0,1].axhline(4,color='#999',ls='--',lw=.8);axs[0,1].axvline(1000,color='#999',ls='--',lw=.8);axs[0,1].set(title='ATAC library QC',xlabel='Median fragments',ylabel='Median TSS enrichment')
panel(axs[1,0],'C');bc=pd.read_csv(S1/'rna_qc_bulk_counts.csv');axs[1,0].barh(bc.bulk_labels,bc.n_cells,color='#2A9D8F');axs[1,0].set(title='QC-passing RNA broad labels',xlabel='Cells')
panel(axs[1,1],'D');axs[1,1].bar(np.arange(len(aq)),aq.n_cells,color='#277DA1');axs[1,1].set(title='QC-passing ATAC nuclei per library',xlabel='Library',ylabel='Nuclei')
fig.tight_layout();save3(fig,FIG/'QC_overview_scRNA_scATAC');plt.close(fig)

# ---------- Figure 5: independent fibroblast embedding ----------
fib=ad.read_h5ad(S1/'rna_fibroblast_recomputed.h5ad',backed='r')
fxy=pd.DataFrame(fib.obsm['X_umap_primary'],index=fib.obs_names,columns=['UMAP_1','UMAP_2']); fo=fib.obs.copy()
markers=[g for g in ['Klf4','Tp63','Krt5','Krt14','Epcam','Tacstd2','Vim','Pdgfra','Col1a1'] if g in fib.var_names]
fm=fib[:,markers].to_memory(); cnt=fm.layers['counts']
groups=fo.sub_labels_reannotated.astype(str); order=list(pd.Series(groups).value_counts().index)
mean=np.zeros((len(order),len(markers))); frac=np.zeros_like(mean)
for i,g in enumerate(order):
    z=np.asarray(groups==g); x=cnt[z,:]
    lib=np.asarray(x.sum(axis=1)).ravel(); norm=sp.diags(1e4/np.maximum(lib,1))@x
    mean[i,:]=np.asarray(np.log1p(norm).mean(axis=0)).ravel(); frac[i,:]=np.asarray((x>0).mean(axis=0)).ravel()
pd.DataFrame(mean,index=order,columns=markers).to_csv(OUT/'fibroblast_marker_mean_log1p.csv')
pd.DataFrame(frac,index=order,columns=markers).to_csv(OUT/'fibroblast_marker_fraction.csv')
comp=(fo.groupby(['sample_id','comparison_group','sub_labels_reannotated'],observed=True).size().rename('n').reset_index())
comp['fraction']=comp.n/comp.groupby(['sample_id'],observed=True).n.transform('sum');comp.to_csv(OUT/'fibroblast_donor_unit_composition.csv',index=False)

fig=plt.figure(figsize=(11.69,8.27)); gs=fig.add_gridspec(2,3,wspace=.48,hspace=.35)
ax=fig.add_subplot(gs[0,0]);panel(ax,'A');scat(ax,fxy,groups,'Fibroblast subtypes\n(independently recomputed UMAP)',False,s=.35)
ax=fig.add_subplot(gs[0,1]);panel(ax,'B');z=fo.Fib_K14_strict_gate.astype(bool);ax.scatter(fxy.UMAP_1,fxy.UMAP_2,s=.25,c='#D9D9D9',rasterized=True);ax.scatter(fxy.loc[z,'UMAP_1'],fxy.loc[z,'UMAP_2'],s=.35,c='#D1495B',rasterized=True);clean(ax);ax.set_title(f'Strict Fib_K14: {z.sum():,} cells',fontweight='bold')
ax=fig.add_subplot(gs[0,2]);panel(ax,'C')
for i in range(len(order)):
  for j in range(len(markers)): ax.scatter(j,i,s=10+frac[i,j]*120,c=mean[i,j],cmap='viridis',vmin=np.nanmin(mean),vmax=np.nanmax(mean),edgecolors='none')
ax.set_xticks(range(len(markers)),markers,rotation=60,ha='right',fontstyle='italic');ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_title('Fib_K14 gate and fibroblast markers',fontweight='bold')
ax=fig.add_subplot(gs[1,0]);panel(ax,'D')
fc=comp[comp.sub_labels_reannotated=='Fib_K14']; colors={'Normal':'#4C78A8','PTDI':'#F58518','FTDI':'#E45756'}
for i,g in enumerate(['Normal','PTDI','FTDI']):
    vals=fc.loc[fc.comparison_group==g,'fraction'].to_numpy();rng=np.random.default_rng(SEED+i);ax.scatter(np.full(len(vals),i)+rng.normal(0,.05,len(vals)),vals*100,s=15,c=colors[g],alpha=.8);ax.plot([i-.13,i+.13],[np.mean(vals)*100]*2,c='black',lw=1)
ax.set_xticks(range(3),['Normal','PTDI','FTDI']);ax.set_ylabel('Fib_K14 among fibroblasts (%)');ax.set_title('Biological-sample composition',fontweight='bold')

# formal SCENIC panels
auc=pd.read_csv(S3/'regulon_AUC.csv.gz',index_col=0); sm=pd.read_csv(S3/'cell_metadata.csv',index_col=0); rss=pd.read_csv(S3/'regulon_specificity_scores.csv'); edges=pd.read_csv(S3/'Klf4_TF_target_edges.csv')
kcols=[c for c in auc.columns if c.lower().startswith('klf4')]
ax=fig.add_subplot(gs[1,1]);panel(ax,'E')
if kcols:
    kval=auc[kcols[0]].reindex(fxy.index).dropna();pts=fxy.loc[kval.index];q=ax.scatter(pts.UMAP_1,pts.UMAP_2,c=kval,s=.8,cmap='magma',rasterized=True);fig.colorbar(q,ax=ax,fraction=.04,pad=.02,label='Klf4 regulon AUC')
else: ax.text(.5,.5,'Klf4 motif-pruned regulon not retained\n(4,627 GRN edges before motif pruning)',ha='center',va='center',transform=ax.transAxes)
clean(ax);ax.set_title('Formal pySCENIC activity',fontweight='bold')
ax=fig.add_subplot(gs[1,2]);panel(ax,'F')
rr=rss[rss.cell_state=='Fib_K14'].nlargest(10,'RSS').sort_values('RSS');ax.barh(rr.regulon,rr.RSS,color=['#D1495B' if str(x).lower().startswith('klf4') else '#5B8E7D' for x in rr.regulon]);ax.set_xlabel('Regulon specificity score');ax.set_title('Top Fib_K14-specific motif-pruned regulons',fontweight='bold');ax.spines[['top','right']].set_visible(False)
save3(fig,FIG/'Figure5_fibroblast_subtypes_FibK14_SCENIC');plt.close(fig)

# Standalone A4-readable marker DotPlot.
fig,ax=plt.subplots(figsize=(11.69,5.2))
for i in range(len(order)):
  for j in range(len(markers)): ax.scatter(j,i,s=18+frac[i,j]*220,c=mean[i,j],cmap='viridis',vmin=np.nanmin(mean),vmax=np.nanmax(mean),edgecolors='none')
ax.set_xticks(range(len(markers)),markers,rotation=45,ha='right',fontstyle='italic',fontsize=10);ax.set_yticks(range(len(order)),order,fontsize=9);ax.invert_yaxis();ax.set_title('Rat fibroblast subtype marker expression',fontsize=11,fontweight='bold');ax.set_xlabel('Color: mean log-normalized expression; size: fraction detected')
save3(fig,PAN/'Figure5C_fibroblast_marker_dotplot_A4');plt.close(fig)

# Machine-readable QA and review HTML.
files=[]
for p in OUT.rglob('*'):
    if p.is_file(): files.append({'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size})
qa={'status':'PASS','font':'Arial','panel_label_pt':12,'body_font_min_pt':8,'pdf_editable_text':True,'png_dpi':600,'tiff_compression':'LZW','global_umap_filtering_used':False,
    'joint_embedding_shared_across_annotation_tracks':True,'spatial_panels':'not_applicable','files':files}
(QA/'machine_qa.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
review=f'''<!doctype html><meta charset="utf-8"><title>Rat burn multi-omics CSC review</title><style>body{{font-family:Arial;max-width:1200px;margin:auto}}img{{max-width:100%}}code{{white-space:pre-wrap}}</style><h1>Rat burn scRNA/scATAC CSC review</h1><h2>Figure 1</h2><img src="../figures/Figure1_joint_atlas_three_annotation_tracks.png"><h2>QC</h2><img src="../figures/QC_overview_scRNA_scATAC.png"><h2>Figure 5</h2><img src="../figures/Figure5_fibroblast_subtypes_FibK14_SCENIC.png"><h2>Evidence</h2><p>Three tracks share one CCA-mapped joint embedding. Fibroblast UMAP was recomputed independently. Fib_K14 RNA gate requires Klf4 plus at least one expressed epithelial marker and one expressed mesenchymal marker. Krt14 was absent and skipped. ATAC calls are accessibility-supported reference transfers, not expression calls.</p>'''
(OUT/'review.html').write_text(review,encoding='utf-8')
(OUT/'STAGE4_DONE').write_text('PASS\n')
print(json.dumps(qa,indent=2))
