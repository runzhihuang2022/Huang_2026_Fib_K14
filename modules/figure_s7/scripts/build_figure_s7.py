from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from scipy.stats import spearmanr
R=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
D=R/'source_data';O=R/'outputs';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6})
fig=plt.figure(figsize=(210/25.4,297/25.4),facecolor='white')
def tx(x,y,s,**kw):return fig.text(x/210,1-y/297,s,va='top',**kw)
def axmm(x,y,w,h):return fig.add_axes([x/210,1-(y+h)/297,w/210,h/297])
def title(letter,x,y,s):tx(x,y,letter,fontsize=12,fontweight='bold');tx(x+6,y+.4,s,fontsize=8,fontweight='bold')
for letter,y,groups in [('A',6,['KLF4-OE','OEC']),('B',48,['IGFL1','NC'])]:
 title(letter,6,y,'')
 for j,t in enumerate([3,5,7,10,12]):tx(31.25+j*38,y+1,f'{t} dpb',fontsize=7,ha='center')
 for k,g in enumerate(groups):
  tx(7,y+11+k*15,g,fontsize=6.5,rotation=90,ha='center',va='center') if False else None
  fig.text(9/210,1-(y+16+k*15)/297,g,rotation=90,fontsize=6.5,va='center',ha='center')
  for j in range(5):
   a=axmm(14+j*38, y+8+k*15,34.5,11.9)
   a.imshow(Image.open(D/f'HE_{letter}_{k*5+j:02}.png'));a.set_axis_off()
 tx(204,y+37,'Scale bars, 1 mm',fontsize=6,ha='right')

title('C',6,92,'Shared and distinct transcriptional responses')
ks={n:pd.read_csv(D/f'{n}_DESeq2_all.csv').set_index('gene_name')for n in ['KLF4','IGFL1']}
genes=['KLF4','KRT14','KRT17','PDGFRA','VIM','EPCAM','TACSTD2','TGFA','HBEGF','EGFR','DUSP6','NFKBIA','MKI67','PCNA','FGF7']
arr=np.array([[ks[n].loc[g,'log2FoldChange'] for n in ks]for g in genes]);a=axmm(27,105,58,65)
im=a.imshow(arr,cmap='RdBu_r',vmin=-5,vmax=5,aspect='auto');a.set_yticks(range(len(genes)),genes,fontstyle='italic');a.set_xticks([0,1],['KLF4-OE','IGFL1']);a.xaxis.tick_top();a.tick_params(axis='both',length=0,pad=3)
for j,n in enumerate(ks):
 for i,g in enumerate(genes):
  row=ks[n].loc[g];stars='***'if row.padj<.001 else '**'if row.padj<.01 else '*'if row.padj<.05 else ''
  a.text(j,i,f'{row.log2FoldChange:+.2f}{stars}',ha='center',va='center',fontsize=6.5,color='white'if abs(row.log2FoldChange)>2.8 else 'black')
cax=axmm(89,105,2.7,65);cb=fig.colorbar(im,cax=cax);cb.set_ticks([-5,0,5]);cb.set_label('log2 fold change',fontsize=6.5,labelpad=3)
tx(8,174,'KLF4-OE: 2 vs 2 libraries; IGFL1: 2 vs 3 libraries.',fontsize=6)
tx(8,178,'* FDR <0.05; ** <0.01; *** <0.001 (gene-wide BH).',fontsize=6)
tx(8,182,'IGFL1 raw counts in the KLF4 experiment: 0 / 0 / 0 / 0.',fontsize=6)

title('D',109,92,'Pathway-level transcriptional responses')
terms=['Cell cycle','NF-kappa B signaling pathway','TNF signaling pathway','PI3K-Akt signaling pathway','MAPK signaling pathway','ErbB signaling pathway','Epithelial Cell Differentiation','Epithelial Cell Migration']
labels=['Cell cycle','NF-κB','TNF','PI3K–AKT','MAPK','ErbB','Epithelial differentiation','Epithelial migration']
a=axmm(150,108,43,60);norm=Normalize(-2.1,2.1);cmap=plt.get_cmap('RdBu_r');stats=[]
for j,n in enumerate(['KLF4','IGFL1']):
 dd=pd.read_csv(D/f'{n}_GSEA_mechanism_focus.csv')
 for i,term in enumerate(terms):
  row=dd[dd.pathway.str.split('::').str[-1].str.replace(r' \(GO:.*','',regex=True)==term].iloc[0]
  a.scatter(j,i,s=35,facecolor=cmap(norm(row.NES))if row.padj<.05 else 'white',edgecolor=cmap(norm(row.NES)),linewidth=1)
  a.text(j+.13,i,f'{row.NES:+.2f}',fontsize=6,va='center')
  stats.append({'experiment':n,'pathway':row.pathway,'NES':row.NES,'FDR':row.padj})
a.set_xlim(-.25,1.65);a.set_ylim(7.6,-.7);a.set_yticks(range(8),labels);a.set_xticks([0,1],['KLF4-OE','IGFL1']);a.xaxis.tick_top();a.tick_params(length=0,pad=3);a.spines[['left','bottom']].set_visible(False)
tx(111,174,'Numbers: NES. Filled symbols: full-library FDR <0.05.',fontsize=6)
tx(111,178,'Open symbols: FDR ≥0.05; GO + KEGG corrected jointly.',fontsize=6)
tx(111,182,'IGFL1: NF-κB FDR=0.0233; cell cycle FDR=0.000409.',fontsize=6)

title('E',6,191,'Patient-level associations within Fib_K14')
p=pd.read_csv(D/'human_scrna_patient_subtype_expression.csv');p=p[p.sub_labels=='Fib_K14'];cs=pd.read_csv(D/'human_patient_axis_module_correlations.csv')
for x,y,xl,yl,xx in [('KLF4_mean','IGFL1_mean','KLF4','IGFL1',18),('IGFL1_mean','wound_keratin_score','IGFL1','Wound-keratin score',65)]:
 row=cs[(cs.subtype=='Fib_K14')&(cs.x==x)&(cs.y==y)].iloc[0];assert len(p)==row.n_patients==20
 assert abs(spearmanr(p[x],p[y]).statistic-row.spearman_rho)<1e-12
 a=axmm(xx,213,33,48);a.scatter(p[x],p[y],s=18,edgecolor='#217a91',facecolor='white',linewidth=.8)
 a.set_xlabel(xl);a.set_ylabel(yl,labelpad=2);a.tick_params(pad=2);a.locator_params(nbins=3)
 a.set_title(f'ρ = {row.spearman_rho:.3f}\nFDR = {row.FDR_BH:.3g}',fontsize=7,pad=6)
tx(8,272,'One point per patient (n=20); samples averaged per patient.',fontsize=6)
tx(8,277,'Gene expression: mean log1p(CP10K).',fontsize=6)
tx(8,282,'Module: KRT6A/B, KRT16/17, LAMC2, LAMB3, TACSTD2.',fontsize=6)
tx(8,287,'Spearman correlations; BH across all tested associations.',fontsize=6)

title('F',109,191,'Ligand and candidate receptor distribution')
d=pd.read_csv(D/'human_scrna_donor_balanced_sources_receivers.csv').set_index('sub_labels')
sel=['KC_Spinous_Mig','KC_Granular','KC_Basal_Mig','Fib_K14','SAC_SG_Progenitor','CD4_T','CD8_T','Mast','cDC2','pDC']
a=axmm(150,205,47,61);yy=np.arange(len(sel))
for k,(g,col)in enumerate([('IGFL1','#329fbb'),('IGFLR1','#ee7960')]):a.barh(yy+(k-.5)*.30,d.loc[sel,g+'_pct'],height=.28,color=col,label=g)
a.set_yticks(yy,[s+f' ({int(d.loc[s,"n_patients"])})' for s in sel],fontsize=6);a.invert_yaxis();a.set_xlim(0,22);a.set_xticks([0,10,20]);a.set_xlabel('Mean detected cells (%)',fontsize=6.5);a.legend(frameon=False,loc='upper center',bbox_to_anchor=(.48,1.16),ncol=2,fontsize=6)
tx(111,277,'Parentheses: patients per cell state; patient-weighted means.',fontsize=6)
tx(111,282,f'Fib_K14: IGFLR1 detected in {d.loc["Fib_K14","IGFLR1_pct"]:.2f}% of cells.',fontsize=6)
tx(111,287,'Detection: >0 counts; ≥20 cells per sample–state group.',fontsize=6)
for ext in ['pdf','svg','png','tiff']:
 fig.savefig(O/f'Figure_S7_V12.{ext}',dpi=600,facecolor='white',**({'pil_kwargs':{'compression':'tiff_lzw'}}if ext=='tiff'else{}))
fig.savefig(O/'Figure_S7_V12_preview.png',dpi=180,facecolor='white')
from matplotlib.transforms import Bbox
for label,(x,y,w,h) in {'A':(5,5,201,39),'B':(5,47,201,39),'C':(5,91,99,96),'D':(108,91,98,96),'E':(5,190,99,102),'F':(108,190,98,102)}.items():
 for ext in ['svg','pdf']:
  fig.savefig(O/f'Figure_S7{label}_V12.{ext}',bbox_inches=Bbox.from_bounds(x/25.4,(297-y-h)/25.4,w/25.4,h/25.4),facecolor='white')
Image.open(O/'Figure_S7_V12.png').resize((2480,3508),Image.Resampling.LANCZOS).save(O/'Figure_S7_V12_Word.png',dpi=(300,300))
pd.DataFrame(stats).to_csv(O/'S7D_displayed_pathways.csv',index=False)
aud={'page_mm':[210,297],'new_panel_labels':['C','D','E','F'],'retained_panel_labels':['A','B'],'font_family':'Arial','vector_font_range_pt':[6,12],'raster_labels':'Original embedded histology scale-bar text preserved; smaller source text is not independently editable. New group/time labels and scale-bar caption are vector text.','histology_adjustment':'Crop full pre-existing frames; uniform proportional scaling; no color/intensity alteration','IGFL1_KLF4_log2FC':ks['IGFL1'].loc['KLF4','log2FoldChange'],'IGFL1_KLF4_FDR':ks['IGFL1'].loc['KLF4','padj'],'FibK14_IGFLR1_pct':d.loc['Fib_K14','IGFLR1_pct'],'patient_correlations_recomputed':True,'statistical_models_rerun_this_command':False,'AI_status':'Editable SVG/PDF provided; native Illustrator AI conversion pending'}
(O/'figure_QA.json').write_text(json.dumps(aud,indent=2),encoding='utf8');plt.close(fig)
print(json.dumps(aud))
