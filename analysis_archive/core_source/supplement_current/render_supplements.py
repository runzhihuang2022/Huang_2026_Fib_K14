"""Matched A-H layout from auditable sources. H gaps are explicit, never zero-filled."""
from pathlib import Path
import json,textwrap,math,re,copy,xml.etree.ElementTree as ET,argparse
import numpy as np,pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from scipy import stats
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
ET.register_namespace('','http://www.w3.org/2000/svg')
ET.register_namespace('xlink','http://www.w3.org/1999/xlink')
cfg=json.loads((ROOT/'source_data/figure_config.json').read_text('utf8'))
IF=json.loads((ROOT/'source_data/IF_source_and_inset_manifest.json').read_text('utf8'))
W,H=210.,297.;MM=1/25.4
mpl.rcParams.update({'font.family':'Arial','font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'svg.fonttype':'none',
 'svg.image_inline':True,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False,
 'axes.linewidth':.7,'lines.linewidth':.8,'figure.facecolor':'white','axes.facecolor':'white'})
np.random.seed(20261016)
panel_boxes={'A':[4,4,54,81],'B':[62,4,144,81],'C':[4,91,24,62],
 'D':[32,91,89,62],'E':[125,91,81,62],'F':[4,158,146,98],
 'G':[154,158,52,98],'H':[4,264,202,29]}
fontrows=[];manifests=[]
label_offsets={
 'KC_Basal_Mig':[-2,-.8],'KC_Spinous':[-1.8,.7],
 'KC_Spinous_Mig':[3,.1],'KC_Spinous_Mat':[.5,1.4],'KC_Granular':[2,.2],
 'SAC_SG_Progenitor':[-2,-1],'SAC_SG_Dark':[.2,1],
 'SAC_HF_IRS':[-2,-2],'SAC_SG_Clear':[2,-1],'SAC_HF_HFSC':[-2.5,1],
 'SAC_HF_DP_DS':[4,1],'SAC_HF_Matrix':[3,3],'SAC_HF_ORS':[4,-.2]}

def setup():
 fig=plt.figure(figsize=(W*MM,H*MM),dpi=160)
 def ax(x,y,w,h,gid=None):
  a=fig.add_axes([x/W,1-(y+h)/H,w/W,h/H]);a.set_gid(gid) if gid else None;return a
 def text(x,y,t,size=7.5,**kw):return fig.text(x/W,1-y/H,t,fontsize=size,**kw)
 for l,(x,y,w,h) in panel_boxes.items():text(x,y+3,l,12,fontweight='bold',va='center')
 return fig,ax,text
def clean(a):
 a.spines[['top','right']].set_visible(False);a.tick_params(length=2,width=.7,pad=1.5)
def hidden(a):
 a.set(xticks=[],yticks=[])
 for s in a.spines.values():s.set_visible(False)
def abbreviation(s,lineage):return s.replace(lineage+'_','',1)
def save(fig,stem):
 def owner(x,y):
  candidates=[(l,b)for l,b in panel_boxes.items() if b[0]-1<=x<=b[0]+b[2]+1 and b[1]-3<=y<=b[1]+b[3]+3]
  if candidates:return candidates[0][0]
  return min(panel_boxes,key=lambda l:(x-panel_boxes[l][0]-panel_boxes[l][2]/2)**2+(y-panel_boxes[l][1]-panel_boxes[l][3]/2)**2)
 for i,a in enumerate(fig.axes):
  if not(a.get_gid() or '').startswith('Panel_'):
   b=a.get_position();a.set_gid(f'Panel_{owner((b.x0+b.width/2)*W,(1-b.y0-b.height/2)*H)}_axis_{i}')
 for i,t in enumerate(fig.texts):
  x,y=t.get_position();t.set_gid(f'Panel_{owner(x*W,(1-y)*H)}_text_{i}')
 fig.canvas.draw()
 for t in fig.findobj(match=mpl.text.Text):
  if not t.get_text():continue
  if t.get_fontsize()<7:t.set_fontsize(7)
  fontrows.append(dict(figure=stem,text=t.get_text(),font=t.get_fontfamily()[0],size=t.get_fontsize()))
 fig.savefig(ROOT/f'outputs/assembled/{stem}.svg',dpi=600)
 svg=ET.parse(ROOT/f'outputs/assembled/{stem}.svg');ns={'s':'http://www.w3.org/2000/svg'}
 for letter in panel_boxes:
  page=copy.deepcopy(svg.getroot());figure=page.find('s:g[@id="figure_1"]',ns)
  for child in list(figure):
   if not child.attrib.get('id','').startswith(f'Panel_{letter}_'):figure.remove(child)
  ET.ElementTree(page).write(ROOT/f'outputs/panels/{stem}_Panel_{letter}.svg',encoding='utf-8',xml_declaration=True)
 fig.savefig(ROOT/f'QC/{stem}_vector_master.pdf',dpi=600)
 # Small proof first; final 600-dpi raster is rendered from Illustrator-saved AI.
 fig.savefig(ROOT/f'QC/{stem}_preview.png',dpi=220)
 plt.close(fig)

parser=argparse.ArgumentParser();parser.add_argument('--lineage',choices=['KC','SAC']);args=parser.parse_args()
for lineage,figure_name in [('KC','FigureS2'),('SAC','FigureS3')]:
 if args.lineage and lineage!=args.lineage:continue
 subtypes=cfg['LINEAGES'][lineage];target=cfg['targets'][lineage];palette=cfg['RAW_COLORS'];program=cfg['marker_programs'][lineage]
 emb=np.load(ROOT/f'source_data/{lineage}_embedding_inputs.npz');xy=emb['xy'];labels=emb['labels'];proxy=emb['neighbour_density_log2FC']
 patient=pd.read_csv(ROOT/f'source_data/{lineage}_patient_proportions.tsv',sep='\t');counts=pd.read_csv(ROOT/f'source_data/{lineage}_sample_counts.tsv',sep='\t')
 da=pd.read_csv(ROOT/f'source_data/{lineage}_recovered_count_comparison.tsv',sep='\t').sort_values('P_sample_nominal',ascending=False)
 fig,ax,text=setup()
 text(30,12,'Keratinocyte subtypes' if lineage=='KC' else 'Skin appendage subtypes',9.5,fontweight='bold',ha='center')
 a=ax(7,17,49,49,'Panel_A_UMAP');order=np.random.default_rng(2026).permutation(len(xy))
 for st in subtypes:
  ii=order[labels[order]==st];a.scatter(xy[ii,0],xy[ii,1],s=.12,c=palette[st],alpha=.8,linewidths=0,rasterized=True)
 ext=np.ptp(xy,axis=0);centers=[]
 for st in subtypes:
  center=np.median(xy[labels==st],axis=0);offset=np.array(label_offsets.get(st,[0,0]))
  pos=center+offset;centers.append(pos)
  if np.linalg.norm(offset):a.plot([center[0],pos[0]],[center[1],pos[1]],color='.4',lw=.45,zorder=2)
  a.text(*pos,abbreviation(st,lineage),fontsize=7,fontweight='bold',ha='center',bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=.1))
 a.set_aspect('equal');hidden(a)
 legend=ax(6,69,51,15,'Panel_A_legend');legend.axis('off')
 handles=[Line2D([0],[0],marker='o',linestyle='',markerfacecolor=palette[s],markeredgecolor='none',markersize=3)for s in subtypes]
 legend.legend(handles,[abbreviation(s,lineage)for s in subtypes],ncol=2 if lineage=='SAC' else 3,loc='upper center',frameon=False,handletextpad=.2,columnspacing=.5,labelspacing=.3)

 # B: identical percent-size / Reds expression scale and block separators.
 blocks=pd.read_csv(ROOT/f'source_data/{lineage}_marker_blocks.tsv',sep='\t');genes=blocks.gene.tolist();owners=blocks.subtype.tolist()
 pct=pd.read_csv(ROOT/f'source_data/{lineage}_dotplot_percent_expressed.tsv',sep='\t',index_col=0).loc[subtypes,genes].to_numpy()
 z=pd.read_csv(ROOT/f'source_data/{lineage}_dotplot_average_expression_zscore.tsv',sep='\t',index_col=0).loc[subtypes,genes].to_numpy()
 b=ax(87,20,119,48,'Panel_B_dotplot');xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(subtypes)))
 dots=b.scatter(xx.ravel(),yy.ravel(),s=(1+45*(pct.ravel()/100)**1.5),c=z.ravel(),cmap='Reds',vmin=-1,vmax=2.5,edgecolors='none')
 b.set(xlim=(-.5,len(genes)-.5),ylim=(len(subtypes)-.5,-.5),yticks=range(len(subtypes)),yticklabels=subtypes,xticks=range(len(genes)))
 b.set_xticklabels(genes,rotation=90,fontstyle='italic');b.tick_params(length=1.5,pad=1)
 for j in range(1,len(owners)):
  if owners[j]!=owners[j-1]:b.axvline(j-.5,color='#EF8F84',lw=.6,ls=(0,(3,2)))
 for j,st in enumerate(owners):b.add_patch(Rectangle((j-.5,-1.01),1,.26,color=palette[st],clip_on=False,lw=0))
 text(87,12,'Percent expressed',7.5)
 for k,val in enumerate([25,50,75,100]):
  lx=117+6.5*k;text(lx,9,str(val),7,ha='center');la=ax(lx-1.3,11,2.6,2.6);la.scatter([0],[0],s=1+45*(val/100)**1.5,c='black');la.set(xlim=(-1,1),ylim=(-1,1));hidden(la)
 text(148,12,'Average expression',7.5)
 ca=ax(178,11,28,1.8,'Panel_B_expression_key');cb=fig.colorbar(dots,cax=ca,orientation='horizontal',ticks=[0,1,2]);cb.outline.set_visible(False);ca.xaxis.set_ticks_position('top');ca.tick_params(labelsize=7,length=1,pad=.5)

 # C: patient-mean composition inside the requested parent lineage.
 c=ax(10,101,17.5,44,'Panel_C_composition');bottom=np.zeros(4)
 for st in subtypes:
  vals=np.array([patient.loc[(patient.sub_labels==st)&(patient.Time_category==m),'proportion'].mean()*100 for m in cfg['MODES']])
  c.bar(range(4),vals,bottom=bottom,color=palette[st],width=.78,edgecolor='white',linewidth=.3);bottom+=vals
 c.set(ylim=(0,100),yticks=[0,25,50,75,100],xticks=range(4),ylabel=f'Proportion of {lineage} cells (%)');c.set_xticklabels(['Mode 1','Mode 2','Mode 3','Mode 4'],rotation=40,ha='right');clean(c)

 # D: original nearest-neighbour proxy is named accurately; counts have same metric in point plot and test.
 text(37,96,'Neighbourhood density change',7,fontweight='bold')
 d=ax(35,101,25,28,'Panel_D_neighbour_proxy');lim=max(2.5,np.quantile(abs(proxy),.995))
 d.scatter(xy[:,0],xy[:,1],c=proxy,cmap='coolwarm',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),s=.2,linewidths=0,rasterized=True);d.set_aspect('equal');hidden(d)
 targetxy=xy[labels==target];lo=np.quantile(targetxy,.01,axis=0);hi=np.quantile(targetxy,.99,axis=0)
 d.add_patch(Rectangle(lo,*(hi-lo),fill=False,color=palette[target],linestyle=(0,(3,2)),lw=.65))
 ca=ax(43,131,15,1.5);cb=fig.colorbar(mpl.cm.ScalarMappable(norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),cmap='coolwarm'),cax=ca,orientation='horizontal',ticks=[-2,0,2]);cb.outline.set_visible(False);ca.tick_params(labelsize=7,length=1,pad=.5)
 text(47,138,'Log2 density ratio',7,ha='center')
 countax=ax(39,142,21,10.5,'Panel_D_target_counts')
 for j,m in enumerate(cfg['MODES']):
  v=counts.loc[counts.Time_category==m,target].to_numpy();sem=v.std(ddof=1)/np.sqrt(len(v))
  countax.bar(j,v.mean(),color=palette[target],alpha=.6,width=.62,edgecolor='black',lw=.35)
  countax.errorbar(j,v.mean(),yerr=sem,c='black',capsize=1,lw=.5)
  countax.scatter(j+np.random.uniform(-.14,.14,len(v)),v,s=3,color='#444444',alpha=.8,zorder=3)
 countax.set(xticks=range(4),xticklabels=['M1','M2','M3','M4']);countax.tick_params(labelsize=7);countax.set_title(target,fontsize=7,pad=2);clean(countax)
 countax.set_ylabel('Cells / sample',fontsize=7,labelpad=1)
 text(94,96,'Exploratory count contrast',7.5,fontweight='bold',ha='center')
 d2=ax(95,101,25,50,'Panel_D_recovered_counts');tickpositions=[]
 for i,(_,r) in enumerate(da.iterrows()):
  st=r.sub_labels;base=counts.loc[counts.Time_category==cfg['MODES'][0],st].mean();v=counts.loc[counts.Time_category.isin(cfg['MODES'][1:]),st].to_numpy()
  xx=np.log2((v+.5)/(base+.5));cc=plt.cm.coolwarm(Normalize(-7,7)(xx));d2.scatter(xx,i+np.random.uniform(-.22,.22,len(v)),s=7,c=cc,edgecolors='none')
  q=r.BH_FDR_within_lineage;qs='<0.001' if q<.001 else f'={q:.3f}';color=palette[target] if st==target else ('#C74849' if q<.05 and r.log2FC_count>0 else '#506CBD' if q<.05 else 'black')
  d2.text(-.08,i,abbreviation(st,lineage)+f' (q{qs})',transform=d2.get_yaxis_transform(),ha='right',va='center',fontsize=7,color=color,fontweight='bold' if st==target else 'normal')
  tickpositions.append(i)
 d2.axvline(0,color='.55',lw=.7,ls='--');d2.set(ylim=(-.6,len(da)-.4),yticks=[],xlabel='Log2FC (M2+3+4 / M1)');d2.xaxis.label.set_size(7);d2.set_xticks([-8,-4,0,4]);clean(d2)

 # E: actual lineage-specific marker enrichment; never copied Fib_K14 terms.
 enrichment=pd.read_csv(ROOT/f'source_data/{lineage}_enrichment_displayed.tsv',sep='\t').head(6)
 e=ax(172,101,34,50,'Panel_E_functions');terms=[]
 for term in enrichment.Term:
  term=re.sub(r'\s*\(GO:\d+\)','',term);terms.append('\n'.join(textwrap.wrap(term.lower(),29)))
 values=-np.log10(enrichment['Adjusted P-value'].clip(lower=np.nextafter(0,1)).to_numpy())
 colors=plt.cm.viridis(np.linspace(.92,.23,len(values)))
 e.barh(np.arange(len(values)),values,color=colors,edgecolor='.3',lw=.4,height=.72)
 e.set(yticks=np.arange(len(values)),yticklabels=terms,xlabel='-log10(adj. P)');e.invert_yaxis();e.tick_params(axis='y',labelsize=7,pad=2);e.xaxis.grid(True,color='.85',lw=.35);clean(e)
 for i,q in enumerate(enrichment['Adjusted P-value']):e.text(values[i]+max(values)*.018,i,'***' if q<.001 else '**' if q<.01 else '*',fontsize=7,va='center')
 e.set_xlim(0,max(values)*1.15);text(166,97,target+' enriched functions',8,fontweight='bold',ha='center')

 # F: all map transforms and samples inherited; coordinate ratio never distorted.
 headers=['Distribution','Co-expression']+program
 colw=28.4
 for j,title in enumerate(headers):text(19+j*colw,169,title,7.5,ha='center',fontstyle='italic' if j>=2 else 'normal')
 spatial_manifest=json.loads((ROOT/'source_data/spatial_source_manifest.json').read_text('utf8'))
 for mode in range(1,5):
  data=np.load(ROOT/f'source_data/spatial_Mode{mode}.npz');pos=data['xy'];expr=data['expression'][:,[list(data['genes']).index(g)for g in program]];annot=data['labels']==target
  positive=(expr>0).all(1);co=np.zeros(len(pos));co[positive]=np.exp(np.log(expr[positive]).mean(1))
  if annot.sum():center=np.median(pos[annot],axis=0)
  elif positive.sum():center=np.median(pos[positive],axis=0)
  else:center=np.median(pos,axis=0)
  full=[pos[:,0].min(),pos[:,0].max(),pos[:,1].min(),pos[:,1].max()];side=.13*(full[1]-full[0]);zoom=[center[0]-side/2,center[0]+side/2,center[1]-side/2,center[1]+side/2]
  y=174+(mode-1)*20.5;text(5,y+9,f'Mode {mode}',7.5,rotation=90,ha='center',va='center')
  for col in range(5):
   vals=co if col==1 else expr[:,col-2] if col>=2 else annot.astype(float)
   lim=np.quantile(vals[vals>0],.99) if (vals>0).any() else 1;scaled=np.clip(vals/max(lim,1e-12),0,1)
   for iszoom,rx,rw in [(False,8+col*colw,20),(True,29+col*colw,7)]:
    a=ax(rx,y+1.5,rw,16 if not iszoom else 7,f'Panel_F_Mode{mode}_{col}_{"zoom" if iszoom else "full"}')
    a.scatter(pos[:,0],pos[:,1],s=.15 if not iszoom else .6,c='#E9E9E9',lw=0,rasterized=True)
    if col==0:
     a.scatter(pos[annot,0],pos[annot,1],s=.5 if not iszoom else 1.7,c='red',lw=0,rasterized=True)
    else:
     use=vals>0;a.scatter(pos[use,0],pos[use,1],s=.28 if not iszoom else 1,c=scaled[use],cmap='viridis',vmin=0,vmax=1,lw=0,rasterized=True)
     # Binned Gaussian-smoothed weighted density, 50/85% of peak, not cell percentages.
     hist,xe,ye=np.histogram2d(pos[use,0],pos[use,1],bins=[100,60],range=[full[:2],full[2:]],weights=vals[use]);smooth=gaussian_filter(hist,1.8)
     if smooth.max()>0:
      a.contour((xe[:-1]+xe[1:])/2,(ye[:-1]+ye[1:])/2,smooth.T,levels=[smooth.max()*.5,smooth.max()*.85],colors=['#F5C242','#F07860'],linewidths=.45)
    limits=zoom if iszoom else full;a.set(xlim=limits[:2],ylim=limits[2:],aspect='equal',xticks=[],yticks=[])
    for s in a.spines.values():s.set_visible(iszoom);s.set_linewidth(.4);s.set_color('.55')
    if not iszoom:a.add_patch(Rectangle((zoom[0],zoom[2]),side,side,fill=False,ec='.35',lw=.45,ls=(0,(2,1))))
    bar_um=200 if iszoom else 1000;length=bar_um/.5;dx=limits[1]-limits[0];dy=limits[3]-limits[2]
    bx,by=limits[0]+dx*.04,limits[2]+dy*.04
    if length<dx*.75:
     a.plot([bx,bx+length],[by,by],color='black',lw=.8)
     if col==0 and mode==1:
      if iszoom:text(rx+rw/2,y+11.5,f'{bar_um} µm',7,ha='center')
      else:a.text(bx,by+dy*.05,f'{bar_um} µm',fontsize=7,va='bottom')
  text(8,y+18.5,spatial_manifest[mode-1]['clinical_time'].replace('_SPTDI2',' SPTDI'),7,color='.35')

 # G: source fluorescence pixels only. Insets explicitly geometric, not counted.
 gcolors=['#00A000','#E03030','#B49900','#1549E8'] if lineage=='KC' else ['black']*4
 gnames=['MMP3','SOX9','KRT6B','DAPI'] if lineage=='KC' else ['SCGB2A2','S100A7','KRT77','DAPI']
 gx=155
 for j,(name,color) in enumerate(zip(gnames,gcolors)):
  text(gx,168+(j//2)*3.5,name,7.5,color=color);gx+=20 if j%2==0 else -20
 for mode in range(1,5):
  rec=next(r for r in IF if r['lineage']==lineage and r['mode']==mode);im=Image.open(rec['path']).convert('RGB');w,h=im.size;roi=rec['roi_native_xyxy'];y=177+(mode-1)*20.5
  a=ax(155,y,38,15,'Panel_G_IF_Mode'+str(mode));a.imshow(im,interpolation='none');hidden(a)
  a.add_patch(Rectangle((roi[0],roi[1]),roi[2]-roi[0],roi[3]-roi[1],fill=False,ec='white',lw=.45))
  a.text(.02,.95,f'Mode {mode}',transform=a.transAxes,color='white',fontsize=7,ha='left',va='top')
  a=ax(194,y,12,12,'Panel_G_IF_zoom_Mode'+str(mode));a.imshow(Image.open(rec['zoom_path']),interpolation='none');hidden(a)
 text(155,258,'Original scale bars retained',7,color='.4')

 # H: no surrogate gate, no false absence/zero for unmeasured proteins.
 text(12,267,'Spatial proteomics - marker coverage incomplete',8,fontweight='bold')
 spx=np.load(ROOT/'source_data/proteomics_available_features.npz');pos=spx['xy'];columns=['Distribution','Co-expression']+program
 for j,name in enumerate(columns):
  x=10+j*39.6;text(x+17,273,name,7.5,ha='center',fontstyle='italic' if j>=2 else 'normal')
  a=ax(x,276,34,14,f'Panel_H_{j}');hidden(a)
  if j<2:
   a.set_facecolor('#F7F7F7');a.text(.5,.5,'Not evaluated\nmissing marker(s)',transform=a.transAxes,ha='center',va='center',fontsize=7,color='.4')
  elif name not in spx.files:
   a.set_facecolor('#F7F7F7');a.text(.5,.5,'Not in the\nprotein library',transform=a.transAxes,ha='center',va='center',fontsize=7,color='.4')
  else:
   vals=spx[name];hi=max(np.quantile(vals,.99),1e-8);a.scatter(pos[:,0],pos[:,1],c=np.clip(vals/hi,0,1),cmap='viridis',vmin=0,vmax=1,s=1,lw=0,rasterized=True);a.set_aspect('equal')
 text(12,294,'Available protein maps are not a three-marker-positive cell-state assignment.',7,color='.4')
 stem=f'{figure_name}_{lineage}_Figure2_A-H_matched_Arial_A4_v01_REVIEW'
 save(fig,stem)
 manifests.append(dict(figure=figure_name,lineage=lineage,target=target,canvas_mm=[210,297],source_object='data/external/local_c75f45c8d1a1',source_annotation='sub_labels preserved',cells=len(labels),panels=panel_boxes,H_status='incomplete source coverage; no triple-marker gate computed',G_scale_status='original bars preserved; inset physical scale unavailable',SAC_channel_status='black names pending confirmation' if lineage=='SAC' else 'supplied channel labels retained',output_stem=stem))
 print(stem,'rendered',flush=True)
fontpath=ROOT/'QC/font_audit_pre_AI.csv'
fonts=pd.DataFrame(fontrows)
if fontpath.exists():
 old=pd.read_csv(fontpath);fonts=pd.concat([old[~old.figure.isin(fonts.figure.unique())],fonts],ignore_index=True)
fonts.to_csv(fontpath,index=False,encoding='utf-8-sig')
manifestpath=ROOT/'source_data/figure_manifests.json'
if manifestpath.exists():
 old=json.loads(manifestpath.read_text('utf8'));manifests=[m for m in old if m['figure'] not in [n['figure'] for n in manifests]]+manifests
manifestpath.write_text(json.dumps(manifests,ensure_ascii=False,indent=2),'utf8')
