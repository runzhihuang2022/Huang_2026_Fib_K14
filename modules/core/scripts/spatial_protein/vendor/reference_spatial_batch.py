"""Read-only human Stereo-seq five-panel plotting; no clustering or inference.

python spatial_batch.py --config config.json --cell-type Fib_K14 \
    --genes PDGFRA KRT14 KRT5 --logic BASE_AND_ANY --output results/run01
python spatial_batch.py --config config.json --batch jobs.csv --output results/batch01
"""
from __future__ import annotations
import argparse, csv, hashlib, html, json, platform, re, shutil, sys, traceback
from datetime import datetime
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from PIL import Image

HERE=Path(__file__).resolve().parent
EXPECTED_SHAPE=(387061,54236)
CANONICAL='01_canonical/stereoseq/stereoseq_11samples_compatible.h5ad'
def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
def dump(path,obj):
    def numpy_scalar(value):
        if isinstance(value,np.generic):return value.item()
        if isinstance(value,np.ndarray):return value.tolist()
        raise TypeError(type(value).__name__)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False,default=numpy_scalar),encoding='utf-8')
def read_col(node):
    if isinstance(node,h5py.Dataset):
        return node.asstr()[:] if h5py.check_string_dtype(node.dtype) else node[:]
    if 'categories' in node and 'codes' in node:
        cats=read_col(node['categories']);codes=node['codes'][:]
        a=np.empty(len(codes),dtype=object);a[:]=None;a[codes>=0]=cats[codes[codes>=0]];return a
    raise ValueError(f'Unsupported metadata encoding: {node.name}')
def rotated(xy,angle):
    t=np.deg2rad(angle);r=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    return np.asarray(xy,dtype=float)@r.T
def full_limits(xy,crop_factor=1):
    lo=xy.min(0);hi=xy.max(0);c=(lo+hi)/2;span=(hi-lo)/crop_factor
    pad=np.maximum(span,1)*.03
    return [float(c[0]-span[0]/2-pad[0]),float(c[0]+span[0]/2+pad[0]),float(c[1]-span[1]/2-pad[1]),float(c[1]+span[1]/2+pad[1])]
def gate(v,logic,threshold=0):
    if logic=='ALL':return (v>threshold).all(1)
    if logic=='BASE_AND_ANY':return (v[:,0]>threshold)&(v[:,1:]>threshold).any(1)
    raise ValueError('logic must be ALL or BASE_AND_ANY')
def coexpression(v,logic):
    if logic=='ALL':return np.cbrt(v[:,0]*v[:,1]*v[:,2])
    if logic=='BASE_AND_ANY':return np.sqrt(v[:,0]*np.maximum(v[:,1],v[:,2]))
    raise ValueError(logic)
def safe_name(s):
    name=re.sub(r'[^A-Za-z0-9_.-]+','_',s).strip('._')
    if not name:raise ValueError('cell_type must contain filename-safe characters')
    return name
def load_jobs(args):
    if args.batch:
        with open(args.batch,encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
        jobs=[{'cell_type':r['cell_type'],'genes':[r['marker1'],r['marker2'],r['marker3']],
               'logic':r['logic'],'color':r.get('color') or '#FF0000'} for r in rows]
    else:
        if not args.cell_type or not args.genes:raise ValueError('Use --batch or --cell-type and exactly three --genes')
        jobs=[{'cell_type':args.cell_type,'genes':args.genes,'logic':args.logic,'color':args.color}]
    if not jobs:raise ValueError('No plotting jobs')
    stems=[]
    for j in jobs:
        if len(j['genes'])!=3 or len(set(j['genes']))!=3:raise ValueError('Exactly three distinct, exact gene symbols required')
        gate(np.ones((1,3)),j['logic']);matplotlib.colors.to_rgba(j['color']);stems.append(safe_name(j['cell_type']))
    if len(stems)!=len(set(stems)):raise ValueError('Duplicate output cell_type names; use separate runs for alternate gene sets')
    return jobs
def audit_route(cfg,out):
    root=Path(cfg['route_root']).resolve()
    with open(root/'00_manifest/ASSET_MANIFEST_SHA256.tsv',encoding='utf-8-sig') as f: manifest={r['relative_path']:r for r in csv.DictReader(f,delimiter='\t')}
    evidence=[]
    for rel in ['dataset.yaml','00_manifest/cohort_metadata.tsv',CANONICAL]:
        p=root/rel
        if rel not in manifest:raise ValueError(f'Unregistered input: {rel}')
        actual=sha256(p);expected=manifest[rel]['sha256']
        if actual.lower()!=expected.lower():raise ValueError(f'SHA-256 mismatch: {p}')
        evidence.append({'path':str(p),'bytes':p.stat().st_size,'sha256':actual})
    # Preserve the registry declaration, not just a path guessed from file names.
    metadata=pd.read_csv(root/'00_manifest/cohort_metadata.tsv',sep='\t',keep_default_na=False)
    stereo=metadata[metadata.Modality=='Stereo-seq'].copy()
    if len(stereo)!=11 or stereo['Sample ID'].duplicated().any():raise ValueError('Stereo-seq whitelist is not 11 unique sections')
    stereo['Cells / bins / spots']=pd.to_numeric(stereo['Cells / bins / spots'])
    if stereo['Cells / bins / spots'].sum()!=EXPECTED_SHAPE[0]:raise ValueError('Cohort count mismatch')
    if out:
        dump(out/'input_SHA256.json',evidence)
        shutil.copy2(root/'dataset.yaml',out/'routed_dataset.yaml')
        metadata.to_csv(out/'cohort_metadata.tsv',sep='\t',index=False)
    return root/CANONICAL,stereo
def read_input(path,genes):
    """Read CSR in row chunks once, densifying only the requested gene columns."""
    with h5py.File(path,'r') as f:
        shape=tuple(f['X'].attrs['shape'])
        if shape!=EXPECTED_SHAPE:raise ValueError(f'Unexpected canonical shape: {shape}')
        obs=f['obs'];index_key=obs.attrs.get('_index','_index')
        var=f['var'];vkey=var.attrs.get('_index','_index');names=read_col(var[vkey]).astype(str)
        if len(np.unique(names))!=len(names):raise ValueError('Duplicate var_names; explicit mapping required')
        lookup={g:i for i,g in enumerate(names)}
        missing=sorted(set(genes)-set(lookup))
        if missing:raise ValueError(f'Missing exact marker symbols (not replaced with zeros): {missing}')
        ids=np.array([lookup[g] for g in genes]);matrix=f['X']
        if matrix.attrs.get('encoding-type')!='csr_matrix':raise ValueError('This verified reader expects canonical CSR X')
        indptr=matrix['indptr'][:];values=np.empty((shape[0],len(ids)),dtype=np.float32)
        for start in range(0,shape[0],8192):
            end=min(shape[0],start+8192);lo,hi=int(indptr[start]),int(indptr[end])
            block=csr_matrix((matrix['data'][lo:hi],matrix['indices'][lo:hi],indptr[start:end+1]-lo),shape=(end-start,shape[1]))
            values[start:end]=block[:,ids].toarray()
        if not np.isfinite(values).all() or (values<0).any():raise ValueError('Selected expression is negative or non-finite')
        samples=read_col(obs['sample_batch_new'])
        if any(x is None for x in samples):raise ValueError('Missing sample_batch_new')
        barcodes=read_col(obs[index_key]).astype(str)
        if len(np.unique(barcodes))!=len(barcodes):raise ValueError('Duplicate bin identifiers')
        xy=np.asarray(f['obsm/spatial'],dtype=float)
        if not np.isfinite(xy).all():raise ValueError('Invalid spatial coordinates')
        annotations={k:read_col(obs[k]) for k in ['sub_labels','bulk_labels'] if k in obs}
        audit={'shape':list(shape),'expression_source':'X, existing processed nonnegative expression; not raw counts',
               'has_log1p_metadata':'log1p' in f['uns'],'layers':list(f['layers']),
               'selected_gene_range':[float(values.min()),float(values.max())],
               'selected_gene_noninteger_fraction':float(np.mean(values[values>0]%1!=0)) if (values>0).any() else 0,
               'annotation_columns':list(annotations)}
    return dict(xy=xy,samples=np.asarray(samples,str),barcodes=barcodes,values=values,genes=genes,annotations=annotations,audit=audit)
def label_vector(cfg,data,jobs):
    mode=cfg['distribution_mode']
    if mode=='marker_positive':return None
    if mode!='annotation':raise ValueError('distribution_mode must be annotation or marker_positive')
    col=cfg.get('annotation_column','sub_labels');external=cfg.get('annotation_csv')
    if external:
        p=Path(external);expected=cfg.get('annotation_sha256')
        if not expected or sha256(p)!=expected:raise ValueError('External annotation CSV requires its exact SHA-256')
        df=pd.read_csv(p,dtype=str,keep_default_na=False)
        if not {'bin_id',col}.issubset(df):raise ValueError('Annotation CSV needs bin_id and annotation column')
        if df.bin_id.duplicated().any() or set(df.bin_id)!=set(data['barcodes']):raise ValueError('Annotation bin_id must match canonical bins one-to-one')
        labels=df.set_index('bin_id').loc[data['barcodes'],col].to_numpy()
    elif col in data['annotations']:labels=np.asarray(data['annotations'][col],str)
    else:raise ValueError(f'{col} is absent. Supply audited annotation CSV; no marker-based fallback is performed.')
    if pd.isna(labels).any() or any(str(x).strip()=='' for x in labels):raise ValueError('Missing annotations')
    for j in jobs:
        if j['cell_type'] not in labels:raise ValueError(f"Exact subtype absent: {j['cell_type']}")
    return labels
def density(xy,v,limits,cfg,seed):
    """Original expression-weighted gaussian_kde, Scott bandwidth divided by 1.5.
    Contours are 50% / 85% of peak density, NOT probability mass or cell fraction.
    """
    x0,x1,y0,y1=limits;sel=(v>0)&(xy[:,0]>=x0)&(xy[:,0]<=x1)&(xy[:,1]>=y0)&(xy[:,1]<=y1)
    ids=np.flatnonzero(sel);info={'positive_bins':int(len(ids)),'status':'not_estimable','approximate':False}
    if len(ids)<6:info['reason']='fewer than six positive bins';return None,info
    cap=int(cfg.get('kde_max_points',0))
    if cap and len(ids)>cap:
        ids=np.random.default_rng(seed).choice(ids,cap,replace=False);info['approximate']=True
    if np.linalg.matrix_rank(np.cov(xy[ids].T))<2:info['reason']='singular spatial covariance';return None,info
    try:
        kde=gaussian_kde(xy[ids].T,weights=v[ids]);kde.set_bandwidth(kde.factor/1.5)
        gx,gy=np.meshgrid(np.linspace(x0,x1,100),np.linspace(y0,y1,100))
        z=kde(np.vstack([gx.ravel(),gy.ravel()])).reshape(gx.shape)
    except (ValueError,np.linalg.LinAlgError) as ex:info['reason']=str(ex);return None,info
    info.update(status='ok',used_bins=int(len(ids)),peak=float(z.max()))
    return (gx,gy,z),info
def style():
    font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
    plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
def anatomy(ax,section):
    for item in section.get('annotation_paths',[]):
        xy=np.asarray(item['xy']);ax.plot(xy[:,0],xy[:,1],color='black',lw=.45,ls=(0,(2,2)) if item.get('dashed') else '-',zorder=9)
    for item in section.get('annotation_labels',[]):
        ax.annotate(item['text'],xy=item['xy'],xytext=item['text_axes'],textcoords='axes fraction',fontsize=7,
                    ha='center',va='center',annotation_clip=False,arrowprops=dict(arrowstyle='-',lw=.4,color='black'),zorder=12)
def scatter_panel(ax,item,row,limits,zoom=False,cfg=None):
    xy=item['xy'];x0,x1,y0,y1=limits
    keep=(xy[:,0]>=x0)&(xy[:,0]<=x1)&(xy[:,1]>=y0)&(xy[:,1]<=y1)
    if row==0:
        target=item['target']
        for m,c,s,z in [(keep&~target,'#E0E0E0',.025 if not zoom else .16,1),(keep&target,item['color'],.14 if not zoom else .6,4)]:
            ax.scatter(xy[m,0],xy[m,1],c=c,s=s,edgecolors='none',rasterized=True,zorder=z)
    else:
        v=item['features'][row-1];ids=np.flatnonzero(keep);ids=ids[np.argsort(v[ids],kind='stable')]
        ax.scatter(xy[ids,0],xy[ids,1],c=v[ids],s=.025 if not zoom else .22,cmap='viridis',vmin=0,vmax=1,edgecolors='none',rasterized=True)
        dd=item['density'][(row,zoom)]
        if dd is not None:
            gx,gy,z=dd;ax.contour(gx,gy,z,levels=[z.max()*.5,z.max()*.85],colors=['#DAA520','#FF4500'],linewidths=.45,zorder=6)
    ax.set(xlim=(x0,x1),ylim=(y0,y1),aspect='equal',xticks=[],yticks=[])
    for spine in ax.spines.values():spine.set_visible(zoom);spine.set_linewidth(.4)
    if not zoom:
        r=item['roi'];ax.add_patch(Rectangle((r[0],r[2]),r[1]-r[0],r[3]-r[2],fill=False,lw=.5,ls='--',zorder=10))
        if row==0:anatomy(ax,item['section'])
    if row==0 and not item['target'].any():ax.text(.5,.5,'No target bins',transform=ax.transAxes,fontsize=7,ha='center')
def save_figure(fig,base,formats,dpi):
    fig.canvas.draw()
    sizes=[t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text()]
    if min(sizes)<7 or max(sizes)>12:raise ValueError(f'Font QA failed: {min(sizes)}, {max(sizes)}')
    base.parent.mkdir(exist_ok=True,parents=True)
    for ext in formats:
        opts={'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
        fig.savefig(base.with_suffix('.'+ext),dpi=dpi,bbox_inches=None,facecolor='white',**opts)
        if ext in ['png','tiff']:
            with Image.open(base.with_suffix('.'+ext)) as im:rgb=im.convert('RGB')
            rgb.save(base.with_suffix('.'+ext),dpi=(dpi,dpi),**({'compression':'tiff_lzw'} if ext=='tiff' else {}))
    plt.close(fig)
    return {'file_stem':base.name,'width_mm':float(fig.get_figwidth()*25.4),'height_mm':float(fig.get_figheight()*25.4),'font_range_pt':[min(sizes),max(sizes)],'dpi':dpi}
def render(items,job,rows,out,cfg,basename):
    n=len(items);W=210;H=17+len(rows)*27+6
    fig=plt.figure(figsize=(W/25.4,H/25.4));left=25;right=184;gap=2;cw=(right-left-(n-1)*gap)/n
    labels=['Distribution','Co-expression']+job['genes'];overview_axes=[];zoom_axes=[]
    for c,item in enumerate(items):
        x=left+c*(cw+gap)
        fig.text((x+cw/2)/W,1-5/H,item['section']['display_label'],ha='center',va='center',fontsize=10,fontweight='bold')
        fig.text((x+cw/2)/W,1-10/H,item['section'].get('timepoint',''),ha='center',va='center',fontsize=7)
        for rr,row in enumerate(rows):
            top=18+rr*27
            if c==0:fig.text(.008,1-(top+12)/H,labels[row],ha='left',va='center',fontsize=8 if row<2 else 9,fontstyle='italic' if row>=2 else 'normal',fontweight='bold' if row<2 else 'normal')
            # 3:1 overview / zoom; equal aspect, no stretch or tight cropping.
            ax=fig.add_axes([x/W,1-(top+23)/H,cw*.73/W,20/H]);az=fig.add_axes([(x+cw*.77)/W,1-(top+23)/H,cw*.23/W,20/H])
            scatter_panel(ax,item,row,item['limits'],cfg=cfg);scatter_panel(az,item,row,item['roi'],zoom=True,cfg=cfg)
            if rr==0:overview_axes.append(ax);zoom_axes.append(az)
    if 0 in rows:
        fig.text(.891,.91,'Target bins',fontsize=7,fontweight='bold')
        fig.legend(handles=[Line2D([],[],marker='o',ls='',color=job['color'],markersize=3,label='Positive' if cfg['distribution_mode']=='marker_positive' else 'Annotated'),Line2D([],[],marker='o',ls='',color='#DADADA',markersize=3,label='Other bins')],loc='upper left',bbox_to_anchor=(.876,.90),frameon=False,fontsize=7,handletextpad=.3)
    if len(rows)>1:
        fig.text(.89,.74,'KDE peak',fontsize=7,fontweight='bold')
        fig.legend(handles=[Line2D([],[],color='#DAA520',lw=.6,label='50%'),Line2D([],[],color='#FF4500',lw=.6,label='85%')],loc='upper left',bbox_to_anchor=(.876,.73),frameon=False,fontsize=7)
        fig.text(.89,.56,'Expression',fontsize=7,fontweight='bold')
        cax=fig.add_axes([.899,.36,.012,.17]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=cax);cb.set_ticks([0,1],labels=['Low','High']);cb.ax.tick_params(labelsize=7,length=0)
        fig.text(.891,.29,'Scale (µm)',fontsize=7,fontweight='bold');fig.text(.891,.26,'Full / zoom',fontsize=7)
        fig.canvas.draw()
        for i,(ax,az) in enumerate(zip(overview_axes,zoom_axes)):
            yy=.215-i*.04;fig.text(.891,yy,f'M{i+1}',fontsize=7)
            for a,xx,value in [(ax,.935,1000),(az,.975,200)]:
                physical=a.get_position().width*value/cfg['um_per_coordinate']/(a.get_xlim()[1]-a.get_xlim()[0])
                fig.add_artist(Line2D([xx,xx+physical],[yy+.004]*2,transform=fig.transFigure,color='black',lw=.8))
        fig.text(.891,.027,'1000 / 200',fontsize=7)
    else:
        if rows[0]>0:
            fig.text(.891,.94,'Expression',fontsize=7,fontweight='bold')
            cax=fig.add_axes([.90,.70,.012,.17]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap='viridis'),cax=cax);cb.set_ticks([0,1],labels=['Low','High']);cb.ax.tick_params(labelsize=7,length=0)
            fig.legend(handles=[Line2D([],[],color='#DAA520',lw=.6,label='50% peak'),Line2D([],[],color='#FF4500',lw=.6,label='85% peak')],loc='upper left',bbox_to_anchor=(.876,.65),frameon=False,fontsize=7)
        fig.text(.891,.38,'1000 / 200 µm',fontsize=7);fig.canvas.draw()
        for i,(ax,az) in enumerate(zip(overview_axes,zoom_axes)):
            yy=.29-i*.075;fig.text(.891,yy,f'M{i+1}',fontsize=7)
            for a,xx,value in [(ax,.935,1000),(az,.975,200)]:
                physical=a.get_position().width*value/cfg['um_per_coordinate']/(a.get_xlim()[1]-a.get_xlim()[0])
                fig.add_artist(Line2D([xx,xx+physical],[yy+.012]*2,transform=fig.transFigure,color='black',lw=.8))
    qa=save_figure(fig,out/basename,cfg['formats'],cfg['dpi']);return qa
def run(args):
    config_path=Path(args.config).resolve();cfg=json.loads(config_path.read_text(encoding='utf-8'))
    if cfg['normalization'] not in ['global_quantile','section_quantile']:raise ValueError('Unknown normalization')
    if cfg.get('positive_threshold',0)<0 or cfg['um_per_coordinate']<=0:raise ValueError('Invalid threshold or scale')
    if cfg['dpi']<600:raise ValueError('Final outputs require at least 600 dpi')
    if not {'pdf','png','tiff'}.issubset(cfg['formats']):raise ValueError('PDF, PNG and TIFF are required')
    if cfg.get('kde_max_points',0)<0:raise ValueError('kde_max_points must be zero or a positive integer')
    jobs=load_jobs(args);out=Path(args.output or HERE/'results'/datetime.now().strftime('%Y%m%d_%H%M%S')).resolve()
    route=Path(cfg['route_root']).resolve()
    if out==route or route in out.parents:raise ValueError('Outputs must be outside the read-only data package')
    out.mkdir(parents=True,exist_ok=False)
    try:
        dump(out/'run_config.json',cfg);dump(out/'jobs.json',jobs)
        shutil.copy2(__file__,out/'spatial_batch.py')
        data_path,cohort=audit_route(cfg,out)
        genes=list(dict.fromkeys(g for j in jobs for g in j['genes']));data=read_input(data_path,genes)
        counts=pd.Series(data['samples']).value_counts()
        if set(counts.index)!=set(cohort['Sample ID']):raise ValueError('Sample whitelist mismatch')
        if any(counts[r['Sample ID']]!=r['Cells / bins / spots'] for _,r in cohort.iterrows()):raise ValueError('Per-section bin counts mismatch')
        dump(out/'expression_audit.json',data['audit']);cohort.to_csv(out/'section_counts.tsv',sep='\t',index=False)
        sections=cfg['sections']
        if not 1<=len(sections)<=4:raise ValueError('Use 1–4 sections per plate; partition the 11 sections into multiple configs')
        if len({s['sample_id'] for s in sections})!=len(sections):raise ValueError('Duplicate sections')
        template=Path(cfg['template_json']);template=template if template.is_absolute() else config_path.parent/template
        if sha256(template)!=cfg['template_sha256']:raise ValueError('Frozen template SHA-256 mismatch')
        frozen=json.loads(template.read_text(encoding='utf-8'));shutil.copy2(template,out/'frozen_template.json')
        for sec in sections:
            sid=sec['sample_id']
            if sid not in counts:raise ValueError(f'Section outside routed whitelist: {sid}')
            actual=cohort.set_index('Sample ID').loc[sid,'Time point']
            if sec['timepoint']!=actual:raise ValueError(f'Wrong timepoint for {sid}')
            if sid not in frozen['sections']:raise ValueError(f'No frozen ROI/annotation geometry for {sid}; provide a verified section template, never borrow another slice')
            sec.update(frozen['sections'][sid]);sec['timepoint']=actual
        lab=label_vector(cfg,data,jobs)
        if args.audit_only:dump(out/'status.json',{'status':'AUDIT_PASS','analysis_started':False});return out
        style();allqa=[];summary=[];density_records=[];failures=[]
        for job in jobs:
            try:
                dest=out/safe_name(job['cell_type']);dest.mkdir();items=[]
                cols=[genes.index(g) for g in job['genes']]
                pooled=data['values'][np.isin(data['samples'],[s['sample_id'] for s in sections])][:,cols]
                global_q=[float(np.quantile(v[v>0],.995)) if (v>0).any() else 1 for v in pooled.T]
                for si,sec in enumerate(sections):
                    m=data['samples']==sec['sample_id'];xy=rotated(data['xy'][m],sec['angle']);v=data['values'][m][:,cols]
                    if not np.allclose(full_limits(xy,sec.get('crop_factor',1)),sec['limits'],atol=1e-4):raise ValueError('Template geometry mismatch: input coordinates or orientation changed')
                    q=global_q if cfg['normalization']=='global_quantile' else [float(np.quantile(a[a>0],.995)) if (a>0).any() else 1 for a in v.T]
                    scaled=np.clip(v/np.asarray(q),0,1);target=gate(v,job['logic'],cfg.get('positive_threshold',0)) if lab is None else lab[m]==job['cell_type']
                    co=coexpression(scaled,job['logic']);features=[co]+[scaled[:,i] for i in range(3)]
                    item=dict(xy=xy,target=target,features=features,color=job['color'],section=sec,limits=sec['limits'],roi=sec['roi'],density={})
                    for row in range(1,5):
                        for zoom,lim in [(False,item['limits']),(True,item['roi'])]:
                            dd,meta=density(xy,features[row-1],lim,cfg,1729+si*20+row*2+zoom)
                            item['density'][(row,zoom)]=dd;density_records.append(dict(cell_type=job['cell_type'],sample_id=sec['sample_id'],row=row,zoom=zoom,**meta))
                    items.append(item)
                    summary.append(dict(cell_type=job['cell_type'],sample_id=sec['sample_id'],display_label=sec['display_label'],timepoint=sec['timepoint'],patient=cohort.set_index('Sample ID').loc[sec['sample_id'],'Patient ID'],bins=len(v),target_bins=int(target.sum()),definition=cfg['distribution_mode'],logic=job['logic'],genes=';'.join(job['genes']),q995=';'.join(map(str,q))))
                    np.savez_compressed(dest/f'section_{si+1}_plot_data.npz',xy=xy,target=target,coexpression=co,scaled_genes=scaled,bin_id=data['barcodes'][m],genes=np.asarray(job['genes']))
                dump(dest/'frozen_ROI_annotations.json',sections)
                allqa.append(render(items,job,list(range(5)),dest,cfg,'Composite_5rows_210mm'))
                for row,name in enumerate(['Distribution','Coexpression']+job['genes']):allqa.append(render(items,job,[row],dest,cfg,f'{row+1:02}_{safe_name(name)}'))
                dump(dest/'STATUS.json',{'status':'PASS','distribution_definition':cfg['distribution_mode'],'coexpression_is_not_lineage_tracing':True})
            except Exception as ex:
                failures.append({'job':job['cell_type'],'error':str(ex),'traceback':traceback.format_exc()})
        pd.DataFrame(summary).to_csv(out/'plot_summary.tsv',sep='\t',index=False)
        dump(out/'density_audit.json',density_records);dump(out/'figure_QA.json',allqa);dump(out/'failures.json',failures)
        dump(out/'environment.json',{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'h5py':h5py.__version__})
        cards=[]
        for job in jobs:
            d=safe_name(job['cell_type']);png=f'{d}/Composite_5rows_210mm.png'
            if not (out/png).exists():continue
            links=[]
            for i,name in enumerate(['Distribution','Coexpression']+job['genes']):
                stem=f'{d}/{i+1:02}_{safe_name(name)}'
                links.append('<p>'+html.escape(name)+' — '+' | '.join(f'<a href="{stem}.{ext}">{ext.upper()}</a>' for ext in ['pdf','png','tiff'])+'</p>')
            cards.append(f'<h2>{html.escape(job["cell_type"])}</h2><p>{html.escape(str(job))}</p><a href="{d}/Composite_5rows_210mm.pdf">PDF</a> <a href="{d}/Composite_5rows_210mm.tiff">TIFF</a><img src="{png}" style="max-width:1400px;width:100%">'+''.join(links))
        (out/'OPEN_ME_FIRST.html').write_text('<meta charset="utf-8"><h1>Five-panel Stereo-seq plots</h1><p>Distribution: '+html.escape(cfg['distribution_mode'])+'. X is processed expression. KDE contours are fractions of peak density, not cell percentages. Marker-positive bins are not independently validated subtype identities. Mode display follows the original morphological panel mapping, not chronological ordering. No IF images are generated.</p>'+''.join(cards),encoding='utf-8')
        dump(out/'status.json',{'status':'PARTIAL_FAILURE' if failures else 'PASS','jobs_requested':len(jobs),'jobs_failed':len(failures),'new_statistical_inference':False})
        if failures:raise RuntimeError(f'{len(failures)} job(s) failed; see failures.json. Successful jobs preserved.')
    except Exception as ex:
        dump(out/'RUN_ERROR.json',{'error':str(ex),'traceback':traceback.format_exc()});raise
    finally:
        entries=[f'{p.relative_to(out).as_posix()}\t{sha256(p)}' for p in sorted(out.rglob('*')) if p.is_file() and p.name!='SHA256.tsv']
        (out/'SHA256.tsv').write_text('path\tsha256\n'+'\n'.join(entries),encoding='utf-8')
    return out
def parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--config',default=str(HERE/'config.json'));p.add_argument('--cell-type');p.add_argument('--genes',nargs=3)
    p.add_argument('--logic',choices=['ALL','BASE_AND_ANY'],default='BASE_AND_ANY');p.add_argument('--color',default='#FF0000')
    p.add_argument('--batch');p.add_argument('--output');p.add_argument('--audit-only',action='store_true');return p
if __name__=='__main__':
    args=parser().parse_args();print(run(args),flush=True)
