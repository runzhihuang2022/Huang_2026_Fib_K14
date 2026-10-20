from pathlib import Path
import json,h5py,numpy as np,pandas as pd,hashlib
from scipy.sparse import csr_matrix
O=Path(__file__).resolve().parents[1]
roots=json.loads((O/'analysis/roots.json').read_text('utf8'))['packages']
modules={
'basal_epithelial':['KRT5','KRT14','KRT15','TP63','ITGA6','ITGB4','COL17A1'],
'epithelial_junction':['EPCAM','CDH1','CLDN1','CLDN4','TJP1','DSG1','DSP','OCLN'],
'epithelial_differentiation':['KRT1','KRT10','IVL','FLG','LOR','GRHL3','ZNF750','OVOL1','IRF6'],
'wound_keratin':['KRT6A','KRT6B','KRT16','KRT17','LAMC2','LAMB3','TACSTD2'],
'proliferation':['MKI67','TOP2A','PCNA','MCM2','MCM5','CDK1','CCNB1','UBE2C'],
'fibroblast_matrix':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','COL6A1','COL6A2'],
'myofibroblast':['ACTA2','TAGLN','MYL9','COL11A1','POSTN','CTHRC1'],
'MAPK_response':['FOS','FOSB','JUN','JUNB','EGR1','DUSP1','DUSP4','DUSP6','SPRY2','ETV4','ETV5'],
'NFKB_response':['NFKBIA','NFKBIZ','TNFAIP3','ICAM1','CXCL1','CXCL2','CXCL8','CCL2','BIRC3','RELB'],
'JAK_STAT_response':['SOCS3','CISH','PIM1','BCL3','IRF1','STAT1','SOCS1'],
'TGFB_response':['SERPINE1','CCN2','SMAD7','PMEPA1','SKIL','TGFBI'],
'WNT_response':['AXIN2','LEF1','TCF7','NKD1','LGR5','MYC','WNT5A'],
'NOTCH_response':['HES1','HEY1','HEY2','NRARP','DTX1'],
}
axis=['KLF4','IGFL1','IGFL2','IGFL3','IGFL4','IGFLR1','TMEM149','IGF1','IGF1R','AREG','HBEGF','TGFA','EGF','EGFR','FGF7','FGF10','FGFR2','HGF','MET','TGFB1','TGFBR1','TGFBR2','CXCL12','CXCR4','TNF','TNFRSF1A','IL6','IL6R','IL6ST','PDGFA','PDGFB','PDGFRA','PDGFRB','AKT1','MTOR','PIK3CA','STAT3','BAX','BCL2','FADD','TRADD','CASP8','SOX9','MMP3','SCGB2A2','S100A7','KRT77','VIM','CDKN1A']
genes=list(dict.fromkeys(axis+sum(modules.values(),[])))
(O/'analysis/defined_modules.json').write_text(json.dumps(dict(modules=modules,axis=axis,genes=genes),indent=2),'utf8')
def dec(x):return np.array([v.decode()if isinstance(v,bytes)else str(v)for v in x])
def readcol(g,k):
 d=g[k]
 if isinstance(d,h5py.Group):
  cats=dec(d['categories'][:]);codes=d['codes'][:];return np.array([cats[c]if c>=0 else ''for c in codes])
 return dec(d[:])
def run(path,did,sample_default=None):
 print('READ',did,flush=True)
 with h5py.File(path,'r')as f:
  v=f['var'];vn=readcol(v,v.attrs.get('_index','_index'));lookup={g.upper():i for i,g in enumerate(vn)}
  indices=[lookup.get(g,-1)for g in genes];present=np.array(indices)>=0;sel=np.array(indices)[present]
  x=f['X'];shape=tuple(x.attrs['shape'])if isinstance(x,h5py.Group)else x.shape
  probe=x['data'][:min(100000,len(x['data']))]if isinstance(x,h5py.Group)else x[:100].ravel()
  countlike=bool(np.all(probe>=0)and np.mean(np.abs(probe-np.round(probe))<1e-6)>.999)
  values=np.zeros((shape[0],len(genes)),np.float32);totals=np.zeros(shape[0]);detected=np.zeros_like(values,dtype=bool)
  ptr=x['indptr'][:] if isinstance(x,h5py.Group) else None
  for start in range(0,shape[0],4000):
   end=min(start+4000,shape[0])
   if ptr is not None:
    i,j=int(ptr[start]),int(ptr[end]);m=csr_matrix((x['data'][i:j],x['indices'][i:j],ptr[start:end+1]-i),shape=(end-start,shape[1]))
   else:m=csr_matrix(x[start:end,:])
   block=m[:,sel].toarray();totals[start:end]=np.asarray(m.sum(1)).ravel();detected[start:end,present]=block>0
   if countlike:block=np.log1p(block/np.maximum(totals[start:end,None],1)*1e4)
   values[start:end,present]=block
  ob=f['obs'];meta=pd.DataFrame({'obs_id':readcol(ob,ob.attrs.get('_index','_index'))})
  for k in ['sample_id','patient','group','time_point','bulk_labels','sub_labels','sample_batch_new','sample_batch','x','y','orig.ident','Time','tissue_group','leiden_res_3.0']:
   if k in ob:meta[k]=readcol(ob,k)
  key=next((k for k in ['sample_id','sample_batch_new','sample_batch']if k in meta),None)
  meta['analysis_sample']=meta[key] if key else sample_default or did
  meta.to_csv(O/f'data/{did}_obs.csv.gz',index=False)
  np.savez_compressed(O/f'data/{did}_panel.npz',values=values,positive=detected,totals=totals,genes=genes,present=present)
  module_matrix={}
  for name,gs in modules.items():
   gi=[genes.index(g)for g in gs if present[genes.index(g)]]
   module_matrix[name]=values[:,gi].mean(1)if gi else np.full(shape[0],np.nan)
  md=pd.DataFrame(module_matrix)
  md.to_csv(O/f'data/{did}_module_scores.csv.gz',index=False)
  pd.DataFrame([dict(dataset=did,gene=g,present=bool(present[z]),source_feature=vn[indices[z]]if present[z]else '',n_positive=int(detected[:,z].sum()),n_observations=shape[0],positive_pct=100*detected[:,z].mean() if present[z]else np.nan,expression_unit='log1p(CP10K)'if countlike else 'canonical_X_no_retransform',probe_integer_fraction=float(np.mean(probe==np.round(probe))))for z,g in enumerate(genes)]).to_csv(O/f'tables/{did}_gene_availability.csv',index=False)
  for gran in ['all','bulk_labels','sub_labels']:
   keys=['analysis_sample']+([gran]if gran!='all'and gran in meta else [])
   if gran!='all'and len(keys)==1:continue
   rows=[]
   for keysval,idx in meta.groupby(keys,sort=False).groups.items():
    keysval=keysval if isinstance(keysval,tuple)else(keysval,);a=np.asarray(idx);r=dict(zip(keys,keysval));r['n_obs']=len(a)
    for k in ['patient','group','time_point']:
     if k in meta:r[k]='|'.join(meta.loc[a,k].unique())
    for z,g in enumerate(genes):r[g+'_mean']=float(values[a,z].mean())if present[z]else np.nan;r[g+'_pct']=100*float(detected[a,z].mean())if present[z]else np.nan
    for nm in modules:r[nm+'_score']=float(md.loc[a,nm].mean())
    rows.append(r)
   pd.DataFrame(rows).to_csv(O/f'tables/{did}_{gran}_sample_expression.csv',index=False)
  print('DONE',did,shape,'countlike',countlike,'axis availability',[(g,present[genes.index(g)],int(detected[:,genes.index(g)].sum()))for g in ['KLF4','IGFL1','IGFLR1']],flush=True)
jobs=[(Path(roots['human_scrna'])/'01_canonical/human_epithelialization_scrna_279305cells.h5ad','human_scrna'),(Path(roots['human_spatial'])/'01_canonical/stereoseq/stereoseq_11samples_compatible.h5ad','human_stereo'),(Path(roots['human_spatial'])/'01_canonical/visium/visium_2sections_4449spots.h5ad','human_visium')]
jobs += [(p,'rat_spatial_'+p.stem)for p in sorted((Path(roots['rat_adult_spatial'])/'01_canonical/sections').glob('*.h5ad'))]
manifest=pd.read_csv(Path(roots['skin_spatial_proteomics'])/'00_manifest/ASSET_MANIFEST_SHA256.tsv',sep='\t')
protein_audit=[]
for row in manifest[manifest.asset_id.isin(['processed_h5ad','final_visualization_h5ad'])].itertuples():
 p=Path(row.source_path);digest=hashlib.sha256(p.read_bytes()).hexdigest();assert p.stat().st_size==row.bytes and digest==row.sha256
 protein_audit.append(dict(asset_id=row.asset_id,path=str(p),sha256=digest,hash_match=True))
 if row.asset_id=='final_visualization_h5ad':jobs.append((p,'spatial_proteomics'))
pd.DataFrame(protein_audit).to_csv(O/'tables/protein_external_source_integrity.csv',index=False)
for p,did in jobs:
 if not(O/f'tables/{did}_all_sample_expression.csv').exists():run(p,did)
