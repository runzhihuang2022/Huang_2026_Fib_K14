#!/usr/bin/env python3
import json, os, gc
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import anndata as ad
import scanpy as sc

SEED = 20260811
np.random.seed(SEED)
SRC = Path('data/external/local_3658084a88f1')
ROOT = Path('data/external/local_6f6697f200d0')
OUT = ROOT/'working'/'stage1_rna'
OUT.mkdir(parents=True, exist_ok=True)

def log(msg):
    print(msg, flush=True)

def bool_col(s):
    return s.astype(str).str.lower().isin(['true','1','yes','doublet'])

def sample_balanced_indices(obs, label_col, cap, always=None):
    keep=[]
    always = np.zeros(len(obs), dtype=bool) if always is None else np.asarray(always)
    keep.extend(np.flatnonzero(always).tolist())
    used=set(keep)
    for (lab,sample), idx in obs.groupby([label_col,'sample_id'], observed=True).indices.items():
        idx=np.asarray(idx, dtype=int)
        idx=np.asarray([i for i in idx if i not in used], dtype=int)
        if len(idx)==0: continue
        n=min(len(idx), max(20, int(np.ceil(cap/max(1, obs[label_col].nunique())))))
        if len(idx)>n: idx=np.random.default_rng(SEED+len(keep)).choice(idx,n,replace=False)
        keep.extend(idx.tolist())
    # fill each label to cap across samples
    cur=pd.Series(obs.iloc[keep][label_col]).value_counts()
    used=set(keep)
    for lab, idx in obs.groupby(label_col, observed=True).indices.items():
        need=max(0, cap-int(cur.get(lab,0)))
        pool=np.asarray([i for i in idx if i not in used],dtype=int)
        if need and len(pool):
            add=np.random.default_rng(SEED+len(keep)).choice(pool,min(need,len(pool)),replace=False)
            keep.extend(add.tolist())
            used.update(add.tolist())
    return np.asarray(sorted(set(keep)),dtype=int)

log('READ_H5AD')
a=ad.read_h5ad(SRC)
o=a.obs.copy()
o['sample_id']=o['sample_id'].astype(str)
o['donor_id']=o['sample_id']
o['comparison_group']=o['Grade'].astype(str)
o['bulk_labels']=o['bulk_labels.x'].astype(str)
o['sub_labels_source']=o['sub_labels'].astype(str)

qc=(o['nFeature_RNA.x'].between(200,7500) & o['nCount_RNA.x'].between(500,100000) &
    (o['mt_percent.x']<20) & (o['HB_percent.x']<10))
if 'is_doublet' in o: qc &= ~bool_col(o['is_doublet'])
o['analysis_qc_pass']=qc.to_numpy()

genes=a.var_names.astype(str)
requested=['Klf4','Tp63','Krt5','Krt14','Epcam','Tacstd2','Vim','Pdgfra','Col1a1']
present=[g for g in requested if g in genes]
absent=[g for g in requested if g not in genes]
epi=[g for g in ['Tp63','Krt5','Krt14','Epcam','Tacstd2'] if g in genes]
mes=[g for g in ['Vim','Pdgfra','Col1a1'] if g in genes]
cols={g:genes.get_loc(g) for g in present}
X=a.X
def anypos(gs):
    y=X[:,[cols[g] for g in gs]]
    return np.asarray((y>0).sum(axis=1)).ravel()>0
klf=np.asarray(X[:,cols['Klf4']].toarray()).ravel()>0 if sp.issparse(X) else X[:,cols['Klf4']]>0
strict=(o['analysis_qc_pass'].to_numpy() & (o['bulk_labels'].to_numpy()=='Fibroblast') & klf & anypos(epi) & anypos(mes))
o['Fib_K14_strict_gate']=strict
o['sub_labels_reannotated']=o['sub_labels_source'].copy()
o.loc[(o['sub_labels_source']=='Fib_K14') & ~strict,'sub_labels_reannotated']='Fib_K14_unresolved'
o.loc[strict,'sub_labels_reannotated']='Fib_K14'
a.obs=o

pd.DataFrame({'marker':requested,'rat_symbol':requested,'present':[x in present for x in requested],
              'action':['used' if x in present else 'skipped_absent' for x in requested]}).to_csv(OUT/'marker_availability.csv',index=False)
o.groupby(['sample_id','comparison_group'],observed=True).agg(n_cells=('sample_id','size'),n_qc=('analysis_qc_pass','sum'),n_strict_FibK14=('Fib_K14_strict_gate','sum')).reset_index().to_csv(OUT/'rna_qc_by_sample.csv',index=False)
o.loc[o.analysis_qc_pass,'bulk_labels'].value_counts().rename_axis('bulk_labels').reset_index(name='n_cells').to_csv(OUT/'rna_qc_bulk_counts.csv',index=False)
o.loc[o.analysis_qc_pass,'sub_labels_reannotated'].value_counts().rename_axis('sub_labels').reset_index(name='n_cells').to_csv(OUT/'rna_qc_subtype_counts.csv',index=False)

log(f'QC_PASS {int(qc.sum())}/{len(qc)} STRICT_FIBK14 {int(strict.sum())}')
aq=a[a.obs.analysis_qc_pass].copy()

# A sample-balanced global atlas is used for newly computed primary and Harmony-sensitivity geometry.
# All QC-passing cells remain in the canonical object and are projected onto this reference later.
idx=sample_balanced_indices(aq.obs,'sub_labels_reannotated',cap=2500,always=aq.obs.Fib_K14_strict_gate.to_numpy())
ref=aq[idx].copy()
log(f'BALANCED_REFERENCE {ref.n_obs}')
ref.layers['counts']=ref.X.copy()
sc.pp.normalize_total(ref,target_sum=1e4)
sc.pp.log1p(ref)
sc.pp.highly_variable_genes(ref,n_top_genes=3000,flavor='seurat',batch_key=None)
rhv=ref[:,ref.var.highly_variable].copy(); sc.pp.scale(rhv,max_value=10)
sc.tl.pca(rhv,n_comps=50,svd_solver='arpack',random_state=SEED)
sc.pp.neighbors(rhv,n_neighbors=30,n_pcs=40,use_rep='X_pca',random_state=SEED,key_added='primary')
sc.tl.umap(rhv,neighbors_key='primary',random_state=SEED)
sc.tl.leiden(rhv,neighbors_key='primary',resolution=0.8,key_added='leiden_primary',random_state=SEED)
ref.obsm['X_pca']=rhv.obsm['X_pca']; ref.obsm['X_umap_primary']=rhv.obsm['X_umap'].copy()
ref.obs['leiden_primary']=rhv.obs['leiden_primary'].astype(str).to_numpy()

try:
    sc.external.pp.harmony_integrate(rhv,key='sample_id',basis='X_pca',adjusted_basis='X_pca_harmony',max_iter_harmony=20)
    sc.pp.neighbors(rhv,n_neighbors=30,use_rep='X_pca_harmony',random_state=SEED,key_added='harmony')
    sc.tl.umap(rhv,neighbors_key='harmony',random_state=SEED)
    sc.tl.leiden(rhv,neighbors_key='harmony',resolution=0.8,key_added='leiden_harmony_sensitivity',random_state=SEED)
    ref.obsm['X_pca_harmony']=rhv.obsm['X_pca_harmony']; ref.obsm['X_umap_harmony_sensitivity']=rhv.obsm['X_umap'].copy()
    ref.obs['leiden_harmony_sensitivity']=rhv.obs['leiden_harmony_sensitivity'].astype(str).to_numpy()
    harmony='PASS'
except Exception as e:
    harmony=f'FAIL: {type(e).__name__}: {e}'
    log(harmony)
ref.obsm['X_umap']=ref.obsm['X_umap_primary'].copy()

# independent fibroblast geometry; never filter the global coordinates for Figure 5
fib=aq[aq.obs.bulk_labels=='Fibroblast'].copy()
fib.layers['counts']=fib.X.copy()
sc.pp.normalize_total(fib,target_sum=1e4)
sc.pp.log1p(fib)
sc.pp.highly_variable_genes(fib,n_top_genes=3000,flavor='seurat',batch_key=None)
fhv=fib[:,fib.var.highly_variable].copy(); sc.pp.scale(fhv,max_value=10)
sc.tl.pca(fhv,n_comps=50,svd_solver='arpack',random_state=SEED)
sc.pp.neighbors(fhv,n_neighbors=30,n_pcs=40,use_rep='X_pca',random_state=SEED,key_added='primary')
sc.tl.umap(fhv,neighbors_key='primary',random_state=SEED)
fib.obsm['X_pca']=fhv.obsm['X_pca']; fib.obsm['X_umap_primary']=fhv.obsm['X_umap'].copy()
for res in [0.3,0.5,0.8,1.0]:
    sc.tl.leiden(fhv,neighbors_key='primary',resolution=res,key_added=f'leiden_r{res}',random_state=SEED)
    fib.obs[f'leiden_r{res}']=fhv.obs[f'leiden_r{res}'].astype(str).to_numpy()
try:
    sc.external.pp.harmony_integrate(fhv,key='sample_id',basis='X_pca',adjusted_basis='X_pca_harmony',max_iter_harmony=20)
    sc.pp.neighbors(fhv,n_neighbors=30,use_rep='X_pca_harmony',random_state=SEED,key_added='harmony')
    sc.tl.umap(fhv,neighbors_key='harmony',random_state=SEED)
    fib.obsm['X_pca_harmony']=fhv.obsm['X_pca_harmony']; fib.obsm['X_umap_harmony_sensitivity']=fhv.obsm['X_umap'].copy()
except Exception as e: log(f'FIB_HARMONY_FAIL {e}')
fib.obsm['X_umap']=fib.obsm['X_umap_primary'].copy()

log('WRITE_OBJECTS')
ref.write_h5ad(OUT/'rna_balanced_reference.h5ad',compression='gzip')
fib.write_h5ad(OUT/'rna_fibroblast_recomputed.h5ad',compression='gzip')

# Seurat/ArchR-compatible sparse exports, raw counts only
for obj,name in [(ref,'integration_reference'),(fib,'fibroblast')]:
    d=OUT/name; d.mkdir(exist_ok=True)
    mmwrite(d/'counts_genes_by_cells.mtx',obj.layers['counts'].T)
    pd.DataFrame({'gene_symbol':obj.var_names}).to_csv(d/'genes.tsv',sep='\t',index=False)
    m=obj.obs.copy(); m.insert(0,'cell_barcode',obj.obs_names); m.to_csv(d/'metadata.csv',index=False)
    pd.DataFrame(obj.obsm['X_pca'],index=obj.obs_names).to_csv(d/'pca.csv.gz',compression='gzip')
    pd.DataFrame(obj.obsm['X_umap_primary'],index=obj.obs_names,columns=['UMAP_1','UMAP_2']).to_csv(d/'umap_primary.csv.gz',compression='gzip')

# SCENIC training matrix: balanced within fibroblast states, all strict Fib_K14 retained up to 5000.
fobs=fib.obs.copy()
always=fobs.Fib_K14_strict_gate.to_numpy()
if always.sum()>5000:
    pos=np.flatnonzero(always); sel=np.random.default_rng(SEED).choice(pos,5000,replace=False); always[:]=False; always[sel]=True
sidx=sample_balanced_indices(fobs,'sub_labels_reannotated',cap=1200,always=always)
sf=fib[sidx].copy()
detect=np.asarray((sf.layers['counts']>0).sum(axis=0)).ravel()
keepg=detect>=max(20,int(0.005*sf.n_obs))
sf=sf[:,keepg]
sd=OUT/'scenic'; sd.mkdir(exist_ok=True)
mmwrite(sd/'rat_counts_cells_by_genes.mtx',sf.layers['counts'])
pd.DataFrame({'cell_id':sf.obs_names}).to_csv(sd/'cells.tsv',sep='\t',index=False)
pd.DataFrame({'rat_symbol':sf.var_names}).to_csv(sd/'rat_genes.tsv',sep='\t',index=False)
sf.obs[['sample_id','comparison_group','bulk_labels','sub_labels_source','sub_labels_reannotated','Fib_K14_strict_gate']].to_csv(sd/'cell_metadata.csv')

summary={'seed':SEED,'source':str(SRC),'n_input':a.n_obs,'n_qc_pass':aq.n_obs,'n_reference':ref.n_obs,
         'n_fibroblast':fib.n_obs,'n_strict_FibK14':int(fib.obs.Fib_K14_strict_gate.sum()),
         'n_scenic_cells':sf.n_obs,'n_scenic_genes':sf.n_vars,'markers_present':present,'markers_absent_skipped':absent,
         'primary_geometry':'uncorrected PCA/UMAP','harmony_sensitivity':harmony,
         'batch_governance':'sample_id equals biological unit; Harmony is sensitivity only'}
(OUT/'stage1_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(OUT/'STAGE1_DONE').write_text('PASS\n')
log(json.dumps(summary,indent=2))
