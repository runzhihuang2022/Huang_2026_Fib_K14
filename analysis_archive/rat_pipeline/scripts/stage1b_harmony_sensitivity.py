#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np,anndata as ad,scanpy as sc,harmonypy as hm
ROOT=Path('data/external/local_6f6697f200d0')
OUT=ROOT/'working'/'stage1_rna'
summary={}
for fn,label in [('rna_balanced_reference.h5ad','global_reference'),('rna_fibroblast_recomputed.h5ad','fibroblast')]:
    p=OUT/fn; a=ad.read_h5ad(p)
    z=hm.run_harmony(np.asarray(a.obsm['X_pca']),a.obs,'sample_id',max_iter_harmony=20,random_state=20260811)
    corr=np.asarray(z.Z_corr)
    if corr.shape[0]!=a.n_obs and corr.shape[1]==a.n_obs: corr=corr.T
    if corr.shape[0]!=a.n_obs: raise ValueError(f'Harmony shape {corr.shape}, expected n_obs={a.n_obs}')
    a.obsm['X_pca_harmony_sensitivity']=corr
    sc.pp.neighbors(a,n_neighbors=30,use_rep='X_pca_harmony_sensitivity',random_state=20260811,key_added='harmony_sensitivity')
    sc.tl.umap(a,neighbors_key='harmony_sensitivity',random_state=20260811)
    a.obsm['X_umap_harmony_sensitivity']=a.obsm['X_umap'].copy()
    sc.tl.leiden(a,neighbors_key='harmony_sensitivity',resolution=0.8,key_added='leiden_harmony_sensitivity',random_state=20260811)
    a.obsm['X_umap']=a.obsm['X_umap_primary'].copy()
    a.write_h5ad(p,compression='gzip')
    summary[label]={'n_cells':a.n_obs,'harmony_shape':list(corr.shape),'status':'PASS'}
(OUT/'harmony_sensitivity_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(OUT/'STAGE1B_DONE').write_text('PASS\n')
print(json.dumps(summary,indent=2))
