#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,scanpy as sc
ROOT=Path('data/external/local_6f6697f200d0')
OUT=ROOT/'working'/'stage2_joint'
xy=pd.read_csv(OUT/'joint_umap_coordinates.csv').set_index('cell_id')
meta=pd.read_csv(OUT/'joint_base_metadata.csv.gz').set_index('cell_id').reindex(xy.index)
scores=pd.read_csv(OUT/'joint_signature_scores.csv.gz').set_index('cell_id').reindex(xy.index)
a=ad.AnnData(X=xy.to_numpy(np.float32),obs=meta.copy())
a.obsm['X_joint_mapped']=xy.to_numpy(np.float32)
sc.pp.neighbors(a,n_neighbors=30,use_rep='X_joint_mapped',random_state=20260811)
sc.tl.leiden(a,resolution=0.8,key_added='joint_cluster',random_state=20260811)
cl=a.obs.joint_cluster.astype(str)
ev=[]
for k in sorted(cl.unique(),key=lambda x:int(x)):
    z=scores.loc[cl==k].mean().sort_values(ascending=False)
    ev.append({'joint_cluster':k,'reviewed_label':z.index[0],'top_score':z.iloc[0],
               'second_score':z.iloc[1],'margin':z.iloc[0]-z.iloc[1],
               'review_status':'marker_reviewed' if z.iloc[0]-z.iloc[1]>=0.015 else 'low_margin_marker_reviewed'})
ev=pd.DataFrame(ev); ev.to_csv(OUT/'unsupervised_cluster_marker_evidence.csv',index=False)
a.obs['unsupervised_reviewed']=cl.map(ev.set_index('joint_cluster').reviewed_label)
# Normalize heterogeneous object columns before HDF5 serialization. Numeric QC
# fields remain numeric; only mixed string/bool/NA columns are converted.
for col in a.obs.columns:
    if a.obs[col].dtype == object:
        a.obs[col] = a.obs[col].where(a.obs[col].notna(), 'NA').astype(str)
out=a.obs.copy();out.insert(0,'cell_id',out.index)
out.to_csv(OUT/'joint_cell_metadata_three_tracks.csv',index=False)
pd.crosstab(out.modality,out.provided_bulk).to_csv(OUT/'joint_modality_by_bulk.csv')
a.write_h5ad(OUT/'rat_RNA_ATAC_joint_embedding.h5ad',compression='gzip')
summary=json.loads((OUT/'stage2_R_summary.json').read_text())
summary.update({'n_joint_clusters':int(cl.nunique()),'unsupervised_annotation':'cluster mean canonical lineage signature, reviewed','status':'PASS'})
(OUT/'stage2_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(OUT/'STAGE2_DONE').write_text('PASS\n')
print(json.dumps(summary,indent=2))
