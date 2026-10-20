#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import loompy

ROOT=Path('data/external/local_6f6697f200d0')
SRC=ROOT/'working'/'stage1_rna'/'scenic'
OUT=ROOT/'working'/'stage3_scenic'; OUT.mkdir(parents=True,exist_ok=True)
HOMO=Path('data/external/local_2f021543f21c')

X=mmread(SRC/'rat_counts_cells_by_genes.mtx').tocsr()
cells=pd.read_csv(SRC/'cells.tsv',sep='\t').cell_id.astype(str).to_numpy()
rat=pd.read_csv(SRC/'rat_genes.tsv',sep='\t').rat_symbol.astype(str).to_numpy()
meta=pd.read_csv(SRC/'cell_metadata.csv',index_col=0).reindex(cells)
h=pd.read_csv(HOMO,sep='\t',header=None,names=['hid','taxid','gene_id','symbol','protein_gi','protein_acc'],dtype=str)
r=h[h.taxid=='10116'][['hid','symbol']].rename(columns={'symbol':'rat_symbol'})
m=h[h.taxid=='10090'][['hid','symbol']].rename(columns={'symbol':'mouse_symbol'})
mp=r.merge(m,on='hid')
# Conservative one-to-one HomoloGene groups only.
nr=mp.groupby('hid').rat_symbol.nunique(); nm=mp.groupby('hid').mouse_symbol.nunique()
good=set(nr[(nr==1)&(nm==1)].index)
mp=mp[mp.hid.isin(good)].drop_duplicates('rat_symbol')
mapdict=dict(zip(mp.rat_symbol,mp.mouse_symbol))
mapped=np.array([mapdict.get(g,'') for g in rat])
keep=mapped!=''
X=X[:,keep]; rat2=rat[keep]; mouse=mapped[keep]

# Sum any duplicate mouse symbols after conservative mapping.
uniq, inv=np.unique(mouse,return_inverse=True)
A=sp.csr_matrix((np.ones(len(inv),dtype=np.float32),(np.arange(len(inv)),inv)),shape=(len(inv),len(uniq)))
Y=(X@A).tocsr()
det=np.asarray((Y>0).sum(axis=0)).ravel()
sel=det>=max(20,int(0.005*Y.shape[0]))
Y=Y[:,sel]; uniq=uniq[sel]

loom=OUT/'fibroblast_mouse_ortholog_counts.loom'
loompy.create(str(loom),Y.T.astype(np.float32),{'Gene':uniq},{'CellID':cells})
pd.DataFrame({'rat_symbol':rat,'mouse_symbol':mapped,'mapped_one_to_one':keep}).to_csv(OUT/'rat_mouse_homologene_audit.csv',index=False)
meta.to_csv(OUT/'cell_metadata.csv')
summary={'n_cells':int(Y.shape[0]),'n_rat_genes_input':int(len(rat)),'n_one_to_one_mapped':int(keep.sum()),
         'n_mouse_genes_after_detection':int(Y.shape[1]),'mapping_source':str(HOMO),
         'mapping_policy':'NCBI HomoloGene taxid 10116 to 10090; one-rat-to-one-mouse groups only'}
(OUT/'prepare_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
