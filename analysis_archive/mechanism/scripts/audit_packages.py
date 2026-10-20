from pathlib import Path
import json,hashlib,time
import yaml,h5py,numpy as np,pandas as pd

OUT=Path(__file__).resolve().parents[1]
BASE=Path('data/external/local_d12ca8436e91')
for p in ['data','h5ad','analysis','tables','figures','reports','logs','scripts','visual_gallery_all']:(OUT/p).mkdir(parents=True,exist_ok=True)
packages={'human_scrna':'人单细胞测序分析','human_spatial':'人时空组学数据分析','rat_adult_scrna':'大鼠单细胞测序分析','rat_adult_spatial':'大鼠时空组学数据分析','rat_adult_scatac':'大鼠scATAC-seq分析','bulk_klf4_oe':'KLF4-OE bulk-RNA-seq分析','bulk_igfl1_treatment':'IGFL1处理 bulk-RNA-seq分析','skin_spatial_proteomics':'空间蛋白组分析'}
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
    return h.hexdigest()
def decode(v):return np.array([x.decode()if isinstance(x,bytes)else str(x)for x in v])
def col(g,key):
    o=g[key]
    if isinstance(o,h5py.Group):
        c=o['codes'][:];cats=decode(o['categories'][:]);return np.array([cats[i]if i>=0 else ''for i in c])
    return decode(o[:])
records=[];objects=[]
for did,name in packages.items():
    root=BASE/name;dataset=yaml.safe_load((root/'dataset.yaml').read_text('utf8'))
    manifest=pd.read_csv(root/'00_manifest/ASSET_MANIFEST_SHA256.tsv',sep='\t')
    cohort=root/'00_manifest/cohort_metadata.tsv'
    rec=dict(dataset_id=did,root=str(root),dataset_yaml=dataset,manifest_rows=len(manifest),cohort_metadata_present=cohort.exists())
    if cohort.exists():
        cohortdf=pd.read_csv(cohort,sep='\t');rec['cohort_rows']=len(cohortdf)
        cohortdf.to_csv(OUT/f'tables/{did}_cohort_metadata_SOURCE.tsv',sep='\t',index=False)
    else:rec['metadata_note']='Package lacks 00_manifest/cohort_metadata.tsv; use registered canonical sample sheet/metadata, with provenance visible. Do not infer missing identifiers.'
    canonical=manifest.loc[manifest.relative_path.astype(str).str.replace('\\','/',regex=False).str.startswith('01_canonical/')] if 'relative_path' in manifest else manifest.iloc[:0]
    # Reference-only protein package has external absolute paths: its manifest is kept as the authority and audited when its object is used.
    for row in canonical.itertuples():
        p=root/str(row.relative_path).replace('\\','/')
        check=dict(dataset_id=did,file=str(p),expected_bytes=int(row.bytes),expected_sha256=row.sha256)
        check['exists']=p.exists();assert check['exists'],str(p)
        assert p.stat().st_size==int(row.bytes),('Size mismatch',p)
        check['sha256']=sha(p);check['hash_match']=check['sha256']==row.sha256
        assert check['hash_match'],('SHA256 mismatch',p)
        if p.suffix=='.h5ad':
            with h5py.File(p,'r')as f:
                x=f['X'];shape=[int(n) for n in (x.attrs['shape'] if isinstance(x,h5py.Group) else x.shape)]
                check.update(shape=shape,X_encoding=str(x.attrs.get('encoding-type','dense')),layers=list(f.get('layers',{})),raw_present='raw'in f,obs_keys=list(f['obs']),var_keys=list(f['var']))
                check['obs_uniques']={k:int(len(np.unique(col(f['obs'],k))))for k in ['sample_id','sample_batch_new','sample','patient','group','time_point','bulk_labels','sub_labels']if k in f['obs']}
                for k in ['sample_id','sample_batch_new','sample']:
                    if k in f['obs']:
                        labels=col(f['obs'],k);check[k+'_counts']=pd.Series(labels).value_counts().to_dict()
                check['var_names_preview']=col(f['var'],f['var'].attrs.get('_index','_index'))[:8].tolist()
        objects.append(check);print(did,p.name,'HASH PASS',check.get('shape',''),flush=True)
    rec['canonical_files_checked']=len(canonical)
    records.append(rec)
    (OUT/'reports/asset_audit.json').write_text(json.dumps(dict(packages=records,objects=objects),ensure_ascii=False,indent=2),'utf8')
(OUT/'analysis/roots.json').write_text(json.dumps(dict(output=str(OUT),packages={k:str(BASE/v)for k,v in packages.items()}),ensure_ascii=False,indent=2),'utf8')
pd.DataFrame([{k:v for k,v in o.items()if not isinstance(v,(list,dict))}for o in objects]).to_csv(OUT/'tables/asset_integrity_audit.tsv',sep='\t',index=False)
print('ALL PACKAGE CANONICAL HASH CHECKS COMPLETE',flush=True)
