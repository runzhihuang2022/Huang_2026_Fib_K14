from pathlib import Path
import json,hashlib,pandas as pd
O=Path(__file__).resolve().parents[1]
cache=Path('data/external/local_955fb33db097')
rows=[];lines=[]
for nm in ['GO_Biological_Process_2023','KEGG_2021_Human']:
 p=cache/(nm+'.gmt.json');d=json.loads(p.read_text('utf8'))
 rows.append(dict(library=nm,source=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),n_sets=len(d),note='Existing version-locked Enrichr library; no assumption of current release'))
 for term,gs in d.items():lines.append('\t'.join([nm+'::'+term,'local_version_locked']+list(dict.fromkeys(gs))))
(O/'data/GO_KEGG_version_locked.gmt').write_text('\n'.join(lines)+'\n','utf8')
pd.DataFrame(rows).to_csv(O/'tables/pathway_gene_set_provenance.csv',index=False)
