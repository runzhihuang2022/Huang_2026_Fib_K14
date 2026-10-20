from pathlib import Path
import json,pandas as pd,hashlib
O=Path(__file__).resolve().parents[1]
roots=json.loads((O/'analysis/roots.json').read_text('utf8'))['packages']
d=pd.read_excel(Path(roots['bulk_igfl1_treatment'])/'01_canonical/IGFL1_raw_counts.xlsx',sheet_name='counts')
assert d.shape==(62700,7) and not d.gene_id.duplicated().any()
d[['gene_id','gene_name']].to_csv(O/'tables/gene_annotation_source.tsv',sep='\t',index=False)
d.set_index('gene_id')[['PBS1','PBS2','PBS3','IGFL-1','IGFL-2']].to_csv(O/'data/IGFL1_counts.csv')
k=pd.read_csv(Path(roots['bulk_klf4_oe'])/'01_canonical/raw_counts_matrix.csv')
assert set(k.gene_id)==set(d.gene_id)
k.to_csv(O/'data/KLF4_counts.csv',index=False)
print('Bulk integer counts and shared gene annotation prepared',flush=True)
