from pathlib import Path
import json,hashlib,html,markdown
import pandas as pd,numpy as np
O=Path(__file__).resolve().parents[1]
report=O/'reports/KLF4_IGFL1_受体与下游机制_整合分析报告.md'
body=markdown.markdown(report.read_text('utf8'),extensions=['tables','fenced_code'])
plots=json.loads((O/'reports/figure_index.json').read_text('utf8'))
gallery=''.join(f'<section><h2>{html.escape(n)}</h2><a href="../figures/{n}.pdf">PDF</a> | <a href="../figures/{n}.svg">SVG</a><img src="../figures/{n}.png" alt="{n}"></section>'for n in plots)
links=''.join(f'<li><a href="../tables/{p.name}">{html.escape(p.name)}</a></li>'for p in sorted((O/'tables').glob('*')))
css='body{font-family:Arial,"Microsoft YaHei",sans-serif;max-width:1180px;margin:30px auto;padding:20px;color:#222;line-height:1.7}h1,h2{color:#216780}h2{margin-top:2em}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:8px;border:1px solid #ddd;text-align:left}th{background:#eef4f7}img{max-width:100%;height:auto}a{color:#126c96}section{padding:10px 0;border-bottom:1px solid #ddd}code{overflow-wrap:anywhere}@media print{body{font-size:10pt;margin:0}h2{break-after:avoid}img{max-height:210mm;object-fit:contain}table{font-size:8pt}section{break-before:page}}'
(O/'reports/整合分析_报告与图表.html').write_text(f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>KLF4–IGFL1 mechanism audit</title><style>{css}</style></head><body><p><a href="#figures">结果图</a> | <a href="#tables">全量分析表</a></p>{body}<h1 id="figures">结果图</h1>{gallery}<h1 id="tables">全量分析表</h1><ul>{links}</ul></body></html>','utf8')
# Quantitative delivery checks independent of plotting code.
checks=[]
for nm,n_samples in [('KLF4',4),('IGFL1',5)]:
 c=pd.read_csv(O/f'data/{nm}_counts.csv');r=pd.read_csv(O/f'tables/{nm}_DESeq2_all.csv');e=pd.read_csv(O/f'tables/{nm}_edgeR_QL_sensitivity.csv');n=pd.read_csv(O/f'tables/{nm}_normalized_counts.csv')
 assert c.shape==(62700,n_samples+1) and not c.gene_id.duplicated().any()
 assert len(r)==int((c.iloc[:,1:].sum(1)>=10).sum()) and len(r)==len(e)==len(n)
 assert r.gene_id.is_unique and all(r.padj.dropna().between(0,1))
 g=pd.read_csv(O/f'tables/{nm}_GSEA_all_GO_KEGG.csv');f=pd.read_csv(O/f'tables/{nm}_GSEA_mechanism_focus.csv')
 assert len(f)>0 and set(f.pathway)<=set(g.pathway)
 checks.append(dict(experiment=nm,n_samples=n_samples,genes_input=len(c),genes_tested=len(r),DESeq2_FDR05=int((r.padj<.05).sum()),GSEA_sets_tested=len(g),GSEA_focus_terms=len(f),status='pass'))
a=pd.read_csv(O/'tables/bulk_axis_raw_counts_including_filtered_genes.csv');assert (a[a.gene_name=='IGFL1'].total_raw_count==0).all()
assert len(pd.read_csv(O/'tables/human_scrna_all_sample_expression.csv'))==38
assert len(pd.read_csv(O/'tables/human_stereo_all_sample_expression.csv'))==11
assert len(pd.read_csv(O/'tables/human_visium_all_sample_expression.csv'))==2
pd.DataFrame(checks).to_csv(O/'reports/quantitative_delivery_QA.csv',index=False)
pd.DataFrame([dict(file=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())for folder in ['tables','figures','reports','scripts']for p in sorted((O/folder).glob('*'))if p.is_file()and p.name!='delivery_manifest_sha256.csv']).to_csv(O/'reports/delivery_manifest_sha256.csv',index=False)
print(pd.DataFrame(checks).to_string(index=False),flush=True)
