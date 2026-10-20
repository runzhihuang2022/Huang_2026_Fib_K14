from pathlib import Path
import json,hashlib,shutil,pandas as pd,numpy as np
from scipy.stats import spearmanr,fisher_exact
O=Path(__file__).resolve().parents[1];roots=json.loads((O/'analysis/roots.json').read_text('utf8'))['packages'];rows=[]
# Integrity check all reused processed ATAC evidence against its package manifest.
root=Path(roots['rat_adult_scatac']);manifest=pd.read_csv(root/'00_manifest/ASSET_MANIFEST_SHA256.tsv',sep='\t')
chosen=manifest[manifest.relative_path.str.contains('16_FibK14_KLF4_regulatory_network') & manifest.relative_path.str.endswith(('.csv','.md'))]
for row in chosen.itertuples():
 p=root/row.relative_path.replace('\\','/');s=hashlib.sha256(p.read_bytes()).hexdigest();assert s==row.sha256 and p.stat().st_size==row.bytes
 rows.append(dict(role='processed_ATAC_evidence',path=str(p),sha256=s,hash_match=True))
 if p.suffix=='.csv':shutil.copy2(p,O/'tables'/('ATAC_ARCHIVED_'+p.name))
# Previously processed, mixed-species SCENIC is an ancillary hypothesis audit, not fresh human-only validation.
strict=Path('data/external/local_da1799c2973e')
for nm in ['IGFL1_KLF4_regulatory_evidence_summary.csv','KLF4_IGFL1_strict_orthology_audit.csv','KLF4_human_cervid_IGFL1_SCENIC_motif_pruned_targets.csv']:
 p=strict/nm;rows.append(dict(role='ancillary_ARCHIVED_not_independent_current_cohort',path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),hash_match='no_registered_manifest_reference_fingerprint_only'));shutil.copy2(p,O/'tables'/('ANCILLARY_'+nm))
animal=Path('data/external/local_84fd7bfad3ee')
for nm in ['healing_per_slice.csv','healing_summary_mean_SD_SEM.csv','healing_Welch_Holm_EXPLORATORY.csv']:
 p=animal/'source_data'/nm;rows.append(dict(role='mouse_HE_REVIEW_slice_not_animal_inference',path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),hash_match='fingerprint_only'));shutil.copy2(p,O/'tables'/('MOUSE_REVIEW_'+nm))
pd.DataFrame(rows).to_csv(O/'tables/reused_processed_evidence_provenance.csv',index=False)
# Cross-perturbation effect concordance: independent samples but unconfirmed cell identity / batch.
k=pd.read_csv(O/'tables/KLF4_DESeq2_all.csv');i=pd.read_csv(O/'tables/IGFL1_DESeq2_all.csv');e1=pd.read_csv(O/'tables/KLF4_edgeR_QL_sensitivity.csv');e2=pd.read_csv(O/'tables/IGFL1_edgeR_QL_sensitivity.csv')
m=k.merge(i,on=['gene_id','gene_name'],suffixes=('_KLF4','_IGFL1'))
m=m.merge(e1[['gene_id','logFC','FDR']].rename(columns={'logFC':'edgeR_logFC_KLF4','FDR':'edgeR_FDR_KLF4'}),on='gene_id').merge(e2[['gene_id','logFC','FDR']].rename(columns={'logFC':'edgeR_logFC_IGFL1','FDR':'edgeR_FDR_IGFL1'}),on='gene_id')
m.to_csv(O/'tables/cross_perturbation_all_genes.csv',index=False)
axis=json.loads((O/'analysis/defined_modules.json').read_text('utf8'))['axis'];m[m.gene_name.isin(axis)].to_csv(O/'tables/cross_perturbation_axis_DESeq2_edgeR.csv',index=False)
z=m.dropna(subset=['log2FoldChange_KLF4','log2FoldChange_IGFL1']);r,p=spearmanr(z.log2FoldChange_KLF4,z.log2FoldChange_IGFL1)
summary=dict(n_shared_tested_genes=len(z),effect_spearman_rho=float(r),effect_spearman_p_descriptive=float(p),same_direction_both_FDR05=int(((z.padj_KLF4<.05)&(z.padj_IGFL1<.05)&(z.log2FoldChange_KLF4*z.log2FoldChange_IGFL1>0)).sum()),opposite_direction_both_FDR05=int(((z.padj_KLF4<.05)&(z.padj_IGFL1<.05)&(z.log2FoldChange_KLF4*z.log2FoldChange_IGFL1<0)).sum()),note='Genes are not independent biological replicates; concordance is descriptive, not mediation or pathway direction')
(O/'reports/cross_perturbation_summary.json').write_text(json.dumps(summary,indent=2),'utf8')
print(summary,flush=True)
print(m[m.gene_name.isin(['KLF4','KRT14','KRT17','PDGFRA','EGFR','HBEGF','TGFA','DUSP6','MKI67'])][['gene_name','log2FoldChange_KLF4','padj_KLF4','edgeR_FDR_KLF4','log2FoldChange_IGFL1','padj_IGFL1','edgeR_FDR_IGFL1']].to_string(index=False),flush=True)
