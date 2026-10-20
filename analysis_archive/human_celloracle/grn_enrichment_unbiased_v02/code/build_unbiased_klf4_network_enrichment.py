#!/usr/bin/env python3
from pathlib import Path
import argparse, json, re
import numpy as np
import pandas as pd
import gseapy as gp

EPITHELIAL_PATTERNS = [
    r"epitheli", r"epiderm", r"keratinocyte", r"keratinization",
    r"cornif", r"skin development", r"epithelial cell migration",
    r"epithelial tube", r"cell-cell junction", r"wound healing",
]

def read_lines(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]

def is_epithelial_term(term: str) -> bool:
    return any(re.search(p, term, flags=re.I) for p in EPITHELIAL_PATTERNS)

def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--background", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args(); args.output.mkdir(parents=True,exist_ok=True)

    edges=pd.read_csv(args.edges,sep="\t")
    sig=edges.loc[edges.FDR < 0.05].copy().sort_values(
        ["evidence_score","coef_abs"],ascending=False
    )
    sig.insert(0,"unbiased_rank",np.arange(1,len(sig)+1))
    sig.to_csv(args.output/"Fib_K14_KLF4_all_FDR05_network_edges.tsv",sep="\t",index=False)
    positive=sig.loc[sig.direction.eq("positive")].copy()
    negative=sig.loc[sig.direction.eq("negative")].copy()
    positive.target.drop_duplicates().to_csv(
        args.output/"Fib_K14_KLF4_all_FDR05_positive_targets.txt",index=False,header=False
    )
    negative.target.drop_duplicates().to_csv(
        args.output/"Fib_K14_KLF4_all_FDR05_negative_targets.txt",index=False,header=False
    )

    background=read_lines(args.background)
    libraries=["GO_Biological_Process_2025","KEGG_2021_Human"]
    enr=gp.enrichr(
        gene_list=positive.target.drop_duplicates().astype(str).tolist(),
        gene_sets=libraries, organism="human", background=background,
        outdir=None, cutoff=1.0, verbose=True,
    )
    ora=enr.results.copy()
    ora["query_genes_n"]=positive.target.nunique()
    ora["background_genes_n"]=len(background)
    ora["minus_log10_FDR"]=-np.log10(ora["Adjusted P-value"].clip(lower=1e-300))
    ora["epithelialization_related"]=ora.Term.map(is_epithelial_term)
    ora=ora.sort_values(["Adjusted P-value","P-value"])
    ora.to_csv(args.output/"Fib_K14_KLF4_unbiased_positive_ORA_all.tsv",sep="\t",index=False)
    ora.loc[ora["Adjusted P-value"] < 0.05].to_csv(
        args.output/"Fib_K14_KLF4_unbiased_positive_ORA_FDR05.tsv",sep="\t",index=False
    )
    epi=ora.loc[ora.epithelialization_related].copy()
    epi.to_csv(args.output/"Fib_K14_KLF4_epithelialization_ORA_all.tsv",sep="\t",index=False)
    epi_sig=epi.loc[epi["Adjusted P-value"] < 0.05].copy()
    epi_sig.to_csv(args.output/"Fib_K14_KLF4_epithelialization_ORA_FDR05.tsv",sep="\t",index=False)

    # The epithelialization-focused target map is derived objectively from the
    # ten best keyword-matched terms. It is kept even when no term passes FDR,
    # with the evidence state retained explicitly.
    epi_focus=epi.head(10).copy()
    epi_focus.to_csv(
        args.output/"Fib_K14_KLF4_epithelialization_ORA_top10.tsv",
        sep="\t",index=False
    )
    rows=[]
    for row in epi_focus.itertuples(index=False):
        genes=[g for g in str(row.Genes).split(";") if g]
        for gene in genes:
            rows.append({
                "target":gene, "Term":row.Term, "Gene_set":row.Gene_set,
                "Adjusted_P_value":np.nan,
                "evidence_state":"",
            })
    mapping=pd.DataFrame(rows)
    if len(mapping):
        # Restore exact adjusted p-values without depending on itertuples field renaming.
        pmap=epi_focus.set_index(["Gene_set","Term"])["Adjusted P-value"].to_dict()
        mapping["Adjusted_P_value"]=[
            pmap[(gs,term)] for gs,term in zip(mapping.Gene_set,mapping.Term)
        ]
        mapping["evidence_state"]=np.where(
            mapping.Adjusted_P_value < 0.05,
            "FDR<0.05",
            "epithelialization-related trend; FDR>=0.05",
        )
        mapping=mapping.merge(
            positive[["target","coef_mean","coef_abs","FDR","evidence_score","strict_filtered_link"]],
            on="target",how="left"
        ).sort_values(["Adjusted_P_value","evidence_score"],ascending=[True,False])
    mapping.to_csv(args.output/"Fib_K14_KLF4_epithelialization_target_term_map.tsv",sep="\t",index=False)

    checkpoint={
        "status":"PASS",
        "species":"human",
        "cluster":"Fib_K14",
        "TF":"KLF4",
        "selection_rule":"all KLF4 edges with within-TF BH FDR < 0.05; no manual gene additions",
        "network_edges":int(len(sig)),
        "positive_targets":int(positive.target.nunique()),
        "negative_targets":int(negative.target.nunique()),
        "strict_CellOracle_targets":int(sig.strict_filtered_link.astype(bool).sum()),
        "ora_query":"all FDR<0.05 positive KLF4 targets",
        "ora_query_n":int(positive.target.nunique()),
        "background_n":len(background),
        "libraries":libraries,
        "ORA_FDR05_terms":int((ora["Adjusted P-value"]<0.05).sum()),
        "epithelialization_keyword_rule":EPITHELIAL_PATTERNS,
        "epithelialization_FDR05_terms":int(len(epi_sig)),
        "epithelialization_focus_rule":"union of genes in ten lowest-FDR keyword-matched epithelialization terms",
        "epithelialization_focus_targets":int(mapping.target.nunique()) if len(mapping) else 0,
        "statistical_scope":"gene-set ORA; not donor-replicated clinical inference",
    }
    (args.output/"unbiased_network_enrichment_checkpoint.json").write_text(
        json.dumps(checkpoint,indent=2),encoding="utf-8"
    )
    print(json.dumps(checkpoint,indent=2))

if __name__=="__main__":
    main()
