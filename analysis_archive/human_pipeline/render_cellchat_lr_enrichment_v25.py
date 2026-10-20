from __future__ import annotations
import argparse, json, re
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MM=1/25.4
def safe(x): return re.sub(r"[^A-Za-z0-9_.-]+","_",str(x)).strip("_")

def main():
    q=argparse.ArgumentParser(); q.add_argument("--tables",type=Path,required=True); q.add_argument("--output",type=Path,required=True); a=q.parse_args()
    mpl.rcParams.update({"font.family":"Arial","font.size":6,"axes.titlesize":8,"axes.labelsize":7,"xtick.labelsize":6,"ytick.labelsize":6,"pdf.fonttype":42,"ps.fonttype":42})
    a.output.mkdir(parents=True,exist_ok=True); tab=a.output.parent.parent/"05_tables"/"cellchat_lr_enrichment"; tab.mkdir(parents=True,exist_ok=True)
    bgfile=a.tables/"CELLCHAT_SIGNALING_BACKGROUND_GENES.tsv"; background=pd.read_csv(bgfile,sep="\t").iloc[:,0].dropna().astype(str).unique().tolist()
    try:
        import gseapy as gp
        sets=gp.get_library(name="GO_Biological_Process_2023",organism="human")
    except Exception as e:
        gp=None; sets=None; library_error=str(e)
    records=[]; all_results=[]
    for fp in sorted(a.tables.glob("*_all_LR.tsv")):
        d=pd.read_csv(fp,sep="\t"); stem=fp.name.removesuffix("_all_LR.tsv")
        genes=sorted(set(d.get("ligand",pd.Series(dtype=str)).dropna().astype(str))|set(d.get("receptor",pd.Series(dtype=str)).dropna().astype(str)))
        pd.Series(genes,name="LR_gene_symbol").to_csv(tab/f"{stem}_all_LR_genes.tsv",sep="\t",index=False)
        res=pd.DataFrame(); status="NOT_ESTIMABLE"; reason="fewer than 3 LR genes"
        if gp is None: reason=library_error
        elif len(genes)>=3:
            try:
                res=gp.enrich(gene_list=genes,gene_sets=sets,background=background,outdir=None,cutoff=1.0).results
                status="PASS"; reason=""
            except Exception as e: reason=str(e)
        if not res.empty:
            res.insert(0,"analysis",stem); res.to_csv(tab/f"{stem}_GO_BP.tsv",sep="\t",index=False); all_results.append(res)
        fig,ax=plt.subplots(figsize=(52.5*MM,52.5*MM))
        if not res.empty:
            x=res.copy(); x["Adjusted P-value"]=pd.to_numeric(x["Adjusted P-value"],errors="coerce"); x=x.sort_values("Adjusted P-value").head(10).iloc[::-1]
            x["score"]=-np.log10(x["Adjusted P-value"].clip(lower=1e-300)); terms=x["Term"].astype(str).str.replace(r" \(GO:\d+\)$","",regex=True).str.slice(0,42)
            ax.barh(range(len(x)),x.score,color="#2C7FB8",height=.72); ax.set_yticks(range(len(x)),terms); ax.set_xlabel(r"$-\log_{10}$(adjusted P)")
        else:
            ax.text(.5,.5,"Enrichment not estimable\n"+reason[:70],ha="center",va="center",transform=ax.transAxes,fontsize=7); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(stem.replace("_"," ")+"\nLR-pair GO enrichment",fontweight="bold",pad=2)
        ax.spines[["top","right"]].set_visible(False); fig.subplots_adjust(left=.42,right=.96,bottom=.17,top=.84)
        out=a.output/f"{stem}_all_LR_GO_BP_barplot"
        fig.savefig(out.with_suffix(".pdf")); fig.savefig(out.with_suffix(".png"),dpi=600); plt.close(fig)
        records.append({"analysis":stem,"status":status,"reason":reason,"n_LR_pairs":len(d),"n_LR_genes":len(genes),"n_terms":len(res)})
    if all_results: pd.concat(all_results,ignore_index=True).to_csv(tab/"ALL_CellChat_LR_GO_BP_enrichment.tsv.gz",sep="\t",index=False,compression="gzip")
    qa={"status":"PASS_WITH_EXPLICIT_NOT_ESTIMABLE","engine":"GO_Biological_Process_2023 over-representation analysis","background":"CellChat signaling database genes","background_n":len(background),"square_mm":[52.5,52.5],"font":"Arial","font_pt_range":[6,8],"records":records}
    (a.output.parent.parent/"10_qa"/"CELLCHAT_LR_ENRICHMENT_AUDIT.json").write_text(json.dumps(qa,indent=2),encoding="utf-8")
if __name__=="__main__": main()
