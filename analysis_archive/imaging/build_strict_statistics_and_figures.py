from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT=Path(r"data/external/local_de73c97db013")
OUT=ROOT/"IF_colocalized_cell_statistics_strict_compartment_20260825"
OLD=ROOT/"IF_colocalized_cell_statistics_all_slides_20260825/04_per_slide/per_slide.tsv"
PY=Path(r"data/external/local_06653b588c43")
HELPER=Path(r"data/external/local_62cfde47fcb0")

for family,ylabel in [("KC","SOX9+ MMP3+ KRT6B+ cell ratio (%)"),("Fib_K14","PDGFRA+ KRT14+ cell ratio (%)")]:
    subprocess.run([str(PY),str(HELPER),"--input",str(OUT/f"04_per_slide/{family}_inference_input.tsv"),
                    "--output-dir",str(OUT/f"05_statistics/{family}"),"--figure-dir",str(OUT/f"06_figures/{family}"),"--ylabel",ylabel],check=True)

groups=[]; tests=[]
for family in ["KC","Fib_K14"]:
    g=pd.read_csv(OUT/f"05_statistics/{family}/group_summary.tsv",sep="\t"); g.insert(0,"marker_family",family); groups.append(g)
    t=pd.read_csv(OUT/f"05_statistics/{family}/statistical_tests.tsv",sep="\t"); t.insert(0,"marker_family",family); tests.append(t)
pd.concat(groups,ignore_index=True).to_csv(OUT/"05_statistics/group_summary.tsv",sep="\t",index=False)
pd.concat(tests,ignore_index=True).to_csv(OUT/"05_statistics/statistical_tests.tsv",sep="\t",index=False)

new=pd.read_csv(OUT/"04_per_slide/per_slide.tsv",sep="\t",dtype={"sample_id":str})
old=pd.read_csv(OLD,sep="\t",dtype={"sample_id":str})[["family","sample_id","numerator_n","denominator_n","ratio_pct"]]
comparison=new.merge(old,on=["family","sample_id"],suffixes=("_strict","_old_centroid_disk"),how="left")
comparison["ratio_change_pct_points"]=comparison.ratio_pct_strict-comparison.ratio_pct_old_centroid_disk
comparison.to_csv(OUT/"05_sensitivity/old_vs_strict_method_comparison.tsv",sep="\t",index=False)

plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"ps.fonttype":42})
fig,axes=plt.subplots(1,2,figsize=(7.0,2.8)); colors={"SPTDI":"#2E9BC2","DPTDI":"#E76552"}; order=["Mode 1","Mode 2","Mode 3"]
for ax,(family,ylabel) in zip(axes,[("KC","SOX9+ MMP3+ KRT6B+\ncell ratio (%)"),("Fib_K14","PDGFRA+ KRT14+\ncell ratio (%)")]):
    include=new["inference_included"].astype(str).str.lower().eq("true")
    d=new[(new.family==family)&include].copy(); rng=np.random.default_rng(20260825)
    for group in ["SPTDI","DPTDI"]:
        means=[]
        for i,mode in enumerate(order):
            y=d[(d.group==group)&(d["mode"]==mode)].ratio_pct.to_numpy();
            ax.scatter(i+rng.uniform(-.045,.045,len(y)),y,s=20,facecolors="white",edgecolors=colors[group],linewidths=.9,zorder=3)
            means.append(np.mean(y) if len(y) else np.nan)
        ax.plot(range(3),means,color=colors[group],marker="o",ms=3.5,lw=1.2,label=group)
    ax.set_xticks(range(3),["Mode\n1","Mode\n2","Mode\n3"]); ax.set_ylabel(ylabel); ax.spines[["top","right"]].set_visible(False); ax.legend(frameon=False,fontsize=7)
fig.tight_layout()
for ext in ["png","pdf","svg"]: fig.savefig(OUT/f"06_figures/KC_Fib_strict_compartment.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight")
plt.close(fig)

(OUT/"05_statistics/analysis_notes.json").write_text(json.dumps({
    "unit_of_analysis":"independent whole slide/subject; cells are nested",
    "thresholds":"pooled robust Otsu per marker/compartment with group labels unused; controls absent",
    "status":"EXPLORATORY",
    "limitations":["no no-primary/isotype/known positive-negative controls","no independent pan-membrane channel","MMP3 is secreted/pericellular","007 and 32-1 group mapping requires confirmation","CUT_007 excluded from inference"]
},ensure_ascii=False,indent=2),encoding="utf-8")
print(comparison[["family","sample_id","ratio_pct_strict","ratio_pct_old_centroid_disk","ratio_change_pct_points"]].to_string(index=False))
