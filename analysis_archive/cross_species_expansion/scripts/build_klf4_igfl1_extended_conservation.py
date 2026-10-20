from __future__ import annotations

import json, shutil
from pathlib import Path
import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from build_fibk14_density_reference_style import sc_specs, marker_call, subset_condition
from build_panelA_annotation_human_rat_flexible_others_planarian_alias import fibk14_annotation_call

mpl.rcParams.update({"font.family":"Arial","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
 "font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":8,"ytick.labelsize":8,
 "pdf.fonttype":42,"ps.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})

ROOT=Path(r"data/external/local_535c3822c456")
FIG,TAB,SCR,REP=[ROOT/x for x in ["figures","tables","scripts","reports"]]
KLF4_TARGETS=["CDKN1A","KRT14","KRT5","KRT17","EPCAM","TACSTD2","CDH1","CLDN1","DSP","JUP","ITGA6","ITGB4","TP63","NOTCH1","IVL"]
MODULES={
 "KLF4 repair proxy":KLF4_TARGETS,
 "ECM niche":["COL1A1","COL1A2","COL3A1","COL6A1","FN1","POSTN","LUM","DCN","MMP2"],
 "Growth-factor niche":["FGF7","FGF10","HGF","IGF1","TGFB1","PDGFA","VEGFA","WNT5A"],
 "Inflammatory niche":["IL6","CXCL8","CXCL2","CCL2","TNF","IL1B","PTGS2"],
 "Epithelial repair":["KRT14","KRT5","KRT17","EPCAM","TACSTD2","TP63","ITGA6","ITGB4"]}

def ensure():
    for p in [FIG,TAB,SCR,REP]:p.mkdir(parents=True,exist_ok=True)

def find_gene(a,g):
    names={str(x).upper():str(x) for x in a.var_names}
    if g.upper() in names:return names[g.upper()],"X"
    if a.raw is not None:
        raw={str(x).upper():str(x) for x in a.raw.var_names}
        if g.upper() in raw:return raw[g.upper()],"raw"
    return "",""

def vec(a,g,layer):
    x=a[:,g].X if layer=="X" else a.raw[:,g].X
    return np.asarray(x.toarray()).ravel() if hasattr(x,"toarray") else np.asarray(x).ravel()

def target_mask(a,s):
    if s.label.startswith(("Human","Rat embryo","Rat adult")):m,_=fibk14_annotation_call(a)
    else:m=marker_call(a,allow_krt5_fallback=False)[0]
    return m & subset_condition(a,s)

def sample_col(a):
    for c in ["sample_id","samples","sample","h5ad_sample","orig.ident","orig.ident.x"]:
        if c in a.obs:return c
    return ""

def condition_short(x):
    return {"Acomys unwounded":"Acomys control"}.get(x,x)

def summarize():
    specs=sc_specs()[:8];rows=[];availability=[];corr=[]
    for s in specs:
        a=ad.read_h5ad(s.path,backed="r");mask=target_mask(a,s);sc=sample_col(a)
        sid=a.obs[sc].astype(str).to_numpy() if sc else np.repeat("pooled",a.n_obs)
        requested=sorted(set(["KLF4","IGFL1"]+sum(MODULES.values(),[])))
        var_lookup={str(x).upper():str(x) for x in a.var_names}
        matched={g:var_lookup[g.upper()] for g in requested if g.upper() in var_lookup}
        genes={}
        if matched:
            ordered=list(dict.fromkeys(matched.values()))
            xm=a[:,ordered].X
            xm=xm.toarray() if hasattr(xm,"toarray") else np.asarray(xm)
            col={s:i for i,s in enumerate(ordered)}
            genes={g:np.asarray(xm[:,col[s]]).ravel() for g,s in matched.items()}
        availability.append({"group":condition_short(s.label),"cohort":s.cohort,"KLF4_available":"KLF4" in genes,
                             "IGFL1_available":"IGFL1" in genes,"KLF4_symbol":find_gene(a,"KLF4")[0],
                             "IGFL1_symbol":find_gene(a,"IGFL1")[0],"source":str(s.path)})
        for sample in np.unique(sid[mask]):
            k=mask&(sid==sample)
            rec={"group":condition_short(s.label),"cohort":s.cohort,"sample":sample,"n_target":int(k.sum())}
            for g in ["KLF4","IGFL1"]:
                if g in genes:
                    rec[f"{g}_mean"]=float(np.mean(genes[g][k]));rec[f"{g}_positive_pct"]=float(np.mean(genes[g][k]>0)*100)
                else:rec[f"{g}_mean"]=np.nan;rec[f"{g}_positive_pct"]=np.nan
            for mod,gl in MODULES.items():
                avail=[genes[g] for g in gl if g in genes]
                rec[mod]=float(np.mean(np.vstack([x[k] for x in avail]))) if avail and k.sum() else np.nan
                rec[f"{mod}_n_genes"]=len(avail)
            rows.append(rec)
        idx=np.where(mask)[0]
        if len(idx)>5000:idx=np.random.default_rng(20260713).choice(idx,5000,replace=False)
        for regulator in ["KLF4","IGFL1"]:
            if regulator not in genes:continue
            for mod,gl in MODULES.items():
                avail=[genes[g] for g in gl if g in genes and g!=regulator]
                if not avail:continue
                score=np.mean(np.vstack([x[idx] for x in avail]),axis=0)
                rho,p=stats.spearmanr(genes[regulator][idx],score)
                corr.append({"group":condition_short(s.label),"regulator":regulator,"module":mod,
                             "rho":rho,"p_value":p,"n_cells":len(idx),"n_module_genes":len(avail)})
        try:a.file.close()
        except Exception:pass
    d=pd.DataFrame(rows);av=pd.DataFrame(availability);co=pd.DataFrame(corr)
    if len(co):co["p_adj_bh"]=multipletests(co.p_value,method="fdr_bh")[1]
    d.to_csv(TAB/"sample_level_expression_and_modules.csv",index=False,encoding="utf-8-sig")
    av.to_csv(TAB/"KLF4_IGFL1_feature_space_availability.csv",index=False,encoding="utf-8-sig")
    co.to_csv(TAB/"regulator_niche_module_correlations.csv",index=False,encoding="utf-8-sig")
    return d,av,co

def effects(d):
    pairs=[("Human","Human wound","Human normal"),("Rat","Rat adult","Rat embryo"),
           ("Reindeer","Reindeer velvet","Reindeer back"),("Acomys","Acomys wound","Acomys control")]
    rows=[]
    for cohort,a,b in pairs:
        for metric in ["KLF4_mean","IGFL1_mean","KLF4 repair proxy","ECM niche","Growth-factor niche","Inflammatory niche","Epithelial repair"]:
            xa=d.loc[d.group.eq(a),metric].dropna().to_numpy();xb=d.loc[d.group.eq(b),metric].dropna().to_numpy()
            if len(xa)<1 or len(xb)<1:rows.append({"cohort":cohort,"contrast":f"{a} vs {b}","metric":metric,"effect":np.nan,"p_value":np.nan,"n_a":len(xa),"n_b":len(xb)});continue
            pooled=np.sqrt(((len(xa)-1)*np.var(xa,ddof=1)+(len(xb)-1)*np.var(xb,ddof=1))/max(len(xa)+len(xb)-2,1)) if len(xa)>1 and len(xb)>1 else np.nan
            effect=(np.mean(xa)-np.mean(xb))/pooled if np.isfinite(pooled) and pooled>0 else np.log2((np.mean(xa)+1e-6)/(np.mean(xb)+1e-6))
            p=stats.mannwhitneyu(xa,xb,alternative="two-sided").pvalue if len(xa)>1 and len(xb)>1 else np.nan
            rows.append({"cohort":cohort,"contrast":f"{a} vs {b}","metric":metric,"effect":effect,"p_value":p,"n_a":len(xa),"n_b":len(xb),"mean_a":np.mean(xa),"mean_b":np.mean(xb)})
    e=pd.DataFrame(rows);valid=e.p_value.notna();e.loc[valid,"p_adj_bh"]=multipletests(e.loc[valid,"p_value"],method="fdr_bh")[1]
    e.to_csv(TAB/"cross_species_effect_sizes.csv",index=False,encoding="utf-8-sig");return e

def plot_expression(d,av):
    order=["Human normal","Human wound","Rat embryo","Rat adult","Reindeer velvet","Reindeer back","Acomys wound","Acomys control"]
    fig,axs=plt.subplots(1,2,figsize=(140/25.4,72/25.4),gridspec_kw={"width_ratios":[1.55,1]},facecolor="white")
    ax=axs[0]
    for yi,g in enumerate(order):
        for xi,gene in enumerate(["KLF4","IGFL1"]):
            sub=d[d.group.eq(g)];vals=sub[f"{gene}_mean"].dropna();pct=sub[f"{gene}_positive_pct"].dropna()
            if vals.empty:ax.scatter(xi,yi,marker="x",s=35,c="#777")
            else:ax.scatter(xi,yi,s=20+2.2*pct.mean(),c=vals.mean(),cmap="viridis",vmin=0,vmax=max(d.KLF4_mean.max(),1),edgecolor="#333",linewidth=.4)
    ax.set_xticks([0,1],["KLF4","IGFL1"]);ax.set_yticks(range(len(order)),order);ax.invert_yaxis();ax.set_title("Expression in Fib_K14",loc="left")
    ax=axs[1];m=av.set_index("group").reindex(order)[["KLF4_available","IGFL1_available"]].astype(float)
    im=ax.imshow(m,cmap=mpl.colors.ListedColormap(["#E5E5E5","#3B8BC2"]),vmin=0,vmax=1,aspect="auto")
    ax.set_xticks([0,1],["KLF4","IGFL1"]);ax.set_yticks(range(len(order)),[]);ax.set_title("Feature availability",loc="left")
    for i in range(len(order)):
        for j in range(2):ax.text(j,i,"yes" if m.iloc[i,j] else "not in space",ha="center",va="center",fontsize=8,color="white" if m.iloc[i,j] else "#555")
    fig.tight_layout()
    for ext in ["png","pdf"]:fig.savefig(FIG/f"01_KLF4_IGFL1_expression_availability.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig)

def plot_effects(e):
    metrics=["KLF4_mean","IGFL1_mean","KLF4 repair proxy","ECM niche","Growth-factor niche","Inflammatory niche","Epithelial repair"]
    cohorts=["Human","Rat","Reindeer","Acomys"]
    mat=e.pivot(index="metric",columns="cohort",values="effect").reindex(index=metrics,columns=cohorts)
    fig,ax=plt.subplots(figsize=(105/25.4,78/25.4),facecolor="white");lim=np.nanmax(np.abs(mat.to_numpy()));lim=max(lim,1)
    im=ax.imshow(mat,cmap="RdBu_r",vmin=-lim,vmax=lim,aspect="auto")
    ax.set_xticks(range(4),cohorts);ax.set_yticks(range(len(metrics)),[x.replace("_mean","") for x in metrics])
    for i in range(len(metrics)):
        for j in range(4):
            v=mat.iloc[i,j];ax.text(j,i,"NA" if not np.isfinite(v) else f"{v:.2f}",ha="center",va="center",fontsize=8)
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.03);cb.set_label("Effect size\n(first - second)")
    ax.set_title("Cross-species injury/regeneration effects",loc="left");fig.tight_layout()
    for ext in ["png","pdf"]:fig.savefig(FIG/f"02_cross_species_effect_size_heatmap.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig)

def plot_correlations(co):
    if co.empty:return
    groups=list(dict.fromkeys(co.group));mods=list(MODULES)
    fig,axs=plt.subplots(2,1,figsize=(150/25.4,92/25.4),sharex=True,facecolor="white")
    for ax,reg in zip(axs,["KLF4","IGFL1"]):
        mat=co[co.regulator.eq(reg)].pivot(index="module",columns="group",values="rho").reindex(index=mods,columns=groups)
        im=ax.imshow(mat,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto");ax.set_yticks(range(len(mods)),mods);ax.set_title(f"{reg}–module association",loc="left")
        for i in range(len(mods)):
            for j in range(len(groups)):
                v=mat.iloc[i,j]
                if np.isfinite(v):ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=7)
    axs[-1].set_xticks(range(len(groups)),groups,rotation=40,ha="right")
    fig.subplots_adjust(left=.18,right=.88,bottom=.27,top=.95,hspace=.30)
    cax=fig.add_axes([.91,.35,.018,.36]);fig.colorbar(im,cax=cax);cax.set_title("ρ",fontsize=8,pad=2)
    for ext in ["png","pdf"]:fig.savefig(FIG/f"03_regulator_niche_module_correlations.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig)

def plot_target_conservation(d):
    order=["Human normal","Human wound","Rat embryo","Rat adult","Reindeer velvet","Reindeer back","Acomys wound","Acomys control"]
    mods=list(MODULES);mat=d.groupby("group")[mods].mean().reindex(order)
    z=(mat-mat.mean(axis=0))/mat.std(axis=0).replace(0,np.nan)
    fig,ax=plt.subplots(figsize=(105/25.4,72/25.4),facecolor="white");im=ax.imshow(z,cmap="RdBu_r",vmin=-2,vmax=2,aspect="auto")
    ax.set_xticks(range(len(mods)),[x.replace(" niche","") for x in mods],rotation=38,ha="right");ax.set_yticks(range(len(order)),order)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.03,label="Column z-score");ax.set_title("Conserved repair and niche programs",loc="left");fig.tight_layout()
    for ext in ["png","pdf"]:fig.savefig(FIG/f"04_KLF4_target_and_niche_program_heatmap.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig)

def main():
    ensure();d,av,co=summarize();e=effects(d);plot_expression(d,av);plot_effects(e);plot_correlations(co);plot_target_conservation(d)
    shutil.copy2(Path(__file__),SCR/Path(__file__).name)
    report={"KLF4_regulon":"curated epithelial-repair target-module proxy; not SCENIC",
            "IGFL1":"direct feature-space evidence only in human and reindeer; absence elsewhere is not proof of biological absence",
            "statistics":"sample-level summaries where replicate sample IDs exist; cell-level correlations are descriptive",
            "outputs":4,"editable_pdf":"Arial Type42"}
    (REP/"report.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps({"output_root":str(ROOT),**report},ensure_ascii=False,indent=2))

if __name__=="__main__":main()
