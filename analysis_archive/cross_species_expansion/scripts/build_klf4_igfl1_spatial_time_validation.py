from __future__ import annotations
import json, shutil
from pathlib import Path
import anndata as ad
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

from build_fibk14_density_reference_style import sc_specs

mpl.rcParams.update({"font.family":"Arial","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
 "font.size":8,"axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":8,"ytick.labelsize":8,
 "pdf.fonttype":42,"ps.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
ROOT=Path(r"data/external/local_535c3822c456")
FIG,TAB,SCR,REP=[ROOT/x for x in ["figures","tables","scripts","reports"]]
RAT=Path(r"data/external/local_fb7c1ee2a2f4")
RS=[("Embryo","E16.5 6H","E16.5_6H_C02938E1_Annotated_Spatial.h5ad"),
    ("P5","P5 3D","P5_3D_B02621C2_Annotated_Spatial.h5ad"),
    ("Adult","Adult 10d","Adult_10d_A02988C5_Annotated_Spatial.h5ad")]

def gene(a,g):
    lookup={str(x).upper():str(x) for x in a.var_names};s=lookup.get(g.upper(),"")
    if not s:return None,""
    x=a[:,s].X;return (np.asarray(x.toarray()).ravel() if hasattr(x,"toarray") else np.asarray(x).ravel()),s

def den(coords,weights,bins=150):
    x,y=coords[:,0],coords[:,1];ix=np.clip(((x-x.min())/(x.max()-x.min()+1e-9)*(bins-1)).astype(int),0,bins-1);iy=np.clip(((y-y.min())/(y.max()-y.min()+1e-9)*(bins-1)).astype(int),0,bins-1)
    h=np.zeros((bins,bins));np.add.at(h,(iy,ix),weights);h=gaussian_filter(h,1.3);v=h[iy,ix];hi=np.quantile(v[v>0],.99) if np.any(v>0) else 1;return np.clip(v/hi,0,1)

def rat_spatial():
    rows=[];fig,axs=plt.subplots(1,4,figsize=(180/25.4,55/25.4),gridspec_kw={"width_ratios":[1,1,1,1.30]},facecolor="white")
    for ax,(stage,label,name) in zip(axs[:3],RS):
        a=ad.read_h5ad(RAT/name,backed="r");key="X_rot_spatial" if "X_rot_spatial" in a.obsm else "spatial";c=np.asarray(a.obsm[key])[:,:2]
        fib=a.obs["sub_labels"].astype(str).to_numpy()=="Fib_K14";v,sym=gene(a,"KLF4")
        if v is None:v=np.zeros(a.n_obs)
        score=den(c,np.where(fib,v,0));idx=np.argsort(score);idx=idx[np.linspace(0,len(idx)-1,min(35000,len(idx)),dtype=int)]
        ax.scatter(c[idx,0],c[idx,1],c=score[idx],s=.12,cmap="magma",vmin=0,vmax=1,linewidths=0);ax.set_title(label);ax.axis("off");ax.set_aspect("equal")
        left,right=np.quantile(c[:,0],[.25,.75]);band=(right-left)*.10;x=c[:,0];front=(np.abs(x-left)<=band)|(np.abs(x-right)<=band);center=(x>left+band)&(x<right-band);reg=np.where(front,"Front",np.where(center,"Center","Peri"))
        for r in ["Center","Front","Peri"]:
            k=fib&(reg==r);rows.append({"stage":stage,"sample":label,"region":r,"n_Fib_K14":int(k.sum()),"KLF4_symbol":sym,"KLF4_mean":float(v[k].mean()) if k.sum() else np.nan,"KLF4_positive_pct":float((v[k]>0).mean()*100) if k.sum() else np.nan,"IGFL1_available":False})
        try:a.file.close()
        except Exception:pass
    d=pd.DataFrame(rows);ax=axs[3];width=.22
    for i,st in enumerate(["Embryo","P5","Adult"]):
        for j,(r,col) in enumerate(zip(["Center","Front","Peri"],["#7A9CC6","#D95F5F","#B5B5B5"])):
            y=d[(d.stage==st)&(d.region==r)].KLF4_mean.iloc[0];ax.bar(i+(j-1)*width,y,width,color=col,edgecolor="#444",linewidth=.4)
    ax.set_xticks(range(3),["Emb.","P5","Adult"]);ax.set_ylabel("KLF4 mean");ax.legend([mpl.patches.Patch(color=x) for x in ["#7A9CC6","#D95F5F","#B5B5B5"]],["Center","Front","Peri"],frameon=False,fontsize=8)
    fig.subplots_adjust(left=.03,right=.99,bottom=.18,top=.88,wspace=.38)
    for ext in ["png","pdf"]:fig.savefig(FIG/f"05_rat_spatial_KLF4_region_validation.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig);d.to_csv(TAB/"rat_spatial_KLF4_region_summary.csv",index=False,encoding="utf-8-sig")

def reindeer_time():
    s=next(x for x in sc_specs() if x.label=="Reindeer velvet");a=ad.read_h5ad(s.path,backed="r");fib=a.obs["Fib_K14"].astype(str).to_numpy()=="Fib_K14";tissue=a.obs["tissue"].astype(str).to_numpy();time=a.obs["timepoint"].astype(str).to_numpy();sample=a.obs["sample"].astype(str).to_numpy()
    genes={g:gene(a,g)[0] for g in ["KLF4","IGFL1"]};rows=[]
    for ti,tlabel in [("antler_velvet","Velvet"),("back_skin","Back")]:
        for tp in ["0","3","7","14"]:
            k=fib&(tissue==ti)&(time==tp)
            for g,v in genes.items():
                rows.append({"tissue":tlabel,"time_day":int(tp),"gene":g,"n_cells":int(k.sum()),"mean_expression":float(v[k].mean()) if k.sum() else np.nan,"positive_pct":float((v[k]>0).mean()*100) if k.sum() else np.nan})
    d=pd.DataFrame(rows);fig,axs=plt.subplots(1,2,figsize=(120/25.4,55/25.4),facecolor="white")
    for ax,g in zip(axs,["KLF4","IGFL1"]):
        q=d[d.gene==g]
        for tissue,col,mark in [("Velvet","#2A9D8F","o"),("Back","#6C757D","s")]:
            z=q[q.tissue==tissue];ax.plot(z.time_day,z.mean_expression,marker=mark,color=col,label=tissue,lw=1.2);ax.scatter(z.time_day,z.mean_expression,s=18+z.positive_pct*1.5,color=col,edgecolor="#333",linewidth=.4)
        ax.set_title(g);ax.set_xlabel("Day");ax.set_ylabel("Mean expression");ax.set_xticks([0,3,7,14]);ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ["png","pdf"]:fig.savefig(FIG/f"06_reindeer_KLF4_IGFL1_timecourse.{ext}",dpi=600 if ext=="png" else None,facecolor="white")
    plt.close(fig);d.to_csv(TAB/"reindeer_KLF4_IGFL1_timecourse.csv",index=False,encoding="utf-8-sig")
    try:a.file.close()
    except Exception:pass

def main():
    for p in [FIG,TAB,SCR,REP]:p.mkdir(parents=True,exist_ok=True)
    rat_spatial();reindeer_time();shutil.copy2(Path(__file__),SCR/Path(__file__).name)
    print(json.dumps({"output_root":str(ROOT),"new_figures":["05 rat spatial KLF4","06 reindeer KLF4/IGFL1 timecourse"]},ensure_ascii=False))
if __name__=="__main__":main()
