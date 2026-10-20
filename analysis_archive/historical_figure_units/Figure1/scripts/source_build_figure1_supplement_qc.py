from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import anndata as ad
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from types import SimpleNamespace
from anndata._io.specs import read_elem

OUT = Path(__file__).resolve().parent
SC = Path(r"data/external/local_decd995edda5")
SPATIAL_ROOT = Path(r"data/external/local_be81937a46ed")
CONFIG = Path(r"data/external/local_02dfb9d3f713")
VISIUM = Path(r"data/external/local_fabb25142153")

SECTIONS = [
    ("5dpb", "BW32_A01597A3_SDSDB_5dpb", 135),
    ("12dpb DPTDI-1", "BW13_1_B1_DSDB_12dpb", 45),
    ("12dpb SPTDI-1", "BW14_1_C1_SSDB_12dpb", 135),
    ("12dpb DPTDI-2", "BW13_A3_DSDB_12dpb", 135),
    ("12dpb SPTDI-2", "BW14_B3_SSDB_12dpb", 135),
    ("19dpb", "BW15D_C6_SDSDB_19dpb", -225),
    ("19dpb p1", "BW15D_1_D1_SDSDB_19dpb", 140),
    ("26dpb p1", "BW81_C02846B6_SDSDB_26dpb_part1", 0),
    ("26dpb p2", "BW81_C02846B6_SDSDB_26dpb_part2", 0),
    ("2mph", "2mph_A03699G6.SCT", 45),
    ("Normal", "NS_C02847B1", 180),
]
FOCUS_STATES = ["KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous", "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular", "SAC_SG_Progenitor", "SAC_SG_Dark"]
GENE_SETS = {
    "Epithelialization": ["KRT14", "KRT5", "KRT10", "KRT1", "IVL", "SBSN", "LOR", "FLG", "AREG", "MMP3", "MMP9", "KRT17", "KRT6A"],
    "Skin development": ["TP63", "WNT10A", "LEF1", "SHH", "BMP2", "EDAR", "SOX9"],
    "Inflammatory response": ["IL1B", "IL6", "TNF", "CXCL8", "CCL2", "CXCL1", "CXCL2"],
    "Vascularization": ["VEGFA", "PECAM1", "CD34", "KDR", "FLT1", "ANGPT1", "VWF"],
    "Matrix remodeling / fibrosis": ["ACTA2", "COL1A1", "COL1A2", "COL3A1", "FN1", "TGFB1", "CTGF"],
}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.6,
})
sns.set_theme(style="white", context="paper")


def load_config():
    spec = importlib.util.spec_from_file_location("raw_annotation_config", CONFIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.raw_annotation_config


def rgba_hex(c):
    return c[:7] if len(c) >= 7 else c


def rotate_xy(x, y, angle):
    x, y = np.asarray(x, float), np.asarray(y, float)
    t = np.deg2rad(angle)
    xc, yc = np.nanmedian(x), np.nanmedian(y)
    xr = (x-xc)*np.cos(t) - (y-yc)*np.sin(t)
    yr = (x-xc)*np.sin(t) + (y-yc)*np.cos(t)
    return xr, yr


def spatial_file(sample_id):
    p = SPATIAL_ROOT / f"{sample_id}_Annotated.h5ad"
    if not p.exists():
        hits = list(SPATIAL_ROOT.rglob(f"{sample_id}_Annotated.h5ad"))
        if not hits:
            raise FileNotFoundError(sample_id)
        p = hits[0]
    return p


def read_obs_obsm(path):
    """Read only metadata/coordinates, bypassing legacy null values in uns/log1p."""
    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        obsm = read_elem(f["obsm"]) if "obsm" in f else {}
    return SimpleNamespace(obs=obs, obsm=obsm, n_obs=len(obs))


def grouped_marker_means(path, gene_indices, labels, states, n_vars, chunk=5000):
    """Stream a CSR H5AD matrix and calculate state means for selected genes."""
    sums = np.zeros((len(states), len(gene_indices)), dtype=np.float64)
    counts = np.array([(labels == s).sum() for s in states], dtype=np.int64)
    state_code = {s:i for i,s in enumerate(states)}
    with h5py.File(path, "r") as f:
        g = f["X"]
        indptr, indices, data = g["indptr"], g["indices"], g["data"]
        for start in range(0, len(labels), chunk):
            stop = min(start + chunk, len(labels))
            p = np.asarray(indptr[start:stop+1], dtype=np.int64)
            lo, hi = int(p[0]), int(p[-1])
            mat = sparse.csr_matrix((np.asarray(data[lo:hi]), np.asarray(indices[lo:hi]), p-lo), shape=(stop-start, n_vars))[:, gene_indices]
            block_labels = labels[start:stop]
            for s in np.unique(block_labels):
                if s in state_code:
                    sums[state_code[s]] += np.asarray(mat[block_labels == s].sum(axis=0)).ravel()
    return sums / np.maximum(counts[:, None], 1)


def get_xy(obs, obsm, angle):
    if "spatial" in obsm:
        xy = np.asarray(obsm["spatial"])
        x, y = xy[:, 0], xy[:, 1]
    else:
        x, y = obs["x"].to_numpy(), obs["y"].to_numpy()
    return rotate_xy(x, y, angle)


def save(fig, stem, dpi=300):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def broad_compartment(state):
    if state.startswith("KC_"): return "Keratinocyte"
    if state.startswith("SAC_"): return "Skin appendage"
    if state.startswith("Fib_"): return "Fibroblast"
    if state.startswith("Endo_") or state == "Pericyte": return "Endothelial / mural"
    if state in {"M1_Macrophage","M2_Macrophage","LAM","Neutrophil","Mast","cDC1","cDC2","pDC","LGHS","B_Mature","B_Plasma","CD4_T","CD8_T","NK"}: return "Immune"
    return "Other skin"


def main():
    cfg = load_config()
    sc = ad.read_h5ad(SC, backed="r")
    obs = sc.obs.copy()
    umap = np.asarray(sc.obsm["X_umap"])

    # Auditable metadata and QC tables.
    sc_rows = []
    for sid, d in obs.groupby("sample_id", observed=True):
        row = {"species":"Human", "modality":"scRNA-seq", "sample_id":str(sid), "section_id":"",
               "patient_id":str(d["patient"].iloc[0]), "injury_depth":str(d["group"].iloc[0]),
               "repair_stage":str(d["Time_category"].iloc[0]), "time_point":str(d["time_point"].iloc[0]),
               "n_observations":len(d)}
        sc_rows.append(row)
    spatial_rows, spatial_qc = [], []
    spatial_cache = {}
    for label, sid, angle in SECTIONS:
        a = read_obs_obsm(spatial_file(sid))
        spatial_cache[sid] = a
        d = a.obs
        spatial_rows.append({"species":"Human", "modality":"Stereo-seq", "sample_id":sid, "section_id":label,
                             "patient_id":"Not available in source object", "injury_depth":"DPTDI" if "DPTDI" in label else ("SPTDI" if "SPTDI" in label else "Not encoded"),
                             "repair_stage":label, "time_point":label.split()[0], "n_observations":a.n_obs})
        for metric in ["nCount_Spatial", "nFeature_Spatial"]:
            x = pd.to_numeric(d[metric], errors="coerce")
            spatial_qc.append({"sample_id":sid,"section_id":label,"metric":metric,"n":len(x),"min":x.min(),"q01":x.quantile(.01),"q25":x.quantile(.25),"median":x.median(),"q75":x.quantile(.75),"q99":x.quantile(.99),"max":x.max()})
    v = read_obs_obsm(VISIUM)
    vobs = v.obs.copy()
    vgroup = "sample_batch_new" if "sample_batch_new" in vobs.columns else ("orig.ident" if "orig.ident" in vobs.columns else None)
    visium_rows=[]
    if vgroup:
        for sid,d in vobs.groupby(vgroup, observed=True):
            visium_rows.append({"species":"Human","modality":"10x Visium","sample_id":str(sid),"section_id":str(sid),"patient_id":"Not available in source object","injury_depth":"Not encoded","repair_stage":"Not encoded","time_point":"Not encoded","n_observations":len(d)})
    pd.DataFrame(sc_rows+spatial_rows+visium_rows).to_csv(OUT/"cohort_metadata.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(spatial_qc).to_csv(OUT/"spatial_qc_summary.csv", index=False, encoding="utf-8-sig")

    qrows=[]
    for sid,d in obs.groupby("sample_id", observed=True):
        for metric in ["nFeature_RNA","nCount_RNA","mt_percent"]:
            x=pd.to_numeric(d[metric],errors="coerce")
            qrows.append({"sample_id":str(sid),"metric":metric,"n":len(x),"min":x.min(),"q01":x.quantile(.01),"q25":x.quantile(.25),"median":x.median(),"q75":x.quantile(.75),"q99":x.quantile(.99),"max":x.max()})
    pd.DataFrame(qrows).to_csv(OUT/"scrna_qc_by_sample.csv",index=False,encoding="utf-8-sig")

    marker_rows=[]
    for state,x in cfg.items():
        marker_rows.append({"cell_state":state,"compartment":broad_compartment(state),"color":rgba_hex(x["color"]),
                            "marker_1":x["markers"][0],"marker_2":x["markers"][1],"marker_3":x["markers"][2],"marker_4":x["markers"][3],
                            "cutoff_1":x["cutoffs"][0],"cutoff_2":x["cutoffs"][1],"cutoff_3":x["cutoffs"][2],"cutoff_4":x["cutoffs"][3]})
    pd.DataFrame(marker_rows).to_csv(OUT/"cell_state_markers.csv",index=False,encoding="utf-8-sig")
    pd.DataFrame([{"process":k,"genes":", ".join(v),"n_genes":len(v)} for k,v in GENE_SETS.items()]).to_csv(OUT/"figure1_gene_sets.csv",index=False,encoding="utf-8-sig")

    # A: object-backed cohort structure; 49-patient umbrella cohort is retained as a declared manuscript claim only.
    fig, ax = plt.subplots(1,3,figsize=(10.5,3.2))
    injury = obs.groupby("group",observed=True)["sample_id"].nunique().sort_values(ascending=False)
    injury.plot.bar(ax=ax[0],color="#3B82F6"); ax[0].set(title="scRNA-seq samples by injury depth",ylabel="Samples",xlabel="")
    stage = obs.groupby("Time_category",observed=True)["sample_id"].nunique().reindex(["0-7dpi","8-14dpi","15-28dpi","1-2mph","Normal"]).dropna()
    stage.plot.bar(ax=ax[1],color="#10B981"); ax[1].set(title="scRNA-seq samples by repair stage",ylabel="Samples",xlabel="")
    mods=pd.Series({"scRNA-seq samples":obs["sample_id"].nunique(),"Stereo-seq sections":len(SECTIONS),"Visium sections":len(visium_rows)})
    mods.plot.bar(ax=ax[2],color=["#4F46E5","#EA580C","#0891B2"]); ax[2].set(title="Figure 1 molecular datasets",ylabel="Count",xlabel="")
    for a in ax:
        a.spines[["top","right"]].set_visible(False); a.tick_params(axis="x",rotation=35)
    fig.suptitle("A  Cohort structure supported by live molecular objects\nDeclared umbrella cohort: 49 patients; molecular object: 23 patient IDs / 38 scRNA samples",fontweight="bold",y=1.05)
    fig.tight_layout(); save(fig,"Figure_S1A_cohort_structure")

    # B: QC before/after (raw server sample optional, current filtered source always included).
    rng=np.random.default_rng(20260721); idx=rng.choice(len(obs),min(50000,len(obs)),replace=False)
    after=obs.iloc[idx][["nFeature_RNA","nCount_RNA","mt_percent"]].copy(); after["QC_stage"]="After filtering"
    raw_csv=OUT/"raw_qc_50k.csv"
    qdat=pd.concat([pd.read_csv(raw_csv),after],ignore_index=True) if raw_csv.exists() else after
    long=qdat.melt(id_vars="QC_stage",var_name="Metric",value_name="Value")
    fig,axs=plt.subplots(1,3,figsize=(9,3.2))
    for a,m in zip(axs,["nFeature_RNA","nCount_RNA","mt_percent"]):
        sns.violinplot(data=long[long.Metric==m],x="QC_stage",y="Value",ax=a,inner="quart",cut=0,palette=["#CBD5E1","#2563EB"][:long.QC_stage.nunique()])
        a.set_title(m); a.set_xlabel(""); a.spines[["top","right"]].set_visible(False); a.tick_params(axis="x",rotation=20)
    note="Before: server pbmc.rds (316,194 cells); After: pbmc_final.h5ad (279,305 cells)" if raw_csv.exists() else "After-filtering object only; raw QC export pending"
    fig.suptitle("B  scRNA-seq quality control\n"+note,fontweight="bold",y=1.03); fig.tight_layout(); save(fig,"Figure_S1B_scrna_QC_before_after")

    # C: UMAP metadata controls.
    idx=rng.choice(len(obs),min(100000,len(obs)),replace=False)
    fig,axs=plt.subplots(1,3,figsize=(11,3.4))
    for a,col,title in zip(axs,["patient","time_point","group"],["Patient ID","Sampling time","Injury depth"]):
        vals=obs.iloc[idx][col].astype(str); cats=sorted(vals.unique()); pal=dict(zip(cats,sns.color_palette("husl",len(cats))))
        for cat in cats:
            mask=(vals.to_numpy()==cat); a.scatter(umap[idx[mask],0],umap[idx[mask],1],s=0.3,c=[pal[cat]],rasterized=True,label=cat,linewidths=0)
        a.set(title=title,xlabel="UMAP1",ylabel="UMAP2"); a.set_xticks([]); a.set_yticks([])
        if len(cats)<=12: a.legend(markerscale=8,fontsize=5,frameon=False,bbox_to_anchor=(1.01,1),loc="upper left")
    fig.suptitle("C  Metadata distribution across the human scRNA-seq atlas",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_S1C_UMAP_metadata")

    # D: 47-state marker heatmap using mean normalized expression and within-gene z-score.
    states=list(cfg.keys()); genes=[]
    for s in states:
        for g in cfg[s]["markers"]:
            if g not in genes and g in sc.var_names: genes.append(g)
    labels=obs["sub_labels"].astype(str).to_numpy(); means=[]
    gene_indices=[int(sc.var_names.get_loc(g)) for g in genes]
    means=grouped_marker_means(SC,gene_indices,labels,states,sc.n_vars)
    z=(means-means.mean(axis=0))/(means.std(axis=0)+1e-8); z=np.clip(z,-2.5,2.5)
    fig,ax=plt.subplots(figsize=(16,9)); sns.heatmap(z,ax=ax,cmap="vlag",center=0,vmin=-2.5,vmax=2.5,xticklabels=genes,yticklabels=states,cbar_kws={"label":"Mean expression z-score"})
    ax.set_title("D  Canonical marker expression across 47 annotated cell states",fontweight="bold"); ax.tick_params(axis="x",labelrotation=90,labelsize=5); ax.tick_params(axis="y",labelsize=6)
    fig.tight_layout(); save(fig,"Figure_S1D_marker_heatmap",dpi=240)

    # E: all 11 sections x two spatial QC metrics.
    fig,axs=plt.subplots(len(SECTIONS),2,figsize=(8,24))
    for i,(label,sid,angle) in enumerate(SECTIONS):
        a=spatial_cache[sid]; x,y=get_xy(a.obs,a.obsm,angle)
        for j,metric in enumerate(["nCount_Spatial","nFeature_Spatial"]):
            val=pd.to_numeric(a.obs[metric],errors="coerce").to_numpy(); vmax=np.nanquantile(val,.99)
            im=axs[i,j].scatter(x,y,c=val,s=0.3,cmap="viridis",vmin=0,vmax=vmax,rasterized=True,linewidths=0)
            axs[i,j].set_aspect("equal"); axs[i,j].axis("off"); axs[i,j].set_title(f"{label} | {metric}",fontsize=7)
            fig.colorbar(im,ax=axs[i,j],fraction=.025,pad=.01)
    fig.suptitle("E  Spatial transcriptomic QC across all 11 registered Figure 1 Stereo-seq sections",fontweight="bold",y=.995)
    fig.tight_layout(rect=[0,0,1,.99]); save(fig,"Figure_S1E_all11_spatial_QC",dpi=250)

    # F/G: selected state transfer scores and categorical validation at 19dpb.
    sid="BW15D_C6_SDSDB_19dpb"; a=spatial_cache[sid]; x,y=get_xy(a.obs,a.obsm,-225)
    fig,axs=plt.subplots(3,3,figsize=(9,8))
    for ax,state in zip(axs.ravel(),FOCUS_STATES):
        col=f"{state}_MergeScore"; val=pd.to_numeric(a.obs[col],errors="coerce").to_numpy(); vmax=np.nanquantile(val,.995)
        im=ax.scatter(x,y,c=val,s=.6,cmap="magma",vmin=0,vmax=vmax,rasterized=True,linewidths=0); ax.axis("off"); ax.set_aspect("equal"); ax.set_title(state,fontsize=7); fig.colorbar(im,ax=ax,fraction=.03,pad=.01)
    fig.suptitle("F  Single-cell-to-spatial mapping scores | representative 19 dpb section",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_S1F_mapping_scores_19dpb",dpi=280)

    fig,axs=plt.subplots(3,3,figsize=(9,8)); labels19=a.obs["sub_labels"].astype(str).to_numpy()
    for ax,state in zip(axs.ravel(),FOCUS_STATES):
        mask=labels19==state; ax.scatter(x,y,c="#E5E7EB",s=.45,rasterized=True,linewidths=0); ax.scatter(x[mask],y[mask],c=rgba_hex(cfg[state]["color"]),s=.9,rasterized=True,linewidths=0)
        ax.axis("off"); ax.set_aspect("equal"); ax.set_title(f"{state} (n={mask.sum():,})",fontsize=7)
    fig.suptitle("G  Spatial validation of epithelial and appendage-associated cell states | 19 dpb",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_S1G_state_validation_19dpb",dpi=280)

    # Compact composite for manuscript assembly.
    panels=[("A",OUT/"Figure_S1A_cohort_structure.png"),("B",OUT/"Figure_S1B_scrna_QC_before_after.png"),("C",OUT/"Figure_S1C_UMAP_metadata.png"),("D",OUT/"Figure_S1D_marker_heatmap.png"),("E",OUT/"Figure_S1E_all11_spatial_QC.png"),("F",OUT/"Figure_S1F_mapping_scores_19dpb.png"),("G",OUT/"Figure_S1G_state_validation_19dpb.png")]
    fig=plt.figure(figsize=(16,20)); gs=fig.add_gridspec(4,2,height_ratios=[1,1.05,1.5,1.1])
    positions=[gs[0,0],gs[0,1],gs[1,0],gs[1,1],gs[2,:],gs[3,0],gs[3,1]]
    for (lab,p),pos in zip(panels,positions):
        ax=fig.add_subplot(pos); ax.imshow(plt.imread(p)); ax.axis("off"); ax.text(-.01,1.02,lab,transform=ax.transAxes,fontsize=15,fontweight="bold",va="bottom")
    fig.suptitle("Figure S1 | Cohort structure, sequencing quality, cell annotation and spatial mapping controls",fontsize=14,fontweight="bold",y=.995)
    fig.tight_layout(rect=[0,0,1,.99]); save(fig,"Figure_S1_QC_composite",dpi=220)

    summary={"scRNA_cells_after_filtering":int(sc.n_obs),"scRNA_samples":int(obs.sample_id.nunique()),"scRNA_patients_in_object":int(obs.patient.nunique()),
             "declared_umbrella_cohort_patients":49,"Stereo_sections":len(SECTIONS),"Stereo_cell_bins":int(sum(x["n_observations"] for x in spatial_rows)),
             "Visium_sections":len(visium_rows),"Visium_spots":int(sum(x["n_observations"] for x in visium_rows)),"cell_states":len(cfg),
             "raw_scRNA_cells_before_filtering":316194}
    (OUT/"verified_counts.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
