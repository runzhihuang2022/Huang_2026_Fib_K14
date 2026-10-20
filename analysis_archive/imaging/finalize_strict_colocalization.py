from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import tifffile
import zarr


PROJECT = Path(r"data/external/local_6b4e1e668702")
BASE = PROJECT / "Figure 3D 细胞计数"
OLD = BASE / "IF_colocalized_cell_statistics_all_slides_20260825"
OUT = BASE / "IF_colocalized_cell_statistics_strict_compartment_20260825"
MAP = BASE / "KC_Fib_SPTDI_DPTDI_whole_slide_revision_20260825/01_source_tables/sample_group_mode_mapping.tsv"
KC_RAW = PROJECT / "KC_Spinous_Mig"
FIB_OME = BASE / "Fib_K14_output/01_ome_tiff"


def sample_id(path: Path, family: str) -> str:
    s = path.name.upper()
    if family == "KC" and s.startswith("CUT-"):
        return "CUT_007"
    for key in ["W14-1", "32-1", "001A", "003", "007", "BS1", "BW3", "BW4", "BW7", "BW8", "BW9"]:
        if key in s:
            return key
    raise ValueError(path)


def norm(x):
    x=np.asarray(x,dtype=np.float32); nz=x[x>0]
    if not len(nz): return np.zeros_like(x)
    lo,hi=np.percentile(nz,[1,99.7]); return np.clip((x-lo)/(hi-lo+1e-9),0,1)**0.72


def overlay_kc(sample: str, cells: pd.DataFrame, dest: Path):
    path = next(p for p in KC_RAW.glob("*.tif") if sample_id(p,"KC")==sample)
    arr=tifffile.memmap(path,series=0); _,h,w=arr.shape; stride=max(1,int(np.ceil(max(h,w)/1400)))
    b=norm(arr[2,::stride,::stride]); g=norm(arr[1,::stride,::stride]); r=norm(arr[3,::stride,::stride])
    rgb=np.stack([r,g,b],axis=-1); x=cells.loc[cells.co_positive,"centroid_x_px"]/stride; y=cells.loc[cells.co_positive,"centroid_y_px"]/stride
    fig,ax=plt.subplots(figsize=(10,6),facecolor="black"); ax.imshow(rgb); ax.scatter(x,y,s=3,facecolors="none",edgecolors="#ffffff",linewidths=.35)
    ax.set_title(f"KC {sample}: strict compartment-positive cells",color="white",fontsize=9); ax.axis("off"); fig.savefig(dest,dpi=180,facecolor="black",bbox_inches="tight"); plt.close(fig)


def overlay_fib(sample: str, cells: pd.DataFrame, dest: Path):
    path = next(p for p in FIB_OME.glob("*.ome.tiff") if sample_id(p,"Fib_K14")==sample)
    store=tifffile.imread(path,series=0,level=6,aszarr=True); arr=zarr.open(store,mode="r"); block=np.asarray(arr[[0,1,4],:,:]); h,w=block.shape[1:]
    with tifffile.TiffFile(path) as tf: h2,w2=tf.series[0].levels[2].shape[1:]
    rgb=np.stack([norm(block[2]),norm(block[1]),norm(block[0])],axis=-1)
    x=cells.loc[cells.co_positive,"centroid_x_px"]*w/w2; y=cells.loc[cells.co_positive,"centroid_y_px"]*h/h2; store.close()
    fig,ax=plt.subplots(figsize=(7,7),facecolor="black"); ax.imshow(rgb); ax.scatter(x,y,s=4,facecolors="none",edgecolors="#ffffff",linewidths=.4)
    ax.set_title(f"Fib_K14 {sample}: strict compartment-positive cells",color="white",fontsize=9); ax.axis("off"); fig.savefig(dest,dpi=180,facecolor="black",bbox_inches="tight"); plt.close(fig)


for rel in ["00_manifest","02_qc/overlays/KC","02_qc/overlays/Fib_K14","04_per_slide","05_statistics/KC","05_statistics/Fib_K14","06_figures","07_excel","08_audit"]:
    (OUT/rel).mkdir(parents=True,exist_ok=True)

mapping=pd.read_csv(MAP,sep="\t",dtype=str).fillna("")
mapping=pd.concat([mapping,pd.DataFrame([{"family":"KC","sample":"CUT_007","group":"DPTDI","day_dpb":"7","mode":"Mode 1","group_mapping_provenance":"cropped derivative of 007; descriptive only"}])],ignore_index=True)
mapping.insert(2,"sample_label",mapping.family+":"+mapping["sample"])
mapping.to_csv(OUT/"00_manifest/sample_metadata.tsv",sep="\t",index=False)

per=pd.read_csv(OUT/"04_per_slide/per_slide_unmapped.tsv",dtype={"sample_id":str},sep="\t")
per=per.merge(mapping[["family","sample","sample_label","group","day_dpb","mode","group_mapping_provenance"]],left_on=["family","sample_id"],right_on=["family","sample"],how="left").drop(columns="sample")
per.loc[per.sample_id.eq("CUT_007"), "analysis_region"]="cropped_field"
per["subject_id"]=per.sample_id.replace({"CUT_007":"007"}); per["inference_included"]=per.sample_id.ne("CUT_007")
per["inference_exclusion_reason"]=np.where(per.inference_included,"","cropped derivative; not an independent whole slide")
per["group_mapping_status"]=np.where(per.sample_id.isin(["007","32-1","CUT_007"]),"REVIEW_REQUIRED","SUPPORTED")
per["segmentation_qc_status"]="REVIEW_REQUIRED"
cols=["family","sample_id","sample_label","subject_id","group","day_dpb","mode","positive_definition","numerator_n","denominator_n","ratio_pct","analysis_region","inference_included","inference_exclusion_reason","segmentation","primary_expansion","threshold_status","group_mapping_status","group_mapping_provenance","segmentation_qc_status"]
per[cols].to_csv(OUT/"04_per_slide/per_slide.tsv",sep="\t",index=False)
for family in ["KC","Fib_K14"]:
    x=per[per.family.eq(family)][cols]
    x.to_csv(OUT/f"04_per_slide/{family}_per_slide.tsv",sep="\t",index=False)
    x[x.inference_included].to_csv(OUT/f"04_per_slide/{family}_inference_input.tsv",sep="\t",index=False)

old_manifest=OLD/"00_manifest/input_manifest.tsv"
if old_manifest.exists(): shutil.copy2(old_manifest,OUT/"00_manifest/input_manifest.tsv")
channel=pd.DataFrame([
    ["KC",0,"SOX9","nuclear"],["KC",1,"MMP3","pericellular; secreted-protein attribution limited"],["KC",2,"DAPI","nuclear segmentation"],["KC",3,"KRT6B","cytoplasmic"],
    ["Fib_K14",0,"DAPI","nuclear segmentation"],["Fib_K14",1,"PDGFRA","outer-cell ring; no independent membrane marker"],["Fib_K14",4,"KRT14","cytoplasmic"]],columns=["family","channel_index_0based","marker","analysis_compartment"])
channel.to_csv(OUT/"00_manifest/channel_map.tsv",sep="\t",index=False)

qc=[]
for family in ["KC","Fib_K14"]:
    for p in sorted((OUT/"03_per_cell"/family).glob("*_per_cell.tsv.gz")):
        sample=p.name.split("_per_cell")[0]; cells=pd.read_csv(p,sep="\t",compression="gzip")
        dest=OUT/"02_qc/overlays"/family/f"{sample}_strict_overlay.png"
        (overlay_kc if family=="KC" else overlay_fib)(sample,cells,dest)
        row=per[(per.family==family)&(per.sample_id==sample)].iloc[0]
        qc.append({"family":family,"sample_id":sample,"denominator_n":row.denominator_n,"numerator_n":row.numerator_n,"overlay":str(dest),"status":"REVIEW_REQUIRED"})
pd.DataFrame(qc).to_csv(OUT/"02_qc/qc_summary.tsv",sep="\t",index=False)

(OUT/"README_请先阅读.md").write_text("""# 严格分区共染重算\n\n本版本不再使用重叠的核中心圆盘。KC使用核内SOX9、互斥细胞质KRT6B和核周细胞质MMP3；Fib使用互斥细胞边缘PDGFRA和细胞质KRT14。主分析扩张为3 µm（Fib）或3原生像素（KC），同时提供2/3/5敏感性分析。\n\n项目缺少无一抗/同型/已知阴阳性对照及独立泛膜通道，因此阈值为按全分析家族、未知分组的稳健Otsu自动阈值，结果必须标记EXPLORATORY。MMP3为分泌蛋白，KC结果应解释为MMP3-associated。CUT_007为裁剪视野，不进入推断统计。\n""",encoding="utf-8")
print(per[cols].to_string(index=False))
