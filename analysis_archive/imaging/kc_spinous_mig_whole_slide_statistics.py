from __future__ import annotations

import csv, gzip, json, re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.stats import kruskal, mannwhitneyu
from skimage.filters import gaussian
from statsmodels.stats.multitest import multipletests

from fib_k14_whole_slide_statistics import TILE, MARGIN, segment_tile, p90_disk, robust_otsu

ROOT = Path(r"data/external/local_6b4e1e668702")
INPUT = ROOT / "KC_Spinous_Mig"
OUT = ROOT / "Figure 3D 细胞计数" / "KC_Spinous_Mig_whole_slide_statistics_native_pixel_v3"
DS = 1
MIN_PX = 8.0
MAX_PX = 500.0
MARKER_RADIUS_PX = 5
BG_SIGMA_PX = 15.0
META = {
    "007": (7, "Mode 1"), "32-1": (6, "Mode 1"), "BS1": (12, "Mode 2"),
    "BW3": (19, "Mode 3"), "BW4": (19, "Mode 3"), "BW7": (22, "Mode 3"),
    "BW8": (22, "Mode 3"), "BW9": (26, "Mode 3"), "W14-1": (14, "Mode 2"),
    "001A": (None, "UNMAPPED"),
}


def sid(path: Path) -> str:
    s = path.name.upper()
    for key in ["W14-1", "32-1", "001A", "007", "BS1", "BW3", "BW4", "BW7", "BW8", "BW9"]:
        if key in s:
            return key
    raise ValueError(path)


def analyze(path: Path, out_path: Path):
    sample = sid(path)
    with tifffile.TiffFile(path) as tf, gzip.open(out_path, "wt", encoding="utf-8", newline="") as fh:
        raw = [pg.asarray(out="memmap") for pg in tf.pages[:4]]
        # ImageJ LUT and morphology audit: C1 red=SOX9(570),
        # C2 green=MMP3(520), C3 blue/nuclear=DAPI, C4 grayscale=KRT6B(690).
        # Reorder to
        # DAPI, MMP3, SOX9, KRT6B for analysis.
        pages = [raw[2], raw[1], raw[0], raw[3]]
        h, w = pages[0].shape
        h2, w2 = (h + DS - 1)//DS, (w + DS - 1)//DS
        wr = csv.writer(fh, delimiter="\t")
        wr.writerow(["sample", "cell_id", "x_native_px", "y_native_px", "area_px2",
                     "MMP3_p90_bgsub", "SOX9_p90_bgsub", "KRT6B_p90_bgsub"])
        n = 0; radius = MARKER_RADIUS_PX; bg = BG_SIGMA_PX
        from skimage.measure import regionprops
        for y0 in range(0, h2, TILE):
            for x0 in range(0, w2, TILE):
                cy1, cx1 = min(h2, y0+TILE), min(w2, x0+TILE)
                ry0, rx0 = max(0, y0-MARGIN), max(0, x0-MARGIN)
                ry1, rx1 = min(h2, cy1+MARGIN), min(w2, cx1+MARGIN)
                block = [np.asarray(a[ry0*DS:min(h, ry1*DS):DS, rx0*DS:min(w, rx1*DS):DS]) for a in pages]
                shape = tuple(min(x.shape[i] for x in block) for i in [0, 1])
                block = [x[:shape[0], :shape[1]] for x in block]
                labels = segment_tile(block[0])
                markers = []
                for a in block[1:4]:
                    a = a.astype(np.float32)
                    markers.append(np.maximum(a-gaussian(a, sigma=bg, preserve_range=True), 0))
                for r in regionprops(labels):
                    if not (MIN_PX <= r.area <= MAX_PX): continue
                    yy, xx = r.centroid; gy, gx = yy+ry0, xx+rx0
                    if not (y0 <= gy < cy1 and x0 <= gx < cx1): continue
                    n += 1
                    vals = [p90_disk(a, yy, xx, radius) for a in markers]
                    wr.writerow([sample, n, f"{gx:.3f}", f"{gy:.3f}", f"{r.area:.3f}"] + [f"{v:.6g}" for v in vals])
    return n


def infer(per: pd.DataFrame, statdir: Path):
    mapped = per[per["mode"].str.startswith("Mode")]
    order = ["Mode 1", "Mode 2", "Mode 3"]
    vals = [mapped.loc[mapped["mode"] == m, "triple_positive_ratio_pct"].to_numpy() for m in order]
    rows = []
    h, p = kruskal(*vals)
    rows.append({"test_family":"omnibus", "contrast":"Mode 1 vs Mode 2 vs Mode 3", "test":"Kruskal-Wallis", "statistic":h, "p_raw":p, "p_BH":p})
    pairs=[]
    for i,a in enumerate(order):
        for b in order[i+1:]:
            x=mapped.loc[mapped["mode"]==a,"triple_positive_ratio_pct"]; y=mapped.loc[mapped["mode"]==b,"triple_positive_ratio_pct"]
            u,p=mannwhitneyu(x,y,alternative="two-sided",method="auto")
            pairs.append({"test_family":"pairwise","contrast":f"{a} vs {b}","test":"Mann-Whitney U","statistic":u,"p_raw":p})
    qs=multipletests([r["p_raw"] for r in pairs],method="fdr_bh")[1]
    for r,q in zip(pairs,qs): r["p_BH"]=q
    pd.DataFrame(rows+pairs).to_csv(statdir/"KC_statistical_tests.tsv",sep="\t",index=False)
    mapped.groupby("mode")["triple_positive_ratio_pct"].agg(["count","mean","std","median","min","max"]).reset_index().to_csv(statdir/"KC_mode_summary.tsv",sep="\t",index=False)


def plot(per: pd.DataFrame, figdir: Path):
    order=["Mode 1","Mode 2","Mode 3"]; df=per[per["mode"].isin(order)]
    plt.rcParams.update({"font.family":"Arial","font.size":8,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
    fig,ax=plt.subplots(figsize=(42/25.4,28/25.4)); rng=np.random.default_rng(20260825)
    for i,m in enumerate(order):
        y=df.loc[df["mode"]==m,"triple_positive_ratio_pct"].to_numpy()
        ax.scatter(i+rng.uniform(-.055,.055,len(y)),y,s=16,facecolors="white",edgecolors="#3198B8",linewidths=.8,zorder=3)
    g=df.groupby("mode")["triple_positive_ratio_pct"].agg(["mean","sem"]).reindex(order)
    ax.plot(range(3),g["mean"],color="#3198B8",marker="o",markersize=3.5,linewidth=1.2)
    ax.fill_between(range(3),(g["mean"]-g["sem"]).to_numpy(),(g["mean"]+g["sem"]).to_numpy(),color="#3198B8",alpha=.13,linewidth=0)
    ax.set_xticks(range(3),["Mode\n1","Mode\n2","Mode\n3"]); ax.set_ylabel("SOX9+ MMP3+\nKRT6B+\ncell ratio (%)",fontsize=8,labelpad=0)
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(width=.8,length=3,labelsize=8)
    fig.subplots_adjust(left=.40,right=.96,bottom=.27,top=.95)
    for ext in ["pdf","svg","png"]: fig.savefig(figdir/f"KC_whole_slide_lineplot.{ext}",dpi=600 if ext=="png" else None)
    plt.close(fig)


def main():
    cell=OUT/"01_per_cell"; stat=OUT/"02_statistics"; fig=OUT/"03_figures"; audit=OUT/"04_audit"
    for d in [cell,stat,fig,audit]: d.mkdir(parents=True,exist_ok=True)
    paths=[p for p in sorted(INPUT.glob("MAX_*.tif")) if not p.name.lower().startswith("cut-")]
    for p in paths:
        op=cell/f"{sid(p)}_whole_slide_cells.tsv.gz"
        if not op.exists(): print("ANALYZE",sid(p),flush=True); analyze(p,op)
    allc=pd.concat([pd.read_csv(p,sep="\t",dtype={"sample":str}) for p in sorted(cell.glob("*.tsv.gz"))],ignore_index=True)
    th={m:robust_otsu(allc[f"{m}_p90_bgsub"].to_numpy()) for m in ["MMP3","SOX9","KRT6B"]}
    for m,t in th.items(): allc[f"{m}_positive"]=allc[f"{m}_p90_bgsub"]>t
    allc["triple_positive"]=allc[[f"{m}_positive" for m in th]].all(axis=1)
    rows=[]
    for s,x in allc.groupby("sample",sort=True):
        day,mode=META[s]
        rows.append({"sample":s,"day_dpb":day,"mode":mode,"nuclei_count":len(x),
                     **{f"{m}_positive_n":int(x[f'{m}_positive'].sum()) for m in th},
                     "triple_positive_n":int(x.triple_positive.sum()),"triple_positive_ratio_pct":float(x.triple_positive.mean()*100)})
    per=pd.DataFrame(rows); per.to_csv(stat/"KC_whole_slide_per_slide.tsv",sep="\t",index=False)
    infer(per,stat); plot(per,fig)
    (audit/"analysis_config.json").write_text(json.dumps({"scope":"whole slide; ROI not used","source_channel_order":{"C1":"SOX9","C2":"MMP3","C3":"DAPI","C4":"KRT6B"},"analysis_order":["DAPI","MMP3","SOX9","KRT6B"],"channel_order_evidence":"ImageJ LUTs red, green, blue, grayscale; C3 punctate nuclear morphology; wavelengths in filename","calibration_status":"CALIBRATION_REQUIRED; source TIFF has no micrometres-per-pixel metadata","analysis_scale":"native pixels","nuclear_area_px2":[MIN_PX,MAX_PX],"marker_radius_px":MARKER_RADIUS_PX,"background_sigma_px":BG_SIGMA_PX,"locked_thresholds":th,"inference_unit":"slide","status":"exploratory; validate automated segmentation"},indent=2),encoding="utf-8")
    print(per.to_string(index=False)); print(th)

if __name__=="__main__": main()
