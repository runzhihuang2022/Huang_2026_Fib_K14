"""Uniform H&E, Arial 8–10 pt, nominal scale bars and missing-aware four-marker gate."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import sys
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import anndata as ad
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from PIL import Image


ROOT=Path(r"data/external/local_e63c701db393")
BASE=ROOT/"12_Fib_K14_KRT14_PDGFRA_same_position_20261014_v01"
INPUT=ROOT/"10_xaxis_reflected_HE_merged_atlas/h5ad/spatial_proteomics.xaxis_reflected_epidermis_up.HE_aligned.h5ad"
HE=ROOT/"13_Fib_K14_HE_background_20261014_v01/assets/HE_epidermis_up_photo_band.png"
META=ROOT/"02_extracted_matrices/protein_metadata.tsv"
AI_IDS=ROOT/"02_extracted_matrices/ai_quant_accessions.tsv"
PPT=Path(r"data/external/local_8853af99994d")
OUT=ROOT/"14_Fib_K14_four_marker_HE_Arial8to10_20261014_v01"
PLOTS,TABLES,AUDIT,ASSETS,CODE,H5=[OUT/s for s in ["plots","tables","audit","assets","code","h5ad"]]
GENES=["KRT14","PDGFRA","KLF4","TP63"]
EXTENT=(-.5,64.5,-.5,38.5)
GRID_UM=100.0
HE_WEIGHT=.72
RED="#ED2024"
ORANGE="#E88718"
GREY="#C9C9C9"
BAR_RECORDS=[]
FIGURE_RECORDS=[]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
    return h.hexdigest()


def dump(path,value):
    def cast(obj):
        if isinstance(obj,np.ndarray):return obj.tolist()
        if isinstance(obj,np.generic):return obj.item()
        if isinstance(obj,Path):return str(obj)
        raise TypeError(type(obj).__name__)
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False,default=cast),encoding="utf-8")


def dense(m):return m.toarray() if hasattr(m,"toarray") else np.asarray(m)


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def inside(xy,b):return (xy[:,0]>=b[0])&(xy[:,0]<=b[1])&(xy[:,1]>=b[2])&(xy[:,1]<=b[3])


def backdrop(ax,he,limits,zoom=False):
    ax.imshow(he,extent=EXTENT,origin="upper",interpolation="bilinear",alpha=1,zorder=0)
    ax.set(xlim=limits[:2],ylim=limits[2:],aspect="equal",xticks=[],yticks=[])
    for spine in ax.spines.values():spine.set_visible(zoom);spine.set_linewidth(.55)


def roi_box(ax,b,label=None):
    ax.add_patch(Rectangle((b[0],b[2]),b[1]-b[0],b[3]-b[2],fill=False,ec="black",lw=.6,ls=(0,(4,2)),zorder=10))
    if label:ax.text(b[0]+.5,b[3]-.5,label,fontsize=8,va="top",bbox=dict(facecolor="white",edgecolor="none",pad=.5),zorder=11)


def scale_bar(ax,limits,zoom,panel_id):
    um=200. if zoom else 1000.
    n_units=um/GRID_UM
    x1=limits[1]-(.8 if zoom else 1.2)
    y=limits[2]+(.7 if zoom else 1.2)
    line=ax.plot([x1-n_units,x1],[y,y],color="black",lw=1.1,zorder=21,
                 path_effects=[pe.Stroke(linewidth=2.8,foreground="white"),pe.Normal()])[0]
    ax.text(x1,y+(.25 if zoom else .55),f"{um:g} µm",ha="right",va="bottom",fontsize=8,
            color="black",fontstyle="normal",bbox=dict(facecolor="white",edgecolor="none",alpha=.85,pad=.4),zorder=22)
    BAR_RECORDS.append(dict(panel_id=panel_id,length_um=um,data_length=n_units,line=line,ax=ax))


def panel(ax,he,xy,target,feature,limits,roi,zoom,panel_id,color=RED,dd=None,missing_text=None,all_points=False):
    backdrop(ax,he,limits,zoom)
    visible=inside(xy,limits)
    if missing_text is not None:
        ax.scatter(*xy[visible].T,c=GREY,s=1.2 if not zoom else 3,alpha=.10,linewidths=0,zorder=1)
        ax.text(.5,.58,missing_text,transform=ax.transAxes,ha="center",va="center",fontsize=8,
                bbox=dict(facecolor="white",alpha=.90,edgecolor="none",pad=2),zorder=15)
    else:
        other=visible&~target
        if not all_points:ax.scatter(*xy[other].T,c=GREY,s=1.4 if not zoom else 3.6,alpha=.10,linewidths=0,zorder=1)
        ids=np.flatnonzero(visible if all_points else visible&target)
        if feature is None:
            ax.scatter(*xy[ids].T,c=color,s=2.8 if not zoom else 7,alpha=.94,linewidths=0,zorder=4)
        else:
            ids=ids[np.argsort(feature[ids],kind="stable")]
            ax.scatter(*xy[ids].T,c=feature[ids],cmap="viridis",vmin=0,vmax=1,s=2.8 if not zoom else 7.2,
                       alpha=.97,linewidths=0,zorder=4)
            if dd is not None:
                gx,gy,z=dd
                ax.contour(gx,gy,z,levels=[.5*z.max(),.85*z.max()],colors=["#DAA520","#FF4500"],linewidths=.5,zorder=6)
    if not zoom:roi_box(ax,roi)
    scale_bar(ax,limits,zoom,panel_id)


def export(fig,base,n_spatial):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    sizes=[t.get_fontsize() for t in texts]
    assert min(sizes)>=8-1e-6 and max(sizes)<=10+1e-6,(min(sizes),max(sizes))
    assert all(t.get_fontfamily()==["Arial"] and t.get_fontstyle()=="normal" for t in texts)
    outside=[]
    for t in texts:
        b=t.get_window_extent(renderer)
        if b.x0<-1 or b.y0<-1 or b.x1>fig.bbox.x1+1 or b.y1>fig.bbox.y1+1:outside.append(t.get_text())
    if outside:raise RuntimeError(f"Text outside canvas: {outside}")
    fig.savefig(base.with_name(base.name+"_preview").with_suffix(".png"),dpi=180,facecolor="white")
    for ext in ["pdf","svg","png","tiff"]:
        opts={"pil_kwargs":{"compression":"tiff_lzw"}} if ext=="tiff" else {}
        fig.savefig(base.with_suffix("."+ext),dpi=600,facecolor="white",**opts)
        if ext in ["png","tiff"]:
            with Image.open(base.with_suffix("."+ext)) as im:rgb=im.convert("RGB")
            rgb.save(base.with_suffix("."+ext),dpi=(600,600),**({"compression":"tiff_lzw"} if ext=="tiff" else {}))
    for record in BAR_RECORDS:
        if record["ax"].figure is fig:
            pts=record["line"].get_transform().transform(record["line"].get_xydata())
            drawn_mm=float(np.linalg.norm(pts[1]-pts[0])/fig.dpi*25.4)
            record["drawn_mm"]=drawn_mm
            record["nominal_data_units_per_drawn_mm"]=record["data_length"]/drawn_mm
    with fitz.open(base.with_suffix(".pdf")) as d:
        p=d[0]
        spans=[s for b in p.get_text("dict")["blocks"] if "lines" in b for line in b["lines"] for s in line["spans"] if s["text"].strip()]
        assert min(s["size"] for s in spans)>=7.99 and max(s["size"] for s in spans)<=10.01
        assert all("Arial" in s["font"] for s in spans)
        assert not any("Italic" in s["font"] or "Oblique" in s["font"] for s in spans)
        imgs=[m for m in p.get_image_info() if min(fitz.Rect(m["bbox"]).width,fitz.Rect(m["bbox"]).height)>15]
        assert len(imgs)==n_spatial,(base.name,len(imgs),n_spatial)
        FIGURE_RECORDS.append(dict(stem=base.name,canvas_mm=[fig.get_figwidth()*25.4,fig.get_figheight()*25.4],
                                   font_range_pt=[min(s["size"] for s in spans),max(s["size"] for s in spans)],
                                   fonts=sorted({s["font"] for s in spans}),he_background_panels=len(imgs),text_within_canvas=True))
    plt.close(fig)


def compact(he,xy,rois,limits,rows,title,subtitle,legend_color,legend_label,counts,base,strict=False):
    W=210.;H=187. if len(rows)==4 else 264.
    fig=plt.figure(figsize=(W/25.4,H/25.4))
    fig.text(.5,1-4.5/H,title,ha="center",va="center",fontsize=10,fontweight="bold")
    fig.text(.5,1-10.7/H,subtitle,ha="center",va="center",fontsize=8)
    full_w,zoom_w,gap,row_h=59.4,19.8,.6,38.5
    full_h=full_w*(limits[3]-limits[2])/(limits[1]-limits[0])
    for ri,r in enumerate(rois):
        x=[22.,104.][ri]
        fig.text((x+full_w/2)/W,1-17/H,f"Full slice: ROI {r['roi_id']}",ha="center",fontsize=8.5,fontweight="bold")
        fig.text((x+full_w+gap+zoom_w/2)/W,1-17/H,f"Zoom {r['roi_id']}",ha="center",fontsize=8.5,fontweight="bold")
        for row,config in enumerate(rows):
            ytop=21.7+row*row_h
            if ri==0:fig.text(1.1/W,1-(ytop+full_h/2)/H,config["label"],ha="left",va="center",fontsize=8.5,
                              fontweight="bold" if row<2 else "normal",fontstyle="normal")
            ax=fig.add_axes([x/W,1-(ytop+full_h)/H,full_w/W,full_h/H])
            az=fig.add_axes([(x+full_w+gap)/W,1-(ytop+(full_h+zoom_w)/2)/H,zoom_w/W,zoom_w/H])
            for a,zoom,lim,name in [(ax,False,limits,"full"),(az,True,r["bounds"],r["roi_id"])]:
                panel(a,he,xy,config["target"],config.get("feature"),lim,r["bounds"],zoom,
                      f"{base.name}/row{row}/{r['roi_id']}/{name}",color=legend_color,
                      dd=config.get("density",{}).get(name),missing_text=config.get("missing"),all_points=config.get("all_points",False))
    lx=187.2/W
    fig.text(lx,.865,"Spot status",fontsize=8,fontweight="bold")
    fig.legend(handles=[Line2D([],[],marker="o",ls="",color=legend_color,markersize=3,label=legend_label),
                        Line2D([],[],marker="o",ls="",color=GREY,markersize=3,label="Other")],
               loc="upper left",bbox_to_anchor=(lx-.015,.85),frameon=False,fontsize=8,handlelength=.5,handletextpad=.3,borderpad=0)
    fig.text(lx,.72,"KDE peak",fontsize=8,fontweight="bold")
    if strict:
        fig.text(lx,.692,"Not shown",fontsize=8)
        fig.text(lx,.67,"Gate pending",fontsize=8)
    else:
        fig.legend(handles=[Line2D([],[],color="#DAA520",lw=.6,label="50%"),Line2D([],[],color="#FF4500",lw=.6,label="85%")],
                   loc="upper left",bbox_to_anchor=(lx-.015,.71),frameon=False,fontsize=8,handlelength=1,borderpad=0,handletextpad=.3)
    fig.text(lx,.575,"Abundance",fontsize=8,fontweight="bold")
    cax=fig.add_axes([lx+.012,.35,.012,.20]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap="viridis"),cax=cax)
    cb.set_ticks([0,1],labels=["Low","High"]);cb.ax.tick_params(labelsize=8,length=0)
    fig.text(lx,.29,"Zoom size",fontsize=8,fontweight="bold")
    fig.text(lx,.265,"1.4 x 1.4 mm",fontsize=8)
    fig.text(lx,.21,f"A: {counts[0]}",fontsize=8)
    fig.text(lx,.184,f"B: {counts[1]}",fontsize=8)
    if strict:
        fig.text(.5,.025,"Orange: KRT14+/PDGFRA+/TP63+ candidates awaiting KLF4. Fib_K14 assignment remains unresolved.",ha="center",fontsize=8)
    else:
        fig.text(.5,.025,"Two-protein reference only. Fib_K14 now requires KRT14, PDGFRA, KLF4 and TP63 to be positive.",ha="center",fontsize=8)
    fig.text(.5,.009,"Nominal scale: 100 µm/grid, inferred from the reported 6.5 x 3.9 mm field and 65 x 39 grids.",ha="center",fontsize=8)
    export(fig,base,len(rows)*4)


def main():
    if (OUT/"STATUS.json").exists() and "--replace-own-draft" not in sys.argv:
        raise FileExistsError("Use a new output version")
    for d in [PLOTS,TABLES,AUDIT,ASSETS,CODE,H5]:d.mkdir(parents=True,exist_ok=True)
    protected=[INPUT,HE,META,AI_IDS,PPT,BASE/"run_config.json",BASE/"tables/all_2029_points_with_two_protein_gate.tsv"]
    source_hashes={str(p):sha(p) for p in protected}
    sp=module("reference",BASE/"code/reference_spatial_batch.py")
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.labelsize":8,
                         "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"figure.titlesize":10,
                         "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","figure.facecolor":"white"})
    # Nominal calibration is derived from the company report, never borrowed from RNA data.
    with ZipFile(PPT) as z:
        title="".join(t.text or "" for t in ET.fromstring(z.read("ppt/slides/slide3.xml")).findall(".//{http://schemas.openxmlformats.org/drawingml/2006/main}t"))
    assert "3.9*6.5mm" in title
    calibration=dict(status="NOMINAL_DERIVED_FROM_REPORTED_FIELD_NOT_INDEPENDENTLY_CALIBRATED",source_pptx=str(PPT),source_slide=3,
                     source_text=title,field_width_um=6500.,field_height_um=3900.,grid_columns=65,grid_rows=39,
                     um_per_grid_x=6500/65,um_per_grid_y=3900/39,full_bar_um=1000,zoom_bar_um=200)
    dump(AUDIT/"scale_calibration_provenance.json",calibration)
    he_original=np.asarray(Image.open(HE).convert("RGB"))
    he=np.round(he_original.astype(float)*HE_WEIGHT+255*(1-HE_WEIGHT)).astype(np.uint8)
    Image.fromarray(he).save(ASSETS/"HE_uniform_background_RGB.png")
    he_settings=dict(source=str(HE),source_sha256=sha(HE),background_asset="assets/HE_uniform_background_RGB.png",
                     precomposed_white_blend_weight=HE_WEIGHT,render_alpha=1.,extent=list(EXTENT),
                     panel_specific_normalization=False,color_adjustment_per_panel=False,orientation="epidermis_up")
    dump(AUDIT/"uniform_HE_display_settings.json",he_settings)

    a=ad.read_h5ad(INPUT)
    xy=np.asarray(a.obsm["spatial"],float)
    assert a.n_obs==2029 and np.array_equal(xy[:,0],np.asarray(a.obsm["spatial_before_xaxis_reflection"])[:,0])
    assert np.allclose(xy[:,1]+np.asarray(a.obsm["spatial_before_xaxis_reflection"])[:,1],38)
    data=pd.read_csv(BASE/"tables/all_2029_points_with_two_protein_gate.tsv",sep="\t")
    assert np.array_equal(xy,data[["x_display","y_display"]].to_numpy())
    rois=json.loads((BASE/"run_config.json").read_text(encoding="utf-8"))["rois"]
    limits=sp.full_limits(xy)
    available=[g for g in GENES if g in a.var_names];missing=[g for g in GENES if g not in a.var_names]
    assert missing==["KLF4"],missing
    raw=dense(a[:,available].layers["raw_intensity"]).astype(float)
    expression=dense(a[:,available].X).astype(float)
    assert np.array_equal(raw>0,expression>0)
    observed_positive=np.all(raw>0,axis=1)
    assert int(observed_positive.sum())==217
    four_status=pd.array([pd.NA if x else False for x in observed_positive],dtype="boolean")
    labels=np.where(observed_positive,"unresolved_missing_KLF4","excluded_by_available_markers")
    raw_map={g:raw[:,i] for i,g in enumerate(available)}
    scaled={"KRT14":data.KRT14_display_scaled.to_numpy(float),"PDGFRA":data.PDGFRA_display_scaled.to_numpy(float)}
    tp63=expression[:,available.index("TP63")];tp63_scale=float(np.quantile(tp63[tp63>0],.995));scaled["TP63"]=np.clip(tp63/tp63_scale,0,1)
    stage=[]
    for g in GENES:
        stage.append(dict(marker=g,available=g in available,n_positive=int((raw_map[g]>0).sum()) if g in available else None,
                          threshold=0,status="measured" if g in available else "absent_from_source_protein_list"))
    pd.DataFrame(stage).to_csv(TABLES/"four_marker_availability.tsv",sep="\t",index=False,na_rep="NA")
    # Verify absence in both the original protein-identification table and AI accessions.
    metadata=pd.read_csv(META,sep="\t");ids=pd.read_csv(AI_IDS,sep="\t")
    hits=metadata.astype(str).apply(lambda c:c.str.contains(r"KLF4|O43474|Krueppel.like factor 4|Kruppel.like factor 4",case=False,regex=True)).any(axis=1)
    accession_hits=ids.astype(str).apply(lambda c:c.str.contains("O43474",case=False,regex=False)).any(axis=1)
    assert not hits.any() and not accession_hits.any()
    dump(AUDIT/"KLF4_source_audit.json",dict(metadata_rows=len(metadata),ai_accession_rows=len(ids),
                                           searched_terms=["KLF4","O43474","Krueppel-like factor 4","Kruppel-like factor 4"],
                                           matching_original_protein_rows=int(hits.sum()),matching_ai_accession_rows=int(accession_hits.sum()),
                                           missing_value_policy="KLF4 is unavailable; no zero substitution or KLF13 substitution"))
    spot=a.obs.copy();spot.insert(0,"spot_id",a.obs_names.astype(str))
    spot["x_display"]=xy[:,0];spot["y_display"]=xy[:,1]
    for g in GENES:spot[f"{g}_input_intensity"]=raw_map.get(g,np.full(a.n_obs,np.nan))
    spot["KLF4_available"]=False
    spot["available_three_positive"]=observed_positive
    spot["Fib_K14"]=four_status
    spot["Fib_K14_status"]=labels
    spot.to_csv(TABLES/"all_spots_four_marker_gate.tsv",sep="\t",index=False,na_rep="NA",float_format="%.10g")
    spot.loc[observed_positive].to_csv(TABLES/"217_three_positive_candidates_pending_KLF4.tsv",sep="\t",index=False,na_rep="NA",float_format="%.10g")
    pd.DataFrame([dict(definition="KRT14>0 AND PDGFRA>0 AND KLF4>0 AND TP63>0",status="NOT_FULLY_EVALUABLE_MISSING_KLF4",
                       n_all=2029,n_strict_Fib_K14=np.nan,n_unresolved_candidates=217,n_excluded_by_available_markers=1812)]).to_csv(TABLES/"strict_Fib_K14_result.tsv",sep="\t",index=False,na_rep="NA")
    # Encode unknown membership explicitly and retain all source variables.
    a.obs["Fib_K14"]=four_status
    a.obs["Fib_K14_four_marker_status"]=pd.Categorical(labels)
    a.obs["Fib_K14_available_three_positive"]=observed_positive
    a.obs["KLF4_positive"]=pd.array([pd.NA]*a.n_obs,dtype="boolean")
    a.uns["Fib_K14_four_marker_gate_20261014"]=dict(required_markers=GENES,available_markers=available,missing_markers=missing,
                                                   definition="KRT14>0 AND PDGFRA>0 AND KLF4>0 AND TP63>0",
                                                   status="NOT_FULLY_EVALUABLE_MISSING_KLF4",n_unresolved=217,n_excluded=1812,
                                                   marker_input="existing raw_intensity layer, AI-reconstructed source abundance",
                                                   missing_KLF4_is_not_zero=True,nullable_membership_encoding="False where a measured marker fails; NA where all measured markers pass")
    h5out=H5/"spatial_proteomics.Fib_K14_four_marker_missing_KLF4_status.h5ad"
    a.write_h5ad(h5out,compression="gzip")

    # Retype and standardize the two previous reference figures.
    two=(data.KRT14_raw_intensity.to_numpy()>0)&(data.PDGFRA_raw_intensity.to_numpy()>0)
    features=[data.coexpression_display_score.to_numpy(float),scaled["KRT14"],scaled["PDGFRA"]]
    density=[]
    for v in features:
        estimates={}
        for name,bounds in [("full",limits)]+[(r["roi_id"],r["bounds"]) for r in rois]:
            estimates[name],_=sp.density(xy[two],v[two],bounds,{"kde_max_points":0},1729)
        density.append(estimates)
    reference_rows=[dict(label="Distribution",target=two)]+[dict(label=name,target=two,feature=v,density=dd) for name,v,dd in zip(["Co-expression","KRT14","PDGFRA"],features,density)]
    refbase=PLOTS/"02_two_protein_reference_HE_Arial8to10_scale"
    compact(he,xy,rois,limits,reference_rows,"Two-protein reference: KRT14 / PDGFRA",
            "KRT14 > 0 AND PDGFRA > 0: 759/2,029 grid points (37.41%)",RED,"Double +",
            [int((two&inside(xy,r["bounds"])).sum()) for r in rois],refbase)

    strict_rows=[dict(label="Distribution",target=observed_positive),
                 dict(label="Co-expression",target=observed_positive,missing="Not evaluable\nKLF4 missing")]
    for g in GENES:
        if g in scaled:strict_rows.append(dict(label=g,target=observed_positive,feature=scaled[g],all_points=True))
        else:strict_rows.append(dict(label=g,target=observed_positive,missing="Unavailable\nin source data"))
    strictbase=PLOTS/"01_Fib_K14_four_marker_AND_HE_Arial8to10_scale"
    compact(he,xy,rois,limits,strict_rows,"Fib_K14: four-protein AND definition",
            "KRT14+/PDGFRA+/KLF4+/TP63+ cannot be fully evaluated: KLF4 unavailable",ORANGE,"Pending KLF4",
            [int((observed_positive&inside(xy,r["bounds"])).sum()) for r in rois],strictbase,strict=True)

    fig=plt.figure(figsize=(210/25.4,98/25.4))
    fig.text(.5,.955,"Two-protein reference localization on uniform H&E",ha="center",fontsize=10,fontweight="bold")
    fig.text(.5,.89,"KRT14+/PDGFRA+: 759/2,029 (37.41%). This reference is distinct from the four-protein Fib_K14 definition.",ha="center",fontsize=8)
    flags=np.c_[data.KRT14_positive.to_numpy(bool),data.PDGFRA_positive.to_numpy(bool)]
    codes=flags[:,0].astype(int)+2*flags[:,1].astype(int)
    colors=["#BBBBBB","#A75DBA","#189FA9",RED];names=["Neither","KRT14 only","PDGFRA only","Double +"]
    for i,title in enumerate(["KRT14","PDGFRA","Detection classes"]):
        ax=fig.add_axes([.035+i*.325,.32,.295,.50]);backdrop(ax,he,limits)
        if i<2:
            vals=scaled[["KRT14","PDGFRA"][i]];order=np.argsort(vals,kind="stable")
            ax.scatter(*xy[order].T,c=vals[order],cmap="viridis",vmin=0,vmax=1,s=3.1,linewidths=0,alpha=.95,zorder=4)
        else:
            for code,color in enumerate(colors):ax.scatter(*xy[codes==code].T,c=color,s=3.1,linewidths=0,alpha=.95,zorder=4)
            for r in rois:roi_box(ax,r["bounds"],r["roi_id"])
        ax.set_title(title,fontsize=9,fontweight="normal",fontstyle="normal",pad=3)
        scale_bar(ax,limits,False,f"overview/{title}")
    cax=fig.add_axes([.075,.18,.52,.025]);cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap="viridis"),cax=cax,orientation="horizontal")
    cb.set_ticks([0,1],labels=["Low","High"]);cb.ax.tick_params(labelsize=8,length=0,pad=2)
    fig.text(.335,.24,"Per-protein relative abundance",ha="center",fontsize=8)
    fig.legend(handles=[Line2D([],[],marker="o",ls="",color=c,markersize=3,label=f"{n}: {int((codes==i).sum())}") for i,(c,n) in enumerate(zip(colors,names))],
               loc="upper left",bbox_to_anchor=(.672,.27),frameon=False,fontsize=8,handlelength=.5,handletextpad=.3,borderpad=0,labelspacing=.3)
    fig.text(.5,.035,"Nominal scale: 100 µm/grid, inferred from the reported 6.5 x 3.9 mm field and 65 x 39 grids.",ha="center",fontsize=8)
    overviewbase=PLOTS/"03_two_protein_overview_HE_Arial8to10_scale"
    export(fig,overviewbase,3)
    with fitz.open() as combined:
        for b in [strictbase,refbase,overviewbase]:
            with fitz.open(b.with_suffix(".pdf")) as doc:combined.insert_pdf(doc)
        combined.save(PLOTS/"Fib_K14_HE_uniform_scale_Arial8to10_combined.pdf",deflate=True)

    # Independent readback and geometry/font audit.
    result=ad.read_h5ad(h5out,backed="r")
    assert result.obs.Fib_K14.isna().sum()==217
    assert result.obs.Fib_K14.eq(False).sum()==1812
    assert result.obs.Fib_K14.eq(True).sum()==0
    assert "KLF4" not in result.var_names
    assert np.array_equal(result.obsm["spatial"],xy)
    assert np.array_equal(dense(result[:,available].X),expression)
    result.file.close()
    clean_bar_records=[{k:v for k,v in r.items() if k not in ["line","ax"]} for r in BAR_RECORDS]
    assert len(clean_bar_records)==43
    # The same data-space lengths must have the same displayed physical length within a panel class.
    for length in [1000.,200.]:
        records=[r for r in clean_bar_records if r["length_um"]==length]
        assert all(r["drawn_mm"]>0 for r in records)
    assert source_hashes=={str(p):sha(p) for p in protected}
    dump(AUDIT/"scale_bar_geometry.json",clean_bar_records)
    dump(AUDIT/"figure_font_background_QA.json",FIGURE_RECORDS)
    dump(AUDIT/"input_SHA256.json",[{"path":p,"sha256":s,"unchanged":True} for p,s in source_hashes.items()])
    summary=dict(status="FIGURES_COMPLETE_GATE_LIMITED_BY_MISSING_KLF4",n_all=2029,n_strict_Fib_K14=None,
                 missing_markers=missing,n_known_three_positive_candidates=217,n_excluded_by_available_markers=1812,
                 n_KRT14_positive=873,n_PDGFRA_positive=1894,n_TP63_positive=672,
                 figure_count=3,combined_pdf_pages=3,HE_spatial_panels=43,scale_bars=43,
                 font="Arial",font_min_pt=8,font_max_pt=10,protein_symbols="upright",uniform_HE_asset=True,
                 calibration_status=calibration["status"],um_per_grid=100.,source_unchanged=True)
    dump(AUDIT/"FINAL_QA.json",summary)
    dump(OUT/"STATUS.json",summary)
    (OUT/"README_ZH.md").write_text(
        "# Fib_K14四蛋白定义与统一HE图版\n\n"
        "用户定义：KRT14>0 且 PDGFRA>0 且 KLF4>0 且 TP63>0。\n\n"
        "## 检测结果\n\n"
        "KRT14、PDGFRA、TP63分别在873、1894、672个spot检出。当前H5AD、原始5852行蛋白鉴定表、"
        "AI accession列表均不存在KLF4/O43474。不能将缺失KLF4替换成0或KLF13。\n"
        "217个spot满足KRT14+/PDGFRA+/TP63+，列为等待KLF4验证的候选。1812个spot至少有一个已测marker不满足>0，"
        "可判为不满足四阳性。完整四阳性数量为NA，不应报告为0。\n"
        "新H5AD的obs['Fib_K14']为nullable Boolean：上述217个为NA，其余1812个为False；不把217个候选命名为已确认Fib_K14。"
        "obs['Fib_K14_four_marker_status']保留原因。未增加伪造的KLF4表达列，也未改动原表达矩阵或已有注释。\n\n"
        "## 图形设置\n\n"
        "所有43个空间面板使用同一HE像素副本、同一方向/范围和同一白色混合权重0.72。"
        "文字全部Arial，8–10pt。标题与面板分组加粗。此处检测对象为蛋白，KRT14、PDGFRA、KLF4、TP63用正体。"
        "采用HGNC关于基因/蛋白符号区分的建议：https://www.genenames.org/about/guidelines/ 。\n"
        "主图保留固定ROI A/B。四蛋白图中的橙色仅表示三阳性候选。四蛋白Co-expression与KLF4面板明确标注无法计算/缺失；"
        "其他单蛋白行显示全部spot的相对丰度。旧两蛋白图保留为reference，不再作为新Fib_K14四阳性结果。\n\n"
        "## 比例尺\n\n"
        "公司PPT第3页记录视场3.9×6.5mm，当前完整坐标网格39×65，因此按视场外缘推算100µm/grid。"
        "全部全景比例尺1000µm，放大图200µm。该比例尺为公司报告视场推算的名义标定，尚无独立扫描仪/采样标定文件，"
        "图注已保留nominal。来源、推导和每根比例尺的实际图面长度见audit/scale_calibration_provenance.json与scale_bar_geometry.json。\n\n"
        "输出PDF/SVG含矢量文字、散点、ROI与比例尺，HE为嵌入图像；PNG/TIFF为600dpi RGB。\n",
        encoding="utf-8")
    cards=[]
    for b,title in [(strictbase,"四蛋白定义：等待KLF4验证"),(refbase,"更新样式的两蛋白参考图"),(overviewbase,"更新样式的两蛋白全景")]:
        rel="plots/"+b.name
        cards.append(f'<h2>{title}</h2><p>'+" · ".join(f'<a href="{rel}.{ext}">{ext.upper()}</a>' for ext in ["pdf","png","svg","tiff"])+f'</p><img src="{rel}_preview.png">')
    (OUT/"OPEN_ME_FIRST.html").write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Fib_K14四蛋白HE图版</title><style>body{font:16px Arial,"Microsoft YaHei",sans-serif;max-width:1450px;margin:24px auto;padding:0 20px}p{line-height:1.7}img{max-width:100%}a{color:#176590}</style><h1>Fib_K14四蛋白AND门控与统一HE背景</h1><p>当前缺少KLF4，严格四阳性总数为NA。217个三阳性spot列为待验证候选。新obs["Fib_K14"]保留NA/False。</p><p>图形全部Arial 8–10pt、统一HE；比例尺按公司报告视场推算100µm/grid，为名义标定。</p><p><a href="plots/Fib_K14_HE_uniform_scale_Arial8to10_combined.pdf">合并PDF</a> · <a href="tables/217_three_positive_candidates_pending_KLF4.tsv">217个候选坐标</a> · <a href="tables/four_marker_availability.tsv">四蛋白可用性</a> · <a href="README_ZH.md">完整说明</a></p>'+''.join(cards)+'</html>',encoding="utf-8")
    shutil.copy2(Path(__file__),CODE/Path(__file__).name)
    shutil.copy2(BASE/"code/reference_spatial_batch.py",CODE/"reference_spatial_batch.py")
    files=[p for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="SHA256.tsv"]
    pd.DataFrame([dict(relative_path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files]).to_csv(OUT/"SHA256.tsv",sep="\t",index=False)
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":main()
