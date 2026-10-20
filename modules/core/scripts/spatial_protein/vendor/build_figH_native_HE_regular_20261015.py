"""Figure H revision: source-native grid-free H&E, regular fonts, first-column scale labels."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
from zipfile import ZipFile

import anndata as ad
import cv2
import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(r"data/external/local_e63c701db393")
OLD = ROOT / "17_FigF_Fib_K14_Xgt2_epithelial_front_20261015_v01"
BASE_CODE = ROOT / "16_FigF_Fib_K14_four_protein_Xgt1_single_row_20261015_v01/code/build_figF_proteomics_four_Xgt1_20261015.py"
HE_OLD = ROOT / "14_Fib_K14_four_marker_HE_Arial8to10_20261014_v01/assets/HE_uniform_background_RGB.png"
HE_UNBLENDED_OLD = ROOT / "13_Fib_K14_HE_background_20261014_v01/assets/HE_epidermis_up_photo_band.png"
PPT = Path(r"data/external/local_8853af99994d")
USER_IMAGE = Path(r"data/external/local_befeb793c30d")
OUT = ROOT / "18_FigH_Fib_K14_Xgt2_HE_clarity_regular_font_20261015_v01"
PLOTS, TABLES, AUDIT, ASSETS, CODE = [OUT / d for d in ("plots", "tables", "audit", "assets", "code")]
HE_SOURCE = ASSETS / "source_ppt_media/image8.png"
PHOTO_CROP = [88, 252, 759, 710]  # x0,y0,x1,y1, no drawn chart/green outer frame
ROI = [33.5,47.5,25.5,39.5]
EXTENT = [-.5,64.5,-.5,38.5]
TITLES = ["Distribution","Co-expression","PDGFRA","KRT14","KRT5","TACSTD2"]
GENES = ["KRT14","PDGFRA","KRT5","TACSTD2"]
BAR_QA, PANEL_QA, TEXT_QA = [], [], []
POINT_AREA_FACTOR = .45


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


def source_identity(st):
    old=np.asarray(Image.open(HE_UNBLENDED_OLD).convert("RGB"))
    native=np.asarray(Image.open(HE_SOURCE).convert("RGB"))
    x0,y0,x1,y1=PHOTO_CROP
    clean=native[y0:y1,x0:x1][::-1,::-1].copy()
    assert clean.shape==(458,671,3)
    # Resizing below is identity/registration QA ONLY, not a new render source.
    comparable=cv2.resize(clean,(950,572),interpolation=cv2.INTER_CUBIC)
    ga=cv2.cvtColor(old,cv2.COLOR_RGB2GRAY);gb=cv2.cvtColor(comparable,cv2.COLOR_RGB2GRAY)
    sift=cv2.SIFT_create()
    k1,d1=sift.detectAndCompute(ga,None);k2,d2=sift.detectAndCompute(gb,None)
    good=[r for r,s in cv2.BFMatcher().knnMatch(d1,d2,k=2) if r.distance<.7*s.distance]
    x=np.float32([k1[r.queryIdx].pt for r in good]);y=np.float32([k2[r.trainIdx].pt for r in good])
    matrix,inliers=cv2.estimateAffinePartial2D(x,y,method=cv2.RANSAC,ransacReprojThreshold=4)
    use=inliers.ravel()==1
    median=float(np.median(np.linalg.norm(x-y,axis=1)[use]))
    assert len(good)>=100 and int(use.sum())>=100 and median<4
    result=dict(source_pptx=str(PPT),slide=3,embedded_media="pptexternal_asset_00f13cb2921e",embedded_size_px=[1065,750],
                photo_crop_xyxy_px=PHOTO_CROP,native_photo_size_px=[671,458],rotation_degrees=180,
                grayscale_correlation_after_QA_only_resize=float(np.corrcoef(ga.ravel(),gb.ravel())[0,1]),
                SIFT_matches=len(good),RANSAC_inliers=int(use.sum()),median_unwarped_landmark_offset_at_950x572_px=median,
                fitted_similarity_for_QA_only=matrix.tolist(),fitted_transform_applied=False,
                interpretation="Same source tissue photograph, extracted before green sampling-grid overlays and before inpainting. Existing coarse spatial footprint/extent retained, not new cell-level registration.",
                processing="Native rectangular photo crop and 180-degree array reorientation only. No sharpening, generative enhancement, super-resolution, inpainting, local recoloring or tissue reconstruction.",
                real_resolution_limit="The 14-grid ROI has approximately 145 x 164 native pixels, including a small blank surface margin. No full-resolution whole-slide scan supplied.")
    st.dump(AUDIT/"HE_native_source_identity_and_resolution.json",result)
    Image.fromarray(clean).save(ASSETS/"HE_native_grid_free_epidermis_up.png")
    return clean,result


def background(ax,he,limits,zoom):
    ax.imshow(he,extent=EXTENT,origin="upper",interpolation="none",resample=False,alpha=1,zorder=0)
    ax.set(xlim=limits[:2],ylim=limits[2:],aspect="equal",xticks=[],yticks=[])
    for s in ax.spines.values():s.set_visible(zoom);s.set_linewidth(.55)


def panel(st,ax,he,xy,target,feature,limits,zoom,kind,dd,pid,factor):
    background(ax,he,limits,zoom)
    keep=st.inside(xy,limits)
    size=(7.2 if zoom else 2.8)*factor**2*POINT_AREA_FACTOR
    other_size=(3.6 if zoom else 1.4)*factor**2*POINT_AREA_FACTOR
    if kind=="distribution":
        ax.scatter(*xy[keep&~target].T,c="#DADADA",s=other_size,alpha=.14,linewidths=0,zorder=1)
        ax.scatter(*xy[keep&target].T,c="#FF0000",s=size,alpha=.95,linewidths=0,zorder=4)
    else:
        ids=np.flatnonzero(keep if kind=="protein" else keep&target)
        ids=ids[np.argsort(feature[ids],kind="stable")]
        ax.scatter(*xy[ids].T,c=feature[ids],cmap="viridis",vmin=0,vmax=1,s=size,alpha=.90,linewidths=0,zorder=4)
        if dd is not None:
            gx,gy,z=dd
            ax.contour(gx,gy,z,levels=[.5*z.max(),.85*z.max()],colors=["#DAA520","#FF4500"],linewidths=.25,zorder=6)
    if not zoom:
        ax.add_patch(Rectangle((ROI[0],ROI[2]),14,14,fill=False,ec="black",lw=.45,ls=(0,(3,1.5)),zorder=10,clip_on=False))
    PANEL_QA.append(dict(panel_id=pid,limits=list(limits),n_visible=int(keep.sum()),n_Fib_K14_visible=int((keep&target).sum()),
                         HE_asset="assets/HE_native_grid_free_epidermis_up.png",HE_alpha=1,HE_white_blend=False,
                         HE_interpolation="none",HE_resample=False,HE_extent=EXTENT,point_area_pt2=size,point_area_factor=.45))


def scale(st,fig,ax,x,width,height,limits,zoom,pid,show_label,scale_y=24.2,label_y=25.):
    um=200. if zoom else 1000.
    units=um/100.
    physical_mm=width*units/(limits[1]-limits[0])
    left=(x+width-physical_mm) if zoom else x
    right=left+physical_mm
    line=Line2D([left/210,right/210],[1-scale_y/height]*2,transform=fig.transFigure,color="black",lw=.8)
    fig.add_artist(line)
    if show_label:
        fig.text((x+width if zoom else x)/210,1-label_y/height,f"{um:g} µm",ha="right" if zoom else "left",va="top",fontsize=8,fontweight="normal")
    st.BAR_RECORDS.append(dict(panel_id=pid,length_um=um,data_length=units,line=line,ax=ax,expected_drawn_mm=physical_mm))
    BAR_QA.append(dict(panel_id=pid,length_um=um,expected_drawn_mm=physical_mm,numeric_label_visible=show_label))


def font_and_overlap_QA(fig,stem):
    fig.canvas.draw()
    all_text=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    assert all(t.get_fontweight() in ["normal",400] for t in all_text)
    assert all(t.get_fontfamily()==["Arial"] and t.get_fontstyle()=="normal" and 8<=t.get_fontsize()<=8.5 for t in all_text)
    boxes=[(t.get_text(),t.get_window_extent(fig.canvas.get_renderer())) for t in all_text]
    overlaps=[]
    for i,(ta,ba) in enumerate(boxes):
        for tb,bb in boxes[i+1:]:
            v=matplotlib.transforms.Bbox.intersection(ba,bb)
            if v is not None and v.width>.5 and v.height>.5:overlaps.append([ta,tb])
    assert not overlaps,overlaps
    TEXT_QA.append(dict(stem=stem,all_visible_text_regular=True,no_text_overlaps=True,
                        numeric_bar_label_counts={v:sum(t.get_text()==v for t in all_text) for v in ["1000 µm","200 µm"]}))


def render(base,st,he,xy,target,features,densities,limits,standalone):
    geo=base.GEOMETRY;W=210.;H=49. if standalone else 30.
    fig=plt.figure(figsize=(W/25.4,H/25.4))
    group=(W-geo["left_mm"]-geo["right_mm"]-5*geo["column_gap_mm"])/6
    zoom_w=(group-geo["full_zoom_gap_mm"])/4;full_w=3*zoom_w
    full_h=full_w*(limits[3]-limits[2])/(limits[1]-limits[0]);factor=full_w/59.4
    fig.text(2.5/W,1-(geo["image_top_mm"]+full_h/2)/H,"Proteomics",rotation=90,fontsize=8,ha="center",va="center",fontweight="normal")
    stem="FigH_Fib_K14_Xgt2_nativeHE_regular_210mm_"+("with_legend" if standalone else "for_assembly")
    for col,title in enumerate(TITLES):
        x=geo["left_mm"]+col*(group+geo["column_gap_mm"])
        fig.text((x+group/2)/W,1-geo["title_y_mm"]/H,title,ha="center",va="center",fontsize=8.5,fontweight="normal",fontstyle="normal")
        ax=fig.add_axes([x/W,1-(geo["image_top_mm"]+full_h)/H,full_w/W,full_h/H])
        zx=x+full_w+geo["full_zoom_gap_mm"];ztop=geo["image_top_mm"]+(full_h-zoom_w)/2
        az=fig.add_axes([zx/W,1-(ztop+zoom_w)/H,zoom_w/W,zoom_w/H])
        kind="distribution" if col==0 else ("coexpression" if col==1 else "protein")
        feature=None if col==0 else features[title]
        for current,zoom,lim,px,pw,suffix in [(ax,False,limits,x,full_w,"full"),(az,True,ROI,zx,zoom_w,"zoom")]:
            pid=f"{stem}/{title}/{suffix}"
            panel(st,current,he,xy,target,feature,lim,zoom,kind,densities.get((title,suffix)),pid,factor)
            scale(st,fig,current,px,pw,H,lim,zoom,pid,col==0)
    if standalone:
        base.draw_legend(fig,H)
        fig.text(.5,1-41/H,"Fib_K14: 657 / 2,029 spots (32.38%); KRT14, PDGFRA, KRT5 and TACSTD2: X > 2.",ha="center",va="center",fontsize=8,fontweight="normal")
        fig.text(.5,1-45.5/H,"Native grid-free HE; regular Arial text. Same front ROI and nominal scale; no synthetic tissue details.",ha="center",va="center",fontsize=8,fontweight="normal")
    font_and_overlap_QA(fig,stem)
    st.export(fig,PLOTS/stem,12)
    return stem


def morphology_reference(st,he,limits):
    """Unoccluded source reference, without inventing a histological boundary/leading edge."""
    W,H=210.,98.
    fig=plt.figure(figsize=(W/25.4,H/25.4))
    fig.text(.5,1-6/H,"Native HE reference: whole section and user-marked epithelialization-front ROI",ha="center",va="center",fontsize=8.5,fontweight="normal")
    full_w=116.;full_h=full_w*(limits[3]-limits[2])/(limits[1]-limits[0])
    ax=fig.add_axes([5/W,1-(17+full_h)/H,full_w/W,full_h/H])
    az=fig.add_axes([132/W,1-(17+66)/H,66/W,66/H])
    background(ax,he,limits,False);background(az,he,ROI,True)
    ax.add_patch(Rectangle((ROI[0],ROI[2]),14,14,fill=False,ec="black",lw=.6,ls=(0,(3,1.5)),clip_on=False))
    fig.text((5+full_w/2)/W,1-13/H,"Whole section",ha="center",fontsize=8.5,fontweight="normal")
    fig.text(165/W,1-13/H,"Same front ROI; no protein dots",ha="center",fontsize=8.5,fontweight="normal")
    # Broad tissue-layer labels are deliberately omitted: no high-resolution scan or pathologist segmentation was supplied.
    fig.text(.5,1-90/H,"Source photo: 671 x 458 pixels. No sharpening or AI reconstruction; diagnostic layer/edge delineation needs the original scan.",ha="center",fontsize=8,fontweight="normal")
    fig.text(.5,1-95/H,"Existing coarse HE-to-grid mapping retained. The marked front is user-specified, not a newly validated pathological boundary.",ha="center",fontsize=8,fontweight="normal")
    stem="HE_native_whole_section_front_ROI_reference"
    font_and_overlap_QA(fig,stem)
    st.export(fig,PLOTS/stem,2)
    return stem


def main():
    if (OUT/"STATUS.json").exists():raise FileExistsError("Completed output exists; use a new version")
    for p in [PLOTS,TABLES,AUDIT,ASSETS,CODE]:p.mkdir(parents=True,exist_ok=True)
    base=module("previous_FigF",BASE_CODE);st=module("HE_style",base.STYLE);ref=module("reference_spatial",base.REFERENCE)
    protected=[base.INPUT,PPT,HE_SOURCE,HE_OLD,HE_UNBLENDED_OLD,BASE_CODE,OLD/"run_config.json",OLD/"tables/all_2029_spots_Fib_K14_Xgt2.tsv"]
    hashes={str(p):st.sha(p) for p in protected}
    with ZipFile(PPT) as z:assert z.read("pptexternal_asset_00f13cb2921e")==HE_SOURCE.read_bytes()
    he,identity=source_identity(st)
    source_table=OLD/"tables/all_2029_spots_Fib_K14_Xgt2.tsv"
    table=pd.read_csv(source_table,sep="\t")
    a=ad.read_h5ad(base.INPUT)
    assert a.shape==(2029,5846) and list(a.obs_names.astype(str))==list(table.spot_id)
    xy=np.asarray(a.obsm["spatial"],float)
    assert np.array_equal(xy,table[["x_display","y_display"]].to_numpy(float))
    values=base.dense(a[:,GENES].X).astype(float)
    target=(values>2).all(axis=1)
    assert np.array_equal(target,table.Fib_K14) and target.sum()==657
    in_roi=st.inside(xy,ROI)
    assert in_roi.sum()==174 and (target&in_roi).sum()==72
    features={g:table[g+"_display_scaled"].to_numpy(float) for g in GENES}
    features["Co-expression"]=table.coexpression_display_score.to_numpy(float)
    limits=ref.full_limits(xy)
    cache,density_audit={},[]
    for feature in ["Co-expression",*base.DISPLAY_GENES]:
        for region,bounds in [("full",limits),("zoom",ROI)]:
            dd,info=ref.density(xy[target],features[feature][target],bounds,{"kde_max_points":0},1729)
            cache[(feature,region)]=dd;density_audit.append(dict(feature=feature,region=region,**info))
    plt.rcParams.update({"font.family":"Arial","font.size":8,"font.weight":"normal","axes.titleweight":"normal","axes.labelweight":"normal",
                         "axes.titlesize":8.5,"axes.labelsize":8,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,
                         "figure.titlesize":8.5,"figure.titleweight":"normal","pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","figure.facecolor":"white"})
    stems=[render(base,st,he,xy,target,features,cache,limits,v) for v in [False,True]]
    reference_stem=morphology_reference(st,he,limits)
    for stem in stems+[reference_stem]:
        with fitz.open(PLOTS/(stem+".pdf")) as doc:
            spans=[s for block in doc[0].get_text("dict")["blocks"] if "lines" in block for ln in block["lines"] for s in ln["spans"] if s["text"].strip()]
            assert all(s["font"]=="ArialMT" for s in spans)
            if stem in stems:
                words=doc[0].get_text()
                assert words.count("1000 µm")==1 and words.count("200 µm")==1
    bars=[{k:v for k,v in q.items() if k not in ("ax","line")} for q in st.BAR_RECORDS]
    assert len(bars)==24 and sum(q["numeric_label_visible"] for q in BAR_QA)==4
    assert all(np.isclose(q["drawn_mm"],q["expected_drawn_mm"]) for q in bars)
    pd.DataFrame(bars).to_csv(AUDIT/"scale_bar_geometry.tsv",sep="\t",index=False)
    st.dump(AUDIT/"scale_labels_first_column_only.json",BAR_QA)
    st.dump(AUDIT/"text_regular_weight_QA.json",TEXT_QA)
    st.dump(AUDIT/"panel_QA.json",PANEL_QA)
    st.dump(AUDIT/"density_audit.json",density_audit)
    st.dump(AUDIT/"source_sha256.json",hashes)
    for name in ["all_2029_spots_Fib_K14_Xgt2.tsv","Fib_K14_657_spots_Xgt2_coordinates.tsv","epithelial_front_174_spots.tsv","epithelial_front_72_Fib_K14_spots.tsv","shared_epithelial_front_ROI.tsv","four_protein_gate_and_fixed_display_scales.tsv"]:
        shutil.copy2(OLD/"tables"/name,TABLES/name)
    shutil.copy2(OLD/"audit/scale_calibration_provenance.json",AUDIT/"scale_calibration_provenance.json")
    settings=dict(source_photo=str(HE_SOURCE),crop_xyxy=PHOTO_CROP,rotation=180,asset="assets/HE_native_grid_free_epidermis_up.png",native_size_px=[671,458],
                  extent=EXTENT,HE_alpha=1,HE_white_blend=False,interpolation="none",resample=False,uniform_all_panels=True,
                  original_previous_alpha_equivalent=.72,change="Replace grid-inpainted and white-mixed image with source-native grid-free photo; no image-detail synthesis",
                  overlay_point_area_factor=.45,overlay_point_alpha=.90,density_linewidth_pt=.25,
                  unprovided_original_scan=True,biological_pixel_resolution_not_increased=True)
    st.dump(AUDIT/"uniform_HE_display_settings.json",settings)
    assert hashes=={str(p):st.sha(p) for p in protected}
    st.dump(AUDIT/"figure_QA.json",dict(figures=st.FIGURE_RECORDS,all_text_regular_weight=True,fonts=["ArialMT"],font_range_pt=[8,8.5],
              assembly_canvas_mm=[210,30],all12_HEnative_panels=True,scale_bars_per_row=12,numeric_scale_labels_per_row=2,
              numeric_labels_only_in_Distribution_column=True,source_molecular_values_unchanged=True,n_Fib_K14=657,
              shared_front_ROI=ROI,ROI_spots=174,ROI_Fib_K14=72,source_HE_same_tissue_verified=True,no_generative_enhancement=True,
              true_optical_resolution_increased=False,original_sources_unchanged=True))
    st.dump(OUT/"run_config.json",dict(previous_output=str(OLD),input_h5ad=str(base.INPUT),gate="KRT14_X>2 AND PDGFRA_X>2 AND KRT5_X>2 AND TACSTD2_X>2",
              n_Fib_K14=657,expression_source="X_total_intensity_normalized_log1p",ROI=ROI,geometry=base.GEOMETRY,
              preserved=["657spotgate","frontROI","expressionvalues","displayscales","HEextent","210x30mm","600dpi"],
              requested_changes=["source-native grid-free undiluted HE","all fonts regular","only Distribution column has numeric full/zoom scale labels"],
              display_settings=settings,no_panel_letter="H label belongs to the user's composite and is not duplicated in the strip"))
    shutil.copy2(Path(__file__),CODE/Path(__file__).name);shutil.copy2(BASE_CODE,CODE/BASE_CODE.name)
    (OUT/"README_结果说明.md").write_text(
        "# H面板修改：原生无网格HE / 常规字重 / 首列数字标尺\n\n"
        "保持四蛋白X均>2、657/2,029个spot及前沿ROI不变；ROI内仍174个spot/72个Fib_K14。未重定门控或归一化。\n\n"
        "主要变更：\n\n"
        "1. HE取自公司报告第3页内嵌image8.png中同一组织的无网格照片，裁剪框[88,252,759,710]，原生671×458像素。"
        "仅裁剪并旋转180度，不做锐化、超分辨率、生成式增强或修补。使用原色alpha=1，无白色淡化，不再使用旧去网格修补图。"
        "现有坐标映射不变：x=-0.5..64.5,y=-0.5..38.5。SIFT对应关系与残差见audit/HE_native_source_identity_and_resolution.json。"
        "散点面积降至旧版45%，表达点alpha=.90，减少对组织结构遮挡；分子值及色阶不变。\n\n"
        "2. 全部Arial常规字重、正体，8–8.5pt；标题Distribution/Co-expression不再粗体。实际PDF字体仅ArialMT，已核查无Bold/Italic。\n\n"
        "3. 仅首列Distribution的全景/放大分别标1000µm/200µm，其余十个图保留对应长度标尺，不标数字。"
        "每行仍有12根几何长度正确的标尺，名义100µm/grid来源不变。\n\n"
        "拼版文件*_for_assembly为210×30mm；带图例版为210×49mm；均有PDF/SVG和600dpi RGB PNG/LZW TIFF。"
        "未重复添加H字母，便于直接放入已有拼图。另附HE_native_whole_section_front_ROI_reference，不加蛋白点供检查原生组织可见结构。\n\n"
        "重要限制：现有照片不是高清病理扫描；前沿14网格窗口约145×164个原生像素，部分为空白表面边缘。"
        "本次改善的是无网格、无修补、无淡化的可辨性，不是增加真实光学分辨率。"
        "若需可靠区分表皮-真皮边界、判断真正上皮化最前沿或细胞细节，请提供同一切片的原始TIFF/SVS/KFB等扫描。"
        "未凭低分辨率背景绘制未经验证的组织分割线或伪造表皮细节。未用其他病例/切片的高分辨率图片替换。\n\n"
        "本数据为公司AI重建空间蛋白网格，同spot四蛋白丰度不等同同一细胞共表达或已证实成纤维细胞谱系。"
        "原图、原数据及旧版输出均保留。\n",encoding="utf-8")
    files=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=st.sha(p)) for p in OUT.rglob("*") if p.is_file()]
    st.dump(OUT/"STATUS.json",dict(status="COMPLETE_DISPLAY_OPTIMIZATION_WITH_SOURCE_RESOLUTION_LIMIT",QA="PASS",n_Fib_K14=657,
             original_scan_needed_for_pathology_grade_details=True,files=files))
    print(json.dumps(dict(output=str(OUT),n_Fib_K14=657,HE_native_px=[671,458],same_source_identity=identity,
                          figures=st.FIGURE_RECORDS,regular_fonts_only=True,scale_numeric_labels_first_column_only=True),ensure_ascii=False))


if __name__=="__main__":main()
