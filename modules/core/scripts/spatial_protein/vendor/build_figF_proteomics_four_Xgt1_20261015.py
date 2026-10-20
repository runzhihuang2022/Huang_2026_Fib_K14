"""One-row, 210-mm Figure F addition using the user-confirmed four-protein X>1 gate."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil

import anndata as ad
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
INPUT = ROOT / "10_xaxis_reflected_HE_merged_atlas/h5ad/spatial_proteomics.xaxis_reflected_epidermis_up.HE_aligned.h5ad"
PREV = ROOT / "14_Fib_K14_four_marker_HE_Arial8to10_20261014_v01"
HE = PREV / "assets/HE_uniform_background_RGB.png"
STYLE = PREV / "code/revise_fib_k14_four_marker_he_20261014.py"
REFERENCE = ROOT / "12_Fib_K14_KRT14_PDGFRA_same_position_20261014_v01/code/reference_spatial_batch.py"
REF_PDF = Path(r"data/external/local_0001cf51a943")
OUT = ROOT / "16_FigF_Fib_K14_four_protein_Xgt1_single_row_20261015_v01"
PLOTS, TABLES, AUDIT, ASSETS, CODE = [OUT / d for d in ("plots", "tables", "audit", "assets", "code")]
GENES = ["KRT14", "PDGFRA", "KRT5", "TACSTD2"]
DISPLAY_GENES = ["PDGFRA", "KRT14", "KRT5", "TACSTD2"]
TITLES = ["Distribution", "Co-expression", *DISPLAY_GENES]
THRESHOLD = 1.0
RED, GREY = "#FF0000", "#DADADA"
DENSITY_COLORS = ["#DAA520", "#FF4500"]
GEOMETRY = dict(width_mm=210., strip_height_mm=30., annotated_height_mm=49.,
                left_mm=5.5, right_mm=3.5, column_gap_mm=1.6, full_zoom_gap_mm=.6,
                full_to_zoom_ratio=3., image_top_mm=8.4, scale_y_mm=24.2,
                scale_label_y_mm=25.0, title_y_mm=4.7)
PANEL_QA = []


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dense(m):
    return m.toarray() if hasattr(m, "toarray") else np.asarray(m)


def select_roi(xy, target, inside):
    """As in the prior analysis: max positive count, purity, centroid distance, then xy."""
    width = 14
    centroid = xy[target].mean(axis=0)
    candidates = []
    for gx in range(65 - width + 1):
        for gy in range(39 - width + 1):
            bounds = [gx - .5, gx + width - .5, gy - .5, gy + width - .5]
            keep = inside(xy, bounds)
            n, total = int((target & keep).sum()), int(keep.sum())
            if total:
                center = np.array([gx + (width - 1) / 2, gy + (width - 1) / 2])
                rank = (-n, -n / total, float(((center - centroid) ** 2).sum()), gx, gy)
                candidates.append((rank, bounds, n, total))
    candidates.sort(key=lambda v: v[0])
    _, bounds, n, total = candidates[0]
    return dict(roi_id="A", bounds=bounds, n_Fib_K14=n, n_spots=total,
                purity=n / total, fraction_of_all_target=n / int(target.sum()),
                selection="14x14 grids; maximize X>1 four-positive count, then purity, then target-centroid distance, then lexicographic xy")


def scatter_panel(st, ax, he, xy, target, values, limits, roi, zoom, kind, dd, panel_id, factor):
    st.backdrop(ax, he, limits, zoom)
    visible = st.inside(xy, limits)
    # Preserve the previous full/zoom dot dimensions relative to each map width.
    s_target = (7.2 if zoom else 2.8) * factor ** 2
    s_other = (3.6 if zoom else 1.4) * factor ** 2
    if kind == "distribution":
        ax.scatter(*xy[visible & ~target].T, c=GREY, s=s_other, alpha=.35, linewidths=0, zorder=1)
        ax.scatter(*xy[visible & target].T, c=RED, s=s_target, alpha=.95, linewidths=0, zorder=4)
    else:
        if kind == "coexpression":
            ax.scatter(*xy[visible & ~target].T, c=GREY, s=s_other, alpha=.16, linewidths=0, zorder=1)
            ids = np.flatnonzero(visible & target)
        else:
            # Individual protein maps show all observed grid values, not just gate-positive spots.
            ids = np.flatnonzero(visible)
        ids = ids[np.argsort(values[ids], kind="stable")]
        ax.scatter(*xy[ids].T, c=values[ids], s=s_target, cmap="viridis", vmin=0, vmax=1,
                   alpha=.97, linewidths=0, zorder=4)
        if dd is not None:
            gx, gy, z = dd
            ax.contour(gx, gy, z, levels=[.50 * z.max(), .85 * z.max()], colors=DENSITY_COLORS,
                       linewidths=.35, zorder=6)
    if not zoom:
        ax.add_patch(Rectangle((roi[0], roi[2]), roi[1]-roi[0], roi[3]-roi[2],
                               fill=False, lw=.45, ec="black", ls=(0, (3, 1.5)), zorder=10))
    PANEL_QA.append(dict(panel_id=panel_id,limits=list(limits),n_visible=int(visible.sum()),
                         n_target_visible=int((visible & target).sum()),map_kind=kind,
                         HE_sha256=st.sha(HE),HE_alpha=1.,HE_extent=list(st.EXTENT),
                         dot_size_pt2=s_target,protein_values="all_spots" if kind=="protein" else "gate_positive_only"))


def add_scale(st, fig, ax, x, width, y, h, limits, zoom, panel_id):
    """Put real-size bars below each tiny map so 8-pt labels never cover H&E."""
    length_um = 200. if zoom else 1000.
    units = length_um / st.GRID_UM
    line_mm = width * units / (limits[1] - limits[0])
    right = x + width
    if zoom:
        a, b = right - line_mm, right
        text_x, align = right, "right"
    else:
        a, b = x, x + line_mm
        text_x, align = x, "left"
    line = Line2D([a / GEOMETRY["width_mm"], b / GEOMETRY["width_mm"]], [1-y/h]*2,
                  transform=fig.transFigure, lw=.8, color="black")
    fig.add_artist(line)
    fig.text(text_x/GEOMETRY["width_mm"],1-GEOMETRY["scale_label_y_mm"]/h,
             f"{length_um:g} µm", ha=align, va="top", fontsize=8)
    st.BAR_RECORDS.append(dict(panel_id=panel_id,length_um=length_um,data_length=units,
                              line=line,ax=ax,expected_drawn_mm=line_mm))


def draw_legend(fig, height):
    W = GEOMETRY["width_mm"]
    y = 34.
    entries = [(15., RED, "Fib_K14"), (45., GREY, "Other spots")]
    for x, color, label in entries:
        fig.add_artist(Line2D([x/W], [1-y/height], transform=fig.transFigure,
                              ls="", marker="o", markersize=3, color=color))
        fig.text((x+2.5)/W,1-y/height,label,va="center",fontsize=8)
    for x, color, label in [(77.,DENSITY_COLORS[0],"Density 50%"),(112.,DENSITY_COLORS[1],"Density 85%")]:
        fig.add_artist(Line2D([x/W,(x+4.)/W],[1-y/height]*2,transform=fig.transFigure,lw=.6,color=color))
        fig.text((x+5.5)/W,1-y/height,label,va="center",fontsize=8)
    fig.text(178/W,1-31.9/height,"Relative abundance",ha="center",fontsize=8)
    cax=fig.add_axes([160/W,1-36.0/height,37/W,2./height])
    cb=fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap="viridis"),cax=cax,orientation="horizontal")
    cb.set_ticks([])
    fig.text(158.3/W,1-35.0/height,"Low",ha="right",va="center",fontsize=8)
    fig.text(198.7/W,1-35.0/height,"High",ha="left",va="center",fontsize=8)


def render(st, he, xy, target, features, density, limits, roi, standalone=False):
    W=GEOMETRY["width_mm"]
    H=GEOMETRY["annotated_height_mm"] if standalone else GEOMETRY["strip_height_mm"]
    fig=plt.figure(figsize=(W/25.4,H/25.4))
    group_w=(W-GEOMETRY["left_mm"]-GEOMETRY["right_mm"]-5*GEOMETRY["column_gap_mm"])/6
    zoom_w=(group_w-GEOMETRY["full_zoom_gap_mm"])/4
    full_w=3*zoom_w
    full_h=full_w*(limits[3]-limits[2])/(limits[1]-limits[0])
    factor=full_w/59.4
    fig.text(2.5/W,1-(GEOMETRY["image_top_mm"]+full_h/2)/H,"Proteomics",fontsize=8,
             ha="center",va="center",rotation=90)
    stem="FigF_Fib_K14_Xgt1_210mm_single_row"+ ("_with_legend" if standalone else "_for_assembly")
    for col,title in enumerate(TITLES):
        x=GEOMETRY["left_mm"]+col*(group_w+GEOMETRY["column_gap_mm"])
        fig.text((x+group_w/2)/W,1-GEOMETRY["title_y_mm"]/H,title,ha="center",va="center",fontsize=8.5,
                 fontweight="bold" if col<2 else "normal",fontstyle="normal")
        ax=fig.add_axes([x/W,1-(GEOMETRY["image_top_mm"]+full_h)/H,full_w/W,full_h/H])
        zx=x+full_w+GEOMETRY["full_zoom_gap_mm"]
        ztop=GEOMETRY["image_top_mm"]+(full_h-zoom_w)/2
        az=fig.add_axes([zx/W,1-(ztop+zoom_w)/H,zoom_w/W,zoom_w/H])
        kind="distribution" if col==0 else ("coexpression" if col==1 else "protein")
        vals=None if col==0 else features[title]
        for current,zoom,lim,px,pw,suffix in [(ax,False,limits,x,full_w,"full"),(az,True,roi,zx,zoom_w,"zoom")]:
            pid=f"{stem}/{title}/{suffix}"
            dd=density.get((title,suffix))
            scatter_panel(st,current,he,xy,target,vals,lim,roi,zoom,kind,dd,pid,factor)
            add_scale(st,fig,current,px,pw,GEOMETRY["scale_y_mm"],H,lim,zoom,pid)
    if standalone:
        draw_legend(fig,H)
        fig.text(.5,1-41./H,f"Fib_K14: {target.sum():,} / {len(target):,} spots ({100*target.mean():.2f}%); all four proteins: X > 1.",
                 ha="center",va="center",fontsize=8)
        fig.text(.5,1-45.5/H,"Nominal scale: 100 µm/grid. Density contours: 50% / 85% of peak, not spot fractions.",
                 ha="center",va="center",fontsize=8)
    fig.canvas.draw()
    text_boxes=[(t.get_text(),t.get_window_extent(fig.canvas.get_renderer())) for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    overlaps=[]
    for i,(ta,ba) in enumerate(text_boxes):
        for tb,bb in text_boxes[i+1:]:
            intersection=matplotlib.transforms.Bbox.intersection(ba,bb)
            if intersection is not None and intersection.width>.5 and intersection.height>.5:
                overlaps.append([ta,tb])
    assert not overlaps,overlaps
    st.export(fig,PLOTS/stem,12)
    return stem,dict(full_width_mm=full_w,full_height_mm=full_h,zoom_width_mm=zoom_w,
                     per_column_width_mm=group_w,reference_dot_dimension_factor=factor,
                     no_text_overlaps=True,all_column_zoom_bounds_identical=True)


def main():
    if (OUT/"STATUS.json").exists():
        raise FileExistsError("Completed output exists; use a new version")
    for p in (PLOTS,TABLES,AUDIT,ASSETS,CODE):p.mkdir(parents=True,exist_ok=True)
    st=load_module("locked_he_style",STYLE)
    ref=load_module("reference_spatial",REFERENCE)
    protected=[INPUT,HE,STYLE,REFERENCE,REF_PDF,PREV/"audit/scale_calibration_provenance.json",PREV/"audit/uniform_HE_display_settings.json"]
    before={str(p):st.sha(p) for p in protected}
    a=ad.read_h5ad(INPUT)
    assert a.shape==(2029,5846) and a.obs_names.is_unique and a.var_names.is_unique
    assert all(g in a.var_names for g in GENES)
    xy=np.asarray(a.obsm["spatial"],float)
    prior=np.asarray(a.obsm["spatial_before_xaxis_reflection"],float)
    assert np.array_equal(xy[:,0],prior[:,0]) and np.allclose(xy[:,1]+prior[:,1],38)
    assert len(np.unique(xy,axis=0))==a.n_obs
    sub=a[:,GENES]
    values=dense(sub.X).astype(float)
    raw=dense(sub.layers["raw_intensity"]).astype(float)
    assert np.isfinite(values).all() and values.min()>=0
    target=ref.gate(values,"ALL",THRESHOLD)
    assert np.array_equal(target,(values[:,0]>1)&(values[:,1]>1)&(values[:,2]>1)&(values[:,3]>1))
    assert target.sum()==659
    assert np.all(target<=np.all(values>.5,axis=1))
    quantiles=np.quantile(values[target],.995,axis=0)
    scaled=np.clip(values/quantiles,0,1)
    coexpression=np.prod(scaled,axis=1)**.25
    features={g:scaled[:,i] for i,g in enumerate(GENES)}
    features["Co-expression"]=coexpression
    limits=ref.full_limits(xy)
    roi=select_roi(xy,target,st.inside)
    selected=st.inside(xy,roi["bounds"])
    densities,density_audit={},[]
    for feature_name in ["Co-expression",*DISPLAY_GENES]:
        for suffix,bounds in [("full",limits),("zoom",roi["bounds"])]:
            dd,info=ref.density(xy[target],features[feature_name][target],bounds,{"kde_max_points":0},1729)
            densities[(feature_name,suffix)]=dd
            density_audit.append(dict(feature=feature_name,region=suffix,**info))
    all_spots=a.obs.copy()
    all_spots.insert(0,"spot_id",a.obs_names.astype(str))
    all_spots["x_display"]=xy[:,0];all_spots["y_display"]=xy[:,1]
    for i,g in enumerate(GENES):
        all_spots[g+"_X"]=values[:,i]
        all_spots[g+"_input_intensity"]=raw[:,i]
        all_spots[g+"_X_gt1"]=values[:,i]>1
        all_spots[g+"_display_scaled"]=scaled[:,i]
    all_spots["Fib_K14"]=target
    all_spots["Fib_K14_definition"]="KRT14_X>1 AND PDGFRA_X>1 AND KRT5_X>1 AND TACSTD2_X>1"
    all_spots["coexpression_display_score"]=coexpression
    all_spots["in_shared_zoom_ROI"]=selected
    all_spots.to_csv(TABLES/"all_2029_spots_Fib_K14_Xgt1.tsv",sep="\t",index=False,float_format="%.10g")
    all_spots.loc[target].to_csv(TABLES/"Fib_K14_659_spots_Xgt1_coordinates.tsv",sep="\t",index=False,float_format="%.10g")
    all_spots.loc[np.all(values>.5,axis=1)&~target].to_csv(TABLES/"excluded_when_raising_0p5_to_1.tsv",sep="\t",index=False,float_format="%.10g")
    pd.DataFrame([dict(threshold=t,comparison="strictly_greater_than",expression_source="X_log1p_total_intensity_normalized",
                       n_Fib_K14=int(ref.gate(values,"ALL",t).sum()),n_all=len(target),percent=100*ref.gate(values,"ALL",t).mean())
                  for t in [0,.5,1,2,3,4,5,6,7,8]]).to_csv(TABLES/"threshold_sensitivity_counts.tsv",sep="\t",index=False)
    pd.DataFrame([dict(marker=g,source="X",gate_threshold=1,n_X_gt1=int((values[:,i]>1).sum()),
                       display_q995_among_659_Fib_K14=quantiles[i],input_accession=str(sub.var.iloc[i]["Protein accession"]))
                  for i,g in enumerate(GENES)]).to_csv(TABLES/"four_protein_gate_and_display_scales.tsv",sep="\t",index=False)
    sub.var.to_csv(TABLES/"four_protein_identification_metadata.tsv",sep="\t")
    pd.DataFrame([{k:v for k,v in roi.items() if k!="bounds"}|dict(zip(["x_min","x_max","y_min","y_max"],roi["bounds"]))]).to_csv(TABLES/"shared_zoom_ROI.tsv",sep="\t",index=False)
    for col in ["leiden_unsupervised","sub_labels","bulk_labels"]:
        if col in all_spots:
            table=all_spots.groupby(col,observed=True).agg(n_spots=("Fib_K14","size"),n_Fib_K14=("Fib_K14","sum")).reset_index()
            table.to_csv(TABLES/f"Fib_K14_by_{col}.tsv",sep="\t",index=False)
    shutil.copy2(HE,ASSETS/HE.name)
    he=np.asarray(Image.open(ASSETS/HE.name).convert("RGB"))
    for item in ("scale_calibration_provenance.json","uniform_HE_display_settings.json"):
        shutil.copy2(PREV/"audit"/item,AUDIT/item)
    plt.rcParams.update({"font.family":"Arial","font.size":8,"axes.titlesize":9,"axes.labelsize":8,
                         "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"figure.titlesize":10,
                         "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none","figure.facecolor":"white"})
    plots=[]
    for standalone in (False,True):
        stem,geometry=render(st,he,xy,target,features,densities,limits,roi["bounds"],standalone)
        plots.append(dict(stem=stem,**geometry))
    bars=[{k:v for k,v in r.items() if k not in ("ax","line")} for r in st.BAR_RECORDS]
    assert len(bars)==24 and all(np.isclose(r["drawn_mm"],r["expected_drawn_mm"],rtol=1e-9) for r in bars)
    pd.DataFrame(bars).to_csv(AUDIT/"scale_bar_geometry.tsv",sep="\t",index=False)
    after={str(p):st.sha(p) for p in protected}
    assert before==after
    st.dump(AUDIT/"source_sha256.json",before)
    st.dump(AUDIT/"density_audit.json",density_audit)
    st.dump(AUDIT/"panel_QA.json",PANEL_QA)
    st.dump(AUDIT/"figure_QA.json",dict(figures=st.FIGURE_RECORDS,geometry=plots,
                original_sources_unchanged=before==after,Arial_8to10pt=True,protein_names_upright=True,
                same_HE_asset_all_panels=True,same_HE_mapping_all_panels=True,source_orientation_unchanged=True,
                scale_bars_validated=len(bars),scale_is_nominal=True,gate_independently_recomputed=True,
                strict_threshold=1,expression_source="X",n_Fib_K14=659,n_all=2029))
    with fitz.open(REF_PDF) as doc:
        ref_size=[float(doc[0].rect.width*25.4/72),float(doc[0].rect.height*25.4/72)]
    config=dict(input=str(INPUT),sample_scope="one company SPX AI-reconstructed spatial-proteomics section; 2029 spatial grid spots",
                marker_rule=all_spots["Fib_K14_definition"].iloc[0],threshold=THRESHOLD,
                gate_expression_source="X: total-intensity normalized and log1p transformed; not 0-1 display scores or raw intensity",
                user_confirmation="User selected X > 0.5, then explicitly raised threshold to 1 on 2026-10-15.",
                proteins=GENES,display_order=TITLES,n_target=int(target.sum()),percent_target=100*target.mean(),
                display_normalization="Divide each protein X by q99.5 among current target spots, clip to [0,1]; display only, not gating",
                coexpression="Fourth root of product of four scaled abundances; only selected target spots plotted in coexpression column",
                individual_featureplots="All 2029 spatial locations plotted",KDE="Expression-weighted KDE of current target spots; Scott bandwidth/1.5; contours 50/85 percent of peak",
                zoom=roi,full_limits=limits,geometry=GEOMETRY,reference_pdf=str(REF_PDF),reference_pdf_mm=ref_size,
                preserved_parameters=["width210mm","Arial","8to10pt","protein_labels_upright","600dpi","viridis","same_HE_asset","epidermis_up","nominal100um_per_grid"],
                layout_adaptation="Original multi-row canvas reduced to one 30mm-high row; standalone legend copy49mm high. Six columns include TACSTD2. No repeated F label.",
                source_preservation="No original source or previous output modified")
    st.dump(OUT/"run_config.json",config)
    shutil.copy2(Path(__file__),CODE/Path(__file__).name)
    shutil.copy2(STYLE,CODE/"locked_he_style.py")
    shutil.copy2(REFERENCE,CODE/"reference_spatial_batch.py")
    (OUT/"README_结果说明.md").write_text(
        "# Figure F 单行空间蛋白组补图\n\n"
        "已按最新要求将 Fib_K14 定义为 KRT14 >1 且 PDGFRA >1 且 KRT5 >1 且 TACSTD2 >1。"
        "阈值作用于用户确认的 X：逐网格总强度标准化后 log1p 矩阵，不是原始质谱强度，也不是 0–1 绘图色阶。"
        f"结果 659/2,029 个 spot（{100*target.mean():.2f}%）；>0.5 时 660 个，只有 1 个因此排除。\n\n"
        "拼版主文件：plots/FigF_Fib_K14_Xgt1_210mm_single_row_for_assembly.pdf（210 × 30 mm）；"
        "附带完整图例版本为 *_with_legend.pdf（210 × 49 mm）。"
        "两版每列图像及文字物理尺寸一致，仅下方图例/说明空间不同；均导出 SVG、600 dpi PNG 与 LZW TIFF。"
        "PDF/SVG 保留矢量点和 Arial 文字。拼接时按 210 mm 原宽放置，避免缩小后低于 8 pt。\n\n"
        "列顺序：Distribution、Co-expression、PDGFRA、KRT14、KRT5、TACSTD2，每列全景+相同 ROI 放大。"
        "无重复 F 字母或新增 Mode 编号。四个蛋白名称正体，分类标题粗体。"
        "单行版 8–8.5 pt；上一版要求允许范围 8–10 pt。统一复用上一版 HE RGB 素材、映射、透明度与表皮向上方向，未再次翻转。\n\n"
        "Distribution 红点为本规则 Fib_K14，灰点为其余 spot。Co-expression 只展示 659 个规则阳性 spot 的四蛋白几何均值。"
        "各单蛋白列展示所有 2,029 个位置的相对丰度；每蛋白按阳性群体 X 的 q99.5 缩放至 0–1，仅用于显色，不改变 >1 筛选。"
        "KDE 仅在规则阳性群体计算，轮廓为各图密度峰值的 50%/85%，不是覆盖 spot 的比例。"
        "ROI 为 14×14 网格窗口，按目标 spot 数最多、纯度、距离及坐标顺序选择；全部列共用。\n\n"
        "比例尺：全景 1000 µm、放大 200 µm，名义 100 µm/grid，继承公司报告 6.5×3.9 mm 与 65×39 网格推算，未经独立物理校准。"
        "为了保留 8 pt 字号而不遮住微小放大图，标尺及标签置于各图下方。\n\n"
        "解释边界：这是单张切片的公司 AI 重建空间蛋白网格，规则命名 Fib_K14 仅描述同网格四蛋白共丰度，"
        "不能直接证明同一个细胞共表达或独立确认成纤维细胞谱系/EMT。"
        "原始鉴定表中 PDGFRA 为 1 条 unique peptide、覆盖率 1.01%，应结合独立验证；详见四蛋白鉴定元数据。"
        "本图未将旧 KLF4/TP63 规则合并进来，也未覆盖旧分析文件。\n\n"
        "tables/ 提供 659 个位置、全部网格、阈值敏感性数量、放大 ROI 与显示缩放因子；audit/ 提供 SHA256、"
        "标尺几何、字体与实际 PDF 输出 QA。\n",encoding="utf-8")
    files=[dict(path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=st.sha(p)) for p in OUT.rglob("*") if p.is_file()]
    assert all(f["bytes"]>0 for f in files)
    st.dump(OUT/"STATUS.json",dict(status="COMPLETE",n_Fib_K14=int(target.sum()),n_all=2029,threshold=1,
                                 gate_scale="X_log1p",QA="PASS",files=files))
    print(json.dumps(dict(output=str(OUT),n_target=int(target.sum()),percent=100*target.mean(),roi=roi,figures=st.FIGURE_RECORDS),ensure_ascii=False))


if __name__=="__main__":main()
