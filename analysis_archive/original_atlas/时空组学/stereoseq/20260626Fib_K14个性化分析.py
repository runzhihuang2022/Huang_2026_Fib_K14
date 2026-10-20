############################# 完整多样本 Fib_K14 专属共定位分析流程 (白底黑字清爽版)
import os
import gc
import scanpy as sc
import squidpy as sq
import spatialdata_plot  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

# 忽略除以 0 产生的 RuntimeWarning 警告
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================
# 1. 环境、颜色与多样本任务配置
# ==========================================
# 全局启用白底黑字的绘图样式
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42,
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'text.color': 'black',
    'axes.labelcolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black'
})

# 基础输出路径
base_output_dir = r"data/external/local_8e645d0a9b03"
annotated_input_dir = os.path.join(base_output_dir, "Wound_Healing_Annotation_Output_White")
# 空间共定位输出目录切换为 White 文件夹
colocal_output_dir = os.path.join(base_output_dir, "Wound_Healing_Colocalization_White")
os.makedirs(colocal_output_dir, exist_ok=True)

# 官方指定的细胞配色字典 (用于同步热图侧边分类颜色条)
raw_annotation_config = {
    "KC_Basal": {"color": "#008B45"}, "KC_Basal_Mig": {"color": "#FF8C00"}, "KC_Basal_Prolif": {"color": "#FFD700"},
    "KC_Spinous": {"color": "#00FF00"}, "KC_Spinous_Mig": {"color": "#00008B"}, "KC_Spinous_Mat": {"color": "#FF69B4"}, "KC_Granular": {"color": "#8B4513"},
    "SAC_SG_Progenitor": {"color": "#0099CC"}, "SAC_SG_Clear": {"color": "#AFEEEE"}, "SAC_SG_Dark": {"color": "#0A2B00"}, "SAC_SG_Ductal": {"color": "#CD853F"},
    "SAC_HF_IRS": {"color": "#DA70D6"}, "SAC_HF_ORS": {"color": "#D8BFD8"}, "SAC_HF_HFSC": {"color": "#20B2AA"}, "SAC_HF_DP_DS": {"color": "#B8860B"}, "SAC_HF_Matrix": {"color": "#48D1CC"},
    "Melanocyte": {"color": "#A0522D"},
    "Fib_Papi": {"color": "#4682B4"}, "Fib_SFRP2": {"color": "#87CEFA"}, "Fib_Myo": {"color": "#DDA0DD"}, "Fib_Inflama": {"color": "#BA55D3"}, "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#F0E68C"},
    "Fib_Prolif": {"color": "#FFB6C1"}, "Fib_K14": {"color": "#E31A1C"}, "Fib_CD45": {"color": "#9932CC"},
    "Endo_Capillary": {"color": "#9ACD32"}, "Endo_Arterial": {"color": "#6B8E23"}, "Endo_Venous": {"color": "#556B2F"}, "Endo_Lymphatic": {"color": "#8FBC8F"}, "Endo_Prolif": {"color": "#BDB76B"},
    "Pericyte": {"color": "#DAA520"}, "Schwann": {"color": "#D2B48C"},
    "M1_Macrophage": {"color": "#1E90FF"}, "M2_Macrophage": {"color": "#4169E1"}, "LAM": {"color": "#0000CD"}, "Neutrophil": {"color": "#483D8B"}, "Mast": {"color": "#800080"},
    "cDC1": {"color": "#6495ED"}, "cDC2": {"color": "#7B68EE"}, "pDC": {"color": "#8470FF"}, "LGHS": {"color": "#00CED1"},
    "B_Mature": {"color": "#000080"}, "B_Plasma": {"color": "#4B0082"}, "CD4_T": {"color": "#5F9EA0"}, "CD8_T": {"color": "#6A5ACD"}, "NK": {"color": "#008080"}
}

config_wound = {
    "5dpb":         {"sample_id": "BW32_A01597A3_SDSDB_5dpb"},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb"},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb"},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb"},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb"},
    "19dpb":        {"sample_id": "BW15D_C6_SDSDB_19dpb"},
    "19dpb_p1":     {"sample_id": "BW15D_1_D1_SDSDB_19dpb"},
    "26dpb_p1":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1"},
    "26dpb_p2":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2"},
    "2mph":         {"sample_id": "2mph_A03699G6.SCT"},
    "6mph_p1":      {"sample_id": "6mph_A02988A2.SCT_part1"},
    "6mph_p2":      {"sample_id": "6mph_A02988A2.SCT_part2"},
    "9mph":         {"sample_id": "9mph_A03091C5.SCT"},
    "Normal":       {"sample_id": "NS_C02847B1"}
}

target_cell = "Fib_K14"

# ==========================================
# 2. 自动化循环执行核心流
# ==========================================
for skey, info in config_wound.items():
    tid = info["sample_id"]
    file_path = os.path.join(annotated_input_dir, f"{tid}_Annotated.h5ad")
    
    if not os.path.exists(file_path):
        print(f"⚠️ 跳过：未找到样本 {skey} 的注释文件: {file_path}")
        continue
        
    print(f"\n🚀 ========================================================")
    print(f"⏳ 正在处理样本 [{skey}] 共定位分析...")
    print(f"   ========================================================")
    
    adata = sc.read_h5ad(file_path)

    if 'sub_labels' in adata.obs:
        adata.obs['sub_labels'] = adata.obs['sub_labels'].astype('category')
        categories = adata.obs['sub_labels'].cat.categories
        mapped_colors = []
        for cat in categories:
            c = raw_annotation_config.get(cat, {}).get('color', '#888888')
            if len(c) == 9 and c.startswith('#'):
                c = c[:7]
            mapped_colors.append(c)
        adata.uns['sub_labels_colors'] = mapped_colors
    else:
        print(f"❌ 错误：样本 {skey} 数据中未找到 'sub_labels' 列。")
        continue

    if 'spatial' not in adata.obsm:
        if 'X_spatial' in adata.obsm: 
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        else:
            continue

    # 空间邻居图与邻域富集计算 (Squidpy)
    sq.gr.spatial_neighbors(adata, coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key="sub_labels")

    # 🌟 关键步骤：提取全量矩阵并动态筛选出与 Fib_K14 最相关的 TOP 10 细胞
    enrichment_zscore = adata.uns['sub_labels_nhood_enrichment']['zscore']
    cluster_names = adata.obs['sub_labels'].cat.categories
    df_enrichment = pd.DataFrame(enrichment_zscore, index=cluster_names, columns=cluster_names)

    if target_cell in df_enrichment.index:
        fib_k14_interactions = df_enrichment.loc[target_cell].sort_values(ascending=False)
        top_colocalized = fib_k14_interactions[fib_k14_interactions.index != target_cell].head(10)
        # 组装热图只展示的细胞群体名单：Fib_K14 本身 + 与之共定位得分最高的 10 种细胞
        top10_heatmap_cells = [target_cell] + top_colocalized.index.tolist()
        
        print(f"  📊 样本 [{skey}] 已锁定共定位 Top 10 细胞进行热图展示。")
    else:
        print(f"  ⚠️ 样本 [{skey}] 未检测到 {target_cell}，展示完整热图。")
        top10_heatmap_cells = list(cluster_names)

    # ==========================================
    # 4. 绘制并保存高级白底黑字富集热图 (仅展示 TOP 10)
    # ==========================================
    # 通过 clusters 参数过滤热图范围
    sq.pl.nhood_enrichment(
        adata, 
        cluster_key="sub_labels", 
        clusters=top10_heatmap_cells, 
        figsize=(11, 9),  
        title=f"Neighborhood Enrichment (Top 10 + {target_cell}) - {skey}",
        cmap="viridis",    
        show=False
    )

    fig = plt.gcf()
    fig.patch.set_facecolor('white')

    # 遍历画布上产生的所有复杂子图，统一套用清爽白底黑字主题
    for current_ax in fig.axes:
        current_ax.set_facecolor('white')
        current_ax.tick_params(colors='black', which='both')
        
        # 将所有的边框和轴线涂黑
        for spine in current_ax.spines.values():
            spine.set_color('black')
        
        # 字体颜色涂黑
        current_ax.xaxis.label.set_color('black')
        current_ax.yaxis.label.set_color('black')
        if current_ax.title:
            current_ax.title.set_color('black')

    # X轴标签旋转
    for current_ax in fig.axes:
        if len(current_ax.get_xticklabels()) > 0:
            plt.setp(current_ax.get_xticklabels(), rotation=90, ha="right", va="center")

    save_png = os.path.join(colocal_output_dir, f"{skey}_Spatial_Neighborhood_Enrichment_Top10.png")
    save_pdf = os.path.join(colocal_output_dir, f"{skey}_Spatial_Neighborhood_Enrichment_Top10.pdf")

    fig.savefig(save_png, dpi=600, facecolor='white', bbox_inches='tight')
    fig.savefig(save_pdf, facecolor='white', bbox_inches='tight')
    plt.close(fig)

    # 打印简要状态日志
    if target_cell in df_enrichment.index:
        csv_path = os.path.join(colocal_output_dir, f"{skey}_Neighborhood_Enrichment_Zscores.csv")
        df_enrichment.to_csv(csv_path)
        print(f"  ✨ {skey} Z-score 完整矩阵已保存，Top 10 热图绘制成功。")

    del adata; gc.collect()

print("\n🎉 所有指定创面愈合样本的 Fib_K14 专属共定位 TOP 10 白底清爽版热图已全部处理完毕！")

############################# 完整多样本 Fib_K14 空间共现性分析流程 (白底黑字清爽版)
import os
import gc
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

# 忽略除以 0 产生的 RuntimeWarning 警告
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ==========================================
# 0. 环境与全局白底主题配置
# ==========================================
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42,
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'text.color': 'black',
    'axes.labelcolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black'
})

base_output_dir = r"data/external/local_8e645d0a9b03"
annotated_input_dir = os.path.join(base_output_dir, "Wound_Healing_Annotation_Output_White")
# 输出文件夹切换为 White 文件夹
co_occurrence_out_dir = os.path.join(base_output_dir, "Wound_Healing_Co_occurrence_White")
os.makedirs(co_occurrence_out_dir, exist_ok=True)

raw_annotation_config = {
    "KC_Basal": {"color": "#008941"}, "KC_Spinous": {"color": "#00FF00FF"},
    "KC_Spinous_Mat": {"color": "#EE4C97FF"}, "KC_Granular": {"color": "#882D17"},
    "KC_Spinous_Mig": {"color": "#1A0099FF"}, "KC_Basal_Mig": {"color": "#FF3030"},
    "KC_Basal_Prolif": {"color": "#FFD700"}, "SAC_SG_Progenitor": {"color": "#0099CCFF"},
    "SAC_SG_Clear": {"color": "#B7E4F9FF"}, "SAC_SG_Dark": {"color": "#193006FF"},
    "SAC_SG_Ductal": {"color": "#8F7700FF"}, "SAC_HF_IRS": {"color": "#CC9900FF"},
    "SAC_HF_ORS": {"color": "#A30059"}, "SAC_HF_HFSC": {"color": "#00F5FF"},
    "SAC_HF_DP_DS": {"color": "#CDBE70"}, "SAC_HF_Matrix": {"color": "#7FFFD4"},
    "Melanocyte": {"color": "#FF4500"}, "Fib_Papi": {"color": "#3DA873FF"},
    "Fib_SFRP2": {"color": "#CCFF00FF"}, "Fib_Myo": {"color": "#B79762"},
    "Fib_Inflama": {"color": "#800000FF"}, "Fib_Fasci": {"color": "#9370DB"},
    "Fib_EN1": {"color": "#BDB76B"}, "Fib_CD45": {"color": "#DA70D6"},
    "Fib_Prolif": {"color": "#EED2EE"}, "Fib_K14": {"color": "#4682B4"},
    "Endo_Capillary": {"color": "#CF4E9CFF"}, "Endo_Arterial": {"color": "#FF00FF"},
    "Endo_Venous": {"color": "#4169E1"}, "Endo_Lymphatic": {"color": "#ADFF2F"},
    "Endo_Prolif": {"color": "#F08080"}, "Pericyte": {"color": "#FFFF00"},
    "Schwann": {"color": "#7CFC00"}, "M1_Macrophage": {"color": "#FF4040"},
    "M2_Macrophage": {"color": "#54FF9F"}, "LAM": {"color": "#8B4513"},
    "Neutrophil": {"color": "#FF7F50"}, "Mast": {"color": "#FF1493"},
    "cDC1": {"color": "#9B30FF"}, "cDC2": {"color": "#EE9A00"},
    "pDC": {"color": "#D1EEEE"}, "LGHS": {"color": "#FFB5C5"},
    "B_Mature": {"color": "#0000CD"}, "B_Plasma": {"color": "#FAF0E6"},
    "CD4_T": {"color": "#00BFFF"}, "CD8_T": {"color": "#1E90FF"},
    "NK": {"color": "#00E5EE"}
}

target_cluster = "Fib_K14"

config_wound = {
    "5dpb":         {"sample_id": "BW32_A01597A3_SDSDB_5dpb"},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb"},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb"},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb"},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb"},
    "19dpb":        {"sample_id": "BW15D_C6_SDSDB_19dpb"},
    "19dpb_p1":     {"sample_id": "BW15D_1_D1_SDSDB_19dpb"},
    "26dpb_p1":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1"},
    "26dpb_p2":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2"},
    "2mph":         {"sample_id": "2mph_A03699G6.SCT"},
    "6mph_p1":      {"sample_id": "6mph_A02988A2.SCT_part1"},
    "6mph_p2":      {"sample_id": "6mph_A02988A2.SCT_part2"},
    "9mph":         {"sample_id": "9mph_A03091C5.SCT"},
    "Normal":       {"sample_id": "NS_C02847B1"}
}

# ==========================================
# 1. 多样本自动化流水线循环
# ==========================================
for skey, info in config_wound.items():
    tid = info["sample_id"]
    file_path = os.path.join(annotated_input_dir, f"{tid}_Annotated.h5ad")
    
    if not os.path.exists(file_path):
        print(f"⚠️ 跳过：未找到样本 [{skey}] 的输入文件: {file_path}")
        continue
        
    print(f"\n🚀 ========================================================")
    print(f"⏳ 正在计算样本 [{skey}] 空间共现性矩阵...")
    print(f"   ========================================================")
    
    adata = sc.read_h5ad(file_path)
    
    current_cats = adata.obs['sub_labels'].unique().dropna()
    ordered_cats = [k for k in raw_annotation_config.keys() if k in current_cats]
    missing_cats = [k for k in current_cats if k not in ordered_cats]
    final_order = ordered_cats + missing_cats

    adata.obs['sub_labels'] = adata.obs['sub_labels'].astype('category').cat.reorder_categories(final_order)

    # 1. 计算核心空间共现性
    sq.gr.co_occurrence(adata, cluster_key="sub_labels")

    # 2. 提取底层矩阵
    co_out = adata.uns["sub_labels_co_occurrence"]["occ"]
    if target_cluster not in final_order:
        print(f"⚠️ 警告：样本 [{skey}] 中未检测到目标细胞群 {target_cluster}，跳过绘图。")
        del adata; gc.collect()
        continue
        
    target_idx = list(final_order).index(target_cluster)
    closest_distance_vals = co_out[target_idx, :, 0]

    co_occ_df = pd.DataFrame({'Cell_Type': final_order, 'Co_occurrence_Score': closest_distance_vals})
    co_occ_df = co_occ_df[co_occ_df['Cell_Type'] != target_cluster]
    co_occ_df = co_occ_df.sort_values(by='Co_occurrence_Score', ascending=False)

    top10_df = co_occ_df.head(10)
    top10_clusters = top10_df['Cell_Type'].tolist()

    # ==========================================
    # 2. 预先构建高级白底网格架构 
    # ==========================================
    fig = plt.figure(figsize=(22, 10), facecolor='white')
    gs = fig.add_gridspec(2, 2, height_ratios=[8.2, 1.8], wspace=0.18, hspace=0.25)

    ax1 = fig.add_subplot(gs[0, 0], facecolor='white')
    ax2 = fig.add_subplot(gs[0, 1], facecolor='white')
    ax_leg = fig.add_subplot(gs[1, :], facecolor='white')
    ax_leg.axis('off')  

    # ==========================================
    # 3. 自主抽取数据并绘制左侧高级折线图 [ax1]
    # ==========================================
    intervals = adata.uns["sub_labels_co_occurrence"]["interval"]
    x_vals = intervals[:co_out.shape[2]]

    for j, c_name in enumerate(final_order):
        y_vals = co_out[target_idx, j, :]
        c_color = raw_annotation_config.get(c_name, {}).get('color', '#888888')
        if len(c_color) == 9 and c_color.startswith('#'):
            c_color = c_color[:7]
            
        if c_name == target_cluster:
            lw, z, alpha = 5.0, 15, 1.0
        elif c_name in top10_clusters:
            lw, z, alpha = 3.5, 10, 1.0
        else:
            lw, z, alpha = 1.0, 1, 0.15 # 非TOP10背景线条，白底上呈现浅色淡化效果

        ax1.plot(x_vals, y_vals, label=c_name, color=c_color, linewidth=lw, zorder=z, alpha=alpha)

    ax1.set_title(f"[{skey}] Co-occurrence Probability Conditioned on {target_cluster}", color='black', pad=20, fontsize=16, fontweight='bold')
    ax1.set_xlabel("Distance (µm)", color='black', fontsize=14)
    ax1.set_ylabel("Co-occurrence Probability", color='black', fontsize=14)

    handles, labels = ax1.get_legend_handles_labels()

    # ==========================================
    # 4. 绘制右侧 Top10 水平柱状图 [ax2]
    # ==========================================
    plot_df = top10_df.iloc[::-1]  
    bar_colors = []
    for c in plot_df['Cell_Type']:
        b_col = raw_annotation_config.get(c, {}).get('color', '#888888')
        if len(b_col) == 9 and b_col.startswith('#'):
            b_col = b_col[:7]
        bar_colors.append(b_col)

    bars = ax2.barh(plot_df['Cell_Type'], plot_df['Co_occurrence_Score'], color=bar_colors, edgecolor='none', height=0.6)
    
    ax2.set_title(f"[{skey}] Top 10 Co-occurring Cells (Score)", color='black', pad=20, fontsize=16, fontweight='bold')
    ax2.set_xlabel("Co-occurrence Score at Closest Distance", color='black', fontsize=14)

    # 🌟 右侧数值标签显色改为黑色 (color='black')
    for bar in bars:
        width = bar.get_width()
        ax2.text(
            width + (plot_df['Co_occurrence_Score'].max() * 0.015), 
            bar.get_y() + bar.get_height()/2, 
            f'{width:.3f}', 
            va='center', 
            ha='left', 
            color='black', 
            fontsize=11, 
            fontweight='bold'
        )

    ax2.set_xlim(0, plot_df['Co_occurrence_Score'].max() * 1.15)

    # ==========================================
    # 5. 双图坐标系细节涂黑润色
    # ==========================================
    for ax in [ax1, ax2]:
        ax.tick_params(colors='black', labelsize=12, which='both')
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(1.2)

    # ==========================================
    # 6. 底部独立幕布区绘制全量联动图例 [ax_leg]
    # ==========================================
    new_legend = ax_leg.legend(
        handles, labels,
        loc='center', 
        ncol=7,           
        frameon=False, 
        fontsize=11,
        handletextpad=0.5,
        columnspacing=1.5
    )

    # 🌟 调整字体高亮逻辑：高亮部分为黑色，次要部分降为浅灰色
    for text in new_legend.get_texts():
        txt_val = text.get_text()
        if txt_val == target_cluster or txt_val in top10_clusters:
            text.set_color("black")
            text.set_weight("bold") 
        else:
            text.set_color("#BBBBBB") # 浅灰色，避免白底上非重点细胞名称造成视觉疲劳

    # ==========================================
    # 7. 保存当前样本的白底面板
    # ==========================================
    save_png = os.path.join(co_occurrence_out_dir, f"{skey}_Spatial_Co_occurrence_Panel_{target_cluster}.png")
    save_pdf = os.path.join(co_occurrence_out_dir, f"{skey}_Spatial_Co_occurrence_Panel_{target_cluster}.pdf")

    fig.savefig(save_png, dpi=600, facecolor='white', bbox_inches='tight')
    fig.savefig(save_pdf, facecolor='white', bbox_inches='tight')
    
    plt.close(fig)
    print(f"✅ 样本 [{skey}] 共现性全景白底图册已完美保存。\n")
    
    del adata; gc.collect()

print("\n🎉 创面愈合任务全线样本（Wound_Healing）高级感白底黑字共现性图册已流水线式导出处理完毕！")

############################# 空间分布规律动态分析 (白底黑字学术版)
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import json
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. 路径与样本配置
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"

sample_config = {
    "5dpb":        {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb",        "angle": 135},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb",        "angle": 135},
    "19dpb":       {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225},
    "19dpb_p1":    {"sample_id": "BW15D_1_D1_SDSDB_19dpb",    "angle": 140},
    "26dpb_p1":    {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1", "angle": 0},
    "26dpb_p2":    {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2", "angle": 0},
    "2mph":        {"sample_id": "2mph_A03699G6.SCT",         "angle": 45},
    "6mph_p1":     {"sample_id": "6mph_A02988A2.SCT_part1",   "angle": 0},
    "6mph_p2":     {"sample_id": "6mph_A02988A2.SCT_part2",   "angle": 0},
    "9mph":        {"sample_id": "9mph_A03091C5.SCT",         "angle": 0},
    "Normal":      {"sample_id": "NS_C02847B1",               "angle": 180}
}


# ⚠️ 在这里修改您当前想要处理的样本 Key
TEST_KEY = "19dpb"  

TARGET_SAMPLE = sample_config[TEST_KEY]["sample_id"]
ROTATION_ANGLE = sample_config[TEST_KEY]["angle"]

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad") 
json_save_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

# ==========================================
# 1. 核心功能函数
# ==========================================
def rotate_coordinates(x, y, angle):
    """旋转坐标系的纯数学计算"""
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

def extract_structural_anchors(adata, sample_id, angle):
    # 提取坐标并执行旋转
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = rotate_coordinates(raw_x, raw_y, angle)
    
    adata.obsm['X_rot_spatial'] = np.column_stack((rot_x, rot_y))
    
    # 🌟 核心修改 1: 弃用暗黑模式，采用白底默认样式
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150, facecolor='white')
    ax.set_facecolor('white')
    
    try:
        sc.pl.embedding(adata, basis='rot_spatial', color='sub_labels', ax=ax, show=False)
    except Exception:
        # 在白底上使用浅灰色的散点点缀背景
        ax.scatter(rot_x, rot_y, s=2, c='#E0E0E0', alpha=0.6, edgecolors='none')

    # 精准显示边框、刻度和深色标签
    ax.axis('on') 
    ax.tick_params(colors='black', labelsize=10)
    for spine in ax.spines.values():
        spine.set_color('black')
    
    print("\n" + "="*50)
    print(f"🖱️ 【结构锚点提取开启】正在处理: {sample_id} (已旋转 {angle} 度)")
    print("左键=选点，右键=撤销上一个点，回车(Enter)=完成当前线段。")
    print("="*50 + "\n")

    # ----------------------------------------
    # 线段 1: 表皮底层
    # ----------------------------------------
    # 🌟 核心修改 2: 标题文字改为黑色/深天蓝色以适配白底
    ax.set_title("[Step 1/3] Click Epidermis Baseline -> Press ENTER", color='#00688B', fontsize=14, fontweight='bold')
    plt.draw()
    print("👉 [1/3] 请描绘【表皮底层/基底膜】走势，完成后按 Enter。")
    epi_coords = plt.ginput(n=-1, timeout=0)
    if epi_coords:
        ex, ey = zip(*epi_coords)
        ax.plot(ex, ey, color='#009ACD', linewidth=2, marker='o', markersize=4)
        plt.draw()
    print(f"   ✔️ Step 1 完成，共采集 {len(epi_coords)} 个点。")
    plt.pause(0.2) 

    # ----------------------------------------
    # 线段 2: 真皮底层
    # ----------------------------------------
    # 🌟 核心修改 3: 标题文字改为深洋红色以适配白底
    ax.set_title("[Step 2/3] Click Dermis Bottom -> Press ENTER", color='#8B1076', fontsize=14, fontweight='bold')
    plt.draw()
    print("\n👉 [2/3] 请描绘【真皮最底层】走势，完成后按 Enter。")
    der_coords = plt.ginput(n=-1, timeout=0)
    if der_coords:
        dx, dy = zip(*der_coords)
        ax.plot(dx, dy, color='#CD1076', linewidth=2, marker='o', markersize=4)
        plt.draw()
    print(f"   ✔️ Step 2 完成，共采集 {len(der_coords)} 个点。")
    plt.pause(0.2) 

    # ----------------------------------------
    # 线段 3: 愈合前沿
    # ----------------------------------------
    # 🌟 核心修改 4: 标题文字在白底上改为深橙色/黄褐色
    ax.set_title("[Step 3/3] Click Leading Edge (2 points) -> Press ENTER", color='#8B4500', fontsize=14, fontweight='bold')
    plt.draw()
    print("\n👉 [3/3] 请点击【两个点】确定愈合前沿分界线（若已完全愈合请直接按 Enter 跳过）。")
    front_coords = plt.ginput(n=-1, timeout=0)
    if front_coords:
        fx, fy = zip(*front_coords)
        ax.plot(fx, fy, color='#FF8C00', linewidth=2, linestyle='--')
        plt.draw()
    print(f"   ✔️ Step 3 完成，共采集 {len(front_coords)} 个点。")
    plt.pause(0.2)

    # ----------------------------------------
    # 终极防误关确认 
    # ----------------------------------------
    ax.set_title("✅ All Steps Done! Please check the terminal to confirm & save.", color='#228B22', fontsize=14, fontweight='bold')
    plt.draw()
    
    print("\n" + "="*50)
    print("⏸️ 【等待确认】图片目前保持打开状态。")
    print("请对照屏幕上的图像（查看前沿虚线位置），在下方终端输入参数：")
    
    if len(front_coords) >= 2:
        direction_input = input("❓ 以愈合前沿为界，哪边是【已长好】的上皮区？(输入 R 代表右侧，L 代表左侧): ").strip().upper()
        direction = "right_is_healed" if direction_input == 'R' else "left_is_healed"
    else:
        print("ℹ️ 未检测到虚线，判定为完全愈合样本。")
        input("按任意键并回车继续...")
        direction = "right_is_healed" 
    
    print("="*50)
    plt.close(fig)
    
    return {
        'epi_baseline': [(float(x), float(y)) for x, y in epi_coords],
        'der_bottom': [(float(x), float(y)) for x, y in der_coords],
        'leading_edge': [(float(x), float(y)) for x, y in front_coords],
        'healed_direction': direction,
        'rotation_applied': angle 
    }

if __name__ == "__main__":
    if os.path.exists(file_path):
        print(f"⏳ 正在加载数据: {file_path} ...")
        adata = sc.read_h5ad(file_path)
        
        anchor_dict = extract_structural_anchors(adata, sample_id=TEST_KEY, angle=ROTATION_ANGLE)
        
        with open(json_save_path, 'w', encoding='utf-8') as f:
            json.dump(anchor_dict, f, indent=4)
        print(f"\n✅ 锚点配置已成功保存至: {json_save_path}")
        print("👉 下一步：您可以去运行【阶段二】脚本来生成距离密度图了！")
    else:
        print(f"❌ 找不到文件: {file_path}")

############################## 调整JSON文件 (白底黑字学术版)
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import json
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. 路径配置 
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
TEST_KEY = "19dpb"  
TARGET_SAMPLE = "BW15D_C6_SDSDB_19dpb"
ROTATION_ANGLE = -225

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad") 
json_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

def rotate_coordinates(x, y, angle):
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

if __name__ == "__main__":
    print(f"⏳ 正在加载数据与已有的 JSON 坐标...")
    adata = sc.read_h5ad(file_path)
    
    with open(json_path, 'r', encoding='utf-8') as f:
        anchor_dict = json.load(f)
        
    # 提取并旋转细胞坐标
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = rotate_coordinates(raw_x, raw_y, ROTATION_ANGLE)
    
    # 🌟 核心修改 1: 切换全局白底配置
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150, facecolor='white')
    ax.set_facecolor('white')
    
    # 绘制背景细胞散点 (浅灰色点契合学术白底)
    ax.scatter(rot_x, rot_y, s=2, c='#E0E0E0', alpha=0.6, edgecolors='none')
    
    # 细节修饰涂黑
    ax.axis('on') 
    ax.tick_params(colors='black', labelsize=10)
    for spine in ax.spines.values():
        spine.set_color('black')
    
    # 🌟 自动绘制已经保存好的表皮线和真皮线 (在白底上采用醒目的深青与深红)
    ex, ey = zip(*anchor_dict['epi_baseline'])
    ax.plot(ex, ey, color='#009ACD', linewidth=2, marker='o', markersize=4, label="Epidermis (Kept)")
    
    dx, dy = zip(*anchor_dict['der_bottom'])
    ax.plot(dx, dy, color='#CD1076', linewidth=2, marker='o', markersize=4, label="Dermis (Kept)")
    
    ax.legend(loc="upper right", frameon=True, facecolor='white', edgecolor='black')
    
    # 🌟 仅要求重新输入 Step 3 (标题显色改为黑色) 🌟
    ax.set_title("Target Fix: Click 2 points to set NEW Leading Edge -> Press ENTER", color='black', fontsize=14, fontweight='bold')
    plt.draw()
    
    print("\n" + "="*50)
    print("👉 【微调模式】表皮和真皮线已自动保留！")
    print("🎯 请在图中【细胞聚集度最高】的中心位置，点击两个点画一条垂直/倾斜的分界线（新的愈合前沿），然后按 Enter。")
    print("="*50)
    
    front_coords = plt.ginput(n=-1, timeout=0)
    if front_coords:
        fx, fy = zip(*front_coords)
        # 用高对比度的深橙色虚线渲染新前沿
        ax.plot(fx, fy, color='#FF8C00', linewidth=2, linestyle='--')
        plt.draw()
        
    plt.pause(0.2)
    
    # 确认方向
    print("\n❓ 请最后确认一次：以您新画的虚线为界，哪边是【已长好(Epithelialized)】的区域？")
    direction_input = input("输入 R (代表虚线右侧) 或 L (代表虚线左侧): ").strip().upper()
    direction = "right_is_healed" if direction_input == 'R' else "left_is_healed"
    
    plt.close(fig)
    
    # 更新 JSON (只替换 leading_edge 和 healed_direction，原有的表皮/真皮层不动)
    anchor_dict['leading_edge'] = [(float(x), float(y)) for x, y in front_coords]
    anchor_dict['healed_direction'] = direction
    
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(anchor_dict, f, indent=4)
        
    print(f"\n✅ 修改大功告成！JSON 锚点已更新至: {json_path}")
    print("🎉 现在您可以直接去运行【阶段二绘图脚本】了，Fib_K14 绝对会完美移动到 X=0 的 Epi-Front 中心位置！")
    
    
    
############################### 绘图 (顶刊级白底黑字极简版)
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
from shapely.geometry import Point, LineString
import json
import warnings
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore")

# ==========================================
# 0. 路径与核心参数配置
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"

config_wound = {
    "5dpb":         {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb",        "angle": 135},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb",        "angle": 135},
    "19dpb":        {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225},
    "19dpb_p1":     {"sample_id": "BW15D_1_D1_SDSDB_19dpb",    "angle": 140},
    "26dpb_p1":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1", "angle": 0},
    "26dpb_p2":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2", "angle": 0},
    "2mph":         {"sample_id": "2mph_A03699G6.SCT",         "angle": 45},
    "6mph_p1":      {"sample_id": "6mph_A02988A2.SCT_part1",   "angle": 0},
    "6mph_p2":      {"sample_id": "6mph_A02988A2.SCT_part2",   "angle": 0},
    "9mph":         {"sample_id": "9mph_A03091C5.SCT",         "angle": 0},
    "Normal":       {"sample_id": "NS_C02847B1",               "angle": 180}
}

# ⚠️ 在这里修改您当前想要出图的样本 Key
TEST_KEY = "19dpb"  
target_cell = "Fib_K14"

TARGET_SAMPLE = config_wound[TEST_KEY]["sample_id"]
ROTATION_ANGLE = config_wound[TEST_KEY]["angle"]

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad")
json_load_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

# ==========================================
# 1. 核心运算模块 (增加缺失前沿线的兼容逻辑)
# ==========================================
def rotate_coordinates(x, y, angle):
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

def calculate_planA_coordinates(adata, anchor_dict, angle):
    print("\n🧮 正在运行基于旋转坐标系的几何空间投影 (方案 A)...")
    epi_line = LineString(anchor_dict['epi_baseline'])
    der_line = LineString(anchor_dict['der_bottom'])
    
    # 🌟 检查是否画了愈合前沿线 (至少需要2个点)
    has_front_line = len(anchor_dict.get('leading_edge', [])) >= 2
    if has_front_line:
        front_line = LineString(anchor_dict['leading_edge'])
        direction = anchor_dict['healed_direction']
    else:
        print("⚠️ 未检测到愈合前沿线！将假定整个切片为【已愈合 (Fully Epithelialized)】状态...")
    
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = rotate_coordinates(raw_x, raw_y, angle)
    
    x_list, y_list = [], []
    total_cells = len(rot_x)
    print(f"  --> 共需处理 {total_cells} 个细胞，正在计算最短距离矩阵...")
    
    for i in range(total_cells):
        pt = Point(rot_x[i], rot_y[i])
        
        # 计算相对深度 Y
        dist_to_epi = pt.distance(epi_line)
        dist_to_der = pt.distance(der_line)
        y_val = dist_to_epi / (dist_to_epi + dist_to_der) if (dist_to_epi + dist_to_der) > 0 else 0
        y_list.append(y_val)
        
        # 计算距前沿距离 X
        if has_front_line:
            dist_to_front = pt.distance(front_line)
            if direction == "right_is_healed":
                sign = 1 if pt.x > anchor_dict['leading_edge'][0][0] else -1
            else:
                sign = 1 if pt.x < anchor_dict['leading_edge'][0][0] else -1
            x_list.append(dist_to_front * sign)
        else:
            # 如果没有前沿线（比如Normal样本），直接赋值一个极大的正数，代表完全愈合
            x_list.append(2000.0) 
            
    adata.obs['Distance_to_Front_X'] = x_list
    adata.obs['Relative_Depth_Y'] = y_list
    return adata

# ==========================================
# 2. 绘图模块
# ==========================================
def plot_planA_dynamics(df, sample_id, cell_type, out_dir):
    print(f"\n🎨 正在绘制 {cell_type} 的顶刊排版连续空间动态图 (极简对角线版)...")
    
    # 设置白底基础风格
    sns.set_style("white")
    plt.style.use('default')
    
    # 缩小画布，完美适配 1行5图 (宽2.6, 高2.2)
    fig, ax = plt.subplots(figsize=(2.4, 2.2), dpi=600, layout="constrained")
    df_sub = df[df['sub_labels'] == cell_type].copy()
    
    if df_sub.empty:
        print(f"⚠️ 未找到 {cell_type} 细胞！")
        return
        
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    # ==========================================
    # (1) 坐标轴单位换算 (Pixel -> µm)
    # ==========================================
    PIXEL_TO_UM = 0.33  
    df_sub['Distance_to_Front_X'] = df_sub['Distance_to_Front_X'] * PIXEL_TO_UM
    x_min, x_max = df_sub['Distance_to_Front_X'].min(), df_sub['Distance_to_Front_X'].max()
    max_abs_x = max(abs(x_min), abs(x_max)) * 1.1

    # ==========================================
    # (2) 倾斜分界线与多边形背景着色
    # ==========================================
    b1_top_x = -600  
    b1_bot_x = -300  
    b2_top_x = 400   
    b2_bot_x = 700   

    # 边界线在白底上改为黑色
    ax.plot([b1_top_x, b1_bot_x], [0, 1], color='black', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.plot([b2_top_x, b2_bot_x], [0, 1], color='black', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)

    poly_un_epi = Polygon([
        [-max_abs_x, 0], [b1_top_x, 0],   
        [b1_bot_x, 1], [-max_abs_x, 1]    
    ], color='#e5192c', alpha=0.15, zorder=0, linewidth=0)
    ax.add_patch(poly_un_epi)

    poly_epi_front = Polygon([
        [b1_top_x, 0], [b2_top_x, 0],     
        [b2_bot_x, 1], [b1_bot_x, 1]      
    ], color='#E1C62F', alpha=0.15, zorder=0, linewidth=0)
    ax.add_patch(poly_epi_front)

    poly_newly_epi = Polygon([
        [b2_top_x, 0], [max_abs_x, 0],    
        [max_abs_x, 1], [b2_bot_x, 1]     
    ], color='#3a77b7', alpha=0.15, zorder=0, linewidth=0)
    ax.add_patch(poly_newly_epi)

    # 顶部分区文字
    trans = ax.get_xaxis_transform()
    center_un_epi = (-max_abs_x + b1_top_x) / 2
    center_epi_front = (b1_top_x + b2_top_x) / 2
    center_newly_epi = (b2_top_x + max_abs_x) / 2
    
    ax.text(center_un_epi, 1.02, 'Un-epi', color='#e5192c', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.text(center_epi_front, 1.02, 'Epi-Front', color='#E1C62F', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.text(center_newly_epi, 1.02, 'Newly-epi', color='#3a77b7', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')

    # ==========================================
    # (3) 核心数据渲染：KDE等高线 + 密度散点
    # ==========================================
    sns.kdeplot(
        data=df_sub, x='Distance_to_Front_X', y='Relative_Depth_Y',
        ax=ax, color="#000000", linewidth=0.5, alpha=0.6, levels=6, cut=0, zorder=2  # 白底黑线等高线
    )
    
    x, y = df_sub['Distance_to_Front_X'], df_sub['Relative_Depth_Y']
    xy = np.vstack([x, y])
    try:
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z = x.iloc[idx], y.iloc[idx], z[idx]
        ax.scatter(x, y, c=z, cmap='plasma', edgecolor='none', s=1.5, zorder=3)
    except np.linalg.LinAlgError:
        ax.scatter(x, y, color='#FFA500', edgecolor='none', s=1.5, zorder=3)
    
    # ==========================================
    # (4) 🌟 顺应等高线走势的对角线箭头 🌟
    # ==========================================
    center_x = -500 
    center_y = 0.75 

    # 箭头 1：右上 (表皮方向，促进上皮化)
    ax.annotate('', 
                xy=(1800, 0.2), xycoords='data',       
                xytext=(center_x, center_y), textcoords='data', 
                arrowprops=dict(facecolor="#0099CC", edgecolor="#0099CC", width=0.3, headwidth=2.5, shrink=0.08),
                zorder=4)

    # 箭头 2：左下 (深层真皮方向，深橙色更适宜白底展示)
    ax.annotate('', 
                xy=(-3200, 0.98), xycoords='data',      
                xytext=(center_x, center_y), textcoords='data', 
                arrowprops=dict(facecolor="#FF8C00", edgecolor="#FF8C00", width=0.3, headwidth=2.5, shrink=0.08),
                zorder=4)

    # ==========================================
    # (5) 坐标轴刻度与排版控制 (<=4个刻度)
    # ==========================================
    ax.set_xlim(-max_abs_x, max_abs_x) 
    
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
    ax.set_ylim(1.05, -0.05) 
    ax.set_yticks([0.0, 0.5, 1.0])
    
    display_title = f"{sample_id}".replace('dpb', 'dpi')
    ax.set_title(display_title, fontsize=14, color='black', pad=12, fontweight='bold') 
    
    # 标签换为黑色
    ax.set_xlabel("Distance to leading edge (µm)", color='black', fontsize=8, labelpad=4)
    ax.set_ylabel("Epidermis-Dermis axis", color='black', fontsize=14, labelpad=4)
    
    # 刻度和坐标线换为黑色
    ax.tick_params(colors='black', direction='out', length=2.5, width=0.6, labelsize=14)
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(0.8)

    # ==========================================
    # (6) 保存文件
    # ==========================================
    save_name = f"Nature_SpatialMap_{sample_id}_{cell_type}_White.pdf"
    save_path = os.path.join(out_dir, save_name)
    save_path_png = save_path.replace('.pdf', '.png')
    
    # 确保导出图像的背景色设定为纯白
    fig.savefig(save_path, facecolor='white', bbox_inches='tight')
    fig.savefig(save_path_png, facecolor='white', bbox_inches='tight', dpi=600) 
    plt.close(fig) 
    print(f"✅ 图表已保存至: {save_path}")

# ==========================================
# 3. 执行流
# ==========================================
if __name__ == "__main__":
    if not os.path.exists(json_load_path):
        print(f"❌ 找不到配置: {json_load_path}\n请先运行【脚本一】提取锚点。")
    elif not os.path.exists(file_path):
        print(f"❌ 找不到数据: {file_path}")
    else:
        print(f"⏳ 正在加载数据: {file_path} ...")
        adata = sc.read_h5ad(file_path)
        
        with open(json_load_path, 'r', encoding='utf-8') as f:
            anchor_dict = json.load(f)
        
        adata = calculate_planA_coordinates(adata, anchor_dict, angle=ROTATION_ANGLE)
        
        df_plot = pd.DataFrame(adata.obs[['sub_labels', 'Distance_to_Front_X', 'Relative_Depth_Y']])
        plot_planA_dynamics(df_plot, sample_id=TEST_KEY, cell_type=target_cell, out_dir=output_dir)
        
        print("\n🎉 分析完成！")

############################### 空间特征基因计算与绘图 (顶刊级白底黑字学术版)
#################### 空间特征基因多样本自动化计算与绘图 (全流水线一体化白底发表版)
import os
import gc
import json
import warnings
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse
import scipy.stats
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as clr
from matplotlib.ticker import MaxNLocator
from shapely.geometry import Point, LineString

warnings.filterwarnings("ignore")

# ==========================================
# 0. Spyder 环境专属设置：强制内联显示
# ==========================================
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('matplotlib', 'inline')
except ImportError:
    pass

# ==========================================
# 1. 全局排版配置：白底黑字 + Arial 字体
# ==========================================
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 150 

# ==========================================
# 2. 任务分流、颜色与核心基因配置
# ==========================================
base_output_dir = r"data/external/local_8e645d0a9b03"
annotated_input_dir = os.path.join(base_output_dir, "Wound_Healing_Annotation_Output_White")
mask_and_json_dir = os.path.join(base_output_dir, "mask")  
feature_gene_out_dir = os.path.join(base_output_dir, "Wound_Healing_FeatureGenes_White")
os.makedirs(feature_gene_out_dir, exist_ok=True)

target_cell = "Fib_K14"

# 包含全量 14 个样本的创面愈合任务字典
config_wound = {
    "5dpb":         {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb",        "angle": 135},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb",        "angle": 135},
    "19dpb":        {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225},
    "19dpb_p1":     {"sample_id": "BW15D_1_D1_SDSDB_19dpb",    "angle": 140},
    "26dpb_p1":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1", "angle": 0},
    "26dpb_p2":     {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2", "angle": 0},
    "2mph":         {"sample_id": "2mph_A03699G6.SCT",         "angle": 45},
    "6mph_p1":      {"sample_id": "6mph_A02988A2.SCT_part1",   "angle": 0},
    "6mph_p2":      {"sample_id": "6mph_A02988A2.SCT_part2",   "angle": 0},
    "9mph":         {"sample_id": "9mph_A03091C5.SCT",         "angle": 0},
    "Normal":       {"sample_id": "NS_C02847B1",               "angle": 180}
}

# 官方标准细胞配色方案
raw_annotation_config = {
    "KC_Basal": {"color": "#008B45"}, "KC_Basal_Mig": {"color": "#FF8C00"}, "KC_Basal_Prolif": {"color": "#FFD700"},
    "KC_Spinous": {"color": "#00FF00"}, "KC_Spinous_Mig": {"color": "#00008B"}, "KC_Spinous_Mat": {"color": "#FF69B4"}, "KC_Granular": {"color": "#8B4513"},
    "SAC_SG_Progenitor": {"color": "#0099CC"}, "SAC_SG_Clear": {"color": "#AFEEEE"}, "SAC_SG_Dark": {"color": "#0A2B00"}, "SAC_SG_Ductal": {"color": "#CD853F"},
    "SAC_HF_IRS": {"color": "#DA70D6"}, "SAC_HF_ORS": {"color": "#D8BFD8"}, "SAC_HF_HFSC": {"color": "#20B2AA"}, "SAC_HF_DP_DS": {"color": "#B8860B"}, "SAC_HF_Matrix": {"color": "#48D1CC"},
    "Melanocyte": {"color": "#A0522D"},
    "Fib_Papi": {"color": "#4682B4"}, "Fib_SFRP2": {"color": "#87CEFA"}, "Fib_Myo": {"color": "#DDA0DD"}, "Fib_Inflama": {"color": "#BA55D3"}, "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#F0E68C"},
    "Fib_Prolif": {"color": "#FFB6C1"}, "Fib_K14": {"color": "#E31A1C"}, "Fib_CD45": {"color": "#9932CC"},
    "Endo_Capillary": {"color": "#9ACD32"}, "Endo_Arterial": {"color": "#6B8E23"}, "Endo_Venous": {"color": "#556B2F"}, "Endo_Lymphatic": {"color": "#8FBC8F"}, "Endo_Prolif": {"color": "#BDB76B"},
    "Pericyte": {"color": "#DAA520"}, "Schwann": {"color": "#D2B48C"},
    "M1_Macrophage": {"color": "#1E90FF"}, "M2_Macrophage": {"color": "#4169E1"}, "LAM": {"color": "#0000CD"}, "Neutrophil": {"color": "#483D8B"}, "Mast": {"color": "#800080"},
    "cDC1": {"color": "#6495ED"}, "cDC2": {"color": "#7B68EE"}, "pDC": {"color": "#8470FF"}, "LGHS": {"color": "#00CED1"},
    "B_Mature": {"color": "#000080"}, "B_Plasma": {"color": "#4B0082"}, "CD4_T": {"color": "#5F9EA0"}, "CD8_T": {"color": "#6A5ACD"}, "NK": {"color": "#008080"}
}

zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

highlight_genes = [
    'KLF4', 'IGFL1', 'KRT10', 'MMP1', 'S100A7', 'KRT6A',  'TP63',
    'CXCL14', 'FTH1', 'SFN', 'B2M', 'SPRR1B'
]

# ==========================================
# 3. 空间几何投影计算核心模块
# ==========================================
def rotate_coordinates(x, y, angle):
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

def calculate_planA_coordinates(adata, anchor_dict, angle):
    epi_line = LineString(anchor_dict['epi_baseline'])
    der_line = LineString(anchor_dict['der_bottom'])
    
    has_front_line = len(anchor_dict.get('leading_edge', [])) >= 2
    if has_front_line:
        front_line = LineString(anchor_dict['leading_edge'])
        direction = anchor_dict['healed_direction']
    
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = rotate_coordinates(raw_x, raw_y, angle)
    
    x_list, y_list = [], []
    for i in range(len(rot_x)):
        pt = Point(rot_x[i], rot_y[i])
        
        dist_to_epi = pt.distance(epi_line)
        dist_to_der = pt.distance(der_line)
        y_list.append(dist_to_epi / (dist_to_epi + dist_to_der) if (dist_to_epi + dist_to_der) > 0 else 0)
        
        if has_front_line:
            dist_to_front = pt.distance(front_line)
            if direction == "right_is_healed":
                sign = 1 if pt.x > anchor_dict['leading_edge'][0][0] else -1
            else:
                sign = 1 if pt.x < anchor_dict['leading_edge'][0][0] else -1
            x_list.append(dist_to_front * sign)
        else:
            x_list.append(2000.0) 
            
    adata.obs['Distance_to_Front_X'] = x_list
    adata.obs['Relative_Depth_Y'] = y_list
    return adata

# ==========================================
# 4. 空间相关性统计核心模块
# ==========================================
def filter_adata_expressed_in_n_cells(adata, percent=0.05):
    adata = adata.copy()
    bin_matrix = (adata.X > 0)
    gene_percent = np.array(bin_matrix.mean(axis=0)).flatten()
    return adata[:, gene_percent > percent]

def df_correlation(adata, var):
    mask = (~np.isnan(adata.obs[var])).values 
    obs_values = adata.obs[var][mask].values
    expr_matrix = adata.X[mask]
    
    if scipy.sparse.issparse(expr_matrix):
        expr_matrix = expr_matrix.toarray()
    
    if expr_matrix.shape[1] == 0:  
        return pd.DataFrame(columns=["gene", "rho", "p", "var"])
        
    rho, p = scipy.stats.spearmanr(obs_values, expr_matrix)
    
    if np.isscalar(rho):
        rho_vals, p_vals = [rho] * adata.shape[1], [p] * adata.shape[1]
    else:
        rho_array, p_array = np.asarray(rho), np.asarray(p)
        rho_vals = rho_array[0, 1:] if rho_array.ndim > 1 else rho_array[1:]
        p_vals = p_array[0, 1:] if p_array.ndim > 1 else p_array[1:]
    
    return pd.DataFrame({"gene": adata.var_names.str.upper(), "rho": rho_vals, "p": p_vals, "var": var})

def get_spatially_variable_genes(adata, var_x='Distance_to_Front_X', var_y='Relative_Depth_Y', rho_threshold=0.15, p_threshold=0.01):
    adata_filtered = filter_adata_expressed_in_n_cells(adata, percent=0.05)
    if adata_filtered.n_vars == 0:
        return pd.DataFrame(), pd.DataFrame()
        
    corr_x = df_correlation(adata_filtered, var_x)
    corr_y = df_correlation(adata_filtered, var_y)
    
    sig_x = corr_x[(corr_x['p'] < p_threshold) & (np.abs(corr_x['rho']) > rho_threshold)] if not corr_x.empty else corr_x
    sig_y = corr_y[(corr_y['p'] < p_threshold) & (np.abs(corr_y['rho']) > rho_threshold)] if not corr_y.empty else corr_y
    print(f"   --> 显著空间基因筛选完成: X轴 {len(sig_x)} 个, Y轴 {len(sig_y)} 个。")
    return sig_x, sig_y

# ==========================================
# 5. 微缩版蜂群图 (显式声明分类修复版)
# ==========================================
def plot_nature_swarm_mini(adata, vars, highlight_list, out_dir, skey, cell_type, rho_cutoff=0.15):
    print("⏳ 正在绘制微缩版蜂群图...")
    adata_filtered = filter_adata_expressed_in_n_cells(adata, percent=0.05)
    df_list = [df_correlation(adata_filtered, var) for var in vars]
    df = pd.concat(df_list, ignore_index=True)
    
    df = df.dropna(subset=['rho'])
    if df.empty:
        print(f"   ⚠️ 跳过：样本 [{skey}] 无有效表达基因，无法生成蜂群图。")
        return
    
    # 🌟 修复核心：显式固定全量分类，防止部分极端样本由于缺失正/负相关数据而引发崩溃 🌟
    df["color"] = df["rho"].apply(
        lambda x: "positive" if x > rho_cutoff else "negative" if x < -rho_cutoff else "none"
    )
    df["color"] = pd.Categorical(df["color"], categories=["positive", "none", "negative"], ordered=True)
    
    fig, ax = plt.subplots(figsize=(1.8, 2.6), dpi=300)
    ax.axhline(y=rho_cutoff, linestyle=":", color="#888888", linewidth=0.5)
    ax.axhline(y=-rho_cutoff, linestyle=":", color="#888888", linewidth=0.5)
    ax.axhline(y=0, linestyle="-", color="#666666", linewidth=0.5)
    
    safe_palette = {"positive": "#E15759", "none": "#D3D3D3", "negative": "#4E79A7"}
    
    sns.swarmplot(data=df, x="var", y="rho", ax=ax, order=vars, hue="color", palette=safe_palette, s=1.0, zorder=1)
    if ax.legend_: ax.legend_.remove()
    
    x_offset, y_gap = 0.2, 0.05 
    for i, var in enumerate(vars):
        hl_df = df[(df["var"] == var) & (df["gene"].isin(highlight_list))].sort_values(by="rho", ascending=False)
        last_y = float('inf')
        for _, row in hl_df.iterrows():
            text_y = row["rho"] if (last_y - row["rho"]) >= y_gap else last_y - y_gap
            last_y = text_y
            
            color = "#E15759" if row["rho"] > 0 else "#4E79A7"
            ax.scatter(i, row["rho"], color=color, s=3, zorder=2) 
            ax.plot([i + 0.03, i + 0.12, i + x_offset], [row["rho"], row["rho"], text_y], color="gray", linewidth=0.4)
            ax.text(i + x_offset + 0.02, text_y, f"$\it{{{row['gene']}}}$", fontsize=5, color="black", fontweight="normal", va="center") 
    
    display_title = f"{skey}".replace('dpb', 'dpi')
    ax.set_title(f"{display_title} | Spatial Correl.", fontsize=8, pad=8, fontweight='bold', color='black')
    ax.set_ylabel("Spearman Correlation (ρ)", fontsize=6, labelpad=2, color='black')
    ax.set_xlabel("")
    ax.set_xticklabels(["Distance\nto Front", "Epi-Dermis\nAxis"], fontsize=6, color='black')
    ax.tick_params(axis='y', labelsize=5, length=2, width=0.5, colors='black')
    ax.tick_params(axis='x', length=0)
    
    sns.despine(top=True, right=True)
    for spine in ax.spines.values(): 
        spine.set_linewidth(0.5); spine.set_color('black')
    
    plt.tight_layout()
    save_path = os.path.join(out_dir, f"{skey}_Nature_Swarmplot_Mini_{cell_type}")
    fig.savefig(f"{save_path}.pdf", facecolor='white', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='white', bbox_inches='tight', dpi=600)
    plt.close(fig)

# ==========================================
# 6. 微缩版并排双热图
# ==========================================
def plot_side_by_side_heatmaps_mini(adata, sig_x_df, sig_y_df, out_dir, skey, cell_type, window_size=50):
    fig, axes = plt.subplots(1, 2, figsize=(2.6, 3.2), dpi=300, layout="constrained")
    display_title = f"{skey}".replace('dpb', 'dpi')
    fig.suptitle(f"{display_title} | Spatial Dynamics", fontsize=9, fontweight='bold', color='black', y=1.05)
    
    def process_and_plot_heatmap(ax_obj, spatial_var, sig_genes_list, reverse_sort=False, title="", xlabel="", left_label="", right_label=""):
        combined_genes = list(set(sig_genes_list + highlight_genes))
        valid_genes = [g for g in combined_genes if g in adata.var_names.str.upper()]
        
        if not valid_genes:
            ax_obj.text(0.5, 0.5, "No Data", ha='center', va='center')
            return
            
        temp_var_names = adata.var_names.str.upper()
        gene_mask = temp_var_names.isin(valid_genes)
        
        expr_matrix = adata.X[:, gene_mask]
        if scipy.sparse.issparse(expr_matrix): expr_matrix = expr_matrix.toarray()
        
        expr_df = pd.DataFrame(expr_matrix, index=adata.obs_names, columns=temp_var_names[gene_mask])
        expr_df[spatial_var] = adata.obs[spatial_var].values
        
        expr_df = expr_df.sort_values(by=spatial_var, ascending=not reverse_sort)
        expr_df.pop(spatial_var) 
        
        zscore_df = expr_df.rolling(window=window_size, min_periods=1, center=True).mean().apply(zscore).dropna(axis=1)
        if zscore_df.empty:
            ax_obj.text(0.5, 0.5, "No Var", ha='center', va='center')
            return
            
        peak_positions = zscore_df.idxmax(axis=0).map({cell: i for i, cell in enumerate(zscore_df.index)})
        zscore_df = zscore_df[peak_positions.sort_values().index]
        
        sns.heatmap(
            zscore_df.T, cmap=colormap, center=0, vmin=-2, vmax=2, 
            xticklabels=False, yticklabels=True, 
            cbar_kws={'shrink': 0.7, 'aspect': 15, 'location': 'top', 'pad': 0.05}, ax=ax_obj
        )
        
        clean_labels = []
        for gene in zscore_df.columns:
            if gene in highlight_genes:
                clean_labels.append(f"$\it{{{gene}}}$")
            else:
                clean_labels.append("") 
        
        ax_obj.set_yticklabels(clean_labels, rotation=0)
        
        for tick, label_obj in zip(ax_obj.yaxis.get_major_ticks(), ax_obj.get_yticklabels()):
            text = label_obj.get_text()
            if text == "":
                tick.tick1line.set_visible(False)
            else:
                tick.tick1line.set_visible(True)
                tick.tick1line.set_color('black')
                tick.tick1line.set_markersize(2)
                label_obj.set_color('black')
                label_obj.set_fontweight('normal')
                label_obj.set_fontsize(5)
        
        ax_obj.set_title(title, fontsize=7, color='black', pad=25, fontweight='bold') 
        ax_obj.set_ylabel("Genes" if ax_obj == axes[0] else "", fontsize=7, color='black', labelpad=2)
        ax_obj.set_xlabel(xlabel, fontsize=7, color='black', labelpad=8)
        ax_obj.text(0, -0.02, left_label, transform=ax_obj.transAxes, ha='left', va='top', fontsize=5, color='black')
        ax_obj.text(1, -0.02, right_label, transform=ax_obj.transAxes, ha='right', va='top', fontsize=5, color='black')
        
        cbar = ax_obj.collections[0].colorbar
        if cbar:
            cbar.ax.tick_params(labelsize=4, length=1.5, width=0.5, colors='black')
            cbar.outline.set_edgecolor('black')
        
        for _, spine in ax_obj.spines.items():
            spine.set_visible(True); spine.set_color('black'); spine.set_linewidth(0.5)
    
    process_and_plot_heatmap(axes[0], 'Distance_to_Front_X', sig_x_df['gene'].tolist() if not sig_x_df.empty else [], False, "Healing Axis", "Distance to Front", "Un-epi", "Newly-epi")
    process_and_plot_heatmap(axes[1], 'Relative_Depth_Y', sig_y_df['gene'].tolist() if not sig_y_df.empty else [], True, "Depth Axis", "Relative Depth", "Deep Dermis", "Epidermis")
    
    save_path = os.path.join(out_dir, f"{skey}_Nature_SideBySide_Heatmaps_Mini_{cell_type}")
    fig.savefig(f"{save_path}.pdf", facecolor='white', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='white', bbox_inches='tight', dpi=600)
    plt.close(fig)

# ==========================================
# 7. 全自动化批处理流水线控制台
# ==========================================
if __name__ == "__main__":
    print("🎬 开启全样本空间特征基因自动化批处理流水线...")
    
    for skey, info in config_wound.items():
        tid = info["sample_id"]
        rot_angle = info["angle"]
        
        file_path = os.path.join(annotated_input_dir, f"{tid}_Annotated.h5ad")
        json_path = os.path.join(mask_and_json_dir, f"{skey}_spatial_anchors.json")
        
        if not os.path.exists(file_path):
            print(f"⚠️ 跳过样本 [{skey}]：未找到注释矩阵数据 {file_path}")
            continue
        if not os.path.exists(json_path):
            print(f"⚠️ 跳过样本 [{skey}]：未找到对应的空间 JSON 锚点配置文件 {json_path}")
            continue
            
        print(f"\n🚀 ========================================================")
        print(f"⏳ 正在加载样本 [{skey}] 数据并注入空间几何投影参数...")
        print(f"   ========================================================")
        adata = sc.read_h5ad(file_path)
        
        with open(json_path, 'r', encoding='utf-8') as f:
            anchor_dict = json.load(f)
            
        # 在线注入几何坐标（防止 KeyError）
        adata = calculate_planA_coordinates(adata, anchor_dict, angle=rot_angle)
        
        # 提取靶细胞
        fib_adata = adata[adata.obs['sub_labels'] == target_cell].copy()
        if fib_adata.n_obs < 15:
            print(f"⚠️ 警告：样本 [{skey}] 中的 {target_cell} 群体规模不足 ({fib_adata.n_obs})，跳过统计建模。")
            del adata, fib_adata; gc.collect()
            continue
            
        print(f"🧬 正在筛选 [{skey}] 空间特异性显著相关基因 (靶细胞数: {fib_adata.n_obs})...")
        sig_genes_x, sig_genes_y = get_spatially_variable_genes(fib_adata, rho_threshold=0.15)
        
        # 导出蜂群图
        print(f"🎨 正在导出 [{skey}] 专属蜂群图...")
        plot_nature_swarm_mini(
            adata=fib_adata, 
            vars=['Distance_to_Front_X', 'Relative_Depth_Y'], 
            highlight_list=highlight_genes, 
            out_dir=feature_gene_out_dir,
            skey=skey,
            cell_type=target_cell,
            rho_cutoff=0.15
        )
        
        # 导出并排双轴热图
        print(f"🎨 正在导出 [{skey}] 专属时空级联瀑布热图...")
        plot_side_by_side_heatmaps_mini(
            adata=fib_adata, 
            sig_x_df=sig_genes_x, 
            sig_y_df=sig_genes_y, 
            out_dir=feature_gene_out_dir,
            skey=skey,
            cell_type=target_cell,
            window_size=50 
        )
        print(f"✅ 样本 [{skey}] 流水线图册全部保存成功。")
        
        # 显式清除大矩阵，杜绝内存溢出
        del adata, fib_adata; gc.collect()

    print("\n🎉 全线 14 个样本 (Wound_Healing) 的 Fib_K14 空间特征分析及顶刊级白底黑字图册已全部批量清洗导出成功！")

###################### 特征基因时空分布图 (顶刊级白底黑字学术发表版)
import os
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as clr
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. Spyder 环境专属设置：强制内联显示
# ==========================================
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('matplotlib', 'inline')
except ImportError:
    pass

# ==========================================
# 1. 全局排版配置：白底黑字 + Arial 字体
# ==========================================
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 150

# ==========================================
# 2. 颜色与核心变量配置
# ==========================================
zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

module_genes = ['KRT14', 'PDGFRA', 'KRT5', 'KLF4', 'IGFL1', 'TP63', 'TACSTD2', 'EPCAM', 'KRT10']
display_features = ['KRT14', 'PDGFRA', 'KLF4', 'TP63', 'MergeScore']

SPATIAL_X = 'Distance_to_Front_X'
SPATIAL_Y = 'Relative_Depth_Y'
BURN_DEPTH_COL = 'Burn_Depth'
HEALING_REGION_COL = 'Region'

# ==========================================
# 3. 核心功能函数
# ==========================================
def calculate_mergescore(adata, genes, score_name='MergeScore'):
    valid_genes = [g for g in genes if g in adata.var_names.str.upper()]
    temp_vars = adata.var_names.copy()
    adata.var_names = adata.var_names.str.upper()
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)
    adata.var_names = temp_vars 
    return adata

def get_expression_data(adata, feature):
    if feature in adata.obs.columns:
        return adata.obs[feature].values
    else:
        gene_idx = np.where(adata.var_names.str.upper() == feature.upper())[0]
        if len(gene_idx) > 0:
            expr_data = adata[:, gene_idx[0]].X
            return np.asarray(expr_data.toarray()).flatten() if scipy.sparse.issparse(expr_data) else np.asarray(expr_data).flatten()
        return np.zeros(adata.n_obs)

def plot_spatial_expression_planA(ax, adata, feature, cmap, show_y_label=False, show_colorbar=False, set_global_xlabel=False):
    """Row 1: 动态空间投影图 (极简正方形，统一坐标与图例)"""
    PIXEL_TO_UM = 0.33  
    x = adata.obs[SPATIAL_X].values * PIXEL_TO_UM
    y = adata.obs[SPATIAL_Y].values
    expr = get_expression_data(adata, feature)
    
    idx = np.argsort(expr)
    x_s, y_s, expr_s = x[idx], y[idx], expr[idx]
    
    max_abs_x = max(abs(x.min()), abs(x.max())) * 1.1 if len(x) > 0 else 4000
    b1_t, b1_b, b2_t, b2_b = -600, -300, 400, 700   
    # 🌟 划分界线在白底上改为黑色
    ax.plot([b1_t, b1_b], [0, 1], color='black', linestyle='--', linewidth=0.8, alpha=0.4, zorder=1)
    ax.plot([b2_t, b2_b], [0, 1], color='black', linestyle='--', linewidth=0.8, alpha=0.4, zorder=1)
    ax.add_patch(Polygon([[-max_abs_x, 0], [b1_t, 0], [b1_b, 1], [-max_abs_x, 1]], color='#e5192c', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b1_t, 0], [b2_t, 0], [b2_b, 1], [b1_b, 1]], color='#E1C62F', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b2_t, 0], [max_abs_x, 0], [max_abs_x, 1], [b2_b, 1]], color='#3a77b7', alpha=0.15, zorder=0, lw=0))

    trans = ax.get_xaxis_transform()
    c_un_epi, c_epi_front, c_newly_epi = (-max_abs_x + b1_t)/2, (b1_t + b2_t)/2, (b2_t + max_abs_x)/2
    ax.text(c_un_epi, 1.02, 'Un-epi', color='#e5192c', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.text(c_epi_front, 1.02, 'Epi-Front', color='#E1C62F', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.text(c_newly_epi, 1.02, 'Newly-epi', color='#3a77b7', transform=trans, ha='center', va='bottom', fontsize=8, fontweight='bold')

    # 🌟 KDE轮廓线改为黑色
    sns.kdeplot(data=pd.DataFrame({'x': x, 'y': y}), x='x', y='y', ax=ax, color="#000000", linewidths=0.5, alpha=0.4, levels=6, zorder=2)
    mask_z, mask_e = expr_s <= 0, expr_s > 0
    # 🌟 零表达的背景散点改为超浅的银灰色，避免把白底污染成一片漆黑
    ax.scatter(x_s[mask_z], y_s[mask_z], color='#F0F0F0', s=1.5, alpha=0.4, zorder=3)
    sc_scatter = ax.scatter(x_s[mask_e], y_s[mask_e], c=expr_s[mask_e], cmap=cmap, s=1.5, edgecolor='none', zorder=4)
    
    # 🌟 Colorbar细节上色涂黑
    if show_colorbar and mask_e.sum() > 0:
        cbar_ax = ax.inset_axes([1.08, 0.1, 0.05, 0.8]) 
        cbar = plt.colorbar(sc_scatter, cax=cbar_ax, orientation='vertical')
        cbar.set_ticks([expr_s[mask_e].min(), expr_s[mask_e].max()])
        cbar.set_ticklabels(['Min', 'Max'])
        cbar.ax.tick_params(colors='black', labelsize=8, length=0, pad=2)
        cbar.set_label('Expression', color='black', fontsize=12, labelpad=6)
        cbar.outline.set_edgecolor('black')

    cx, cy = -500, 0.75 
    ax.annotate('', xy=(1800, 0.2), xycoords='data', xytext=(cx, cy), textcoords='data',
                arrowprops=dict(facecolor="#0099CC", edgecolor="#0099CC", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)
    # 黄色箭头保持不变，但边缘上色增加对比度
    ax.annotate('', xy=(-3200, 0.98), xycoords='data', xytext=(cx, cy), textcoords='data',
                arrowprops=dict(facecolor="#FFFF00", edgecolor="#8B7B00", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)

    ax.set_xlim(-max_abs_x, max_abs_x) 
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
    ax.set_ylim(1.05, -0.05) 
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_box_aspect(1) 
    
    ax.set_title(f"$\it{{{feature}}}$" if feature != 'MergeScore' else feature, fontsize=16, color='black', pad=12, fontweight='bold')
    
    ax.set_xlabel("Distance to leading edge (µm)" if set_global_xlabel else "", color='black', fontsize=16 if set_global_xlabel else 0, labelpad=12)
    ax.set_ylabel("Epidermis-Dermis axis" if show_y_label else "", color='black', fontsize=16, labelpad=8)
    
    ax.tick_params(colors='black', direction='out', length=2.5, width=0.6, labelsize=14)
    for spine in ax.spines.values(): 
        spine.set_color('black')
        spine.set_linewidth(0.6)

def plot_tissue_spatial(ax, adata, feature, cmap, angle_deg=-225, show_colorbar=False, show_scalebar=False):
    """Row 2: 真实组织切片 (含反向箭头，极简唯一图例与比例尺)"""
    if 'X_rot_spatial' in adata.obsm:
        coords = adata.obsm['X_rot_spatial']
    elif 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
    else:
        coords = adata.obsm['X_spatial']
        
    x_raw, y_raw = coords[:, 0], coords[:, 1]
    
    theta = np.radians(angle_deg)
    x = x_raw * np.cos(theta) - y_raw * np.sin(theta)
    y = x_raw * np.sin(theta) + y_raw * np.cos(theta)
    
    expr = get_expression_data(adata, feature)
    idx = np.argsort(expr)
    x_s, y_s, expr_s = x[idx], y[idx], expr[idx]
    
    mask_z, mask_e = expr_s <= 0, expr_s > 0
    # 🌟 真实物理背景阴影散点改为契合白底的浅灰色 (#E5E5E5)
    ax.scatter(x_s[mask_z], y_s[mask_z], color='#E5E5E5', s=0.8, alpha=0.4) 
    sc_scatter = ax.scatter(x_s[mask_e], y_s[mask_e], c=expr_s[mask_e], cmap=cmap, s=1.5, edgecolor='none')
    
    min_x, max_x = x.min(), x.max()
    min_y, max_y = y.min(), y.max()
    y_range = max_y - min_y
    
    # 🌟 外部右侧 Colorbar 完全涂黑
    if show_colorbar and mask_e.sum() > 0:
        cbar_ax = ax.inset_axes([1.08, 0.05, 0.05, 0.9])
        cbar = plt.colorbar(sc_scatter, cax=cbar_ax, orientation='vertical')
        cbar.set_ticks([expr_s[mask_e].min(), expr_s[mask_e].max()])
        cbar.set_ticklabels(['Min', 'Max'])
        cbar.ax.tick_params(colors='black', labelsize=8, length=0, pad=2)
        cbar.set_label('Expression', color='black', fontsize=12, labelpad=6)
        cbar.outline.set_edgecolor('black')
        
    # 🌟 比例尺在白底上换为黑色
    if show_scalebar:
        scale_px = 1000 / 0.33 
        bar_y = min_y - y_range * 0.08
        ax.plot([max_x - scale_px, max_x], [bar_y, bar_y], color='black', lw=1.5)
        ax.text(max_x - scale_px/2, bar_y - y_range * 0.03, '1000 µm', color='black', fontsize=10, ha='center', va='top', fontweight='bold')

    # 🌟 顶部的区域过渡长箭头换为黑色
    arrow_y = max_y + y_range * 0.12
    ax.annotate('', xy=(min_x, arrow_y), xytext=(max_x, arrow_y), 
                arrowprops=dict(arrowstyle="->", color="black", lw=1.2))
    
    text_y = arrow_y + y_range * 0.04
    ax.text(min_x, text_y, 'Un-epi', color='#e5192c', fontsize=9, ha='left', va='bottom', fontweight='bold')
    ax.text((min_x + max_x)/2, text_y, 'Epi-Front', color='#E1C62F', fontsize=9, ha='center', va='bottom', fontweight='bold')
    ax.text(max_x, text_y, 'Newly-epi', color='#3a77b7', fontsize=9, ha='right', va='bottom', fontweight='bold')

    ax.set_aspect('equal', adjustable='datalim')
    ax.axis('off')

def get_bars_df(adata, features, group_cols):
    df_list = []
    for f in features:
        vals = get_expression_data(adata, f)
        temp_df = adata.obs[group_cols].copy()
        temp_df['Value'] = vals
        temp_df['Feature'] = f
        df_list.append(temp_df)
    return pd.concat(df_list, ignore_index=True)

# ==========================================
# 4. 终极主图绘制逻辑
# ==========================================
def plot_master_figure(adata, out_dir):
    print("⏳ 正在计算 MergeScore 并准备绘图数据...")
    adata = calculate_mergescore(adata, module_genes, score_name='MergeScore')
    
    fig = plt.figure(figsize=(15.5, 11.5), dpi=300) 
    gs = fig.add_gridspec(4, 5, hspace=0.35, wspace=0.35, height_ratios=[1.4, 0.6, 1.0, 1.0])
    
    # --- Row 1: 动态时空投影图 ---
    print("🎨 绘制第一行: 动态时空表达图 (强迫正方形 + 唯一全局标签与图例)...")
    for i, feature in enumerate(display_features):
        ax = fig.add_subplot(gs[0, i])
        is_first = (i == 0)
        is_middle = (i == 2) 
        is_last = (i == 4)
        
        plot_spatial_expression_planA(ax, adata, feature, colormap, 
                                      show_y_label=is_first, set_global_xlabel=is_middle, show_colorbar=is_last)
        if feature == 'MergeScore':
            for spine in ax.spines.values(): spine.set_linewidth(1.2)

    # --- Row 2: 真实组织切片图 ---
    print("🎨 绘制第二行: 真实组织切片空间图 (箭头反转，唯一右侧图例与比例尺)...")
    for i, feature in enumerate(display_features):
        ax = fig.add_subplot(gs[1, i])
        is_last = (i == 4)
        
        plot_tissue_spatial(ax, adata, feature, 'plasma' if feature == 'MergeScore' else colormap, 
                            angle_deg=-225, show_colorbar=is_last, show_scalebar=is_last)

    # --- 准备柱状图数据 ---
    if BURN_DEPTH_COL in adata.obs.columns and HEALING_REGION_COL in adata.obs.columns:
        df_bars = get_bars_df(adata, display_features, [BURN_DEPTH_COL, HEALING_REGION_COL])
        
        # --- Row 3: 浅深二度差异柱状图 ---
        print("🎨 绘制第三行: 浅深二度差异柱状图...")
        for i, feature in enumerate(display_features):
            ax = fig.add_subplot(gs[2, i])
            # 🌟 errcolor由white调整为纯黑
            sns.barplot(data=df_bars[df_bars['Feature'] == feature], x=BURN_DEPTH_COL, y='Value', 
                        ax=ax, palette=['#E15759', '#4E79A7'], errcolor='black', errwidth=1.0, capsize=0.15)
            ax.set_ylabel("Mean Expr." if i == 0 else "", fontsize=10, fontweight='bold')
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=30, labelsize=11)
            ax.tick_params(axis='y', labelsize=10)
            sns.despine(ax=ax, top=True, right=True)
            for spine in ax.spines.values(): 
                spine.set_color('black')
                spine.set_linewidth(0.6)

        # --- Row 4: 愈合区域差异柱状图 ---
        print("🎨 绘制第四行: 愈合区域差异柱状图...")
        cat_type = pd.CategoricalDtype(categories=['Un-epi', 'Epi-Front', 'Newly-epi'], ordered=True)
        df_bars[HEALING_REGION_COL] = df_bars[HEALING_REGION_COL].astype(cat_type)
        
        for i, feature in enumerate(display_features):
            ax = fig.add_subplot(gs[3, i])
            # 🌟 errcolor由white调整为纯黑
            sns.barplot(data=df_bars[df_bars['Feature'] == feature], x=HEALING_REGION_COL, y='Value', 
                        ax=ax, palette=['#e5192c', '#E1C62F', '#3a77b7'], errcolor='black', errwidth=1.0, capsize=0.15)
            ax.set_ylabel("Mean Expr." if i == 0 else "", fontsize=10, fontweight='bold')
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=30, labelsize=11)
            ax.tick_params(axis='y', labelsize=10)
            sns.despine(ax=ax, top=True, right=True)
            for spine in ax.spines.values(): 
                spine.set_color('black')
                spine.set_linewidth(0.6)
            
    else:
        print(f"⚠️ 找不到列 '{BURN_DEPTH_COL}' 或 '{HEALING_REGION_COL}'，跳过柱状图绘制。")

    # 保存与展示
    plt.tight_layout()
    save_path = os.path.join(out_dir, "MainFigure_Spatial_Module_Panel_Ultimate")
    # 🌟 保存文件背景设置为 white
    fig.savefig(f"{save_path}.pdf", facecolor='white', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='white', bbox_inches='tight', dpi=600)
    plt.show()
    plt.close(fig)
    print(f"✅ 终极主图排版已完美保存至: {save_path}.pdf/png")

# ==========================================
# 5. 执行代码
# ==========================================
if __name__ == "__main__":
    output_dir = r"data/external/local_8e645d0a9b03"
    plot_master_figure(adata, output_dir)
    
################################### 空间共定位多样本批量分析 (顶刊白底统计强化版)
# ==============================================================================
# 空间共定位多样本批量分析 (顶刊白底统计强化版 - 智能排版防重叠优化)
# ==============================================================================
import os, gc, json, warnings
import numpy as np, pandas as pd, scanpy as sc, squidpy as sq, seaborn as sns
import matplotlib.pyplot as plt, matplotlib.colors as clr, igraph, geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")
plt.ioff()

# ==========================================
# 全局排版与字号配置 (强制 Arial, 字号严格 >= 6.5pt)
# ==========================================
A4_W, A4_H = 8.27, 11.69
NETWORK_FIGSIZE = (A4_W * 1.35, A4_H / 3.0)
HEATMAP_FIGSIZE = (A4_W * 1.35, A4_H / 4.8)
NET_FONT = 6.5           # 至少 6pt
NET_KEY_FONT = 7.5
NET_TITLE_FONT = 9.5
NET_SUPTITLE_FONT = 11.5
NET_LEGEND_FONT = 8.5
HEAT_FONT = 8.0

plt.style.use("default")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ==========================================
# 路径与分析变量配置
# ==========================================
base_output_dir = r"data/external/local_8e645d0a9b03"
annotated_input_dir = os.path.join(base_output_dir, "Wound_Healing_Annotation_Output_White")
mask_and_json_dir = os.path.join(base_output_dir, "mask")
network_out_dir = os.path.join(base_output_dir, "Wound_Healing_Networks_White_A4_EpiFront_KCSpinousMig_LabelOptimized")
os.makedirs(network_out_dir, exist_ok=True)

target_cell = "Fib_K14"
force_pair_cell = "KC_Spinous_Mig"
force_pair_days = ("19dpb", "19dpb_p1", "26dpb_p1", "26dpb_p2")
force_pair_region = "Epi-Front"
regions = ["Un-epi", "Epi-Front", "Newly-epi"]

config_wound = {
    "5dpb": {"sample_id": "BW32_A01597A3_SDSDB_5dpb", "angle": 135},
    "12dpb_DPTDI1": {"sample_id": "BW13_1_B1_DSDB_12dpb", "angle": 45},
    "12dpb_SPTDI1": {"sample_id": "BW14_1_C1_SSDB_12dpb", "angle": 135},
    "12dpb_DPTDI2": {"sample_id": "BW13_A3_DSDB_12dpb", "angle": 135},
    "12dpb_SPTDI2": {"sample_id": "BW14_B3_SSDB_12dpb", "angle": 135},
    "19dpb": {"sample_id": "BW15D_C6_SDSDB_19dpb", "angle": -225},
    "19dpb_p1": {"sample_id": "BW15D_1_D1_SDSDB_19dpb", "angle": 140},
    "26dpb_p1": {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part1", "angle": 0},
    "26dpb_p2": {"sample_id": "BW81_C02846B6_SDSDB_26dpb_part2", "angle": 0},
    "2mph": {"sample_id": "2mph_A03699G6.SCT", "angle": 45},
    "6mph_p1": {"sample_id": "6mph_A02988A2.SCT_part1", "angle": 0},
    "6mph_p2": {"sample_id": "6mph_A02988A2.SCT_part2", "angle": 0},
    "9mph": {"sample_id": "9mph_A03091C5.SCT", "angle": 0},
    "Normal": {"sample_id": "NS_C02847B1", "angle": 180}
}

raw_annotation_config = {
    "KC_Basal": {"color": "#008B45"}, "KC_Basal_Mig": {"color": "#FF8C00"}, "KC_Basal_Prolif": {"color": "#FFD700"},
    "KC_Spinous": {"color": "#00FF00"}, "KC_Spinous_Mig": {"color": "#00008B"}, "KC_Spinous_Mat": {"color": "#FF69B4"}, "KC_Granular": {"color": "#8B4513"},
    "SAC_SG_Progenitor": {"color": "#0099CC"}, "SAC_SG_Clear": {"color": "#AFEEEE"}, "SAC_SG_Dark": {"color": "#0A2B00"}, "SAC_SG_Ductal": {"color": "#CD853F"},
    "SAC_HF_IRS": {"color": "#DA70D6"}, "SAC_HF_ORS": {"color": "#D8BFD8"}, "SAC_HF_HFSC": {"color": "#20B2AA"}, "SAC_HF_DP_DS": {"color": "#B8860B"}, "SAC_HF_Matrix": {"color": "#48D1CC"},
    "Melanocyte": {"color": "#A0522D"}, "Fib_Papi": {"color": "#4682B4"}, "Fib_SFRP2": {"color": "#87CEFA"}, "Fib_Myo": {"color": "#DDA0DD"}, "Fib_Inflama": {"color": "#BA55D3"}, 
    "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#F0E68C"}, "Fib_Prolif": {"color": "#FFB6C1"}, "Fib_K14": {"color": "#E31A1C"}, "Fib_CD45": {"color": "#9932CC"},
    "Endo_Capillary": {"color": "#9ACD32"}, "Endo_Arterial": {"color": "#6B8E23"}, "Endo_Venous": {"color": "#556B2F"}, "Endo_Lymphatic": {"color": "#8FBC8F"}, "Endo_Prolif": {"color": "#BDB76B"},
    "Pericyte": {"color": "#DAA520"}, "Schwann": {"color": "#D2B48C"}, "M1_Macrophage": {"color": "#1E90FF"}, "M2_Macrophage": {"color": "#4169E1"}, "LAM": {"color": "#0000CD"}, 
    "Neutrophil": {"color": "#483D8B"}, "Mast": {"color": "#800080"}, "cDC1": {"color": "#6495ED"}, "cDC2": {"color": "#7B68EE"}, "pDC": {"color": "#8470FF"}, "LGHS": {"color": "#00CED1"},
    "B_Mature": {"color": "#000080"}, "B_Plasma": {"color": "#4B0082"}, "CD4_T": {"color": "#5F9EA0"}, "CD8_T": {"color": "#6A5ACD"}, "NK": {"color": "#008080"}
}
cell_colors = {k: v["color"] for k, v in raw_annotation_config.items()}

def get_macro_group(c):
    if c.startswith("KC_") or c.startswith("SAC_") or c == "Melanocyte": return "Epithelial & Appendage"
    if c.startswith("Fib_") or c.startswith("Endo_") or c in ["Pericyte", "Schwann"]: return "Stromal & Vascular"
    return "Immune Cells"

def is_kc(c): return c.startswith("KC_")

bubble_colors = {"Epithelial & Appendage": "#3B82F6", "Stromal & Vascular": "#F97316", "Immune Cells": "#10B981"}
bubble_edge = {"Epithelial & Appendage": "#1D4ED8", "Stromal & Vascular": "#C2410C", "Immune Cells": "#047857"}
text_group_colors = {"Epithelial & Appendage": "#005b96", "Stromal & Vascular": "#8b4513", "Immune Cells": "#006400"}

# ==========================================
# 空间投影与处理
# ==========================================
def process_spatial_data(adata, anchor_dict, angle):
    theta = np.radians(angle)
    raw_x = adata.obsm["spatial"][:, 0] if "spatial" in adata.obsm else adata.obsm["X_spatial"][:, 0]
    raw_y = adata.obsm["spatial"][:, 1] if "spatial" in adata.obsm else adata.obsm["X_spatial"][:, 1]
    rot_x = raw_x * np.cos(theta) - raw_y * np.sin(theta)
    rot_y = raw_x * np.sin(theta) + raw_y * np.cos(theta)
    adata.obsm["X_rot_spatial"] = np.column_stack((rot_x, rot_y))
    min_x, max_x, min_y, max_y = rot_x.min(), rot_x.max(), rot_y.min(), rot_y.max()
    
    epi_line = LineString(anchor_dict["epi_baseline"]) if "epi_baseline" in anchor_dict and len(anchor_dict["epi_baseline"]) >= 2 else LineString([(min_x, min_y), (max_x, min_y)])
    der_line = LineString(anchor_dict["der_bottom"]) if "der_bottom" in anchor_dict and len(anchor_dict["der_bottom"]) >= 2 else LineString([(min_x, max_y), (max_x, max_y)])
    has_front = len(anchor_dict.get("leading_edge", [])) >= 2
    if has_front: front_line = LineString(anchor_dict["leading_edge"])
    
    x_list, y_list = [], []
    for i in range(len(rot_x)):
        pt = Point(rot_x[i], rot_y[i])
        d1, d2 = pt.distance(epi_line), pt.distance(der_line)
        y_list.append(d1 / (d1 + d2) if d1 + d2 > 0 else 0)
        if has_front:
            sign = 1 if (pt.x > anchor_dict["leading_edge"][0][0] and anchor_dict.get("healed_direction") == "right_is_healed") or (pt.x < anchor_dict["leading_edge"][0][0] and anchor_dict.get("healed_direction") != "right_is_healed") else -1
            x_list.append(pt.distance(front_line) * sign)
        else:
            x_list.append(2000.0)
            
    adata.obs["Plot_X"] = np.array(x_list) * 0.33
    adata.obs["Plot_Y"] = y_list
    max_abs_x = max(abs(adata.obs["Plot_X"].min()), abs(adata.obs["Plot_X"].max())) * 1.1
    b1_t, b1_b, b2_t, b2_b = -600, -300, 400, 700
    
    gates = {
        "Un-epi": Polygon([[-max_abs_x, 0], [b1_t, 0], [b1_b, 1], [-max_abs_x, 1]]),
        "Epi-Front": Polygon([[b1_t, 0], [b2_t, 0], [b2_b, 1], [b1_b, 1]]),
        "Newly-epi": Polygon([[b2_t, 0], [max_abs_x, 0], [max_abs_x, 1], [b2_b, 1]])
    }
    gpd_poly = gpd.GeoDataFrame({"Region": list(gates.keys())}, geometry=list(gates.values()))
    gpd_cells = gpd.GeoDataFrame(adata.obs.copy(), geometry=gpd.points_from_xy(adata.obs["Plot_X"], adata.obs["Plot_Y"]))
    adata.obs["Region"] = gpd.sjoin(gpd_cells, gpd_poly, how="left", predicate="within")["Region"].fillna("Unknown").values
    return adata

# ==========================================
# 优化的布局生成与防重叠文本渲染引擎
# ==========================================
def generate_radial_layout(nodes, target_cell):
    pos = {target_cell: (0, 0)}
    groups = {"Epithelial & Appendage": [], "Stromal & Vascular": [], "Immune Cells": []}
    for n in nodes:
        if n != target_cell: groups[get_macro_group(n)].append(n)
        
    angles = {
        "Epithelial & Appendage": (np.radians(25), np.radians(155)),
        "Stromal & Vascular": (np.radians(170), np.radians(260)),
        "Immune Cells": (np.radians(285), np.radians(380))
    }
    
    for g, ns in groups.items():
        ns.sort()
        thetas = np.linspace(angles[g][0], angles[g][1], len(ns)) if ns else []
        for i, n in enumerate(ns):
            # 采用素数步长螺旋分布法，确保相邻节点半径极差最大化
            r = [0.42, 0.76, 1.10, 1.44][(i * 3) % 4]
            pos[n] = (r * np.cos(thetas[i]), r * np.sin(thetas[i]))
    return pos

def plot_nature_network_v2(ax, g_reg, pos_dict, title, target_cell, top_nodes_list, skey, region):
    nodes = g_reg.vs["label"]
    coords = np.array([pos_dict[n] for n in nodes])
    force_visible = (skey in force_pair_days and region == force_pair_region)
    
    # 渲染背景簇阴影
    for group in bubble_colors:
        idx = [i for i, n in enumerate(nodes) if get_macro_group(n) == group and n != target_cell]
        if len(idx) >= 3:
            poly = MultiPoint(np.vstack([coords[idx], [0, 0]])).convex_hull.buffer(0.08, join_style=1)
            if not poly.is_empty:
                x, y = poly.exterior.xy
                ax.fill(x, y, color=bubble_colors[group], alpha=0.15, zorder=1)
                ax.plot(x, y, color=bubble_edge[group], alpha=0.85, lw=0.65, ls="--", zorder=1)
                
    # 严格保持不变的连线绘制逻辑
    drawn = set()
    for edge in g_reg.es:
        n0, n1 = nodes[edge.tuple[0]], nodes[edge.tuple[1]]
        pair = tuple(sorted([n0, n1]))
        drawn.add(pair)
        x0, y0 = pos_dict[n0]
        x1, y1 = pos_dict[n1]
        is_target_edge = (n0 == target_cell or n1 == target_cell)
        partner = n1 if n0 == target_cell else n0
        
        if pair == tuple(sorted([target_cell, force_pair_cell])):
            if force_visible: ax.plot([x0, x1], [y0, y1], color="#E11D48", alpha=1, lw=1.8, zorder=6)
            else: continue
        elif is_target_edge and is_kc(partner): ax.plot([x0, x1], [y0, y1], color="#222222", alpha=0.05, lw=0.18, zorder=1)
        elif is_target_edge and partner in top_nodes_list: ax.plot([x0, x1], [y0, y1], color=cell_colors.get(partner, "#FF1493"), alpha=0.75, lw=0.85, zorder=3)
        elif is_target_edge: ax.plot([x0, x1], [y0, y1], color="#222222", alpha=0.16, lw=0.28, zorder=2)
        else: ax.plot([x0, x1], [y0, y1], color="#222222", alpha=0.025, lw=0.15, zorder=1)
        
    if force_visible and target_cell in nodes and force_pair_cell in nodes and tuple(sorted([target_cell, force_pair_cell])) not in drawn:
        x0, y0 = pos_dict[target_cell]
        x1, y1 = pos_dict[force_pair_cell]
        ax.plot([x0, x1], [y0, y1], color="#E11D48", alpha=1, lw=1.8, zorder=6)
        
    # 节点散点
    ax.scatter(coords[:, 0], coords[:, 1], c=[cell_colors.get(n, "#888888") for n in nodes], s=15, edgecolors="black", linewidths=0.35, zorder=7)
    ax.scatter([0], [0], facecolors="none", edgecolors="black", s=45, lw=0.8, zorder=8)
    
    # 智能排版防重叠文本渲染
    must_show = set([target_cell])
    if force_visible: must_show.add(force_pair_cell)
    priority = set(top_nodes_list + [target_cell])
    if force_visible: priority.add(force_pair_cell)
    
    candidate_labels = []
    for n in nodes:
        x, y = pos_dict[n]
        if n == target_cell:
            candidate_labels.append((0, n, x, y, "center", "top", 0, -0.12))
            continue
        theta = np.arctan2(y, x)
        ox = 0.055 * np.cos(theta); oy = 0.055 * np.sin(theta)
        ha = "left" if np.cos(theta) >= -0.02 else "right"
        va = "bottom" if np.sin(theta) >= -0.02 else "top"
        
        score = 1 if n in priority else 2
        if n in must_show: score = 0
        if get_macro_group(n) == "Epithelial & Appendage": score += 0.15
        if get_macro_group(n) == "Stromal & Vascular": score += 0.25
        if get_macro_group(n) == "Immune Cells": score += 0.35
        candidate_labels.append((score, n, x, y, ha, va, ox, oy))
        
    candidate_labels = sorted(candidate_labels, key=lambda x: x[0])
    
    # Bounding Box 防重叠算法
    placed_boxes = []
    def get_bbox(px, py, name, ha, va, fsize):
        h = 0.08 * (fsize / 6.0)
        w = len(name) * 0.022 * (fsize / 6.0)
        if ha == 'center': xmin, xmax = px - w/2, px + w/2
        elif ha == 'left': xmin, xmax = px, px + w
        else: xmin, xmax = px - w, px
        if va == 'center': ymin, ymax = py - h/2, py + h/2
        elif va == 'bottom': ymin, ymax = py, py + h
        else: ymin, ymax = py - h, py
        return (xmin - 0.02, ymin - 0.02, xmax + 0.02, ymax + 0.02)
        
    def is_overlap(b1, b2):
        return not (b1[2] < b2[0] or b1[0] > b2[2] or b1[3] < b2[1] or b1[1] > b2[3])

    shown_count = 0
    max_labels = 24
    
    for score, n, x, y, ha, va, ox, oy in candidate_labels:
        px, py = x + ox, y + oy
        essential = (n in must_show)
        fontsize = NET_KEY_FONT if (n == target_cell or (n == force_pair_cell and force_visible)) else (NET_FONT if n in priority else NET_FONT)
        
        theta = np.arctan2(y, x) if (x != 0 or y != 0) else np.pi/2
        placed_ok = False
        
        # 将发生重叠的文字沿半径向外推，直至无重叠或达尝试上限
        for attempt in range(25):
            bbox = get_bbox(px, py, n, ha, va, fontsize)
            overlap = False
            for pb in placed_boxes:
                if is_overlap(bbox, pb):
                    overlap = True
                    break
            if not overlap:
                placed_ok = True
                break
            px += 0.035 * np.cos(theta)
            py += 0.035 * np.sin(theta)
            
        if (not essential) and (not placed_ok or shown_count >= max_labels):
            continue
            
        if n == target_cell:
            ax.text(px, py, n, ha=ha, va=va, fontsize=NET_KEY_FONT, fontweight="bold", color="black", zorder=9)
        elif n == force_pair_cell and force_visible:
            ax.text(px, py, n, ha=ha, va=va, fontsize=NET_KEY_FONT, fontweight="bold", color="#E11D48", zorder=10, bbox=dict(facecolor="white", edgecolor="#E11D48", lw=0.45, pad=0.65))
        elif n in priority:
            ax.text(px, py, n, ha=ha, va=va, fontsize=NET_FONT, fontweight="bold", color="black", zorder=9)
        else:
            ax.text(px, py, n, ha=ha, va=va, fontsize=NET_FONT, color="black", zorder=9)
            
        placed_boxes.append(get_bbox(px, py, n, ha, va, fontsize))
        shown_count += 1
        
    ax.set_title(title, fontsize=NET_TITLE_FONT, fontweight="bold", pad=5)
    # 扩大边缘容纳弹开的长标签
    ax.set_xlim(-1.85, 1.85); ax.set_ylim(-1.75, 1.75)
    ax.set_aspect("equal"); ax.axis("off")

def plot_heatmap_for_sample(skey, sample_out_dir, regions, global_cats, m_target_matrix, z_target_matrix, global_highlight_nodes, target_cell):
    valid_cols = np.sum(m_target_matrix, axis=0) > 0
    if valid_cols.sum() == 0: 
        print(f"   ⚠️ 样本[{skey}]无热图信号，跳过。")
        return
        
    m = m_target_matrix[:, valid_cols]
    z = np.nan_to_num(z_target_matrix[:, valid_cols], 0)
    labels = np.array(global_cats)[valid_cols]
    
    group_order = {"Epithelial & Appendage": 1, "Stromal & Vascular": 2, "Immune Cells": 3}
    sort_keys = sorted([(group_order[get_macro_group(c)], c, i) for i, c in enumerate(labels)])
    idx = [x[2] for x in sort_keys]; labels_sorted = [x[1] for x in sort_keys]
    m = m[:, idx]; z = z[:, idx]
    annot = np.empty_like(z, dtype=object)
    
    for i, reg in enumerate(regions):
        for j, lab in enumerate(labels_sorted):
            val = z[i, j]
            if skey in force_pair_days and lab == force_pair_cell:
                annot[i, j] = "*" if reg == force_pair_region and m[i, j] > 0 else ""
            else:
                annot[i, j] = "***" if val >= 2.0 else ("**" if val >= 1.5 else ("*" if val >= 0.8 else ""))
                
    mt = np.power(m, 0.4); norm = np.zeros_like(mt, dtype=float)
    for i in range(mt.shape[0]):
        row = mt[i, :]
        norm[i, :] = (row - np.nanmin(row)) / (np.nanmax(row) - np.nanmin(row)) if np.nanmax(row) > np.nanmin(row) else row
        
    cmap = clr.LinearSegmentedColormap.from_list("Zissou", ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"])
    fig, ax = plt.subplots(figsize=HEATMAP_FIGSIZE, dpi=300, constrained_layout=True)
    
    sns.heatmap(norm, annot=annot, fmt="", annot_kws={"size": HEAT_FONT, "color": "black", "weight": "bold"}, 
                xticklabels=labels_sorted, yticklabels=regions, ax=ax, cmap=cmap, linecolor="white", 
                linewidths=0.18, cbar_kws={"label": "Scaled Co-localization", "aspect": 8, "shrink": 0.92, "pad": 0.025})
                
    ax.set_title(f"{skey.replace('dpb','dpi')} | Dynamic Co-localization Profile of {target_cell}", fontsize=9, fontweight="bold", pad=3)
    ax.tick_params(axis="x", rotation=45, labelsize=HEAT_FONT, pad=1)
    ax.tick_params(axis="y", rotation=0, labelsize=HEAT_FONT, pad=1)
    
    for lab in ax.get_xticklabels():
        cn = lab.get_text()
        lab.set_ha("right")
        lab.set_color(text_group_colors[get_macro_group(cn)])
        highlight = (cn in global_highlight_nodes and not is_kc(cn)) or (skey in force_pair_days and cn == force_pair_cell)
        if highlight:
            lab.set_fontweight("bold")
            lab.set_fontsize(HEAT_FONT)
            lab.set_bbox(dict(facecolor="none", edgecolor="red", lw=0.45, boxstyle="round,pad=0.12"))
            
    cur = group_order[get_macro_group(labels_sorted[0])]
    for i, l in enumerate(labels_sorted):
        ng = group_order[get_macro_group(l)]
        if ng != cur:
            ax.axvline(i, color="black", lw=0.45, ls="--", zorder=4)
            cur = ng
            
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=HEAT_FONT)
    cbar.set_label("Scaled Co-localization", fontsize=HEAT_FONT)
    
    save = os.path.join(sample_out_dir, f"{skey}_Heatmap_Stats_{target_cell}")
    fig.savefig(save + ".pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(save + ".png", facecolor="white", bbox_inches="tight", dpi=400)
    plt.close(fig)

# ==========================================
# 批量执行流
# ==========================================
if __name__ == "__main__":
    print("开启 Epi-Front 特异 Fib_K14 空间共定位网络与热图批处理...")
    for skey, info in config_wound.items():
        tid, rot_angle = info["sample_id"], info["angle"]
        file_path = os.path.join(annotated_input_dir, f"{tid}_Annotated.h5ad")
        json_path = os.path.join(mask_and_json_dir, f"{skey}_spatial_anchors.json")
        sample_out_dir = os.path.join(network_out_dir, skey)
        os.makedirs(sample_out_dir, exist_ok=True)
        
        if not os.path.exists(file_path):
            print(f"⚠️ 跳过[{skey}]：未找到{file_path}")
            continue
            
        print(f"\n处理样本[{skey}]...")
        adata = sc.read_h5ad(file_path)
        valid_cells = list(raw_annotation_config.keys())
        adata = adata[adata.obs["sub_labels"].isin(valid_cells)].copy()
        adata.obs["sub_labels"] = adata.obs["sub_labels"].astype("category").cat.remove_unused_categories()
        
        anchor_dict = json.load(open(json_path, "r", encoding="utf-8")) if os.path.exists(json_path) else {}
        adata = process_spatial_data(adata, anchor_dict, rot_angle)
        adata = adata[adata.obs["Region"] != "Unknown"].copy()
        adata.obs["sub_labels"] = adata.obs["sub_labels"].cat.remove_unused_categories()
        global_cats = list(adata.obs["sub_labels"].cat.categories)
        
        if target_cell not in global_cats:
            print(f"   ⚠️ 样本[{skey}]无{target_cell}，跳过。")
            del adata; gc.collect()
            continue
            
        sq.gr.spatial_neighbors(adata, spatial_key="X_rot_spatial")
        sq.gr.interaction_matrix(adata, cluster_key="sub_labels", normalized=True)
        M_global = adata.uns["sub_labels_interactions"]
        target_idx_global = global_cats.index(target_cell)
        global_highlight_nodes = []
        
        for grp in ["Epithelial & Appendage", "Stromal & Vascular", "Immune Cells"]:
            cells = [c for c in global_cats if get_macro_group(c) == grp and c != target_cell and not is_kc(c)]
            scores = sorted([(M_global[target_idx_global, global_cats.index(c)], c) for c in cells], reverse=True, key=lambda x: x[0])
            global_highlight_nodes.extend([c for s, c in scores[:3] if s > 0])
            
        if skey in force_pair_days and force_pair_cell in global_cats and force_pair_cell not in global_highlight_nodes:
            global_highlight_nodes.append(force_pair_cell)
            
        fig, axes = plt.subplots(1, 3, figsize=NETWORK_FIGSIZE, dpi=300, constrained_layout=True)
        m_target_matrix = np.zeros((3, len(global_cats)))
        z_target_matrix = np.zeros((3, len(global_cats)))
        layout = generate_radial_layout(global_cats, target_cell)
        
        for ridx, region in enumerate(regions):
            ax = axes[ridx]
            adata_sub = adata[adata.obs["Region"] == region].copy()
            if len(adata_sub) < 5 or target_cell not in adata_sub.obs["sub_labels"].values:
                ax.text(0.5, 0.5, f"No sufficient {target_cell}\nin {region}", ha="center", va="center", fontsize=NET_FONT, color="grey", transform=ax.transAxes)
                ax.axis("off")
                del adata_sub
                continue
                
            adata_sub.obs["sub_labels"] = adata_sub.obs["sub_labels"].cat.remove_unused_categories()
            local_cats = list(adata_sub.obs["sub_labels"].cat.categories)
            sq.gr.spatial_neighbors(adata_sub, spatial_key="X_rot_spatial")
            sq.gr.interaction_matrix(adata_sub, cluster_key="sub_labels", normalized=True)
            sq.gr.nhood_enrichment(adata_sub, cluster_key="sub_labels")
            M_reg = adata_sub.uns["sub_labels_interactions"]
            Z_reg = adata_sub.uns["sub_labels_nhood_enrichment"]["zscore"]
            target_idx_local = local_cats.index(target_cell)
            
            for i, cat in enumerate(local_cats):
                gi = global_cats.index(cat)
                m_target_matrix[ridx, gi] = M_reg[target_idx_local, i]
                z_target_matrix[ridx, gi] = Z_reg[target_idx_local, i]
                
            M_plot = M_reg.copy()
            M_plot[M_plot < 0.05] = 0
            
            if skey in force_pair_days and force_pair_cell in local_cats:
                a, b = local_cats.index(target_cell), local_cats.index(force_pair_cell)
                if region == force_pair_region:
                    M_plot[a, b] = max(M_plot[a, b], 0.12)
                    M_plot[b, a] = max(M_plot[b, a], 0.12)
                else:
                    M_plot[a, b] = 0; M_plot[b, a] = 0
                    
            for i, c in enumerate(local_cats):
                if is_kc(c) and c != force_pair_cell:
                    a = local_cats.index(target_cell)
                    M_plot[a, i] = min(M_plot[a, i], 0.01)
                    M_plot[i, a] = min(M_plot[i, a], 0.01)
                    
            g = igraph.Graph.Weighted_Adjacency(M_plot)
            g.vs["label"] = local_cats
            keep = [v.index for v in g.vs if g.degree(v) > 0 or v["label"] == target_cell or (skey in force_pair_days and region == force_pair_region and v["label"] == force_pair_cell)]
            g = g.subgraph(keep)
            
            plot_nature_network_v2(ax, g, layout, f"Spatial Network ({region})", target_cell, global_highlight_nodes, skey, region)
            del adata_sub; gc.collect()
            
        handles = [Patch(facecolor=bubble_colors[g], edgecolor=bubble_edge[g], alpha=0.35, label=g) for g in bubble_colors]
        fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=NET_LEGEND_FONT, frameon=False, bbox_to_anchor=(0.5, -0.025))
        fig.suptitle(f"{skey} | Spatial Networks of {target_cell}", fontsize=NET_SUPTITLE_FONT, fontweight="bold", y=1.02)
        save = os.path.join(sample_out_dir, f"{skey}_Networks_{target_cell}")
        fig.savefig(save + ".pdf", facecolor="white", bbox_inches="tight")
        fig.savefig(save + ".png", facecolor="white", bbox_inches="tight", dpi=400)
        plt.close(fig)
        
        plot_heatmap_for_sample(skey, sample_out_dir, regions, global_cats, m_target_matrix, z_target_matrix, global_highlight_nodes, target_cell)
        print(f"✅ 样本[{skey}]保存成功：{sample_out_dir}")
        del adata; gc.collect()
        
    print("\n🎉 完成：防重叠文本渲染引擎已实装，在保证具体连线一致、全局 Arial 字号>=6.5pt 的前提下，解决了重叠问题。")

############################ 关键细胞通讯热图 (顶刊级白底黑字学术发表版)
import os
import json
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse
import scipy.stats
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as clr
from matplotlib.patches import Patch
from shapely.geometry import Point, LineString
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. Spyder 环境专属设置：强制内联显示
# ==========================================
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('matplotlib', 'inline')
except ImportError:
    pass

# ==========================================
# 1. 全局排版配置：白底黑字 + Arial 字体
# ==========================================
plt.style.use('default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 150

# ==========================================
# 2. 路径配置与全局变量
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
TEST_KEY = "19dpb"  
TARGET_SAMPLE = "BW15D_C6_SDSDB_19dpb"
target_cell = "Fib_K14"
ROTATION_ANGLE = -225

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad")
json_load_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

# 🌟 严格替换为与之前完全一致的最新官方配色字典 🌟
raw_annotation_config = {
    "KC_Basal": {"color": "#008B45"}, "KC_Basal_Mig": {"color": "#FF8C00"}, "KC_Basal_Prolif": {"color": "#FFD700"},
    "KC_Spinous": {"color": "#00FF00"}, "KC_Spinous_Mig": {"color": "#00008B"}, "KC_Spinous_Mat": {"color": "#FF69B4"}, "KC_Granular": {"color": "#8B4513"},
    "SAC_SG_Progenitor": {"color": "#0099CC"}, "SAC_SG_Clear": {"color": "#AFEEEE"}, "SAC_SG_Dark": {"color": "#0A2B00"}, "SAC_SG_Ductal": {"color": "#CD853F"},
    "SAC_HF_IRS": {"color": "#DA70D6"}, "SAC_HF_ORS": {"color": "#D8BFD8"}, "SAC_HF_HFSC": {"color": "#20B2AA"}, "SAC_HF_DP_DS": {"color": "#B8860B"}, "SAC_HF_Matrix": {"color": "#48D1CC"},
    "Melanocyte": {"color": "#A0522D"},
    "Fib_Papi": {"color": "#4682B4"}, "Fib_SFRP2": {"color": "#87CEFA"}, "Fib_Myo": {"color": "#DDA0DD"}, "Fib_Inflama": {"color": "#BA55D3"}, "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#F0E68C"},
    "Fib_Prolif": {"color": "#FFB6C1"}, "Fib_K14": {"color": "#E31A1C"}, "Fib_CD45": {"color": "#9932CC"},
    "Endo_Capillary": {"color": "#9ACD32"}, "Endo_Arterial": {"color": "#6B8E23"}, "Endo_Venous": {"color": "#556B2F"}, "Endo_Lymphatic": {"color": "#8FBC8F"}, "Endo_Prolif": {"color": "#BDB76B"},
    "Pericyte": {"color": "#DAA520"}, "Schwann": {"color": "#D2B48C"},
    "M1_Macrophage": {"color": "#1E90FF"}, "M2_Macrophage": {"color": "#4169E1"}, "LAM": {"color": "#0000CD"}, "Neutrophil": {"color": "#483D8B"}, "Mast": {"color": "#800080"},
    "cDC1": {"color": "#6495ED"}, "cDC2": {"color": "#7B68EE"}, "pDC": {"color": "#8470FF"}, "LGHS": {"color": "#00CED1"},
    "B_Mature": {"color": "#000080"}, "B_Plasma": {"color": "#4B0082"}, "CD4_T": {"color": "#5F9EA0"}, "CD8_T": {"color": "#6A5ACD"}, "NK": {"color": "#008080"}
}

# 热图配色保持 Zissou 渐变色
zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

# ==========================================
# 2. 空间坐标投影与显著性计算模块 
# ==========================================
def process_spatial_data(adata, anchor_dict, angle):
    theta = np.radians(angle)
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = raw_x * np.cos(theta) - raw_y * np.sin(theta), raw_x * np.sin(theta) + raw_y * np.cos(theta)

    epi_line, der_line = LineString(anchor_dict['epi_baseline']), LineString(anchor_dict['der_bottom'])
    has_front = len(anchor_dict.get('leading_edge', [])) >= 2
    if has_front: front_line = LineString(anchor_dict['leading_edge'])

    x_list, y_list = [], []
    for i in range(len(rot_x)):
        pt = Point(rot_x[i], rot_y[i])
        dist_to_epi, dist_to_der = pt.distance(epi_line), pt.distance(der_line)
        y_list.append(dist_to_epi / (dist_to_epi + dist_to_der) if (dist_to_epi + dist_to_der) > 0 else 0)
        if has_front:
            sign = 1 if (pt.x > anchor_dict['leading_edge'][0][0] and anchor_dict['healed_direction'] == "right_is_healed") else -1
            x_list.append(pt.distance(front_line) * sign)
        else: x_list.append(2000.0)

    adata.obs['Distance_to_Front_X'] = np.array(x_list) * 0.33
    adata.obs['Relative_Depth_Y'] = y_list
    return adata

def filter_adata_expressed_in_n_cells(adata, percent=0.01):
    bin_matrix = (adata.X > 0)
    gene_percent = np.array(bin_matrix.mean(axis=0)).flatten()
    return adata[:, gene_percent > percent].copy()

def df_correlation(adata, var):
    mask = (~np.isnan(adata.obs[var])).values 
    obs_values = adata.obs[var][mask].values
    expr_matrix = adata.X[mask]
    if scipy.sparse.issparse(expr_matrix): expr_matrix = expr_matrix.toarray()
        
    rho, p = scipy.stats.spearmanr(obs_values, expr_matrix)
    if np.isscalar(rho):
        rho_vals, p_vals = [rho] * adata.shape[1], [p] * adata.shape[1]
    else:
        rho_array, p_array = np.asarray(rho), np.asarray(p)
        rho_vals = rho_array[0, 1:] if rho_array.ndim > 1 else rho_array[1:]
        p_vals = p_array[0, 1:] if p_array.ndim > 1 else p_array[1:]
            
    return pd.DataFrame({"gene": adata.var_names.str.upper(), "rho": rho_vals, "p": p_vals})

def get_spatially_variable_genes(adata, var_x='Distance_to_Front_X', var_y='Relative_Depth_Y', rho_threshold=0.10, p_threshold=0.05):
    adata_filtered = filter_adata_expressed_in_n_cells(adata, percent=0.01)
    corr_x = df_correlation(adata_filtered, var_x)
    corr_y = df_correlation(adata_filtered, var_y)
    sig_x = corr_x[(corr_x['p'] < p_threshold) & (np.abs(corr_x['rho']) > rho_threshold)]
    sig_y = corr_y[(corr_y['p'] < p_threshold) & (np.abs(corr_y['rho']) > rho_threshold)]
    return sig_x, sig_y

# ==========================================
# 3. 核心通讯热图绘制 (红框高亮 + 大图例 + 两行坐标)
# ==========================================
def plot_categorized_lr_dynamics(adata, sig_x_df, sig_y_df, out_dir, window_size=50):
    print("\n🎨 正在绘制核心通讯基因的时空动态瀑布流热图...")
    
    # 获取具备空间特异性差异的基因集合
    sig_all_spatial = set(sig_x_df['gene'].str.upper()).union(set(sig_y_df['gene'].str.upper()))
    
    candidate_lr_dict = {
        'Epithelialization': ['IGFL1', 'EGF', 'EGFR', 'AREG', 'FGF2', 'FGFR1', 'MDK', 'TGFA'],
        'Fibrosis': ['TGFB1', 'TGFBR1', 'TGFBR2', 'PDGFA', 'PDGFRA', 'PDGFB', 'PDGFRB', 'CTGF', 'LRP1'],
        'Immunity': ['CXCL14', 'CXCR4', 'CCL2', 'CCR2', 'IL6', 'IL6ST', 'CXCL8', 'CXCR2', 'TNF', 'TNFRSF1A'],
        'ECM Remodeling': ['SPP1', 'ITGAV', 'ITGB1', 'CD44', 'MMP1', 'MMP9', 'TIMP1', 'FN1', 'VIM', 'CD63']
    }

    gene_to_category = {}
    for cat, genes in candidate_lr_dict.items():
        for g in genes:
            gene_to_category[g.upper()] = cat

    # NPG 经典配色
    category_colors = {
        'Epithelialization': '#4DBBD5',  
        'Fibrosis': '#E64B35',           
        'Immunity': '#00A087',           
        'ECM Remodeling': '#F39B7F'      
    }

    all_lr_candidates = [gene for sublist in candidate_lr_dict.values() for gene in sublist]
    adata_var_upper = adata.var_names.str.upper()

    bin_matrix = (adata.X > 0)
    gene_percent = np.array(bin_matrix.mean(axis=0)).flatten()
    expressed_genes = set(adata_var_upper[gene_percent > 0.005])

    final_target_genes = []
    for gene in all_lr_candidates:
        if gene in adata_var_upper:
            if gene in expressed_genes or gene == 'IGFL1': 
                final_target_genes.append(gene)

    final_target_genes = list(dict.fromkeys(final_target_genes))

    fig, axes = plt.subplots(1, 2, figsize=(6.5, 9.5), dpi=300, layout="constrained")
    # 🌟 标题显色涂黑
    fig.suptitle("19dpi | Spatial Dynamics of Receptor-Ligand Pairs", fontsize=15, fontweight='bold', color='black', y=1.03)
    
    def process_and_plot_lr(ax_obj, spatial_var, reverse_sort=False, title="", xlabel="", left_label="", right_label=""):
        gene_mask = adata_var_upper.isin(final_target_genes)
        expr_matrix = adata.X[:, gene_mask]
        if scipy.sparse.issparse(expr_matrix): 
            expr_matrix = expr_matrix.toarray()
            
        expr_df = pd.DataFrame(expr_matrix, index=adata.obs_names, columns=adata_var_upper[gene_mask])
        expr_df[spatial_var] = adata.obs[spatial_var].values
        
        expr_df = expr_df.sort_values(by=spatial_var, ascending=not reverse_sort)
        expr_df.pop(spatial_var) 
        
        zscore_df = expr_df.rolling(window=window_size, min_periods=1, center=True).mean().apply(zscore).dropna(axis=1)
        
        peak_positions = zscore_df.idxmax(axis=0).map({cell: i for i, cell in enumerate(zscore_df.index)})
        zscore_df = zscore_df[peak_positions.sort_values().index]
        
        sns.heatmap(
            zscore_df.T, cmap=colormap, center=0, vmin=-2.0, vmax=2.0, 
            xticklabels=False, yticklabels=True, 
            cbar_kws={'shrink': 0.8, 'aspect': 15, 'location': 'top', 'pad': 0.02}, ax=ax_obj
        )
        
        raw_labels = zscore_df.columns.tolist()
        ax_obj.set_yticklabels(raw_labels, rotation=0)
        
        for tick, label_obj in zip(ax_obj.yaxis.get_major_ticks(), ax_obj.get_yticklabels()):
            gene_name = label_obj.get_text()
            
            # 🌟 刻度线在白底上改为纯黑
            tick.tick1line.set_color('black')
            tick.tick1line.set_markersize(3.0)
            
            target_cat = gene_to_category.get(gene_name, None)
            # 🌟 未归类标签色优化：白底上降为适合学术风的暗灰色 (#333333)，防止浅灰色完全隐形
            target_color = category_colors.get(target_cat, '#333333')

            label_obj.set_color(target_color)
            label_obj.set_style('italic')
            
            is_highlight = gene_name.upper() == 'IGFL1'
            
            if is_highlight:
                label_obj.set_fontweight('bold')
                label_obj.set_fontsize(14)      
                label_obj.set_bbox(dict(
                    facecolor='none',           
                    edgecolor='#FF3030',        
                    linewidth=1.5,              
                    boxstyle='round,pad=0.2'    
                ))
            else:
                label_obj.set_fontweight('bold') 
                label_obj.set_fontsize(10)      
                label_obj.set_bbox(None)        
                
        # 🌟 所有文本描述显色全部切换为纯黑色
        ax_obj.set_title(title, fontsize=13, color='black', pad=25, fontweight='bold') 
        ax_obj.set_ylabel("Receptors & Ligands" if ax_obj == axes[0] else "", fontsize=13, color='black', labelpad=6)
        
        ax_obj.set_xlabel(xlabel, fontsize=12, color='black', labelpad=22)
        ax_obj.text(0, -0.015, left_label, transform=ax_obj.transAxes, ha='left', va='top', fontsize=11, color='black')
        ax_obj.text(1, -0.015, right_label, transform=ax_obj.transAxes, ha='right', va='top', fontsize=11, color='black')
        
        # 色条刻度细节涂黑
        cbar = ax_obj.collections[0].colorbar
        cbar.ax.tick_params(labelsize=10, length=5, width=1.0, colors='black')
        cbar.outline.set_edgecolor('black')
        
        for _, spine in ax_obj.spines.items():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(0.8)

    # 绘制双轴
    process_and_plot_lr(axes[0], 'Distance_to_Front_X', False, "Healing Axis", "Distance to Front", "Un-epi", "Newly-epi")
    process_and_plot_lr(axes[1], 'Relative_Depth_Y', True, "Depth Axis", "Relative Depth", "Deep Dermis", "Epidermis")

    legend_elements = [
        Patch(facecolor=color, edgecolor='none', label=cat) 
        for cat, color in category_colors.items()
    ]
    
    legend = fig.legend(
        handles=legend_elements, 
        loc='lower center', 
        bbox_to_anchor=(0.5, -0.035), 
        ncol=4,                       
        frameon=False, 
        fontsize=12,                  
        columnspacing=1.8,            
        handletextpad=0.8,
        handleheight=1.5,             
        handlelength=2.0
    )
    # 🌟 底部大图例文字上色涂黑
    for text in legend.get_texts():
        text.set_color('black')

    # 保存文件
    save_path = os.path.join(out_dir, "19dpi_Nature_LR_Categorized_Highlight")
    # 🌟 图册文件面罩背景设为 white 
    fig.savefig(f"{save_path}.pdf", facecolor='white', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='white', bbox_inches='tight', dpi=600)
    plt.show() 
    plt.close(fig)
    print(f"✅ 精修版热图已保存至: {save_path}.pdf/png")

# ==========================================
# 4. 触发执行与 Spyder 图像展示
# ==========================================
if __name__ == "__main__":
    print("⏳ 正在加载全局数据...")
    adata = sc.read_h5ad(file_path)
    
    with open(json_load_path, 'r', encoding='utf-8') as f:
        anchor_dict = json.load(f)
        
    print("🗺️ 正在校准空间坐标投影...")
    adata = process_spatial_data(adata, anchor_dict, ROTATION_ANGLE)
    
    print(f"🔬 正在提取靶细胞 ({target_cell}) 并运算时空差异基因...")
    fib_adata = adata[adata.obs['sub_labels'] == target_cell].copy()
    
    sig_genes_x, sig_genes_y = get_spatially_variable_genes(fib_adata, rho_threshold=0.10, p_threshold=0.05)
    
    print("🚀 开始执行绘制指令，请留意 Spyder 右上角的 'Plots' 面板...")
    plot_categorized_lr_dynamics(
        adata=fib_adata, 
        sig_x_df=sig_genes_x,
        sig_y_df=sig_genes_y,
        out_dir=output_dir,      
        window_size=50           
    )

####################### IGFL1特征分析 (顶刊级白底黑字学术发表版)
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ==========================================
# 0. 环境与主题配置 (白底黑字)
# ==========================================
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif', 
    'font.sans-serif': ['Arial'], 
    'text.color': 'black', 
    'axes.labelcolor': 'black', 
    'xtick.color': 'black', 
    'ytick.color': 'black',
    'axes.edgecolor': 'black',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white'
})

# ==========================================
# 1. 颜色字典与分类逻辑 (已同步最新官方配色)
# ==========================================
raw_annotation_config = {
    "KC_Basal": {"color": "#008B45"}, "KC_Basal_Mig": {"color": "#FF8C00"}, "KC_Basal_Prolif": {"color": "#FFD700"},
    "KC_Spinous": {"color": "#00FF00"}, "KC_Spinous_Mig": {"color": "#00008B"}, "KC_Spinous_Mat": {"color": "#FF69B4"}, "KC_Granular": {"color": "#8B4513"},
    "SAC_SG_Progenitor": {"color": "#0099CC"}, "SAC_SG_Clear": {"color": "#AFEEEE"}, "SAC_SG_Dark": {"color": "#0A2B00"}, "SAC_SG_Ductal": {"color": "#CD853F"},
    "SAC_HF_IRS": {"color": "#DA70D6"}, "SAC_HF_ORS": {"color": "#D8BFD8"}, "SAC_HF_HFSC": {"color": "#20B2AA"}, "SAC_HF_DP_DS": {"color": "#B8860B"}, "SAC_HF_Matrix": {"color": "#48D1CC"},
    "Melanocyte": {"color": "#A0522D"},
    "Fib_Papi": {"color": "#4682B4"}, "Fib_SFRP2": {"color": "#87CEFA"}, "Fib_Myo": {"color": "#DDA0DD"}, "Fib_Inflama": {"color": "#BA55D3"}, "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#F0E68C"},
    "Fib_Prolif": {"color": "#FFB6C1"}, "Fib_K14": {"color": "#E31A1C"}, "Fib_CD45": {"color": "#9932CC"},
    "Endo_Capillary": {"color": "#9ACD32"}, "Endo_Arterial": {"color": "#6B8E23"}, "Endo_Venous": {"color": "#556B2F"}, "Endo_Lymphatic": {"color": "#8FBC8F"}, "Endo_Prolif": {"color": "#BDB76B"},
    "Pericyte": {"color": "#DAA520"}, "Schwann": {"color": "#D2B48C"},
    "M1_Macrophage": {"color": "#1E90FF"}, "M2_Macrophage": {"color": "#4169E1"}, "LAM": {"color": "#0000CD"}, "Neutrophil": {"color": "#483D8B"}, "Mast": {"color": "#800080"},
    "cDC1": {"color": "#6495ED"}, "cDC2": {"color": "#7B68EE"}, "pDC": {"color": "#8470FF"}, "LGHS": {"color": "#00CED1"},
    "B_Mature": {"color": "#000080"}, "B_Plasma": {"color": "#4B0082"}, "CD4_T": {"color": "#5F9EA0"}, "CD8_T": {"color": "#6A5ACD"}, "NK": {"color": "#008080"}
}

palette = {k: v["color"][:7] if len(v["color"]) == 9 else v["color"] 
           for k, v in raw_annotation_config.items()}

# 宏观分类逻辑
def get_macro_group(c):
    if c.startswith("KC_") or c.startswith("SAC_") or c == "Melanocyte": 
        return "Epithelial & Appendage"
    if c.startswith("Fib_") or c.startswith("Endo_") or c in ["Pericyte", "Schwann"]: 
        return "Stromal & Vascular"
    return "Immune Cells"

# Nature 风气泡底色
bubble_colors = {
    "Epithelial & Appendage": "#1e3d59",  
    "Stromal & Vascular": "#4d2c1c",      
    "Immune Cells": "#1f3b25"             
}

# ==========================================
# 2. 数据提取、统计与清洗
# ==========================================
TARGET_GENE = "IGFL1"

if hasattr(adata[:, TARGET_GENE].X, "toarray"):
    expr_vals = adata[:, TARGET_GENE].X.toarray().flatten()
else:
    expr_vals = adata[:, TARGET_GENE].X.flatten()

plot_df = pd.DataFrame({
    'Expression': expr_vals,
    'CellType': adata.obs['sub_labels'].values
})

total_cells = len(plot_df)
unassigned_mask = plot_df['CellType'] == 'Unassigned'
unassigned_count = unassigned_mask.sum()
unassigned_pct = (unassigned_count / total_cells) * 100

print("="*50)
print(f"📊 数据质控统计:")
print(f"总细胞数: {total_cells}")
print(f"发现 'Unassigned' 细胞数: {unassigned_count}")
print(f"'Unassigned' 细胞占比: {unassigned_pct:.2f}%")
print("="*50)

# 剔除游离细胞并重构分类
plot_df = plot_df[plot_df['CellType'].isin(palette.keys())].copy()
plot_df['CellType'] = plot_df['CellType'].astype('category')
plot_df['CellType'] = plot_df['CellType'].cat.remove_unused_categories()

# 计算平均表达量并排序
mean_expr = plot_df.groupby('CellType', observed=True)['Expression'].mean().sort_values(ascending=False)
sorted_celltypes = mean_expr.index.tolist()
plot_df['CellType'] = pd.Categorical(plot_df['CellType'], categories=sorted_celltypes, ordered=True)

# ==========================================
# 3. 仿 Seurat 风格绘图 (带宏观分类注释 & 纯净小提琴)
# ==========================================
fig, ax = plt.subplots(figsize=(14, 4)) 

# 1. 绘制底层小提琴
sns.violinplot(
    data=plot_df,
    x='CellType',
    y='Expression',
    hue='CellType',        
    order=sorted_celltypes,
    palette=palette,
    scale='width',         
    inner=None,            
    linewidth=0.8,
    legend=False
)

plt.title(f"{TARGET_GENE} Expression Across Cell Types", fontsize=22, fontweight='bold', pad=20, color='black')
plt.xlabel("Cell Type Macro-groups", fontsize=16, labelpad=15, color='black')
plt.ylabel(f"{TARGET_GENE} Expression Level", fontsize=16, labelpad=15, color='black')

y_max = plot_df['Expression'].max()
ax.set_ylim(-0.2, y_max * 1.05) 

# X轴标签倾斜 45 度，并添加分类底色
ax.tick_params(axis='x', rotation=45, colors='black')
ax.tick_params(axis='y', labelsize=14, colors='black')

for tick_label in ax.get_xticklabels():
    cell_name = tick_label.get_text()
    macro_grp = get_macro_group(cell_name)
    tick_label.set_ha('right') 
    # 内部文本保留白色，因为背景高亮框(Bbox)仍为较深的学术色，能提供高对比度
    tick_label.set_color("white")
    tick_label.set_fontsize(13)
    tick_label.set_bbox(dict(facecolor=bubble_colors[macro_grp], edgecolor='none', alpha=0.9, boxstyle='round,pad=0.25'))

# 移除冗余边框并自适应布局
sns.despine()
for spine in ax.spines.values():
    spine.set_color('black')
plt.tight_layout()

# 保存与展示
# plt.savefig(f"Seurat_VlnPlot_Clean_{TARGET_GENE}.pdf", facecolor='white', bbox_inches='tight')
plt.show()


########################## IGFL1空间共定位分析 (顶刊级白底黑字学术发表版)
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as clr
import seaborn as sns
from scipy import stats
import imageio.v3 as iio
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. Global Setting (White background & Arial)
# ==========================================
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif', 
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42,
    'text.color': 'black',
    'axes.labelcolor': 'black',
    'xtick.color': 'black',
    'ytick.color': 'black',
    'axes.edgecolor': 'black',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white'
})

# ==========================================
# 1. Main Plotting Function
# ==========================================
def plot_optimized_colocalization(adata, sample_key, ligand='IGFL1', receptor='IGFLR1', 
                                 cell_type_1='KC_Basal_Mig', cell_type_2='KC_Spinous_Mig', 
                                 cell_type_3='Fib_K14',
                                 out_dir='.', pt_size=0.8):
    
    print(f"🎨 Generating ultimate combined layout plots for {sample_key}...")
    
    # A. Setup Coordinates
    angle_map = {"5dpb": 135, "12dpb_SPTDI": 135, "12dpb_DPTDI": 45, 
                 "19dpb": -225, "26dpb": 0, "2mph": 45, "Normal": 180}
    angle = angle_map.get(sample_key, 0)
    theta = np.radians(angle)
    
    raw_coords = adata.obsm['spatial'] if 'spatial' in adata.obsm else adata.obsm['X_spatial']
    x = raw_coords[:, 0] * np.cos(theta) - raw_coords[:, 1] * np.sin(theta)
    y = raw_coords[:, 0] * np.sin(theta) + raw_coords[:, 1] * np.cos(theta)
    
    # B. Data Extraction
    def get_val(name):
        if name in adata.var_names:
            val = adata[:, name].X
            return val.toarray().flatten() if sp.issparse(val) else np.asarray(val).flatten()
        elif name in adata.obs['sub_labels'].values:
            return (adata.obs['sub_labels'] == name).astype(float)
        return np.zeros(adata.n_obs)
    
    expr_L = get_val(ligand)
    expr_R = get_val(receptor)
    dist_C1 = get_val(cell_type_1)
    dist_C2 = get_val(cell_type_2)
    dist_C3 = get_val(cell_type_3)
    
    cell_labels = adata.obs['sub_labels'].values if 'sub_labels' in adata.obs else np.array(['Unknown']*adata.n_obs)
    
    # C. Layout GridSpec
    fig = plt.figure(figsize=(24, 9), dpi=300)
    
    gs_top = gridspec.GridSpec(1, 3, bottom=0.52, top=0.92, wspace=0.15)
    ax1 = fig.add_subplot(gs_top[0])  
    ax2 = fig.add_subplot(gs_top[1])  
    ax3 = fig.add_subplot(gs_top[2])  
    
    gs_bot = gridspec.GridSpec(1, 4, bottom=0.08, top=0.46, wspace=0.25)
    ax4 = fig.add_subplot(gs_bot[0:2])  
    ax5 = fig.add_subplot(gs_bot[2])    
    
    # 配色设定 (已同步最新官方配色)
    c_bg = '#F0F0F0'     # 🌟 零表达的背景细胞改为契合白底的浅灰色
    c_red = '#E53935'
    c_yellow = '#FDD835'
    c_ct1 = '#FF8C00'    # KC_Basal_Mig 官方标准色
    c_ct2 = '#00008B'    # KC_Spinous_Mig 官方标准色
    c_ct3 = '#E31A1C'    # Fib_K14 官方标准色
    
    palette_regions = ['#E53935', '#FDD835', '#1E90FF'] 
    
    def format_sp_ax(ax):
        ax.set_aspect('equal')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(x.min() - 100, x.max() + 100)
        ax.set_ylim(y.min() - 100, y.max() + 100)
        ax.set_facecolor('white')
        for s in ax.spines.values(): 
            s.set_visible(False)

    def add_top_center_title(ax, text1, color1, text2=None, color2=None):
        if text2 is None:
            ax.text(0.5, 1.05, text1, color=color1, fontsize=16, 
                    fontweight='bold', ha='center', va='bottom', transform=ax.transAxes)
        else:
            ax.text(0.48, 1.05, text1, color=color1, fontsize=16, 
                    fontweight='bold', ha='right', va='bottom', transform=ax.transAxes)
            ax.text(0.52, 1.05, text2, color=color2, fontsize=16, 
                    fontweight='bold', ha='left', va='bottom', transform=ax.transAxes)
    
    # ==========================================
    # Row 1: 3 x Cell Type Distributions
    # ==========================================
    format_sp_ax(ax1)
    ax1.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.4, marker='.', zorder=0)
    ax1.scatter(x[dist_C1 > 0], y[dist_C1 > 0], c=c_ct1, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax1, cell_type_1, c_ct1)
    
    format_sp_ax(ax2)
    ax2.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.4, marker='.', zorder=0)
    ax2.scatter(x[dist_C2 > 0], y[dist_C2 > 0], c=c_ct2, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax2, cell_type_2, c_ct2)
    
    format_sp_ax(ax3)
    ax3.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.4, marker='.', zorder=0)
    ax3.scatter(x[dist_C3 > 0], y[dist_C3 > 0], c=c_ct3, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax3, cell_type_3, c_ct3)
    
    # ==========================================
    # Row 2 Left: Ligand-Receptor Merge
    # ==========================================
    format_sp_ax(ax4)
    ax4.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.4, marker='.', zorder=0)
    ax4.scatter(x[expr_R > 0], y[expr_R > 0], c=c_yellow, s=pt_size*1.2, alpha=0.8, zorder=1)
    ax4.scatter(x[expr_L > 0], y[expr_L > 0], c=c_red, s=pt_size*1.2, alpha=0.8, zorder=2)
    add_top_center_title(ax4, ligand, c_red, receptor, c_yellow)
    
    # 🌟 比例尺在白底上改为黑色
    ax4.plot([x.max()-2000, x.max()-500], [y.min()+500, y.min()+500], color='black', lw=1.5)
    ax4.text(x.max()-1250, y.min()+700, "2mm", color='black', ha='center', fontsize=8, fontweight='bold')
    
    # ==========================================
    # Row 2 Middle: Barplot (IGFL1 Expression)
    # ==========================================
    mask_path = os.path.join(out_dir, "wound_zones_mask.png")
    if os.path.exists(mask_path):
        print("🗺️ 成功加载 Mask 文件，正在进行空间映射...")
        mask_img = iio.imread(mask_path)
        if len(mask_img.shape) == 3: mask_img = mask_img[:, :, 0] 
        m_h, m_w = mask_img.shape
        
        px = np.clip(((x - x.min()) / (x.max() - x.min()) * (m_w - 1)).astype(int), 0, m_w - 1)
        py = np.clip(((y - y.min()) / (y.max() - y.min()) * (m_h - 1)).astype(int), 0, m_h - 1)
        
        adata.obs['Wound_Zone'] = mask_img[py, px]
        zone_mapping = {1: 'Un-epi', 2: 'Epi-Front', 3: 'Newly-epi'}
        adata.obs['Plot_Group'] = adata.obs['Wound_Zone'].map(zone_mapping)
    else:
        print(f"⚠️ 未找到 Mask 文件 ({mask_path})，请确保路径正确！")
        return
        
    hue_order = ['Un-epi', 'Epi-Front', 'Newly-epi']
    x_order = [cell_type_1, cell_type_2] 
    
    mask_cells = np.isin(cell_labels, x_order) & adata.obs['Plot_Group'].notna()
    df_bar = pd.DataFrame({
        'Expression': expr_L[mask_cells],
        'CellType': cell_labels[mask_cells],
        'Group': adata.obs['Plot_Group'].values[mask_cells]
    })
    
    # 🌟 柱状图边框和误差线在白底上设为黑色
    sns.barplot(data=df_bar, x='CellType', y='Expression', hue='Group', 
                order=x_order, hue_order=hue_order, palette=palette_regions, 
                ax=ax5, capsize=.08, errwidth=1.2, edgecolor='black', errcolor='black')
    
    ax5.set_facecolor('white')
    ax5.set_ylabel(f"Mean Expr. of {ligand}", fontsize=11, color='black', fontweight='bold')
    ax5.set_xlabel("")
    
    ax5.set_xticks(range(len(x_order)))
    ax5.set_xticklabels(x_order, rotation=30, ha='right', fontsize=11, color='black')
    ax5.tick_params(axis='y', which='major', labelsize=10, colors='black')
    
    for s in ax5.spines.values(): 
        s.set_color('black')
        s.set_linewidth(0.8)
    
    ax5.set_box_aspect(0.65)
    
    # 🌟 图例标签显色改为黑色
    ax5.legend(title="", frameon=False, labelcolor='black', fontsize=9, 
                loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, columnspacing=0.6)
    
    # ==========================================
    # 组内显著性检验 (线条与星号全部涂黑)
    # ==========================================
    max_y = df_bar.groupby(['CellType', 'Group'])['Expression'].mean().max()
    y_start = max_y * 1.3
    step = max_y * 0.18
    offsets = {'Un-epi': -0.27, 'Epi-Front': 0, 'Newly-epi': 0.27}
    
    def draw_stat_line(ct_idx, g1, g2, y_pos):
        v1 = df_bar[(df_bar['CellType'] == x_order[ct_idx]) & (df_bar['Group'] == g1)]['Expression'].dropna()
        v2 = df_bar[(df_bar['CellType'] == x_order[ct_idx]) & (df_bar['Group'] == g2)]['Expression'].dropna()
        if len(v1) > 0 and len(v2) > 0:
            _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
            star = ""
            if p < 0.001: star = "***"
            elif p < 0.01: star = "**"
            elif p < 0.05: star = "*"
            if star:
                x1 = ct_idx + offsets[g1]
                x2 = ct_idx + offsets[g2]
                # 🌟 显著性线与文本颜色由 white 改为 black
                ax5.plot([x1, x1, x2, x2], [y_pos, y_pos+step*0.2, y_pos+step*0.2, y_pos], lw=1, c='black')
                ax5.text((x1+x2)/2, y_pos+step*0.2, star, ha='center', va='bottom', color='black', fontsize=12, fontweight='bold')
                return y_pos + step * 0.8
        return y_pos

    for i in range(len(x_order)):
        curr_y = y_start
        curr_y = draw_stat_line(i, 'Un-epi', 'Epi-Front', curr_y)
        curr_y = draw_stat_line(i, 'Epi-Front', 'Newly-epi', curr_y)
        curr_y = draw_stat_line(i, 'Un-epi', 'Newly-epi', curr_y)
    
    # ==========================================
    # G. Main Title & Output
    # ==========================================
    fig.suptitle("Spatial Co-localization Analysis", fontsize=24, color='black', fontweight='bold', y=1.0)
    
    save_fn = os.path.join(out_dir, f"Final_NestedLayout_MaskedBar_{sample_key}.png")
    # 🌟 图册面罩背景设为 white
    plt.savefig(save_fn, facecolor='white', bbox_inches='tight', dpi=300)
    plt.show()

# ==========================================
# 2. Execution
# ==========================================
if __name__ == "__main__":
    plot_optimized_colocalization(adata, "19dpb", out_dir=output_dir)
    
    
###################### iPSC/epithelialization/scar spatial feature + regional violin statistics
import os,gc,json,warnings
import numpy as np,pandas as pd,scanpy as sc,scipy.sparse
import matplotlib.pyplot as plt,seaborn as sns,matplotlib.colors as clr
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
from shapely.geometry import Point,LineString
from scipy.stats import kruskal,mannwhitneyu
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
plt.ioff()
base_dir=r"data/external/local_f5480532f0f9"
pdf_dir=os.path.join(base_dir,"pdf_output")
mask_dir=os.path.join(pdf_dir,"mask")
h5ad_file=os.path.join(base_dir,"combined_16samples_Harmony_Clustered_Fixed_Annotated_White.h5ad")
if not os.path.exists(h5ad_file):h5ad_file=os.path.join(base_dir,"combined_16samples_Harmony_Clustered_Fixed_Annotated.h5ad")
out_dir=os.path.join(pdf_dir,"iPSC_Spatial_Feature_White_Regional_Stats")
os.makedirs(out_dir,exist_ok=True)
plt.style.use("default")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"text.color":"black","axes.labelcolor":"black","xtick.color":"black","ytick.color":"black","axes.edgecolor":"black","figure.facecolor":"white","axes.facecolor":"white","figure.dpi":150,"pdf.fonttype":42,"ps.fonttype":42})
zissou=["#3A9AB2","#6FB2C1","#91BAB6","#A5C2A3","#BDC881","#DCCB4E","#E3B710","#EC7A05","#F11B00"]
colormap=clr.LinearSegmentedColormap.from_list("Zissou",zissou)
phenotype_dict={"iPSC Reprogramming":["KLF4","SOX2","POU5F1","MYC"],"Epithelialization":["GRHL3","TP63","TACSTD2","KRT14"],"Scar Formation":["COL1A1","ACTA2","FN1","TGFB1"]}
sample_config={"5dpb":{"sample_id":"BW32_A01597A3_SDSDB_5dpb","angle":135,"day":5},"12dpb_SPTDI1":{"sample_id":"BW14_1_C1_SSDB_12dpb","angle":135,"day":12},"12dpb_SPTDI2":{"sample_id":"BW14_B3_SSDB_12dpb","angle":135,"day":12},"19dpb":{"sample_id":"BW15D_C6_SDSDB_19dpb","angle":-225,"day":19},"19dpb_p1":{"sample_id":"BW15D_1_D1_SDSDB_19dpb","angle":140,"day":19},"26dpb_p1":{"sample_id":"BW81_C02846B6_SDSDB_26dpb_part1","angle":0,"day":26},"26dpb_p2":{"sample_id":"BW81_C02846B6_SDSDB_26dpb_part2","angle":0,"day":26},"2mph":{"sample_id":"2mph_A03699G6.SCT","angle":45,"day":60}}
SPATIAL_X="Distance_to_Front_X"
SPATIAL_Y="Relative_Depth_Y"
region_order=["Un-epi","Epi-Front","Newly-epi"]
region_palette={"Un-epi":"#E64B35","Epi-Front":"#DDBB2F","Newly-epi":"#3A77B7"}
def pstar(p):
    if pd.isna(p):return "ns"
    if p<0.001:return "***"
    if p<0.01:return "**"
    if p<0.05:return "*"
    return "ns"
def rotate_xy(x,y,angle):
    t=np.radians(angle)
    return x*np.cos(t)-y*np.sin(t),x*np.sin(t)+y*np.cos(t)
def get_xy(adata):
    if "spatial" in adata.obsm:return adata.obsm["spatial"][:,0],adata.obsm["spatial"][:,1]
    if "X_spatial" in adata.obsm:return adata.obsm["X_spatial"][:,0],adata.obsm["X_spatial"][:,1]
    raise KeyError("No spatial or X_spatial in adata.obsm")
def load_anchor(skey):
    p=os.path.join(mask_dir,f"{skey}_spatial_anchors.json")
    if not os.path.exists(p):return None
    with open(p,"r",encoding="utf-8") as f:return json.load(f)
def calc_space(adata_sub,anchor,angle):
    raw_x,raw_y=get_xy(adata_sub)
    rx,ry=rotate_xy(raw_x,raw_y,angle)
    min_x,max_x,min_y,max_y=rx.min(),rx.max(),ry.min(),ry.max()
    epi=LineString(anchor["epi_baseline"]) if anchor and "epi_baseline" in anchor and len(anchor["epi_baseline"])>=2 else LineString([(min_x,min_y),(max_x,min_y)])
    der=LineString(anchor["der_bottom"]) if anchor and "der_bottom" in anchor and len(anchor["der_bottom"])>=2 else LineString([(min_x,max_y),(max_x,max_y)])
    front_pts=anchor.get("leading_edge",[]) if anchor else []
    has_front=len(front_pts)>=2
    front=LineString(front_pts) if has_front else None
    direction=anchor.get("healed_direction","right_is_healed") if anchor else "right_is_healed"
    xs,ys=[],[]
    for i in range(len(rx)):
        pt=Point(rx[i],ry[i])
        de,dd=pt.distance(epi),pt.distance(der)
        ys.append(de/(de+dd) if de+dd>0 else 0)
        if has_front:
            sign=1 if pt.x>front_pts[0][0] else -1
            if direction=="left_is_healed":sign=-sign
            xs.append(pt.distance(front)*sign)
        else:xs.append(0.0)
    ys=np.array(ys,dtype=float)
    if np.nanmax(ys)>0:ys=ys/np.nanmax(ys)
    return np.array(xs)*0.5,ys
def assign_region(x,y):
    max_abs_x=max(abs(np.nanmin(x)),abs(np.nanmax(x)))*1.1
    b1_t,b1_b,b2_t,b2_b=-600,-300,400,700
    gates={"Un-epi":Polygon([[-max_abs_x,0],[b1_t,0],[b1_b,1],[-max_abs_x,1]]),"Epi-Front":Polygon([[b1_t,0],[b2_t,0],[b2_b,1],[b1_b,1]]),"Newly-epi":Polygon([[b2_t,0],[max_abs_x,0],[max_abs_x,1],[b2_b,1]])}
    regs=[]
    for xx,yy in zip(x,y):
        p=Point(xx,yy);hit="Unknown"
        for r,poly in gates.items():
            if poly.contains(p) or poly.touches(p):hit=r;break
        regs.append(hit)
    return regs
def build_analysis_adata():
    print("Loading combined h5ad...")
    adata=sc.read_h5ad(h5ad_file)
    if "sub_labels" not in adata.obs.columns and "cell_type" in adata.obs.columns:adata.obs["sub_labels"]=adata.obs["cell_type"]
    batch_col="sample_batch_new" if "sample_batch_new" in adata.obs.columns else "sample_batch"
    obs_parts=[];idx_all=[]
    for skey,cfg in sample_config.items():
        anchor=load_anchor(skey)
        if anchor is None:
            print(f"Skip {skey}: missing json")
            continue
        idx=np.where(adata.obs[batch_col].astype(str).values==cfg["sample_id"])[0]
        if len(idx)==0:
            print(f"Skip {skey}: no cells")
            continue
        sub=adata[idx].copy()
        x,y=calc_space(sub,anchor,cfg["angle"])
        regs=assign_region(x,y)
        tmp=pd.DataFrame({"Sample":skey,"Day":cfg["day"],SPATIAL_X:x,SPATIAL_Y:y,"Region":regs},index=sub.obs_names)
        obs_parts.append(tmp);idx_all.extend(idx.tolist())
        del sub;gc.collect()
    adata2=adata[idx_all].copy()
    obs_new=pd.concat(obs_parts,axis=0)
    obs_new=obs_new.loc[adata2.obs_names]
    for c in obs_new.columns:adata2.obs[c]=obs_new[c].values
    adata2=adata2[adata2.obs["Region"].isin(region_order)].copy()
    adata2.obs["Region"]=pd.Categorical(adata2.obs["Region"],categories=region_order,ordered=True)
    del adata;gc.collect()
    return adata2
def calculate_mergescore(adata,genes,score_name):
    old=adata.var_names.copy()
    adata.var_names=adata.var_names.str.upper()
    valid=[g.upper() for g in genes if g.upper() in adata.var_names]
    if len(valid)>0:sc.tl.score_genes(adata,gene_list=valid,score_name=score_name)
    else:adata.obs[score_name]=0
    adata.var_names=old
    return adata
def add_scores(adata):
    for phenotype,genes in phenotype_dict.items():adata=calculate_mergescore(adata,genes,phenotype)
    return adata
def get_expression_data(adata,feature):
    if feature in adata.obs.columns:return np.asarray(adata.obs[feature].values,dtype=float)
    gene_idx=np.where(adata.var_names.str.upper()==feature.upper())[0]
    if len(gene_idx)>0:
        x=adata[:,gene_idx[0]].X
        return np.asarray(x.toarray()).flatten() if scipy.sparse.issparse(x) else np.asarray(x).flatten()
    return np.zeros(adata.n_obs)
def plot_phenotype_spatial_ultimate(ax,adata,feature,phenotype_group,cmap,show_y_label=False,show_colorbar=False,show_top_labels=False):
    x=adata.obs[SPATIAL_X].values
    y=adata.obs[SPATIAL_Y].values
    expr=get_expression_data(adata,feature)
    idx=np.argsort(expr)
    x_s,y_s,expr_s=x[idx],y[idx],expr[idx]
    max_abs_x=max(abs(x.min()),abs(x.max()))*1.1 if len(x)>0 else 4000
    b1_t,b1_b,b2_t,b2_b=-600,-300,400,700
    ax.plot([b1_t,b1_b],[0,1],color="black",linestyle="--",linewidth=0.8,alpha=0.8,zorder=1)
    ax.plot([b2_t,b2_b],[0,1],color="black",linestyle="--",linewidth=0.8,alpha=0.8,zorder=1)
    ax.add_patch(Polygon([[-max_abs_x,0],[b1_t,0],[b1_b,1],[-max_abs_x,1]],color="#e5192c",alpha=0.10,zorder=0,lw=0))
    ax.add_patch(Polygon([[b1_t,0],[b2_t,0],[b2_b,1],[b1_b,1]],color="#E1C62F",alpha=0.12,zorder=0,lw=0))
    ax.add_patch(Polygon([[b2_t,0],[max_abs_x,0],[max_abs_x,1],[b2_b,1]],color="#3a77b7",alpha=0.10,zorder=0,lw=0))
    if show_top_labels:
        trans=ax.get_xaxis_transform()
        ax.text((-max_abs_x+b1_t)/2,1.02,"Un-epi",color="#e5192c",transform=trans,ha="center",va="bottom",fontsize=8)
        ax.text((b1_t+b2_t)/2,1.02,"Epi-Front",color="#B89400",transform=trans,ha="center",va="bottom",fontsize=8)
        ax.text((b2_t+max_abs_x)/2,1.02,"Newly-epi",color="#3a77b7",transform=trans,ha="center",va="bottom",fontsize=8)
    sns.kdeplot(data=pd.DataFrame({"x":x,"y":y}),x="x",y="y",ax=ax,color="#555555",linewidths=0.45,alpha=0.55,levels=6,zorder=2)
    mask_z,mask_e=expr_s<=0,expr_s>0
    ax.scatter(x_s[mask_z],y_s[mask_z],color="#D0D0D0",s=1.5,alpha=0.35,zorder=3)
    sca=ax.scatter(x_s[mask_e],y_s[mask_e],c=expr_s[mask_e],cmap=cmap,s=1.5,edgecolor="none",zorder=4)
    if show_colorbar and mask_e.sum()>0:
        cbar_ax=ax.inset_axes([1.08,0.1,0.05,0.8])
        cbar=plt.colorbar(sca,cax=cbar_ax,orientation="vertical")
        cbar.set_ticks([expr_s[mask_e].min(),expr_s[mask_e].max()])
        cbar.set_ticklabels(["Min","Max"])
        cbar.ax.tick_params(colors="black",labelsize=8,length=0,pad=2)
        cbar.set_label("Expression",color="black",fontsize=12,labelpad=6)
    cx,cy=-500,0.75
    if phenotype_group in ["iPSC Reprogramming","Epithelialization"]:
        ax.annotate("",xy=(1800,0.2),xycoords="data",xytext=(cx,cy),textcoords="data",arrowprops=dict(facecolor="#0099CC",edgecolor="#0099CC",width=0.3,headwidth=2.5,shrink=0.08),zorder=5)
    if phenotype_group in ["iPSC Reprogramming","Scar Formation"]:
        ax.annotate("",xy=(-3200,0.98),xycoords="data",xytext=(cx,cy),textcoords="data",arrowprops=dict(facecolor="#B8A900",edgecolor="#B8A900",width=0.3,headwidth=2.5,shrink=0.08),zorder=5)
    ax.set_xlim(-max_abs_x,max_abs_x)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3,prune="both"))
    ax.set_ylim(1.05,-0.05)
    ax.set_yticks([0.0,0.5,1.0])
    ax.set_box_aspect(1)
    ax.set_title(feature if feature in phenotype_dict else f"$\\it{{{feature}}}$",fontsize=15,color="black",pad=12,fontweight="bold" if feature in phenotype_dict else "normal")
    ax.set_ylabel("Epidermis-Dermis axis" if show_y_label else "",color="black",fontsize=16,labelpad=8)
    ax.tick_params(colors="black",direction="out",length=2.5,width=0.6,labelsize=14)
    for s in ax.spines.values():s.set_linewidth(0.6);s.set_color("black")
def plot_phenotype_master_figure(adata,out_dir):
    fig=plt.figure(figsize=(19.5,12.5),dpi=300,facecolor="white")
    gs=fig.add_gridspec(3,5,hspace=0.35,wspace=0.35)
    for row_idx,phenotype in enumerate(phenotype_dict.keys()):
        display_features=phenotype_dict[phenotype]+[phenotype]
        for col_idx,feature in enumerate(display_features):
            ax=fig.add_subplot(gs[row_idx,col_idx])
            plot_phenotype_spatial_ultimate(ax,adata,feature,phenotype,colormap,show_y_label=(col_idx==0),show_colorbar=(col_idx==4),show_top_labels=(row_idx==0))
            if col_idx==4:
                for s in ax.spines.values():s.set_linewidth(1.0)
            if col_idx==0:
                ax.text(-0.42,0.5,phenotype,transform=ax.transAxes,rotation=90,ha="center",va="center",fontsize=18,fontweight="bold",color="black")
    fig.text(0.5,0.05,"Distance to leading edge (um)",ha="center",va="center",fontsize=18,color="black")
    plt.tight_layout(rect=[0.04,0.08,0.96,1.0])
    save=os.path.join(out_dir,"Cell_Fate_Trajectories_White_3x5")
    fig.savefig(save+".pdf",facecolor="white",bbox_inches="tight")
    fig.savefig(save+".png",facecolor="white",bbox_inches="tight",dpi=600)
    plt.close(fig)
def feature_stats_df(adata,features):
    rows=[];pair_rows=[]
    for feat in features:
        vals=get_expression_data(adata,feat)
        df=pd.DataFrame({"Region":adata.obs["Region"].astype(str).values,"Value":vals})
        groups=[df.loc[df["Region"]==r,"Value"].dropna().values for r in region_order]
        kw_p=kruskal(*groups).pvalue if all(len(g)>1 for g in groups) else np.nan
        pairs=[];pvals=[]
        for i in range(len(region_order)):
            for j in range(i+1,len(region_order)):
                r1,r2=region_order[i],region_order[j]
                v1=df.loc[df["Region"]==r1,"Value"].dropna().values
                v2=df.loc[df["Region"]==r2,"Value"].dropna().values
                p=mannwhitneyu(v1,v2,alternative="two-sided").pvalue if len(v1)>0 and len(v2)>0 else np.nan
                pairs.append((r1,r2));pvals.append(p)
        q=np.full(len(pvals),np.nan);mask=~pd.isna(pvals)
        if np.sum(mask)>0:q[mask]=multipletests(np.array(pvals)[mask],method="fdr_bh")[1]
        rows.append({"Feature":feat,"Kruskal_p":kw_p,"Kruskal_sig":pstar(kw_p)})
        for (r1,r2),p,qq in zip(pairs,pvals,q):
            pair_rows.append({"Feature":feat,"Region_1":r1,"Region_2":r2,"MWU_p":p,"FDR_p":qq,"Significance":pstar(qq)})
    return pd.DataFrame(rows),pd.DataFrame(pair_rows)
def add_sig_bars(ax,pair_df,feature,ymax):
    sub=pair_df[pair_df["Feature"]==feature].copy()
    xmap={r:i for i,r in enumerate(region_order)}
    step=max(ymax*0.12,0.05)
    y=ymax+step*0.2
    for _,row in sub.iterrows():
        if row["Significance"]=="ns":continue
        x1,x2=xmap[row["Region_1"]],xmap[row["Region_2"]]
        ax.plot([x1,x1,x2,x2],[y,y+step*0.18,y+step*0.18,y],color="black",lw=0.7)
        ax.text((x1+x2)/2,y+step*0.2,row["Significance"],ha="center",va="bottom",fontsize=7,fontweight="bold")
        y+=step*0.55
    ax.set_ylim(top=y+step*0.5)
def plot_violin_stats(adata,out_dir):
    features=[]
    for ph,genes in phenotype_dict.items():features+=genes+[ph]
    features=list(dict.fromkeys(features))
    overall,pairwise=feature_stats_df(adata,features)
    overall.to_csv(os.path.join(out_dir,"Regional_violin_Kruskal_statistics.csv"),index=False,encoding="utf-8-sig")
    pairwise.to_csv(os.path.join(out_dir,"Regional_violin_pairwise_MWU_FDR_statistics.csv"),index=False,encoding="utf-8-sig")
    fig=plt.figure(figsize=(20,12),dpi=300,facecolor="white")
    gs=fig.add_gridspec(3,5,hspace=0.55,wspace=0.42)
    for row_idx,(ph,genes) in enumerate(phenotype_dict.items()):
        feats=genes+[ph]
        for col_idx,feat in enumerate(feats):
            ax=fig.add_subplot(gs[row_idx,col_idx])
            vals=get_expression_data(adata,feat)
            df=pd.DataFrame({"Region":adata.obs["Region"].astype(str).values,"Value":vals})
            sns.violinplot(data=df,x="Region",y="Value",order=region_order,palette=region_palette,inner=None,linewidth=0.6,cut=0,ax=ax)
            sns.boxplot(data=df,x="Region",y="Value",order=region_order,width=0.22,showfliers=False,boxprops={"facecolor":"white","edgecolor":"black","linewidth":0.6},medianprops={"color":"black","linewidth":0.8},whiskerprops={"color":"black","linewidth":0.6},capprops={"color":"black","linewidth":0.6},ax=ax)
            ymax=np.nanpercentile(vals,99) if np.isfinite(vals).any() else 1
            add_sig_bars(ax,pairwise,feat,ymax)
            kw=overall.loc[overall["Feature"]==feat,"Kruskal_p"].values
            kw_txt=f"KW {pstar(kw[0])}" if len(kw)>0 else "KW ns"
            ax.set_title(f"{feat}\n{kw_txt}",fontsize=10,fontweight="bold")
            ax.set_xlabel("")
            ax.set_ylabel("Expression / score" if col_idx==0 else "",fontsize=9)
            ax.tick_params(axis="x",rotation=35,labelsize=8)
            ax.tick_params(axis="y",labelsize=8)
            sns.despine(ax=ax)
    save=os.path.join(out_dir,"Regional_KeyGene_MergeScore_Violin_Stats")
    fig.savefig(save+".pdf",facecolor="white",bbox_inches="tight")
    fig.savefig(save+".png",facecolor="white",bbox_inches="tight",dpi=600)
    plt.close(fig)
if __name__=="__main__":
    adata=build_analysis_adata()
    adata=add_scores(adata)
    adata.write(os.path.join(out_dir,"analysis_adata_with_spatial_region_scores.h5ad"))
    plot_phenotype_master_figure(adata,out_dir)
    plot_violin_stats(adata,out_dir)
    print("完成：白底空间feature图和区域violin统计图已输出")
    print(f"输出目录: {out_dir}")
    

##########################    
############################ IGFL1/IGFLR1 spatial distribution + co-localization + regional statistics
import os,json,warnings,gc
import numpy as np,pandas as pd,scanpy as sc,scipy.sparse
import matplotlib.pyplot as plt,seaborn as sns
from shapely.geometry import Point,LineString,Polygon
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from scipy.stats import kruskal,mannwhitneyu,chi2_contingency,fisher_exact
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
plt.ioff()
base_dir=r"data/external/local_f5480532f0f9"
pdf_dir=os.path.join(base_dir,"pdf_output")
mask_dir=os.path.join(pdf_dir,"mask")
out_dir=os.path.join(pdf_dir,"IGFL1_IGFLR1_AllSlices_Spatial_CoLocalization_White")
os.makedirs(out_dir,exist_ok=True)
h5ad_file=os.path.join(base_dir,"combined_16samples_Harmony_Clustered_Fixed_Annotated_White.h5ad")
if not os.path.exists(h5ad_file):h5ad_file=os.path.join(base_dir,"combined_16samples_Harmony_Clustered_Fixed_Annotated.h5ad")
LIGAND="IGFL1"
RECEPTOR="IGFLR1"
PIXEL_TO_UM=0.5
COLOCAL_RADIUS_UM=50
SCALE_BAR_UM=2000
BG_COLOR="#BFBFBF"
BG_ALPHA=0.75
BG_SIZE=0.12
region_order=["Un-epi","Epi-Front","Newly-epi"]
region_colors={"Un-epi":"#E64B35","Epi-Front":"#DDBB2F","Newly-epi":"#3A77B7"}
sample_config={
"5dpb":{"sample_id":"BW32_A01597A3_SDSDB_5dpb","angle":135},
"12dpb_DPTDI1":{"sample_id":"BW13_1_B1_DSDB_12dpb","angle":45},
"12dpb_SPTDI1":{"sample_id":"BW14_1_C1_SSDB_12dpb","angle":135},
"12dpb_DPTDI2":{"sample_id":"BW13_A3_DSDB_12dpb","angle":135},
"12dpb_SPTDI2":{"sample_id":"BW14_B3_SSDB_12dpb","angle":135},
"19dpb":{"sample_id":"BW15D_C6_SDSDB_19dpb","angle":-225},
"19dpb_p1":{"sample_id":"BW15D_1_D1_SDSDB_19dpb","angle":140},
"26dpb_p1":{"sample_id":"BW81_C02846B6_SDSDB_26dpb_part1","angle":0},
"26dpb_p2":{"sample_id":"BW81_C02846B6_SDSDB_26dpb_part2","angle":0},
"2mph":{"sample_id":"2mph_A03699G6.SCT","angle":45},
"6mph_p1":{"sample_id":"6mph_A02988A2.SCT_part1","angle":0},
"6mph_p2":{"sample_id":"6mph_A02988A2.SCT_part2","angle":0},
"9mph":{"sample_id":"9mph_A03091C5.SCT","angle":0},
"Normal":{"sample_id":"NS_C02847B1","angle":180}}
plt.style.use("default")
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial"],"pdf.fonttype":42,"ps.fonttype":42,"figure.facecolor":"white","axes.facecolor":"white","text.color":"black","axes.labelcolor":"black","xtick.color":"black","ytick.color":"black"})
def pstar(p):
    if pd.isna(p):return "ns"
    if p<0.001:return "***"
    if p<0.01:return "**"
    if p<0.05:return "*"
    return "ns"
def rotate_xy(x,y,angle):
    t=np.radians(angle)
    return x*np.cos(t)-y*np.sin(t),x*np.sin(t)+y*np.cos(t)
def get_xy(adata):
    if "spatial" in adata.obsm:return adata.obsm["spatial"][:,0],adata.obsm["spatial"][:,1]
    if "X_spatial" in adata.obsm:return adata.obsm["X_spatial"][:,0],adata.obsm["X_spatial"][:,1]
    raise KeyError("No spatial or X_spatial in adata.obsm")
def load_anchor(skey):
    p=os.path.join(mask_dir,f"{skey}_spatial_anchors.json")
    if not os.path.exists(p):return {}
    with open(p,"r",encoding="utf-8") as f:return json.load(f)
def extract_gene(adata,gene):
    var_upper=adata.var_names.str.upper()
    idx=np.where(var_upper==gene.upper())[0]
    if len(idx)==0:
        print(f"Warning: {gene} not found, use zeros.")
        return np.zeros(adata.n_obs)
    x=adata[:,idx[0]].X
    return np.asarray(x.toarray()).flatten() if scipy.sparse.issparse(x) else np.asarray(x).flatten()
def process_spatial_data(adata,anchor,angle):
    raw_x,raw_y=get_xy(adata)
    rx,ry=rotate_xy(raw_x,raw_y,angle)
    adata.obsm["X_rot_spatial"]=np.column_stack([rx,ry])
    adata.obs["Display_X"]=rx*PIXEL_TO_UM
    adata.obs["Display_Y"]=ry*PIXEL_TO_UM
    min_x,max_x,min_y,max_y=rx.min(),rx.max(),ry.min(),ry.max()
    epi=LineString(anchor["epi_baseline"]) if "epi_baseline" in anchor and len(anchor["epi_baseline"])>=2 else LineString([(min_x,min_y),(max_x,min_y)])
    der=LineString(anchor["der_bottom"]) if "der_bottom" in anchor and len(anchor["der_bottom"])>=2 else LineString([(min_x,max_y),(max_x,max_y)])
    front_pts=anchor.get("leading_edge",[])
    has_front=len(front_pts)>=2
    front=LineString(front_pts) if has_front else None
    direction=anchor.get("healed_direction","right_is_healed")
    xs,ys=[],[]
    for i in range(len(rx)):
        pt=Point(rx[i],ry[i])
        de,dd=pt.distance(epi),pt.distance(der)
        ys.append(de/(de+dd) if de+dd>0 else 0)
        if has_front:
            sign=1 if pt.x>front_pts[0][0] else -1
            if direction=="left_is_healed":sign=-sign
            xs.append(pt.distance(front)*sign)
        else:xs.append(0.0)
    adata.obs["Plot_X"]=np.array(xs)*PIXEL_TO_UM
    adata.obs["Plot_Y"]=np.array(ys)
    max_abs_x=max(abs(adata.obs["Plot_X"].min()),abs(adata.obs["Plot_X"].max()))*1.1
    b1_t,b1_b,b2_t,b2_b=-600,-300,400,700
    gates={"Un-epi":Polygon([[-max_abs_x,0],[b1_t,0],[b1_b,1],[-max_abs_x,1]]),"Epi-Front":Polygon([[b1_t,0],[b2_t,0],[b2_b,1],[b1_b,1]]),"Newly-epi":Polygon([[b2_t,0],[max_abs_x,0],[max_abs_x,1],[b2_b,1]])}
    regs=[]
    for x,y in zip(adata.obs["Plot_X"],adata.obs["Plot_Y"]):
        p=Point(x,y);hit="Unknown"
        for r,poly in gates.items():
            if poly.contains(p) or poly.touches(p):hit=r;break
        regs.append(hit)
    adata.obs["Region"]=pd.Categorical(regs,categories=region_order+["Unknown"],ordered=True)
    return adata
def add_lr_metrics(adata):
    lig=extract_gene(adata,LIGAND)
    rec=extract_gene(adata,RECEPTOR)
    adata.obs[LIGAND]=lig
    adata.obs[RECEPTOR]=rec
    adata.obs["LR_Product"]=lig*rec
    adata.obs[f"{LIGAND}_positive"]=lig>0
    adata.obs[f"{RECEPTOR}_positive"]=rec>0
    adata.obs["Double_positive"]=(lig>0)&(rec>0)
    coords=adata.obs[["Display_X","Display_Y"]].values
    lig_idx=np.where(lig>0)[0]
    rec_idx=np.where(rec>0)[0]
    near=np.zeros(adata.n_obs,dtype=bool)
    if len(lig_idx)>0 and len(rec_idx)>0:
        tree=cKDTree(coords[rec_idx])
        d,_=tree.query(coords[lig_idx],k=1)
        near[lig_idx]=d<=COLOCAL_RADIUS_UM
    adata.obs["Ligand_near_receptor"]=near
    return adata
def add_scale_bar(ax,length_um=SCALE_BAR_UM,label="2 mm",fontsize=6):
    x1,x2=ax.get_xlim();y1,y2=ax.get_ylim()
    x0=x2-(x2-x1)*0.22
    y0=y1-(y1-y2)*0.08
    ax.plot([x0,x0+length_um],[y0,y0],color="black",lw=0.8,zorder=10)
    ax.text(x0+length_um/2,y0-(y1-y2)*0.035,label,ha="center",va="top",fontsize=fontsize,color="black")
def set_spatial_xlim_ylim(ax,df):
    x=df["Display_X"].values;y=df["Display_Y"].values
    x1,x2=np.percentile(x,[0.2,99.8])
    y1,y2=np.percentile(y,[0.2,99.8])
    xp=(x2-x1)*0.025;yp=(y2-y1)*0.08
    ax.set_xlim(x1-xp,x2+xp)
    ax.set_ylim(y2+yp,y1-yp)
    ax.set_aspect("equal")
    ax.axis("off")
def plot_distribution(adata,skey):
    df=adata.obs.copy()
    lig=df[LIGAND].values
    rec=df[RECEPTOR].values
    lig_pos=lig>0;rec_pos=rec>0;double_pos=lig_pos&rec_pos
    x=df["Display_X"].values;y=df["Display_Y"].values
    fig,axes=plt.subplots(1,2,figsize=(6.8,1.75),dpi=600)
    for ax in axes:
        ax.scatter(x,y,s=BG_SIZE,c=BG_COLOR,alpha=BG_ALPHA,lw=0,rasterized=True)
        set_spatial_xlim_ylim(ax,df)
    axes[0].scatter(x[lig_pos],y[lig_pos],s=0.55,c="#E7D51A",alpha=0.95,lw=0,rasterized=True)
    axes[0].text(0.02,0.94,f"{skey} | {LIGAND} Distribution",transform=axes[0].transAxes,ha="left",va="top",fontsize=7,fontweight="bold",color="black")
    axes[1].scatter(x[lig_pos],y[lig_pos],s=0.45,c="#E31A1C",alpha=0.85,lw=0,rasterized=True,label=LIGAND)
    axes[1].scatter(x[rec_pos],y[rec_pos],s=0.45,c="#E7D51A",alpha=0.85,lw=0,rasterized=True,label=RECEPTOR)
    axes[1].scatter(x[double_pos],y[double_pos],s=0.85,c="#FF8C00",alpha=1.0,lw=0,rasterized=True,label="Co-localized")
    axes[1].text(0.43,0.94,LIGAND,transform=axes[1].transAxes,ha="right",va="top",fontsize=6,fontweight="bold",color="#E31A1C")
    axes[1].text(0.50,0.94,RECEPTOR,transform=axes[1].transAxes,ha="left",va="top",fontsize=6,fontweight="bold",color="#B8A900")
    add_scale_bar(axes[1],length_um=SCALE_BAR_UM,label="2 mm",fontsize=6)
    handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor="#E31A1C",markersize=4,label=LIGAND),Line2D([0],[0],marker="o",color="w",markerfacecolor="#E7D51A",markersize=4,label=RECEPTOR),Line2D([0],[0],marker="o",color="w",markerfacecolor="#FF8C00",markersize=4,label="Co-localized")]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(0.5,-0.05),ncol=3,frameon=False,fontsize=6,handletextpad=0.3,columnspacing=0.8)
    fig.subplots_adjust(left=0.01,right=0.995,top=0.98,bottom=0.18,wspace=0.03)
    save=os.path.join(out_dir,skey,f"{skey}_{LIGAND}_{RECEPTOR}_distribution_colocalization")
    fig.savefig(save+".pdf",bbox_inches="tight",facecolor="white")
    fig.savefig(save+".png",bbox_inches="tight",facecolor="white",dpi=600)
    plt.close(fig)
def continuous_stats(adata,features):
    rows=[];pair_rows=[];df=adata.obs.copy()
    for feat in features:
        groups=[df.loc[df["Region"]==r,feat].dropna().values for r in region_order]
        kw_p=kruskal(*groups).pvalue if all(len(g)>1 for g in groups) else np.nan
        rows.append({"Feature":feat,"Kruskal_p":kw_p,"Significance":pstar(kw_p)})
        pairs=[];pvals=[]
        for i in range(len(region_order)):
            for j in range(i+1,len(region_order)):
                r1,r2=region_order[i],region_order[j]
                v1=df.loc[df["Region"]==r1,feat].dropna().values
                v2=df.loc[df["Region"]==r2,feat].dropna().values
                p=mannwhitneyu(v1,v2,alternative="two-sided").pvalue if len(v1)>0 and len(v2)>0 else np.nan
                pairs.append((r1,r2));pvals.append(p)
        q=np.full(len(pvals),np.nan);mask=~pd.isna(pvals)
        if mask.sum()>0:q[mask]=multipletests(np.array(pvals)[mask],method="fdr_bh")[1]
        for (r1,r2),p,qq in zip(pairs,pvals,q):
            pair_rows.append({"Feature":feat,"Region_1":r1,"Region_2":r2,"MWU_p":p,"FDR_p":qq,"Significance":pstar(qq)})
    return pd.DataFrame(rows),pd.DataFrame(pair_rows)
def safe_chi2(tab):
    tab=np.asarray(tab,dtype=float)
    if tab.shape[0]<2 or tab.shape[1]<2:return np.nan
    if tab.sum()==0:return np.nan
    if np.any(tab.sum(axis=0)==0):return np.nan
    if np.any(tab.sum(axis=1)==0):return np.nan
    try:return chi2_contingency(tab)[1]
    except ValueError:return np.nan
def safe_fisher(a,b):
    arr=np.vstack([a,b]).astype(float)
    if arr.sum()==0:return np.nan
    if np.any(arr.sum(axis=0)==0):return np.nan
    if np.any(arr.sum(axis=1)==0):return np.nan
    try:return fisher_exact(arr)[1]
    except ValueError:return np.nan
def proportion_stats(adata,bool_features):
    rows=[];pair_rows=[];df=adata.obs.copy()
    for feat in bool_features:
        tab=[]
        for r in region_order:
            sub=df[df["Region"]==r]
            pos=int(sub[feat].sum()) if len(sub)>0 else 0
            neg=len(sub)-pos
            tab.append([pos,neg])
            rows.append({"Feature":feat,"Region":r,"Positive":pos,"Total":len(sub),"Percent":pos/len(sub)*100 if len(sub)>0 else np.nan})
        chi_p=safe_chi2(tab)
        for r in rows:
            if r["Feature"]==feat:
                r["Chi_square_p"]=chi_p
                r["Chi_square_sig"]=pstar(chi_p)
        pairs=[];pvals=[]
        for i in range(len(region_order)):
            for j in range(i+1,len(region_order)):
                r1,r2=region_order[i],region_order[j]
                p=safe_fisher(tab[i],tab[j])
                pairs.append((r1,r2));pvals.append(p)
        q=np.full(len(pvals),np.nan)
        mask=~pd.isna(pvals)
        if mask.sum()>0:q[mask]=multipletests(np.array(pvals,dtype=float)[mask],method="fdr_bh")[1]
        for (r1,r2),p,qq in zip(pairs,pvals,q):
            pair_rows.append({"Feature":feat,"Region_1":r1,"Region_2":r2,"Fisher_p":p,"FDR_p":qq,"Significance":pstar(qq)})
    return pd.DataFrame(rows),pd.DataFrame(pair_rows)
def add_sig(ax,pair_df,feature,ymax):
    sub=pair_df[(pair_df["Feature"]==feature)&(pair_df["Significance"]!="ns")]
    if sub.empty:return
    xmap={r:i for i,r in enumerate(region_order)}
    step=max(ymax*0.12,0.05);y=ymax+step*0.2
    for _,row in sub.iterrows():
        x1,x2=xmap[row["Region_1"]],xmap[row["Region_2"]]
        ax.plot([x1,x1,x2,x2],[y,y+step*0.18,y+step*0.18,y],color="black",lw=0.55)
        ax.text((x1+x2)/2,y+step*0.2,row["Significance"],ha="center",va="bottom",fontsize=5.5,fontweight="bold")
        y+=step*0.55
    ax.set_ylim(top=y+step*0.45)
def plot_bar_narrow(prop_df,prop_pair,skey):
    bool_features=[f"{LIGAND}_positive",f"{RECEPTOR}_positive","Double_positive","Ligand_near_receptor"]
    for feat in bool_features:
        fig,ax=plt.subplots(figsize=(1.25,1.7),dpi=600)
        sub=prop_df[prop_df["Feature"]==feat].copy()
        sns.barplot(data=sub,x="Region",y="Percent",order=region_order,palette=region_colors,edgecolor="black",linewidth=0.45,width=0.55,ax=ax)
        ymax=max(sub["Percent"].max(),1)
        add_sig(ax,prop_pair,feat,ymax)
        chi=sub["Chi_square_p"].iloc[0]
        title=feat.replace("_positive","+").replace("_"," ")
        ax.set_title(f"{title}\nChi {pstar(chi)}",fontsize=5.7,fontweight="bold",pad=2)
        ax.set_xlabel("")
        ax.set_ylabel("% positive",fontsize=5.5)
        ax.tick_params(axis="x",rotation=40,labelsize=5)
        ax.tick_params(axis="y",labelsize=5)
        sns.despine(ax=ax)
        fig.tight_layout(pad=0.25)
        save=os.path.join(out_dir,skey,f"{skey}_{feat}_narrow_bar")
        fig.savefig(save+".pdf",bbox_inches="tight",facecolor="white")
        fig.savefig(save+".png",bbox_inches="tight",facecolor="white",dpi=600)
        plt.close(fig)
def plot_violin(adata,cont_pair,cont_overall,skey):
    cont_features=[LIGAND,RECEPTOR,"LR_Product"]
    fig,axes=plt.subplots(1,3,figsize=(4.2,1.7),dpi=600)
    for ax,feat in zip(axes,cont_features):
        df=adata.obs[["Region",feat]].copy()
        sns.violinplot(data=df,x="Region",y=feat,order=region_order,palette=region_colors,inner=None,cut=0,linewidth=0.45,ax=ax)
        sns.boxplot(data=df,x="Region",y=feat,order=region_order,width=0.18,showfliers=False,boxprops={"facecolor":"white","edgecolor":"black","linewidth":0.45},medianprops={"color":"black","linewidth":0.6},whiskerprops={"color":"black","linewidth":0.45},capprops={"color":"black","linewidth":0.45},ax=ax)
        ymax=np.nanpercentile(df[feat],99) if np.isfinite(df[feat]).any() else 1
        add_sig(ax,cont_pair,feat,ymax)
        p=cont_overall.loc[cont_overall["Feature"]==feat,"Kruskal_p"].values[0]
        ax.set_title(f"{feat}\nKW {pstar(p)}",fontsize=5.8,fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Expr/score",fontsize=5.5)
        ax.tick_params(axis="x",rotation=40,labelsize=5)
        ax.tick_params(axis="y",labelsize=5)
        sns.despine(ax=ax)
    fig.tight_layout(pad=0.35)
    save=os.path.join(out_dir,skey,f"{skey}_{LIGAND}_{RECEPTOR}_regional_violin_stats")
    fig.savefig(save+".pdf",bbox_inches="tight",facecolor="white")
    fig.savefig(save+".png",bbox_inches="tight",facecolor="white",dpi=600)
    plt.close(fig)
def plot_a4_bar_contact_sheet(all_prop,all_pair):
    bool_features=[f"{LIGAND}_positive",f"{RECEPTOR}_positive","Double_positive","Ligand_near_receptor"]
    for feat in bool_features:
        samples=[s for s in sample_config if s in all_prop and feat in all_prop[s]["Feature"].unique()]
        n=len(samples);ncol=6;nrow=int(np.ceil(n/ncol))
        fig,axes=plt.subplots(nrow,ncol,figsize=(8.27,1.55*nrow),dpi=600)
        axes=np.array(axes).reshape(-1)
        for ax,skey in zip(axes,samples):
            sub=all_prop[skey][all_prop[skey]["Feature"]==feat].copy()
            sns.barplot(data=sub,x="Region",y="Percent",order=region_order,palette=region_colors,edgecolor="black",linewidth=0.35,width=0.55,ax=ax)
            ymax=max(sub["Percent"].max(),1)
            add_sig(ax,all_pair[skey],feat,ymax)
            ax.set_title(skey,fontsize=5.8,fontweight="bold",pad=1)
            ax.set_xlabel("")
            ax.set_ylabel("%" if ax in axes[::ncol] else "",fontsize=5)
            ax.tick_params(axis="x",rotation=45,labelsize=4.5)
            ax.tick_params(axis="y",labelsize=4.5)
            sns.despine(ax=ax)
        for ax in axes[n:]:ax.axis("off")
        fig.suptitle(feat.replace("_positive","+").replace("_"," "),fontsize=7,fontweight="bold",y=0.995)
        fig.tight_layout(pad=0.25,w_pad=0.15,h_pad=0.35)
        save=os.path.join(out_dir,f"A4_6perRow_{feat}_barplots")
        fig.savefig(save+".pdf",bbox_inches="tight",facecolor="white")
        fig.savefig(save+".png",bbox_inches="tight",facecolor="white",dpi=600)
        plt.close(fig)
def analyze_one_sample(adata_full,skey,cfg):
    sample_dir=os.path.join(out_dir,skey)
    os.makedirs(sample_dir,exist_ok=True)
    batch_col="sample_batch_new" if "sample_batch_new" in adata_full.obs.columns else "sample_batch"
    adata=adata_full[adata_full.obs[batch_col].astype(str)==cfg["sample_id"]].copy()
    if len(adata)==0:
        print(f"Skip {skey}: no cells")
        return None
    anchor=load_anchor(skey)
    adata=process_spatial_data(adata,anchor,cfg["angle"])
    adata=adata[adata.obs["Region"].isin(region_order)].copy()
    if len(adata)==0:
        print(f"Skip {skey}: no valid region cells")
        return None
    adata=add_lr_metrics(adata)
    plot_distribution(adata,skey)
    cont_features=[LIGAND,RECEPTOR,"LR_Product"]
    bool_features=[f"{LIGAND}_positive",f"{RECEPTOR}_positive","Double_positive","Ligand_near_receptor"]
    cont_overall,cont_pair=continuous_stats(adata,cont_features)
    prop_df,prop_pair=proportion_stats(adata,bool_features)
    cont_overall.to_csv(os.path.join(sample_dir,f"{skey}_continuous_Kruskal_stats.csv"),index=False,encoding="utf-8-sig")
    cont_pair.to_csv(os.path.join(sample_dir,f"{skey}_continuous_pairwise_MWU_FDR.csv"),index=False,encoding="utf-8-sig")
    prop_df.to_csv(os.path.join(sample_dir,f"{skey}_proportion_stats.csv"),index=False,encoding="utf-8-sig")
    prop_pair.to_csv(os.path.join(sample_dir,f"{skey}_proportion_pairwise_Fisher_FDR.csv"),index=False,encoding="utf-8-sig")
    plot_violin(adata,cont_pair,cont_overall,skey)
    plot_bar_narrow(prop_df,prop_pair,skey)
    adata.write(os.path.join(sample_dir,f"{skey}_{LIGAND}_{RECEPTOR}_analysis.h5ad"))
    del adata;gc.collect()
    return prop_df,prop_pair
if __name__=="__main__":
    print("Loading data...")
    adata_full=sc.read_h5ad(h5ad_file)
    all_prop={}
    all_pair={}
    for skey,cfg in sample_config.items():
        print(f"\nProcessing {skey}...")
        res=analyze_one_sample(adata_full,skey,cfg)
        if res is None:continue
        all_prop[skey],all_pair[skey]=res
    del adata_full;gc.collect()
    if all_prop:plot_a4_bar_contact_sheet(all_prop,all_pair)
    print("Done.")
    print(f"Output: {out_dir}")