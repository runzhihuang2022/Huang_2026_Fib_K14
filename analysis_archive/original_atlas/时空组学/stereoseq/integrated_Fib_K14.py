# ==========================================
# [阶段一：铁证保存版] 批量交互式垂直分割与 Mask 生成器 
# ==========================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import gc

plt.ioff() 

# 📁 路径配置
data_path = r"data/external/local_611e695933ec"
output_dir = r"data/external/local_8e645d0a9b03"
stitched_dir = r"data/external/local_92743553c2f2" 
os.makedirs(output_dir, exist_ok=True)

# 🌟 样本配置字典 (细分了 Mixed_2 和 Mixed_3)
sample_config = {
    "5dpb":        {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135,  "type": "Mixed_2"},
    "12dpb_SPTDI": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135,  "type": "Pure_SPTDI"},
    "12dpb_DPTDI": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45,   "type": "Pure_DPTDI"},
    "19dpb":       {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225, "type": "Mixed_2"}, 
    "26dpb":       {"sample_id": "BW81_C02846B6_SDSDB_26dpb", "angle": 0,    "type": "Mixed_2"},   
    "2mph":        {"sample_id": "2mph_A03699G6.SCT",         "angle": 45,   "type": "Mixed_3"}, 
    "Normal":      {"sample_id": "NS_C02847B1",               "angle": 180,  "type": "Pure_Normal"}  
}

print("⏳ 正在读取母矩阵基石信息...")
adata_full = sc.read_h5ad(data_path, backed='r') 

def rotate_coordinates(x, y, angle):
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

for TEST_KEY, config in sample_config.items():
    TARGET_SAMPLE = config["sample_id"]
    ROTATION_ANGLE = config["angle"]
    SAMPLE_TYPE = config["type"]
    
    mask_csv_path = os.path.join(output_dir, f"{TEST_KEY}_Spot_Region_Annotations.csv")
    mask_png_path = os.path.join(output_dir, f"{TEST_KEY}_Tissue_Mask.png")
    
    # 智能跳过机制
    if os.path.exists(mask_csv_path) and os.path.exists(mask_png_path):
        print(f"\n⏭️ 检测到 {TEST_KEY} 的CSV与Mask图均已存在，自动跳过。")
        continue

    print(f"\n📂 正在处理样本: {TEST_KEY} ...")
    
    # 🌟 智能分流读取机制
    stitched_file = os.path.join(stitched_dir, f"{TARGET_SAMPLE}_Stitched.h5ad")
    if os.path.exists(stitched_file):
        print(f"✨ 发现修补记录！正在加载拼接版样本: {stitched_file}")
        adata = sc.read_h5ad(stitched_file)
        is_stitched = True
    else:
        print("🔗 正在从母矩阵中提取原始样本...")
        adata = adata_full[adata_full.obs['sample_batch_new'] == TARGET_SAMPLE].to_memory()
        is_stitched = False
    
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    x_coords, y_coords = rotate_coordinates(raw_x, raw_y, ROTATION_ANGLE)
    
    mask_dict = {"SPTDI": np.zeros_like(x_coords, dtype=bool), 
                 "DPTDI": np.zeros_like(x_coords, dtype=bool),
                 "Normal": np.zeros_like(x_coords, dtype=bool)}
                 
    split_lines_to_draw = [] # 专门记录您划下的每一刀，用于最后画 Mask 图
    
    if "Mixed" in SAMPLE_TYPE:
        fig_split, ax_split = plt.subplots(figsize=(10, 10), facecolor='black')
        pt_color = '#C8A2C8' if is_stitched else '#444444'
        ax_split.scatter(x_coords, y_coords, c=pt_color, s=3.0, edgecolors='none')
        ax_split.axis('equal'); ax_split.axis('off')
        
        if SAMPLE_TYPE == "Mixed_2":
            ax_split.set_title(f"{TEST_KEY} - Line Splitter\nClick 1 point in gap to draw line", color='white')
            plt.draw()
            print("\n👉 动作 1: 在组织缝隙处【点一下】，生成垂直切割线。")
            pt_line = plt.ginput(1, timeout=-1)
            
            if len(pt_line) == 1:
                split_x_coord = pt_line[0][0]
                split_lines_to_draw.append(split_x_coord)
                ax_split.axvline(x=split_x_coord, color='white', linestyle='--', linewidth=2)
                plt.draw()
                
                print("👉 动作 2: 在【浅二度 SPTDI】上随便【点一下】。")
                pt_sptdi = plt.ginput(1, timeout=-1)
                if len(pt_sptdi) == 1:
                    xC = pt_sptdi[0][0]
                    if xC < split_x_coord:
                        mask_dict["SPTDI"], mask_dict["DPTDI"] = x_coords < split_x_coord, x_coords >= split_x_coord
                    else:
                        mask_dict["SPTDI"], mask_dict["DPTDI"] = x_coords >= split_x_coord, x_coords < split_x_coord
                        
        elif SAMPLE_TYPE == "Mixed_3":
            ax_split.set_title(f"{TEST_KEY} - 3 Regions Splitter\nClick 2 points in gaps to draw 2 lines", color='white')
            plt.draw()
            print("\n👉 动作 1: 有三个区域，请在两个缝隙处【各点一下（共2下）】，生成2条垂直切割线。")
            pt_lines = plt.ginput(2, timeout=-1)
            
            if len(pt_lines) == 2:
                splits = sorted([pt[0] for pt in pt_lines])
                split_lines_to_draw.extend(splits)
                ax_split.axvline(x=splits[0], color='white', linestyle='--', linewidth=2)
                ax_split.axvline(x=splits[1], color='white', linestyle='--', linewidth=2)
                plt.draw()
                
                print("👉 动作 2: 划线完毕！现在请依次点击：1.【DPTDI区】 2.【SPTDI区】 3.【Normal区】")
                pt_dptdi = plt.ginput(1, timeout=-1)
                pt_sptdi = plt.ginput(1, timeout=-1)
                pt_normal = plt.ginput(1, timeout=-1)
                
                def get_zone_mask(x_click):
                    if x_click < splits[0]: return x_coords < splits[0]
                    elif x_click < splits[1]: return (x_coords >= splits[0]) & (x_coords < splits[1])
                    else: return x_coords >= splits[1]
                
                if len(pt_normal) == 1:
                    mask_dict["DPTDI"] = get_zone_mask(pt_dptdi[0][0])
                    mask_dict["SPTDI"] = get_zone_mask(pt_sptdi[0][0])
                    mask_dict["Normal"] = get_zone_mask(pt_normal[0][0])
                    
        plt.close(fig_split)
        
    elif SAMPLE_TYPE == "Pure_SPTDI": mask_dict["SPTDI"] = np.ones_like(x_coords, dtype=bool)
    elif SAMPLE_TYPE == "Pure_DPTDI": mask_dict["DPTDI"] = np.ones_like(x_coords, dtype=bool)
    elif SAMPLE_TYPE == "Pure_Normal": mask_dict["Normal"] = np.ones_like(x_coords, dtype=bool)

    spot_region_labels = np.array(['Unassigned'] * len(x_coords), dtype=object)
    for region_name, mask_bool in mask_dict.items():
        if np.any(mask_bool): spot_region_labels[mask_bool] = region_name

    # ==========================================
    # 🌟 新增：保存可视化的可重复性 Mask 图像
    # ==========================================
    fig_mask, ax_mask = plt.subplots(figsize=(10, 10), facecolor='white')
    ax_mask.scatter(x_coords, y_coords, c='lightgray', s=3.0, edgecolors='none', label='Unassigned')

    color_map = {"SPTDI": '#3a77b7', "DPTDI": '#e5192c', "Normal": '#2ca02c'}
    
    for region_name in ["SPTDI", "DPTDI", "Normal"]:
        mask_bool = mask_dict[region_name]
        if np.any(mask_bool):
            ax_mask.scatter(x_coords[mask_bool], y_coords[mask_bool], c=color_map[region_name], s=3.0, label=region_name, alpha=0.8)

    # 在白底 Mask 图上原样画出您的黑虚线切割痕迹
    data_y_min, data_y_max = y_coords.min(), y_coords.max()
    for sp_x in split_lines_to_draw:
        ax_mask.plot([sp_x, sp_x], [data_y_min, data_y_max], color='black', linestyle='--', linewidth=2.5)

    ax_mask.set_title(f"Tissue Regions Demarcation Mask\n({TEST_KEY})", fontsize=20, pad=15)
    ax_mask.axis('equal'); ax_mask.axis('off')
    ax_mask.legend(markerscale=5, loc='upper right', frameon=True, edgecolor='black')

    fig_mask.savefig(mask_png_path, dpi=600, bbox_inches='tight')
    plt.close(fig_mask)
    print(f"📷 可视化 Mask 铁证图已生成: {mask_png_path}")

    # 保存 CSV
    pd.DataFrame({"X_coord": x_coords, "Y_coord": y_coords, "Assigned_Region": spot_region_labels}).to_csv(mask_csv_path, index=False)
    print(f"✅ {TEST_KEY} 分割结果及数字矩阵已安全保存！")
    
    del adata
    gc.collect()
    plt.close('all')


# ==========================================
# [阶段二] 批量全自动渲染 FeaturePlot 与统计 (智能多区分界版)
# ==========================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scanpy as sc
import gc

plt.ioff()

data_path = r"data/external/local_611e695933ec"
output_dir = r"data/external/local_8e645d0a9b03"
stitched_dir = r"data/external/local_92743553c2f2"

# 🌟 样本配置字典 (同步更新为您的最新设置)
sample_config = {
    "5dpb":        {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_SPTDI": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_DPTDI": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45,   "cutoffs": [0, 0]},
    "19dpb":       {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225, "cutoffs": [0, 0]}, 
    "26dpb":       {"sample_id": "BW81_C02846B6_SDSDB_26dpb", "angle": 0,    "cutoffs": [0, 0]},   
    "2mph":        {"sample_id": "2mph_A03699G6.SCT",         "angle": 45,   "cutoffs": [0, 0]},
    "Normal":      {"sample_id": "NS_C02847B1",               "angle": 180,  "cutoffs": [0, 0]}  
}

focus_cell_type = "Fib_K14"
focus_markers = ["KRT14", "PDGFRA"]
TARGET_COLOR = "#2A6EBBFF" 
bg_spot_size = 1.5           
target_spot_size = 2.5       

print("\n🚀 正在启动全自动渲染与统计流水线...")
adata_full = sc.read_h5ad(data_path, backed='r')

def rotate_coordinates(x, y, angle):
    theta = np.radians(angle)
    return x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta)

for TEST_KEY, config in sample_config.items():
    TARGET_SAMPLE = config["sample_id"]
    ROTATION_ANGLE = config["angle"]
    focus_cutoffs = config["cutoffs"]
    
    mask_csv_path = os.path.join(output_dir, f"{TEST_KEY}_Spot_Region_Annotations.csv")
    if not os.path.exists(mask_csv_path):
        print(f"⚠️ 跳过 {TEST_KEY}: 未找到分区文件。")
        continue

    print(f"\n🎨 正在极速渲染: {TEST_KEY} ...")
    
    # 读取分区 CSV
    region_df = pd.read_csv(mask_csv_path)
    spot_region_labels = region_df["Assigned_Region"].values
    
    # 分流读取
    stitched_file = os.path.join(stitched_dir, f"{TARGET_SAMPLE}_Stitched.h5ad")
    if os.path.exists(stitched_file):
        adata = sc.read_h5ad(stitched_file)
    else:
        adata = adata_full[adata_full.obs['sample_batch_new'] == TARGET_SAMPLE].to_memory()
        
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    x_coords, y_coords = rotate_coordinates(raw_x, raw_y, ROTATION_ANGLE)
    
    # 基因打捞逻辑
    try:
        var_names_pool = adata.raw.var_names; use_raw = True
    except AttributeError:
        var_names_pool = adata.var_names; use_raw = False

    local_mask = np.ones((adata.n_obs,), dtype=bool)
    for g, cut in zip(focus_markers, focus_cutoffs):
        if g in var_names_pool:
            col_data = adata.raw[:, g].X if use_raw else adata[:, g].X
            if sp.issparse(col_data): col_data = col_data.toarray()
            local_mask &= (col_data.flatten() > cut)
        else:
            local_mask &= False
    is_target = local_mask

    # ================= 🌟 核心修改：多区分界线智能检测 =================
    # 获取当前样本中实际存在的所有区域
    all_regions = ["DPTDI", "SPTDI", "Normal"]
    regions_present = [r for r in all_regions if np.any(spot_region_labels == r)]
    
    split_lines = []
    if len(regions_present) >= 2:
        # 计算每个区域的 X 轴中心点，并按位置从左到右排序
        region_centers = sorted([(r, x_coords[spot_region_labels == r].mean()) for r in regions_present], key=lambda x: x[1])
        # 在相邻区域的交界处计算分界坐标
        for i in range(len(region_centers) - 1):
            left_r = region_centers[i][0]
            right_r = region_centers[i+1][0]
            # 取左区域的最大 X 和右区域的最小 X 的平均值
            line_x = (x_coords[spot_region_labels == left_r].max() + x_coords[spot_region_labels == right_r].min()) / 2
            split_lines.append(line_x)
    # =================================================================

    fig_final, ax_final = plt.subplots(figsize=(10, 10), facecolor='black')
    ax_final.scatter(x_coords[~is_target], y_coords[~is_target], c='#333333', s=bg_spot_size, edgecolors='none', alpha=0.8)
    ax_final.scatter(x_coords[is_target], y_coords[is_target], c=TARGET_COLOR, s=target_spot_size, edgecolors='none')

    data_y_min, data_y_max = y_coords.min(), y_coords.max()
    data_x_min, data_x_max = x_coords.min(), x_coords.max()
    y_range = data_y_max - data_y_min
    ax_final.set_ylim(data_y_min - y_range * 0.15, data_y_max + y_range * 0.05)
    ax_final.set_xlim(data_x_min - (data_x_max - data_x_min) * 0.05, data_x_max + (data_x_max - data_x_min) * 0.05)

    # 🌟 核心修改：绘制所有检测到的分界线并打印所有区域标签
    for line_x in split_lines:
        ax_final.plot([line_x, line_x], [data_y_min, data_y_max], color='white', linestyle='--', linewidth=2)
        
    text_y_pos = data_y_min + y_range * 0.01 
    for r_name in regions_present:
        cx = x_coords[spot_region_labels == r_name].mean()
        # 【修改点 1】：fontweight='bold' 改为 'normal'
        ax_final.text(cx, text_y_pos, r_name, color='white', fontsize=26, fontweight='normal', ha='center', va='bottom')

    # 比例尺与标题
    sb_x_start = data_x_max - (data_x_max - data_x_min) * 0.15 
    sb_y_pos = data_y_min - y_range * 0.05 
    ax_final.plot([sb_x_start, sb_x_start + (data_x_max - data_x_min) * 0.15], [sb_y_pos, sb_y_pos], color='white', linewidth=4)
    # 【修改点 2】：fontweight='bold' 改为 'normal'
    ax_final.text(sb_x_start + (data_x_max - data_x_min) * 0.075, sb_y_pos + y_range * 0.015, "1 mm", color='white', fontsize=18, fontweight='normal', ha='center', va='bottom')
    
    # 【修改点 3】：替换 dpb 为 dpi，并将 fontweight='bold' 改为 'normal'
    display_title = TEST_KEY.replace('dpb', 'dpi')
    ax_final.set_title(display_title, color='white', fontsize=32, fontweight='normal', pad=20)
    
    ax_final.set_aspect('equal'); ax_final.axis('off')

    # 注意这里保存的文件名仍然用 TEST_KEY (即包含 dpb)，以保证后面的拼接代码能顺利读取
    fig_final.savefig(os.path.join(output_dir, f"{TEST_KEY}_{focus_cell_type}_BlackStyle_FeaturePlot.png"), dpi=600, facecolor='black', bbox_inches='tight', pad_inches=0.1)
    fig_final.savefig(os.path.join(output_dir, f"{TEST_KEY}_{focus_cell_type}_BlackStyle_FeaturePlot.pdf"), dpi=600, facecolor='black', bbox_inches='tight', pad_inches=0.1)
    plt.close(fig_final)

    # 统计导出
    stats_list = []
    for region_name in all_regions:
        region_mask = (spot_region_labels == region_name)
        total_spots = np.sum(region_mask)
        if total_spots > 0:
            fib_spots = np.sum(region_mask & is_target)
            stats_list.append({
                "Timepoint": TEST_KEY, "Sample": TARGET_SAMPLE, "Region": region_name,
                "Total_Spots": total_spots, "Fib_K14_Count": fib_spots, "Fib_K14_Proportion(%)": round((fib_spots/total_spots)*100, 4)
            })
    
    if stats_list:
        pd.DataFrame(stats_list).to_csv(os.path.join(output_dir, f"{TEST_KEY}_Fib_K14_Statistics.csv"), mode='w', header=True, index=False, encoding='utf-8-sig')

    del adata
    gc.collect()

print("\n🎉 批量分析已全部结束！2mph 的 Normal 分界线已正确展示，字体已取消加粗，且图片标题中已应用 dpi 替换。")

# ==========================================
# [终极拼图器] 单列垂直无缝拼接 (附带智能去黑边收缩间距功能)
# ==========================================
import os
from PIL import Image, ImageChops

# 📁 路径配置
output_dir = r"data/external/local_8e645d0a9b03"
focus_cell_type = "Fib_K14"

# 🌟 设置拼图顺序 (保持 dpb，因为您说现在的图片名仍为 dpb，这样才能读取到文件)
plot_order = [
    "Normal", 
    "5dpb", 
    "12dpb_SPTDI", 
    "12dpb_DPTDI",  
    "19dpb",  
    "26dpb", 
    "2mph"   # 👈 确保这里是 2mph
]

# ✂️ 智能垂直裁边函数：只裁上下多余黑边，保留左右全宽
def trim_vertical(im, padding=80):
    bg = Image.new(im.mode, im.size, (0, 0, 0)) # 纯黑背景
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox() # 获取非黑像素的边界框
    
    if bbox:
        upper = max(0, bbox[1] - padding)
        lower = min(im.size[1], bbox[3] + padding)
        return im.crop((0, upper, im.size[0], lower))
    return im

print("⏳ 正在加载图像并智能裁切上下冗余黑边...")
cropped_images = []
for key in plot_order:
    img_path = os.path.join(output_dir, f"{key}_{focus_cell_type}_BlackStyle_FeaturePlot.png")
    
    if os.path.exists(img_path):
        img = Image.open(img_path)
        trimmed_img = trim_vertical(img, padding=50) 
        cropped_images.append(trimmed_img)
        print(f"   ✂️ 成功压缩: {key} (原高度: {img.size[1]}px -> 新高度: {trimmed_img.size[1]}px)")
    else:
        raise FileNotFoundError(f"❌ 找不到图片: {img_path}")

# 📐 动态计算最终画板尺寸
canvas_width = cropped_images[0].size[0]
canvas_height = sum(img.size[1] for img in cropped_images)

print(f"\n🖼️ 正在生成紧凑版 {canvas_width} x {canvas_height} 的垂直单列长图 ...")
canvas = Image.new('RGB', (canvas_width, canvas_height), color='black')

# 🧩 紧凑排布粘贴
y_offset = 0
for idx, img in enumerate(cropped_images):
    canvas.paste(img, (0, y_offset))
    y_offset += img.size[1] 
    
# 💾 保存最终拼图 (🌟 将最终输出的文件名加上 _dpi 标识)
final_png_path = os.path.join(output_dir, f"Combined_AllSamples_{focus_cell_type}_Vertical_Compact_dpi.png")
final_pdf_path = os.path.join(output_dir, f"Combined_AllSamples_{focus_cell_type}_Vertical_Compact_dpi.pdf")

# 1. 保存高清 PNG
canvas.save(final_png_path, format="PNG", optimize=True)
print(f"✅ 高清 PNG 已保存: {final_png_path}")

# 2. 保存为 PDF
if canvas.mode == 'RGBA':
    pdf_canvas = canvas.convert('RGB')
else:
    pdf_canvas = canvas

pdf_canvas.save(final_pdf_path, format="PDF", resolution=600.0)
print(f"✅ 矢量 PDF 已保存: {final_pdf_path}")

print("\n🎉 所有格式拼接大功告成！")



###########################细胞比例差异
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import imageio.v3 as iio
from matplotlib.patches import Patch

# ==========================================
# 1. 环境配置与路径
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
h5ad_19d_path = os.path.join(output_dir, "BW15D_C6_SDSDB_19dpb_Annotated.h5ad")
mask_path = os.path.join(output_dir, "wound_zones_mask.png")

raw_time_points = ["5dpb", "12dpb_SPTDI", "12dpb_DPTDI", "19dpb", "26dpb", "2mph"]
final_time_points = ["5dpb", "12dpb", "19dpb", "26dpb", "2mph"]
zone_labels = ['Un-epithelialized', 'Epithelializing', 'Newly Epithelialized']

# ==========================================
# 2. 数据准备：子图 A (时间轴对比)
# ==========================================
all_dfs = []
for tp in raw_time_points:
    csv_path = os.path.join(output_dir, f"{tp}_Fib_K14_Statistics.csv")
    if os.path.exists(csv_path):
        df_t = pd.read_csv(csv_path)
        if "12dpb" in tp: df_t['Timepoint'] = "12dpb"
        all_dfs.append(df_t)

if not all_dfs:
    raise ValueError(f"❌ 路径下未找到 CSV 文件，请确认：{output_dir}")

df_all = pd.concat(all_dfs)
time_plot_data = {"Time": final_time_points, "SPTDI": [], "DPTDI": []}
for tp in final_time_points:
    subset = df_all[df_all['Timepoint'] == tp]
    s_val = subset.loc[subset['Region'] == 'SPTDI', 'Fib_K14_Proportion(%)'].mean()
    d_val = subset.loc[subset['Region'] == 'DPTDI', 'Fib_K14_Proportion(%)'].mean()
    time_plot_data["SPTDI"].append(s_val if pd.notna(s_val) else 0)
    time_plot_data["DPTDI"].append(d_val if pd.notna(d_val) else 0)

# ==========================================
# 3. 数据准备：子图 B (19dpb 比例修正版)
# ==========================================
print(f"⏳ 正在读取 19dpb 数据并修正比例计算逻辑...")
adata_19d = sc.read_h5ad(h5ad_19d_path)
mask_img = iio.imread(mask_path)

rotation_angle = -225
raw_x = adata_19d.obsm['spatial'][:, 0] if 'spatial' in adata_19d.obsm else adata_19d.obsm['X_spatial'][:, 0]
raw_y = adata_19d.obsm['spatial'][:, 1] if 'spatial' in adata_19d.obsm else adata_19d.obsm['X_spatial'][:, 1]
theta = np.radians(rotation_angle)
x_c = raw_x * np.cos(theta) - raw_y * np.sin(theta)
y_c = raw_x * np.sin(theta) + raw_y * np.cos(theta)

m_h, m_w = mask_img.shape
px = np.clip(((x_c - x_c.min()) / (x_c.max() - x_c.min()) * (m_w - 1)).astype(int), 0, m_w - 1)
py = np.clip(((y_c - y_c.min()) / (y_c.max() - y_c.min()) * (m_h - 1)).astype(int), 0, m_h - 1)

adata_19d.obs['Wound_Zone'] = mask_img[py, px]
zone_mapping = {1: 'Un-epithelialized', 2: 'Epithelializing', 3: 'Newly Epithelialized'}

df_all_zones = adata_19d.obs[adata_19d.obs['Wound_Zone'].isin([1, 2, 3])].copy()
df_all_zones['Zone_Name'] = df_all_zones['Wound_Zone'].map(zone_mapping)

total_counts = df_all_zones['Zone_Name'].value_counts()
fib_counts = df_all_zones[df_all_zones['sub_labels'] == "Fib_K14"]['Zone_Name'].value_counts()

zone_counts = (fib_counts / total_counts * 100).reindex(zone_labels).fillna(0)
zone_vals = zone_counts.values

# ==========================================
# 4. 绘图：2:1长宽比、全局大字体、文字旋转与图例靠右
# ==========================================
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

fig = plt.figure(figsize=(18, 6), facecolor='black')

# 🌟 修改点 1: 调整 width_ratios 为 5:2，进一步缩小右边图的宽度
# 稍微加大了 wspace，防止两张图贴得太近
gs = fig.add_gridspec(1, 2, wspace=0.3, width_ratios=[5, 2])

colors = {'DPTDI': '#e5192c', 'SPTDI': '#E1C62F', 'Un-epi': '#e5192c', 'Epi-Front': '#E1C62F', 'Newly-epi': '#3a77b7'}
BAR_WIDTH = 0.35

TITLE_SIZE = 32   
LABEL_SIZE = 28   
TICK_SIZE = 24    
ANNOT_SIZE = 24   
LEGEND_SIZE = 28  

def format_ax(ax, title):
    ax.set_facecolor('black')
    ax.set_title(title, color='white', fontsize=TITLE_SIZE, pad=30)
    ax.tick_params(axis='both', colors='white', labelsize=TICK_SIZE)
    ax.set_ylabel('Proportion (%)', color='white', fontsize=LABEL_SIZE)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('white')
        ax.spines[s].set_linewidth(2.0)
    for s in ['top', 'right']: ax.spines[s].set_visible(False)

# --- 子图 A ---
ax1 = fig.add_subplot(gs[0])
x1 = np.arange(len(final_time_points))
r1 = ax1.bar(x1 - BAR_WIDTH/2, time_plot_data["SPTDI"], BAR_WIDTH, color=colors['SPTDI'], alpha=0.9)
r2 = ax1.bar(x1 + BAR_WIDTH/2, time_plot_data["DPTDI"], BAR_WIDTH, color=colors['DPTDI'], alpha=0.9, edgecolor='white', linewidth=0.5)

display_time_points = [tp.replace('dpb', 'dpi') for tp in final_time_points]
ax1.set_xticks(x1)
ax1.set_xticklabels(display_time_points, color='white', rotation=45, ha='right')

format_ax(ax1, 'Temporal Dynamics of Fib_K14')
ax1.set_ylim(0, max(max(time_plot_data["SPTDI"]), max(time_plot_data["DPTDI"])) * 1.3) 

# 独立图例 A
legend_A = [Patch(facecolor=colors['SPTDI'], label='SPTDI'), 
            Patch(facecolor=colors['DPTDI'], label='DPTDI')]
leg1 = ax1.legend(handles=legend_A, loc='upper right', fontsize=LEGEND_SIZE, frameon=True, facecolor='#111111', edgecolor='white')
for text in leg1.get_texts(): text.set_color("white")

# --- 子图 B ---
ax2 = fig.add_subplot(gs[1])
x2 = np.arange(len(zone_labels))
r3 = ax2.bar(x2, zone_vals, width=BAR_WIDTH, color=[colors['Un-epi'], colors['Epi-Front'], colors['Newly-epi']], alpha=0.9, edgecolor='white', linewidth=0.5)

ax2.set_xticks(x2)
ax2.set_xticklabels(zone_labels, color='white', rotation=45, ha='right')

# 🌟 修改点 2: 标题中加入 \n 进行强制换行
format_ax(ax2, 'Fib_K14 Across\nWound Zones (19dpi)')
ax2.set_ylim(0, max(zone_vals) * 1.3)

# 🌟 修改点 3: 已彻底删除右边图 (ax2) 的图例生成代码

# ==========================================
# 5. 数值标注与输出
# ==========================================
def add_labels(ax, rects):
    for r in rects:
        h = r.get_height()
        if h > 0: ax.annotate(f'{h:.2f}', xy=(r.get_x() + r.get_width()/2, h), xytext=(0, 6), textcoords="offset points", ha='center', va='bottom', color='white', fontsize=ANNOT_SIZE)

add_labels(ax1, r1); add_labels(ax1, r2); add_labels(ax2, r3)

plt.tight_layout()

# 1. 保存为高清 PNG
save_path_png = os.path.join(output_dir, "Fib_K14_Consolidated_Row_Final_Aligned_2to1.png")
fig.savefig(save_path_png, dpi=600, facecolor='black', bbox_inches='tight')

# 2. 保存为矢量 PDF
save_path_pdf = os.path.join(output_dir, "Fib_K14_Consolidated_Row_Final_Aligned_2to1.pdf")
fig.savefig(save_path_pdf, facecolor='black', bbox_inches='tight')

plt.show()

print(f"✅ 更新完成！\n1. 右图已按 5:2 比例缩小宽度\n2. 右图标题已分为两行\n3. 右图图例已彻底删除\n已保存至：{save_path_png}")

#############################完整共定位分析流程
import os
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
# 1. 环境与参数配置
# ==========================================
# 设置黑底白字的绘图样式
plt.style.use('dark_background')

# 设置工作目录并加载数据
file_path = r"data/external/local_20a7f572f234"
output_dir = os.path.dirname(file_path)

print(f"⏳ 正在读取数据: {os.path.basename(file_path)} ...")
adata = sc.read_h5ad(file_path)

# ==========================================
# 2. 数据检查与预处理
# ==========================================
if 'sub_labels' in adata.obs:
    adata.obs['sub_labels'] = adata.obs['sub_labels'].astype('category')
else:
    raise ValueError("❌ 数据中未找到 'sub_labels' 列。")

if 'spatial' not in adata.obsm:
    if 'X_spatial' in adata.obsm: 
        adata.obsm['spatial'] = adata.obsm['X_spatial']
    else:
        raise ValueError("❌ 未在 obsm 中找到空间坐标 ('spatial' 或 'X_spatial')。")

# ==========================================
# 3. 空间共定位分析 (Squidpy)
# ==========================================
print("🌐 正在构建空间邻居图并计算邻域富集...")
sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="sub_labels")

# ==========================================
# 4. 绘制并保存邻域富集热图 (修复版)
# ==========================================
print("🎨 正在绘制黑底白字富集热图...")

# 🌟 关键修复：不要提前创建 fig 和 ax。让 Squidpy 自己生成完整的复杂画布，这样才能自带侧边分类颜色条。
# 绘制复杂热图（Squidpy 会自动接管布局）
sq.pl.nhood_enrichment(
    adata, 
    cluster_key="sub_labels", 
    figsize=(12, 10),  
    title="Neighborhood Enrichment (sub_labels)",
    cmap="viridis",    
)

fig = plt.gcf()
fig.patch.set_facecolor('black')

# 🌟 遍历该画布上产生的所有子图，统一套用黑底白字主题
for current_ax in fig.axes:
    current_ax.set_facecolor('black')
    current_ax.tick_params(colors='white')
    
    # 将所有的边框涂白
    for spine in current_ax.spines.values():
        spine.set_color('white')
    
    # 字体涂白
    current_ax.xaxis.label.set_color('white')
    current_ax.yaxis.label.set_color('white')
    if current_ax.title:
        current_ax.title.set_color('white')

# 专门提取主热图旋转 X 轴标签，防止细胞名称挤在一起
for current_ax in fig.axes:
    if len(current_ax.get_xticklabels()) > 0:
        # 提示：旋转 90 度时，ha="right" 结合垂直居中通常对齐效果最好
        plt.setp(current_ax.get_xticklabels(), rotation=90, ha="right", va="center")

# ❌ 删除了 plt.tight_layout()，避免与 clustermap 布局引擎冲突

# 直接保存即可，bbox_inches='tight' 会负责完美裁边
save_png = os.path.join(output_dir, "Spatial_Neighborhood_Enrichment.png")
save_pdf = os.path.join(output_dir, "Spatial_Neighborhood_Enrichment.pdf")

fig.savefig(save_png, dpi=600, facecolor='black', bbox_inches='tight')
fig.savefig(save_pdf, facecolor='black', bbox_inches='tight')

plt.show()
print(f"✅ 热图已成功保存至:\nPNG: {save_png}\nPDF: {save_pdf}")

# ==========================================
# 5. 提取并打印 Fib_K14 共定位数据
# ==========================================
print("\n📊 正在提取 Fib_K14 共定位得分...")
enrichment_zscore = adata.uns['sub_labels_nhood_enrichment']['zscore']
cluster_names = adata.obs['sub_labels'].cat.categories

df_enrichment = pd.DataFrame(enrichment_zscore, index=cluster_names, columns=cluster_names)

csv_path = os.path.join(output_dir, "Neighborhood_Enrichment_Zscores.csv")
df_enrichment.to_csv(csv_path)
print(f"✅ 完整 Z-score 矩阵已保存至: {csv_path}")

target_cell = "Fib_K14"
if target_cell in df_enrichment.index:
    fib_k14_interactions = df_enrichment.loc[target_cell].sort_values(ascending=False)
    
    top_colocalized = fib_k14_interactions[fib_k14_interactions.index != target_cell].head(5)
    
    print(f"\n🏆 与 {target_cell} 共定位最显著的前 5 种细胞 (Z-score):")
    print("-" * 40)
    for cell_type, score in top_colocalized.items():
        print(f"{cell_type:<20}: {score:.2f}")
    print("-" * 40)
else:
    print(f"⚠️ 警告: '{target_cell}' 未在 'sub_labels' 中找到。")
    

######################################Fib_K14共线性分析
import os
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ==========================================
# 0. 环境与配置
# ==========================================
plt.style.use('dark_background')

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
output_dir = r"data/external/local_8e645d0a9b03"

# ==========================================
# 1. 重排分类与提取 Top 10
# ==========================================
current_cats = adata.obs['sub_labels'].unique().dropna()
ordered_cats = [k for k in raw_annotation_config.keys() if k in current_cats]
missing_cats = [k for k in current_cats if k not in ordered_cats]
final_order = ordered_cats + missing_cats

adata.obs['sub_labels'] = adata.obs['sub_labels'].astype('category').cat.reorder_categories(final_order)
hex_colors = [raw_annotation_config.get(c, {}).get('color', '#888888') for c in final_order]
adata.uns['sub_labels_colors'] = hex_colors

print("⏳ 正在计算空间共现性 (Co-occurrence)...")
sq.gr.co_occurrence(adata, cluster_key="sub_labels")

print(f"📊 正在锁定与 {target_cluster} 空间共现最显著的 Top 10 细胞...")
co_out = adata.uns["sub_labels_co_occurrence"]["occ"]
target_idx = list(final_order).index(target_cluster)
closest_distance_vals = co_out[target_idx, :, 0]

co_occ_df = pd.DataFrame({'Cell_Type': final_order, 'Co_occurrence_Score': closest_distance_vals})
co_occ_df = co_occ_df[co_occ_df['Cell_Type'] != target_cluster]
co_occ_df = co_occ_df.sort_values(by='Co_occurrence_Score', ascending=False)

top10_df = co_occ_df.head(10)
top5_clusters = top10_df.head(10)['Cell_Type'].tolist()

# ==========================================
# 2. 劫持画板并创建上下“分层拼图”结构
# ==========================================
print(f"🎨 正在绘制双栏视图及独立图例区...")

# 让 Squidpy 自己画出原始折线图
sq.pl.co_occurrence(adata, cluster_key="sub_labels", clusters=target_cluster)

fig = plt.gcf()
ax1 = plt.gca()

# 🌟 核心拼图架构：设置高度稍高一点 (16x11)，容纳底部图例
fig.set_size_inches(12, 5)
fig.patch.set_facecolor('black')

# 划分网格：2 行 2 列。第一行高(占7份)画图，第二行矮(占3份)画图例
# hspace=0.2 确保上下两层之间有物理隔离的空白带
gs = fig.add_gridspec(2, 2, height_ratios=[8, 2], wspace=0.03, hspace=0.02)

# 将原折线图锁定在 [第1行, 第1列]
ax1.set_subplotspec(gs[0, 0])
ax1.set_facecolor('black')
ax1.set_title(f"Co-occurrence probability conditioned on {target_cluster}", color='white', pad=20, fontsize=16)

# 渲染高亮线条
for line in ax1.lines:
    label = line.get_label()
    if label == target_cluster:
        line.set_linewidth(5); line.set_zorder(15)
    elif label in top5_clusters:
        line.set_linewidth(3.5); line.set_alpha(1.0); line.set_zorder(10)
    else:
        line.set_linewidth(1.0); line.set_alpha(0.15); line.set_zorder(1)

# 提取并移除折线图自带的图例
handles, labels = [], []
legend = ax1.get_legend()
if legend:
    handles = legend.legend_handles
    labels = [t.get_text() for t in legend.get_texts()]
    legend.remove()

# ==========================================
# 3. 绘制右侧 Top10 柱状图 [第1行, 第2列]
# ==========================================
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_facecolor('black')

plot_df = top10_df.iloc[::-1]
bar_colors = [raw_annotation_config.get(c, {}).get('color', '#888888') for c in plot_df['Cell_Type']]

bars = ax2.barh(plot_df['Cell_Type'], plot_df['Co_occurrence_Score'], color=bar_colors, edgecolor='none')

ax2.set_title(f"Top 10 Co-occurring Cells (Score)", color='white', pad=20, fontsize=16)
ax2.set_xlabel("Co-occurrence Score at closest distance", color='white', fontsize=16)

for bar in bars:
    width = bar.get_width()
    ax2.text(width + 0.02, bar.get_y() + bar.get_height()/2, f'{width:.3f}', 
             va='center', ha='left', color='white', fontsize=12, fontweight='bold')

ax2.set_xlim(0, plot_df['Co_occurrence_Score'].max() * 1.15)

for ax in [ax1, ax2]:
    ax.tick_params(colors='white', labelsize=12)
    for spine in ax.spines.values():
        spine.set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')

# ==========================================
# 4. 在底部的独立画板绘制图例 [第2行, 横跨两列]
# ==========================================
# 🌟 创建一个贯穿底部两列的专属透明坐标系
ax_leg = fig.add_subplot(gs[1, :])
ax_leg.axis('off') # 隐藏边框和刻度，让它变成一块纯黑的幕布

# 在这块专属幕布的正中心 (loc='center') 绘制图例
# 改用 ncol=6 增加一点列间距，防止文字太挤
new_legend = ax_leg.legend(
    handles, labels,
    loc='center', 
    ncol=6, 
    frameon=False, 
    fontsize=11
)

# 字体联动高亮逻辑
for text in new_legend.get_texts():
    if text.get_text() in top5_clusters or text.get_text() == target_cluster:
        text.set_color("white")
        text.set_weight("bold") 
    else:
        text.set_color("#888888")

# ==========================================
# 5. 保存
# ==========================================
save_png = os.path.join(output_dir, f"Spatial_Co_occurrence_Panel_{target_cluster}.png")
save_pdf = os.path.join(output_dir, f"Spatial_Co_occurrence_Panel_{target_cluster}.pdf")

fig.savefig(save_png, dpi=600, facecolor='black', bbox_inches='tight')
fig.savefig(save_pdf, facecolor='black', bbox_inches='tight')

plt.show()
print(f"✅ 图表已更新！(图例已置于独立画布区，永不重叠)\n保存至:\n{save_png}")


#################################空间分布规律动态分析
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
    "5dpb":        {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_SPTDI": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_DPTDI": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45,   "cutoffs": [0, 0]},
    "19dpb":       {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225, "cutoffs": [0, 0]}, 
    "26dpb":       {"sample_id": "BW81_C02846B6_SDSDB_26dpb", "angle": 0,    "cutoffs": [0, 0]},   
    "2mph":        {"sample_id": "2mph_A03699G6.SCT",         "angle": 45,   "cutoffs": [0, 0]},
    "Normal":      {"sample_id": "NS_C02847B1",               "angle": 180,  "cutoffs": [0, 0]}  
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
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
    
    try:
        sc.pl.embedding(adata, basis='rot_spatial', color='sub_labels', ax=ax, show=False)
    except Exception:
        ax.scatter(rot_x, rot_y, s=2, c='#888888', alpha=0.5)

    ax.axis('on') 
    
    print("\n" + "="*50)
    print(f"🖱️ 【结构锚点提取开启】正在处理: {sample_id} (已旋转 {angle} 度)")
    print("左键=选点，右键=撤销上一个点，回车(Enter)=完成当前线段。")
    print("="*50 + "\n")

    # ----------------------------------------
    # 线段 1: 表皮底层
    # ----------------------------------------
    ax.set_title("[Step 1/3] Click Epidermis Baseline -> Press ENTER", color='cyan', fontsize=14, fontweight='bold')
    plt.draw()
    print("👉 [1/3] 请描绘【表皮底层/基底膜】走势，完成后按 Enter。")
    epi_coords = plt.ginput(n=-1, timeout=0)
    if epi_coords:
        ex, ey = zip(*epi_coords)
        ax.plot(ex, ey, color='cyan', linewidth=2, marker='o', markersize=4)
        plt.draw()
    print(f"   ✔️ Step 1 完成，共采集 {len(epi_coords)} 个点。")
    plt.pause(0.2) # 缓冲延迟，防止回车键连击穿透

    # ----------------------------------------
    # 线段 2: 真皮底层
    # ----------------------------------------
    ax.set_title("[Step 2/3] Click Dermis Bottom -> Press ENTER", color='magenta', fontsize=14, fontweight='bold')
    plt.draw()
    print("\n👉 [2/3] 请描绘【真皮最底层】走势，完成后按 Enter。")
    der_coords = plt.ginput(n=-1, timeout=0)
    if der_coords:
        dx, dy = zip(*der_coords)
        ax.plot(dx, dy, color='magenta', linewidth=2, marker='o', markersize=4)
        plt.draw()
    print(f"   ✔️ Step 2 完成，共采集 {len(der_coords)} 个点。")
    plt.pause(0.2) # 缓冲延迟，防止回车键连击穿透

    # ----------------------------------------
    # 线段 3: 愈合前沿
    # ----------------------------------------
    ax.set_title("[Step 3/3] Click Leading Edge (2 points) -> Press ENTER", color='yellow', fontsize=14, fontweight='bold')
    plt.draw()
    print("\n👉 [3/3] 请点击【两个点】确定愈合前沿分界线（若已完全愈合请直接按 Enter 跳过）。")
    front_coords = plt.ginput(n=-1, timeout=0)
    if front_coords:
        fx, fy = zip(*front_coords)
        ax.plot(fx, fy, color='yellow', linewidth=2, linestyle='--')
        plt.draw()
    print(f"   ✔️ Step 3 完成，共采集 {len(front_coords)} 个点。")
    plt.pause(0.2)

    # ----------------------------------------
    # 终极防误关确认 (图形仍保持打开状态)
    # ----------------------------------------
    ax.set_title("✅ All Steps Done! Please check the terminal to confirm & save.", color='lime', fontsize=14, fontweight='bold')
    plt.draw()
    
    print("\n" + "="*50)
    print("⏸️ 【等待确认】图片目前保持打开状态。")
    print("请对照屏幕上的图像（查看黄线位置），在下方终端输入参数：")
    
    if len(front_coords) >= 2:
        direction_input = input("❓ 以愈合前沿(黄线)为界，哪边是【已长好】的上皮区？(输入 R 代表右侧，L 代表左侧): ").strip().upper()
        direction = "right_is_healed" if direction_input == 'R' else "left_is_healed"
    else:
        print("ℹ️ 未检测到黄线，判定为完全愈合样本。")
        input("按任意键并回车继续...")
        direction = "right_is_healed" # 占位符，第二步会自动处理全愈合情况
    
    print("="*50)
    
    # 只有当终端交互彻底完成后，才关闭图形
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
        
##############################调整JSON文件
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import json
import warnings

warnings.filterwarnings("ignore")

# ==========================================
# 0. 路径配置 (保持与您之前完全一致)
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
    
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 6), dpi=150)
    
    # 绘制细胞散点
    ax.scatter(rot_x, rot_y, s=2, c='#888888', alpha=0.5)
    ax.axis('on') 
    
    # 🌟 自动绘制您已经画好的 Step 1 (表皮) 和 Step 2 (真皮) 🌟
    ex, ey = zip(*anchor_dict['epi_baseline'])
    ax.plot(ex, ey, color='cyan', linewidth=2, marker='o', markersize=4, label="Epidermis (Kept)")
    
    dx, dy = zip(*anchor_dict['der_bottom'])
    ax.plot(dx, dy, color='magenta', linewidth=2, marker='o', markersize=4, label="Dermis (Kept)")
    
    # 🌟 仅要求重新输入 Step 3 🌟
    ax.set_title("Target Fix: Click 2 points to set NEW Leading Edge -> Press ENTER", color='yellow', fontsize=14, fontweight='bold')
    plt.draw()
    
    print("\n" + "="*50)
    print("👉 【微调模式】表皮和真皮线已自动保留！")
    print("🎯 请在图中【细胞聚集度最高】的中心位置，点击两个点画一条垂直/倾斜的分界线（新的愈合前沿），然后按 Enter。")
    print("="*50)
    
    front_coords = plt.ginput(n=-1, timeout=0)
    if front_coords:
        fx, fy = zip(*front_coords)
        ax.plot(fx, fy, color='yellow', linewidth=2, linestyle='--')
        plt.draw()
        
    plt.pause(0.2)
    
    # 确认方向
    print("\n❓ 请最后确认一次：以您新画的黄线为界，哪边是【已长好(Epithelialized)】的区域？")
    direction_input = input("输入 R (代表黄线右侧) 或 L (代表黄线左侧): ").strip().upper()
    direction = "right_is_healed" if direction_input == 'R' else "left_is_healed"
    
    plt.close(fig)
    
    # 更新 JSON (只替换 front_coords 和 direction，其他原封不动)
    anchor_dict['leading_edge'] = [(float(x), float(y)) for x, y in front_coords]
    anchor_dict['healed_direction'] = direction
    
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(anchor_dict, f, indent=4)
        
    print(f"\n✅ 修改大功告成！JSON 锚点已更新至: {json_path}")
    print("🎉 现在您可以直接去运行【阶段二绘图脚本】了，Fib_K14 绝对会完美移动到 X=0 的 Epi-Front 中心位置！")


###############################绘图
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

warnings.filterwarnings("ignore")

# ==========================================
# 0. 路径与核心参数配置
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"

sample_config = {
    "5dpb":        {"sample_id": "BW32_A01597A3_SDSDB_5dpb",  "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_SPTDI": {"sample_id": "BW14_1_C1_SSDB_12dpb",      "angle": 135,  "cutoffs": [0, 0]},
    "12dpb_DPTDI": {"sample_id": "BW13_1_B1_DSDB_12dpb",      "angle": 45,   "cutoffs": [0, 0]},
    "19dpb":       {"sample_id": "BW15D_C6_SDSDB_19dpb",      "angle": -225, "cutoffs": [0, 0]}, 
    "26dpb":       {"sample_id": "BW81_C02846B6_SDSDB_26dpb", "angle": 0,    "cutoffs": [0, 0]},   
    "2mph":        {"sample_id": "2mph_A03699G6.SCT",         "angle": 45,   "cutoffs": [0, 0]},
    "Normal":      {"sample_id": "NS_C02847B1",               "angle": 180,  "cutoffs": [0, 0]}  
}

# ⚠️ 在这里修改您当前想要出图的样本 Key
TEST_KEY = "19dpb"  
target_cell = "Fib_K14"

TARGET_SAMPLE = sample_config[TEST_KEY]["sample_id"]
ROTATION_ANGLE = sample_config[TEST_KEY]["angle"]

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
    
    # 🌟 修复点：检查是否画了愈合前沿线 (至少需要2个点)
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
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator  # 🌟 新增导入，用于控制刻度数量

def plot_planA_dynamics(df, sample_id, cell_type, out_dir):
    print(f"\n🎨 正在绘制 {cell_type} 的顶刊排版连续空间动态图 (极简对角线版)...")
    sns.set_style("ticks")
    
    # (3) 进一步缩小画布，完美适配 1行5图 (宽2.6, 高2.2)
    fig, ax = plt.subplots(figsize=(2.4, 2.2), dpi=600, layout="constrained")
    df_sub = df[df['sub_labels'] == cell_type].copy()
    
    if df_sub.empty:
        print(f"⚠️ 未找到 {cell_type} 细胞！")
        return
        
    ax.set_facecolor('black')
    fig.patch.set_facecolor('black')

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

    ax.plot([b1_top_x, b1_bot_x], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.plot([b2_top_x, b2_bot_x], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)

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

    # 顶部分区文字 (字号缩小适应小图)
    trans = ax.get_xaxis_transform()
    center_un_epi = (-max_abs_x + b1_top_x) / 2
    center_epi_front = (b1_top_x + b2_top_x) / 2
    center_newly_epi = (b2_top_x + max_abs_x) / 2
    
    ax.text(center_un_epi, 1.02, 'Un-epi', color='#e5192c', transform=trans, ha='center', va='bottom', fontsize=8)
    ax.text(center_epi_front, 1.02, 'Epi-Front', color='#E1C62F', transform=trans, ha='center', va='bottom', fontsize=8)
    ax.text(center_newly_epi, 1.02, 'Newly-epi', color='#3a77b7', transform=trans, ha='center', va='bottom', fontsize=8)

    # ==========================================
    # (3) 核心数据渲染：KDE等高线 + 密度散点
    # ==========================================
    sns.kdeplot(
        data=df_sub, x='Distance_to_Front_X', y='Relative_Depth_Y',
        ax=ax, color="#FFFFFF", linewidth=0.5, alpha=0.6, levels=6, cut=0, zorder=2
    )
    
    x, y = df_sub['Distance_to_Front_X'], df_sub['Relative_Depth_Y']
    xy = np.vstack([x, y])
    try:
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z = x.iloc[idx], y.iloc[idx], z[idx]
        # 散点进一步缩小为 s=1.5
        ax.scatter(x, y, c=z, cmap='plasma', edgecolor='none', s=1.5, zorder=3)
    except np.linalg.LinAlgError:
        ax.scatter(x, y, color='#FFA500', edgecolor='none', s=1.5, zorder=3)
    
    # ==========================================
    # (4) 🌟 顺应等高线走势的对角线箭头 🌟
    # ==========================================
    # 箭头起点 (密度最高的核心区坐标)
    center_x = -500 
    center_y = 0.75 

    # 箭头 1：右上 (表皮方向，促进上皮化)
    ax.annotate('', 
                xy=(1800, 0.2), xycoords='data',       # 箭头尖端
                xytext=(center_x, center_y), textcoords='data', # 箭头尾巴
                arrowprops=dict(facecolor="#0099CCFF", edgecolor="#0099CCFF", width=0.3, headwidth=2.5, shrink=0.08),
                zorder=4)

    # 箭头 2：左下 (深层真皮方向，普通纤维化)
    ax.annotate('', 
                xy=(-3200, 0.98), xycoords='data',      # 箭头尖端
                xytext=(center_x, center_y), textcoords='data', # 箭头尾巴
                arrowprops=dict(facecolor="#FFFF00", edgecolor="#FFFF00", width=0.3, headwidth=2.5, shrink=0.08),
                zorder=4)

    # ==========================================
    # (5) 🌟 坐标轴刻度与排版控制 (<=4个刻度) 🌟
    # ==========================================
    ax.set_xlim(-max_abs_x, max_abs_x) 
    
    # X 轴智能限制最多 4 个刻度 (nbins=3 通常会生成 3-4 个主刻度，如 -4000, 0, 4000)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
    
    # Y 轴硬编码为极简的 3 个刻度
    ax.set_ylim(1.05, -0.05) 
    ax.set_yticks([0.0, 0.5, 1.0])
    
    display_title = f"{sample_id}".replace('dpb', 'dpi')
    ax.set_title(display_title, fontsize=14, color='white', pad=12, fontweight='normal') 
    
    # 轴标签缩小
    ax.set_xlabel("Distance to leading edge (µm)", color='white', fontsize=8, labelpad=4)
    ax.set_ylabel("Epidermis-Dermis axis", color='white', fontsize=14, labelpad=4)
    
    # 刻度数字缩小
    ax.tick_params(colors='white', direction='out', length=2.5, width=0.6, labelsize=14)
    for spine in ax.spines.values():
        spine.set_color('white')
        spine.set_linewidth(0.6)

    # ==========================================
    # (6) 保存文件
    # ==========================================
    save_name = f"Nature_SpatialMap_{sample_id}_{cell_type}.pdf"
    save_path = os.path.join(out_dir, save_name)
    save_path_png = save_path.replace('.pdf', '.png')
    
    fig.savefig(save_path, facecolor='black', bbox_inches='tight')
    fig.savefig(save_path_png, facecolor='black', bbox_inches='tight', dpi=600) 
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
        
        # 传入 rotation_angle 确保跟之前画线时用的是相同的旋转系
        adata = calculate_planA_coordinates(adata, anchor_dict, angle=ROTATION_ANGLE)
        
        df_plot = pd.DataFrame(adata.obs[['sub_labels', 'Distance_to_Front_X', 'Relative_Depth_Y']])
        plot_planA_dynamics(df_plot, sample_id=TEST_KEY, cell_type=target_cell, out_dir=output_dir)
        
        print("\n🎉 分析完成！")

####################空间特征基因计算
import os
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse
import scipy.stats
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as clr
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
# 1. 全局排版配置：黑底白字 + Arial 字体
# ==========================================
plt.style.use('dark_background')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'white'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['figure.dpi'] = 150 

# ==========================================
# 2. 颜色、路径与核心基因配置
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

# 🌟 终极精简版标签列表：严格限制在 10 个以内，均等代表不同空间分布 🌟
highlight_genes = [
    'KLF4', 'IGFL1', 'KRT10', 'MMP1', 'S100A8', 
    'CXCL14', 'FTH1', 'SFN', 'B2M', 'SPRR1B'
]

# ==========================================
# 3. 核心计算工具函数
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
    corr_x = df_correlation(adata_filtered, var_x)
    corr_y = df_correlation(adata_filtered, var_y)
    sig_x = corr_x[(corr_x['p'] < p_threshold) & (np.abs(corr_x['rho']) > rho_threshold)]
    sig_y = corr_y[(corr_y['p'] < p_threshold) & (np.abs(corr_y['rho']) > rho_threshold)]
    print(f"✅ 扫描完成：X轴显著基因 {len(sig_x)} 个，Y轴显著基因 {len(sig_y)} 个。")
    return sig_x, sig_y

# ==========================================
# 4. 微缩版蜂群图 (极简去噪版)
# ==========================================
def plot_nature_swarm_mini(adata, vars, highlight_list, out_dir, rho_cutoff=0.15):
    print("\n⏳ 正在绘制微缩版 19dpi 蜂群图...")
    adata_filtered = filter_adata_expressed_in_n_cells(adata, percent=0.05)
    df = pd.concat([df_correlation(adata_filtered, var) for var in vars], ignore_index=True)
    
    df["color"] = df["rho"].apply(
        lambda x: "positive" if x > rho_cutoff else "negative" if x < -rho_cutoff else "none"
    ).astype("category").cat.reorder_categories(["positive", "none", "negative"], ordered=True)

    fig, ax = plt.subplots(figsize=(1.8, 2.6), dpi=300)
    ax.axhline(y=rho_cutoff, linestyle=":", color="#AAAAAA", linewidth=0.5)
    ax.axhline(y=-rho_cutoff, linestyle=":", color="#AAAAAA", linewidth=0.5)
    ax.axhline(y=0, linestyle="-", color="#DDDDDD", linewidth=0.5)

    sns.swarmplot(data=df, x="var", y="rho", ax=ax, order=vars, hue="color", palette=["#E15759", "#E8E8E8", "#4E79A7"], s=1.0, zorder=1)
    if ax.legend_: ax.legend_.remove()

    x_offset, y_gap = 0.2, 0.05 # 因为标签变少，稍微放宽一点防重叠间距
    for i, var in enumerate(vars):
        hl_df = df[(df["var"] == var) & (df["gene"].isin(highlight_list))].sort_values(by="rho", ascending=False)
        last_y = float('inf')
        for _, row in hl_df.iterrows():
            text_y = row["rho"] if (last_y - row["rho"]) >= y_gap else last_y - y_gap
            last_y = text_y
            
            color = "#E15759" if row["rho"] > 0 else "#4E79A7"
            ax.scatter(i, row["rho"], color=color, s=3, zorder=2) 
            ax.plot([i + 0.03, i + 0.12, i + x_offset], [row["rho"], row["rho"], text_y], color="gray", linewidth=0.4)
            # 🌟 移除了验证基因的高亮，全部使用标准的纯白斜体 🌟
            ax.text(i + x_offset + 0.02, text_y, f"$\it{{{row['gene']}}}$", fontsize=5, color="white", fontweight="normal", va="center") 

    ax.set_title("19dpi | Spatial Correl.", fontsize=8, pad=8, fontweight='bold')
    ax.set_ylabel("Spearman Correlation (ρ)", fontsize=6, labelpad=2)
    ax.set_xlabel("")
    ax.set_xticklabels(["Distance\nto Front", "Epi-Dermis\nAxis"], fontsize=6)
    ax.tick_params(axis='y', labelsize=5, length=2, width=0.5)
    ax.tick_params(axis='x', length=0)
    
    sns.despine(top=True, right=True)
    for spine in ax.spines.values(): spine.set_linewidth(0.5)
    
    plt.tight_layout()
    save_path = os.path.join(out_dir, "19dpi_Nature_Swarmplot_Mini_Fib_K14")
    fig.savefig(f"{save_path}.pdf", facecolor='black', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='black', bbox_inches='tight', dpi=600)
    plt.show()
    plt.close(fig)

# ==========================================
# 5. 微缩版并排双热图 (全图渲染，最多标注 10 基因)
# ==========================================
def plot_side_by_side_heatmaps_mini(adata, sig_x_df, sig_y_df, out_dir, window_size=50):
    print("\n🎨 正在绘制微缩版 19dpi 并排瀑布流热图...")
    
    fig, axes = plt.subplots(1, 2, figsize=(2.6, 3.2), dpi=300, layout="constrained")
    fig.suptitle("19dpi | Spatial Dynamics", fontsize=9, fontweight='bold', color='white', y=1.05)
    
    def process_and_plot_heatmap(ax_obj, spatial_var, sig_genes_list, reverse_sort=False, title="", xlabel="", left_label="", right_label=""):
        # 合并出完整的画图基因库（所有显著基因 + 那 10 个指定高亮基因）
        combined_genes = list(set(sig_genes_list + highlight_genes))
        
        valid_genes = [g for g in combined_genes if g in adata.var_names.str.upper()]
        temp_var_names = adata.var_names.str.upper()
        gene_mask = temp_var_names.isin(valid_genes)
        
        expr_matrix = adata.X[:, gene_mask]
        if scipy.sparse.issparse(expr_matrix): expr_matrix = expr_matrix.toarray()
            
        expr_df = pd.DataFrame(expr_matrix, index=adata.obs_names, columns=temp_var_names[gene_mask])
        expr_df[spatial_var] = adata.obs[spatial_var].values
        
        expr_df = expr_df.sort_values(by=spatial_var, ascending=not reverse_sort)
        expr_df.pop(spatial_var) 
        
        zscore_df = expr_df.rolling(window=window_size, min_periods=1, center=True).mean().apply(zscore).dropna(axis=1)
        
        peak_positions = zscore_df.idxmax(axis=0).map({cell: i for i, cell in enumerate(zscore_df.index)})
        zscore_df = zscore_df[peak_positions.sort_values().index]
        
        # 绘制热图本体（数据量没有减少，瀑布流依然平滑）
        sns.heatmap(
            zscore_df.T, cmap=colormap, center=0, vmin=-2, vmax=2, 
            xticklabels=False, yticklabels=True, 
            cbar_kws={'shrink': 0.7, 'aspect': 15, 'location': 'top', 'pad': 0.05}, ax=ax_obj
        )
        
        # 🌟 核心去噪：只在这 10 个基因的位置生成白字标签 🌟
        clean_labels = []
        for gene in zscore_df.columns:
            if gene in highlight_genes:
                clean_labels.append(f"$\it{{{gene}}}$")
            else:
                clean_labels.append("") 
                
        ax_obj.set_yticklabels(clean_labels, rotation=0)
        
        # 统一格式化标签
        for tick, label_obj in zip(ax_obj.yaxis.get_major_ticks(), ax_obj.get_yticklabels()):
            text = label_obj.get_text()
            if text == "":
                tick.tick1line.set_visible(False)
            else:
                tick.tick1line.set_visible(True)
                tick.tick1line.set_color('white')
                tick.tick1line.set_markersize(2)
                # 统一纯白常规字体
                label_obj.set_color('white')
                label_obj.set_fontweight('normal')
                label_obj.set_fontsize(5)
        
        # 字体与标题美化
        ax_obj.set_title(title, fontsize=7, color='white', pad=25) 
        ax_obj.set_ylabel("Genes" if ax_obj == axes[0] else "", fontsize=7, color='white', labelpad=2)
        ax_obj.set_xlabel(xlabel, fontsize=7, color='white', labelpad=8)
        ax_obj.text(0, -0.02, left_label, transform=ax_obj.transAxes, ha='left', va='top', fontsize=5, color='white')
        ax_obj.text(1, -0.02, right_label, transform=ax_obj.transAxes, ha='right', va='top', fontsize=5, color='white')
        
        cbar = ax_obj.collections[0].colorbar
        cbar.ax.tick_params(labelsize=4, length=1.5, width=0.5, colors='white')
        
        for _, spine in ax_obj.spines.items():
            spine.set_visible(True)
            spine.set_color('white')
            spine.set_linewidth(0.5)

    process_and_plot_heatmap(
        ax_obj=axes[0], 
        spatial_var='Distance_to_Front_X', 
        sig_genes_list=sig_x_df['gene'].tolist(),
        reverse_sort=False, 
        title=f"Healing Axis", 
        xlabel="Distance to Front", 
        left_label="Un-epi", right_label="Newly-epi"
    )

    process_and_plot_heatmap(
        ax_obj=axes[1], 
        spatial_var='Relative_Depth_Y', 
        sig_genes_list=sig_y_df['gene'].tolist(),
        reverse_sort=True, 
        title=f"Depth Axis", 
        xlabel="Relative Depth", 
        left_label="Deep Dermis", right_label="Epidermis"
    )

    save_path = os.path.join(out_dir, "19dpi_Nature_SideBySide_Heatmaps_Mini")
    fig.savefig(f"{save_path}.pdf", facecolor='black', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='black', bbox_inches='tight', dpi=600)
    plt.show() 
    plt.close(fig)
    print(f"✅ 并排微缩热图已保存至: {save_path}.pdf/png")

# ==========================================
# 6. 执行流
# ==========================================
if __name__ == "__main__":
    # ⚠️ 请确保您的 adata 已经加载，且过滤出了 Fib_K14 细胞
    # fib_adata = adata[adata.obs['sub_labels'] == 'Fib_K14'].copy()
    
    sig_genes_x, sig_genes_y = get_spatially_variable_genes(fib_adata, rho_threshold=0.15)
    
    # 1. 绘制极窄微缩版蜂群图 
    plot_nature_swarm_mini(
        adata=fib_adata, 
        vars=['Distance_to_Front_X', 'Relative_Depth_Y'], 
        highlight_list=highlight_genes, 
        out_dir=output_dir,
        rho_cutoff=0.15
    )
    
    # 2. 绘制极窄微缩版并排双热图
    plot_side_by_side_heatmaps_mini(
        adata=fib_adata, 
        sig_x_df=sig_genes_x, 
        sig_y_df=sig_genes_y, 
        out_dir=output_dir,
        window_size=50 
    )


######################特征基因时空分布图
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
# 1. 全局排版配置：黑底白字 + Arial 字体
# ==========================================
plt.style.use('dark_background')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'white'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['axes.facecolor'] = 'black'
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
    ax.plot([b1_t, b1_b], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.plot([b2_t, b2_b], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.add_patch(Polygon([[-max_abs_x, 0], [b1_t, 0], [b1_b, 1], [-max_abs_x, 1]], color='#e5192c', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b1_t, 0], [b2_t, 0], [b2_b, 1], [b1_b, 1]], color='#E1C62F', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b2_t, 0], [max_abs_x, 0], [max_abs_x, 1], [b2_b, 1]], color='#3a77b7', alpha=0.15, zorder=0, lw=0))

    trans = ax.get_xaxis_transform()
    c_un_epi, c_epi_front, c_newly_epi = (-max_abs_x + b1_t)/2, (b1_t + b2_t)/2, (b2_t + max_abs_x)/2
    ax.text(c_un_epi, 1.02, 'Un-epi', color='#e5192c', transform=trans, ha='center', va='bottom', fontsize=8)
    ax.text(c_epi_front, 1.02, 'Epi-Front', color='#E1C62F', transform=trans, ha='center', va='bottom', fontsize=8)
    ax.text(c_newly_epi, 1.02, 'Newly-epi', color='#3a77b7', transform=trans, ha='center', va='bottom', fontsize=8)

    sns.kdeplot(data=pd.DataFrame({'x': x, 'y': y}), x='x', y='y', ax=ax, color="#FFFFFF", linewidths=0.5, alpha=0.6, levels=6, zorder=2)
    mask_z, mask_e = expr_s <= 0, expr_s > 0
    ax.scatter(x_s[mask_z], y_s[mask_z], color='#333333', s=1.5, alpha=0.3, zorder=3)
    sc_scatter = ax.scatter(x_s[mask_e], y_s[mask_e], c=expr_s[mask_e], cmap=cmap, s=1.5, edgecolor='none', zorder=4)
    
    # 🌟 统一 Colorbar：仅在要求时，绘制于图像外部的最右侧 🌟
    if show_colorbar and mask_e.sum() > 0:
        cbar_ax = ax.inset_axes([1.08, 0.1, 0.05, 0.8]) # 外部垂直悬浮
        cbar = plt.colorbar(sc_scatter, cax=cbar_ax, orientation='vertical')
        cbar.set_ticks([expr_s[mask_e].min(), expr_s[mask_e].max()])
        cbar.set_ticklabels(['Min', 'Max'])
        cbar.ax.tick_params(colors='white', labelsize=8, length=0, pad=2)
        cbar.set_label('Expression', color='white', fontsize=12, labelpad=6)

    cx, cy = -500, 0.75 
    ax.annotate('', xy=(1800, 0.2), xycoords='data', xytext=(cx, cy), textcoords='data',
                arrowprops=dict(facecolor="#0099CC", edgecolor="#0099CC", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)
    ax.annotate('', xy=(-3200, 0.98), xycoords='data', xytext=(cx, cy), textcoords='data',
                arrowprops=dict(facecolor="#FFFF00", edgecolor="#FFFF00", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)

    ax.set_xlim(-max_abs_x, max_abs_x) 
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
    ax.set_ylim(1.05, -0.05) 
    ax.set_yticks([0.0, 0.5, 1.0])
    
    # 🌟 强制正方形比例 (Square Aspect Ratio) 🌟
    ax.set_box_aspect(1) 
    
    ax.set_title(f"$\it{{{feature}}}$" if feature != 'MergeScore' else feature, fontsize=16, color='white', pad=12, fontweight='normal')
    
    # 🌟 全局统一 X/Y 轴标签设定 🌟
    ax.set_xlabel("Distance to leading edge (µm)" if set_global_xlabel else "", color='white', fontsize=16 if set_global_xlabel else 0, labelpad=12)
    ax.set_ylabel("Epidermis-Dermis axis" if show_y_label else "", color='white', fontsize=16, labelpad=8)
    
    ax.tick_params(colors='white', direction='out', length=2.5, width=0.6, labelsize=14)
    for spine in ax.spines.values(): spine.set_linewidth(0.6)

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
    ax.scatter(x_s[mask_z], y_s[mask_z], color='#222222', s=0.8, alpha=0.3) 
    sc_scatter = ax.scatter(x_s[mask_e], y_s[mask_e], c=expr_s[mask_e], cmap=cmap, s=1.5, edgecolor='none')
    
    min_x, max_x = x.min(), x.max()
    min_y, max_y = y.min(), y.max()
    y_range = max_y - min_y
    
    # 🌟 唯一的外部右侧 Colorbar 🌟
    if show_colorbar and mask_e.sum() > 0:
        cbar_ax = ax.inset_axes([1.08, 0.05, 0.05, 0.9])
        cbar = plt.colorbar(sc_scatter, cax=cbar_ax, orientation='vertical')
        cbar.set_ticks([expr_s[mask_e].min(), expr_s[mask_e].max()])
        cbar.set_ticklabels(['Min', 'Max'])
        cbar.ax.tick_params(colors='white', labelsize=8, length=0, pad=2)
        cbar.set_label('Expression', color='white', fontsize=12, labelpad=6)
        
    # 🌟 唯一的右侧 Scalebar 🌟
    if show_scalebar:
        scale_px = 1000 / 0.33 
        bar_y = min_y - y_range * 0.08
        ax.plot([max_x - scale_px, max_x], [bar_y, bar_y], color='white', lw=1.5)
        ax.text(max_x - scale_px/2, bar_y - y_range * 0.03, '1000 µm', color='white', fontsize=10, ha='center', va='top')

    # 🌟 核心：愈合方向的箭头从已愈合方向 (max_x) 指向 未愈合方向 (min_x) 🌟
    arrow_y = max_y + y_range * 0.12
    ax.annotate('', xy=(min_x, arrow_y), xytext=(max_x, arrow_y), 
                arrowprops=dict(arrowstyle="->", color="white", lw=1.2))
    
    text_y = arrow_y + y_range * 0.04
    ax.text(min_x, text_y, 'Un-epi', color='#e5192c', fontsize=9, ha='left', va='bottom')
    ax.text((min_x + max_x)/2, text_y, 'Epi-Front', color='#E1C62F', fontsize=9, ha='center', va='bottom')
    ax.text(max_x, text_y, 'Newly-epi', color='#3a77b7', fontsize=9, ha='right', va='bottom')

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
    
    fig = plt.figure(figsize=(15.5, 11.5), dpi=300) # 再次微调宽度以容纳右侧额外的 Colorbar
    
    gs = fig.add_gridspec(4, 5, hspace=0.35, wspace=0.35, height_ratios=[1.4, 0.6, 1.0, 1.0])
    
    # --- Row 1: 动态时空投影图 ---
    print("🎨 绘制第一行: 动态时空表达图 (强迫正方形 + 唯一全局标签与图例)...")
    for i, feature in enumerate(display_features):
        ax = fig.add_subplot(gs[0, i])
        
        # 仅第 1 个有 Y 轴标签，仅第 3 个有全行公用的 X 轴标签，仅第 5 个有 Colorbar
        is_first = (i == 0)
        is_middle = (i == 2) 
        is_last = (i == 4)
        
        plot_spatial_expression_planA(ax, adata, feature, colormap, 
                                      show_y_label=is_first, set_global_xlabel=is_middle, show_colorbar=is_last)
        if feature == 'MergeScore':
            for spine in ax.spines.values(): spine.set_linewidth(1.0)

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
            sns.barplot(data=df_bars[df_bars['Feature'] == feature], x=BURN_DEPTH_COL, y='Value', 
                        ax=ax, palette=['#E15759', '#4E79A7'], errcolor='white', errwidth=1.0, capsize=0.15)
            ax.set_ylabel("Mean Expr." if i == 0 else "", fontsize=10)
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=30, labelsize=11)
            ax.tick_params(axis='y', labelsize=10)
            sns.despine(ax=ax, top=True, right=True)
            for spine in ax.spines.values(): spine.set_linewidth(0.6)

        # --- Row 4: 愈合区域差异柱状图 ---
        print("🎨 绘制第四行: 愈合区域差异柱状图...")
        cat_type = pd.CategoricalDtype(categories=['Un-epi', 'Epi-Front', 'Newly-epi'], ordered=True)
        df_bars[HEALING_REGION_COL] = df_bars[HEALING_REGION_COL].astype(cat_type)
        
        for i, feature in enumerate(display_features):
            ax = fig.add_subplot(gs[3, i])
            sns.barplot(data=df_bars[df_bars['Feature'] == feature], x=HEALING_REGION_COL, y='Value', 
                        ax=ax, palette=['#e5192c', '#E1C62F', '#3a77b7'], errcolor='white', errwidth=1.0, capsize=0.15)
            ax.set_ylabel("Mean Expr." if i == 0 else "", fontsize=10)
            ax.set_xlabel("")
            ax.tick_params(axis='x', rotation=30, labelsize=11)
            ax.tick_params(axis='y', labelsize=10)
            sns.despine(ax=ax, top=True, right=True)
            for spine in ax.spines.values(): spine.set_linewidth(0.6)
            
    else:
        print(f"⚠️ 找不到列 '{BURN_DEPTH_COL}' 或 '{HEALING_REGION_COL}'，跳过柱状图绘制。")

    # 保存与展示
    plt.tight_layout()
    save_path = os.path.join(out_dir, "MainFigure_Spatial_Module_Panel_Ultimate")
    fig.savefig(f"{save_path}.pdf", facecolor='black', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='black', bbox_inches='tight', dpi=600)
    plt.show()
    plt.close(fig)
    print(f"✅ 终极主图排版已完美保存至: {save_path}.pdf/png")

# ==========================================
# 5. 执行代码
# ==========================================
if __name__ == "__main__":
    output_dir = r"data/external/local_8e645d0a9b03"
    plot_master_figure(adata, output_dir)


######################iPSC分析
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
# 1. 全局排版配置：黑底白字 + Arial 字体
# ==========================================
plt.style.use('dark_background')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['text.color'] = 'white'
plt.rcParams['axes.labelcolor'] = 'white'
plt.rcParams['xtick.color'] = 'white'
plt.rcParams['ytick.color'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['figure.facecolor'] = 'black'
plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['figure.dpi'] = 150

# ==========================================
# 2. 核心基因库配置 (对应三大命运轨迹)
# ==========================================
zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

# 🌟 核心更新：IGFL1 替换为 MYC，为后续 5 列布局做铺垫 🌟
phenotype_dict = {
    "iPSC Reprogramming": ['KLF4', 'SOX2', 'POU5F1', 'MYC'], 
    "Epithelialization": ['GRHL3', 'TP63', 'TACSTD2', 'KRT14'],
    "Scar Formation": ['COL1A1', 'ACTA2', 'FN1', 'TGFB1']
}

SPATIAL_X = 'Distance_to_Front_X'
SPATIAL_Y = 'Relative_Depth_Y'

# ==========================================
# 3. 核心功能函数
# ==========================================
def calculate_mergescore(adata, genes, score_name='MergeScore'):
    """计算模块综合评分"""
    valid_genes = [g for g in genes if g in adata.var_names.str.upper()]
    temp_vars = adata.var_names.copy()
    adata.var_names = adata.var_names.str.upper()
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)
    adata.var_names = temp_vars 
    return adata

def get_expression_data(adata, feature):
    """安全提取表达量/模块评分"""
    if feature in adata.obs.columns:
        return adata.obs[feature].values
    else:
        gene_idx = np.where(adata.var_names.str.upper() == feature.upper())[0]
        if len(gene_idx) > 0:
            expr_data = adata[:, gene_idx[0]].X
            return np.asarray(expr_data.toarray()).flatten() if scipy.sparse.issparse(expr_data) else np.asarray(expr_data).flatten()
        return np.zeros(adata.n_obs)

def plot_phenotype_spatial_ultimate(ax, adata, feature, phenotype_group, cmap, 
                                    show_y_label=False, show_colorbar=False, show_top_labels=False):
    """100% 克隆定稿主图参数的轨迹绘制函数"""
    PIXEL_TO_UM = 0.33  
    x = adata.obs[SPATIAL_X].values * PIXEL_TO_UM
    y = adata.obs[SPATIAL_Y].values
    expr = get_expression_data(adata, feature)
    
    idx = np.argsort(expr)
    x_s, y_s, expr_s = x[idx], y[idx], expr[idx]
    
    # 1. 倾斜分界线与多边形背景着色 (linewidth=0.8, alpha=0.8 完全同步)
    max_abs_x = max(abs(x.min()), abs(x.max())) * 1.1 if len(x) > 0 else 4000
    b1_t, b1_b, b2_t, b2_b = -600, -300, 400, 700   
    ax.plot([b1_t, b1_b], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.plot([b2_t, b2_b], [0, 1], color='white', linestyle='--', linewidth=0.8, alpha=0.8, zorder=1)
    ax.add_patch(Polygon([[-max_abs_x, 0], [b1_t, 0], [b1_b, 1], [-max_abs_x, 1]], color='#e5192c', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b1_t, 0], [b2_t, 0], [b2_b, 1], [b1_b, 1]], color='#E1C62F', alpha=0.15, zorder=0, lw=0))
    ax.add_patch(Polygon([[b2_t, 0], [max_abs_x, 0], [max_abs_x, 1], [b2_b, 1]], color='#3a77b7', alpha=0.15, zorder=0, lw=0))

    # 2. 顶部分区文字 (仅在第一行显示，fontsize=8 极简排版)
    if show_top_labels:
        trans = ax.get_xaxis_transform()
        c_un_epi, c_epi_front, c_newly_epi = (-max_abs_x + b1_t)/2, (b1_t + b2_t)/2, (b2_t + max_abs_x)/2
        ax.text(c_un_epi, 1.02, 'Un-epi', color='#e5192c', transform=trans, ha='center', va='bottom', fontsize=8)
        ax.text(c_epi_front, 1.02, 'Epi-Front', color='#E1C62F', transform=trans, ha='center', va='bottom', fontsize=8)
        ax.text(c_newly_epi, 1.02, 'Newly-epi', color='#3a77b7', transform=trans, ha='center', va='bottom', fontsize=8)

    # 3. KDE 等高线 (levels=6) 与散点 (s=1.5)
    sns.kdeplot(data=pd.DataFrame({'x': x, 'y': y}), x='x', y='y', ax=ax, color="#FFFFFF", linewidths=0.5, alpha=0.6, levels=6, zorder=2)
    mask_z, mask_e = expr_s <= 0, expr_s > 0
    ax.scatter(x_s[mask_z], y_s[mask_z], color='#333333', s=1.5, alpha=0.3, zorder=3)
    sc_scatter = ax.scatter(x_s[mask_e], y_s[mask_e], c=expr_s[mask_e], cmap=cmap, s=1.5, edgecolor='none', zorder=4)
    
    # 🌟 4. 唯一的外部右侧 Colorbar 🌟
    if show_colorbar and mask_e.sum() > 0:
        cbar_ax = ax.inset_axes([1.08, 0.1, 0.05, 0.8])
        cbar = plt.colorbar(sc_scatter, cax=cbar_ax, orientation='vertical')
        cbar.set_ticks([expr_s[mask_e].min(), expr_s[mask_e].max()])
        cbar.set_ticklabels(['Min', 'Max'])
        cbar.ax.tick_params(colors='white', labelsize=8, length=0, pad=2)
        cbar.set_label('Expression', color='white', fontsize=12, labelpad=6)

    # 🌟 5. 命运抉择：按表型智能配置箭头指向 🌟
    cx, cy = -500, 0.75 
    if phenotype_group in ["iPSC Reprogramming", "Epithelialization"]:
        # 命运A：向上、向右 (向表皮化与前沿推进)
        ax.annotate('', xy=(1800, 0.2), xycoords='data', xytext=(cx, cy), textcoords='data',
                    arrowprops=dict(facecolor="#0099CC", edgecolor="#0099CC", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)
    
    if phenotype_group in ["iPSC Reprogramming", "Scar Formation"]:
        # 命运B：向下、向左 (向深层纤维化与未愈合区退避)
        ax.annotate('', xy=(-3200, 0.98), xycoords='data', xytext=(cx, cy), textcoords='data',
                    arrowprops=dict(facecolor="#FFFF00", edgecolor="#FFFF00", width=0.3, headwidth=2.5, shrink=0.08), zorder=5)

    # 6. 坐标轴与排版控制
    ax.set_xlim(-max_abs_x, max_abs_x) 
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune='both'))
    ax.set_ylim(1.05, -0.05) 
    ax.set_yticks([0.0, 0.5, 1.0])
    
    # 🌟 核心：强制正方形比例 (Square Aspect Ratio) 🌟
    ax.set_box_aspect(1) 
    
    # 🌟 动态标题逻辑：如果是模块评分列，则取消斜体并加粗 🌟
    if feature in phenotype_dict.keys():
        ax.set_title(feature, fontsize=15, color='white', pad=12, fontweight='bold')
    else:
        ax.set_title(f"$\it{{{feature}}}$", fontsize=16, color='white', pad=12, fontweight='normal')
    
    # Y轴标签 (仅最左侧显示)
    ax.set_ylabel("Epidermis-Dermis axis" if show_y_label else "", color='white', fontsize=16, labelpad=8)
    
    # 刻度细节同步
    ax.tick_params(colors='white', direction='out', length=2.5, width=0.6, labelsize=14)
    for spine in ax.spines.values(): spine.set_linewidth(0.6)

# ==========================================
# 4. 终极 3x5 黄金矩阵组图逻辑
# ==========================================
def plot_phenotype_master_figure(adata, out_dir):
    print("⏳ 正在计算表型对应基因集的评分 (Module Scores)...")
    
    # 🌟 核心步骤：将之前的 MergeScore 改为各自独立的表型评分 🌟
    for phenotype, genes in phenotype_dict.items():
        adata = calculate_mergescore(adata, genes, score_name=phenotype)
    
    # 🌟 画布适度加宽，完美容纳 3x5 结构与侧边图例 🌟
    fig = plt.figure(figsize=(19.5, 12.5), dpi=300)
    gs = fig.add_gridspec(3, 5, hspace=0.35, wspace=0.35)
    
    phenotypes = list(phenotype_dict.keys())
    
    for row_idx, phenotype in enumerate(phenotypes):
        # 展示特征：4 个代表基因 + 1 个基因集综合评分 (直接拿表型当名字)
        display_features = phenotype_dict[phenotype] + [phenotype]
        
        print(f"🎨 绘制第 {row_idx+1} 行: {phenotype} 轨迹...")
        
        for col_idx, feature in enumerate(display_features):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            
            # 排版开关：仅最左侧显 Y 轴，最右侧显 Colorbar，最上层显顶部分区
            is_first_col = (col_idx == 0)
            is_last_col = (col_idx == 4)
            is_first_row = (row_idx == 0)
            
            plot_phenotype_spatial_ultimate(
                ax, adata, feature, phenotype, colormap, 
                show_y_label=is_first_col, 
                show_colorbar=is_last_col, 
                show_top_labels=is_first_row
            )
            
            # 为最后一张综合评分图加亮外边框，彰显这是该行的 Merge Score
            if is_last_col:
                for spine in ax.spines.values(): spine.set_linewidth(1.0)
            
            # 行标签
            if is_first_col:
                ax.text(-0.42, 0.5, phenotype, transform=ax.transAxes, rotation=90, 
                        ha='center', va='center', fontsize=18, fontweight='bold', color='#DDDDDD')

    # 🌟 完美居中的全局 X 轴标题 🌟
    fig.text(0.5, 0.05, "Distance to leading edge (µm)", ha='center', va='center', fontsize=18, color='white')

    # 防止左右文字与图例被裁剪
    plt.tight_layout(rect=[0.04, 0.08, 0.96, 1.0]) 
    
    save_path = os.path.join(out_dir, "Cell_Fate_Trajectories_Ultimate_3x5")
    fig.savefig(f"{save_path}.pdf", facecolor='black', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='black', bbox_inches='tight', dpi=600)
    plt.show()
    plt.close(fig)
    print(f"✅ 细胞命运抉择轨迹图 (终极 3x5 版) 已完美保存至: {save_path}.pdf/png")

# ==========================================
# 5. 执行代码
# ==========================================
if __name__ == "__main__":
    output_dir = r"data/external/local_8e645d0a9b03"
    
    # ⚠️ 强烈建议您使用目标亚群 (Fib_K14) 数据框执行，效果最震撼
    # fib_adata = adata[adata.obs['sub_labels'] == 'Fib_K14'].copy()
    # plot_phenotype_master_figure(fib_adata, output_dir)
    
    plot_phenotype_master_figure(adata, output_dir)


###################################空间共定位分析
import os
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import igraph
import random
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint
from sklearn.preprocessing import MinMaxScaler
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# ==========================================
# 0. 核心配置与字典解析
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
TEST_KEY = "19dpb"  
target_cell = "Fib_K14"
ROTATION_ANGLE = -225
TARGET_SAMPLE = "BW15D_C6_SDSDB_19dpb"

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad")
json_load_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

# 引入精准字典
raw_annotation_config = {
    "KC_Basal": {"color": "#008941"}, "KC_Spinous": {"color": "#00FF00FF"}, "KC_Spinous_Mat": {"color": "#EE4C97FF"},
    "KC_Granular": {"color": "#882D17"}, "KC_Spinous_Mig": {"color": "#1A0099FF"}, "KC_Basal_Mig": {"color": "#FF3030"},
    "KC_Basal_Prolif": {"color": "#FFD700"}, "SAC_SG_Progenitor": {"color": "#0099CCFF"}, "SAC_SG_Clear": {"color": "#B7E4F9FF"},
    "SAC_SG_Dark": {"color": "#193006FF"}, "SAC_SG_Ductal": {"color": "#8F7700FF"}, "SAC_HF_IRS": {"color": "#CC9900FF"},
    "SAC_HF_ORS": {"color": "#A30059"}, "SAC_HF_HFSC": {"color": "#00F5FF"}, "SAC_HF_DP_DS": {"color": "#CDBE70"},
    "SAC_HF_Matrix": {"color": "#7FFFD4"}, "Melanocyte": {"color": "#FF4500"}, "Fib_Papi": {"color": "#3DA873FF"},
    "Fib_SFRP2": {"color": "#CCFF00FF"}, "Fib_Myo": {"color": "#B79762"}, "Fib_Inflama": {"color": "#800000FF"},
    "Fib_Fasci": {"color": "#9370DB"}, "Fib_EN1": {"color": "#BDB76B"}, "Fib_CD45": {"color": "#DA70D6"},
    "Fib_Prolif": {"color": "#EED2EE"}, "Fib_K14": {"color": "#4682B4"}, "Endo_Capillary": {"color": "#CF4E9CFF"},
    "Endo_Arterial": {"color": "#FF00FF"}, "Endo_Venous": {"color": "#4169E1"}, "Endo_Lymphatic": {"color": "#ADFF2F"},
    "Endo_Prolif": {"color": "#F08080"}, "Pericyte": {"color": "#FFFF00"}, "Schwann": {"color": "#7CFC00"},
    "M1_Macrophage": {"color": "#FF4040"}, "M2_Macrophage": {"color": "#54FF9F"}, "LAM": {"color": "#8B4513"},
    "Neutrophil": {"color": "#FF7F50"}, "Mast": {"color": "#FF1493"}, "cDC1": {"color": "#9B30FF"},
    "cDC2": {"color": "#EE9A00"}, "pDC": {"color": "#D1EEEE"}, "LGHS": {"color": "#FFB5C5"},
    "B_Mature": {"color": "#0000CD"}, "B_Plasma": {"color": "#FAF0E6"}, "CD4_T": {"color": "#00BFFF"},
    "CD8_T": {"color": "#1E90FF"}, "NK": {"color": "#00E5EE"}
}

cell_colors = {k: v["color"].replace("FF", "") if len(v["color"])==9 else v["color"] for k, v in raw_annotation_config.items()}

def get_macro_group(c):
    if c.startswith("KC_") or c.startswith("SAC_") or c == "Melanocyte": 
        return "Epithelial & Appendage"
    if c.startswith("Fib_") or c.startswith("Endo_") or c in ["Pericyte", "Schwann"]: 
        return "Stromal & Vascular"
    return "Immune Cells"

bubble_colors = {
    "Epithelial & Appendage": "#1e3d59",  
    "Stromal & Vascular": "#4d2c1c",      
    "Immune Cells": "#1f3b25"             
}

plt.style.use('dark_background')
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'text.color': 'white', 
                     'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white'})

# ==========================================
# 1. 空间投影计算与划分模块
# ==========================================
def process_spatial_data(adata, anchor_dict, angle):
    theta = np.radians(angle)
    raw_x = adata.obsm['spatial'][:, 0] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 0]
    raw_y = adata.obsm['spatial'][:, 1] if 'spatial' in adata.obsm else adata.obsm['X_spatial'][:, 1]
    rot_x, rot_y = raw_x * np.cos(theta) - raw_y * np.sin(theta), raw_x * np.sin(theta) + raw_y * np.cos(theta)
    adata.obsm['X_rot_spatial'] = np.column_stack((rot_x, rot_y))
    
    epi_line, der_line = LineString(anchor_dict['epi_baseline']), LineString(anchor_dict['der_bottom'])
    has_front = len(anchor_dict.get('leading_edge', [])) >= 2
    if has_front: front_line = LineString(anchor_dict['leading_edge'])

    x_list, y_list = [], []
    for i in range(len(rot_x)):
        pt = Point(rot_x[i], rot_y[i])
        dist_to_epi, dist_to_der = pt.distance(epi_line), pt.distance(der_line)
        y_list.append(dist_to_epi / (dist_to_epi + dist_to_der) if (dist_to_epi + dist_to_der) > 0 else 0)
        if has_front:
            sign = 1 if (pt.x > anchor_dict['leading_edge'][0][0] and anchor_dict['healed_direction'] == "right_is_healed") or \
                        (pt.x < anchor_dict['leading_edge'][0][0] and anchor_dict['healed_direction'] != "right_is_healed") else -1
            x_list.append(pt.distance(front_line) * sign)
        else: x_list.append(2000.0)
            
    adata.obs['Plot_X'] = np.array(x_list) * 0.33
    adata.obs['Plot_Y'] = y_list
    
    max_abs_x = max(abs(adata.obs['Plot_X'].min()), abs(adata.obs['Plot_X'].max())) * 1.1
    b1_t, b1_b, b2_t, b2_b = -600, -300, 400, 700   
    gates = {
        "Un-epi": Polygon([[-max_abs_x, 0], [b1_t, 0], [b1_b, 1], [-max_abs_x, 1]]),
        "Epi-Front": Polygon([[b1_t, 0], [b2_t, 0], [b2_b, 1], [b1_b, 1]]),
        "Newly-epi": Polygon([[b2_t, 0], [max_abs_x, 0], [max_abs_x, 1], [b2_b, 1]])
    }
    gpd_poly = gpd.GeoDataFrame({"gates": gpd.GeoSeries(gates)}, geometry="gates")
    gpd_cells = gpd.GeoDataFrame({"cells": gpd.GeoSeries.from_xy(adata.obs['Plot_X'], adata.obs['Plot_Y'])}, geometry="cells")
    adata.obs["Region"] = gpd.sjoin(gpd_cells, gpd_poly, how="left")["index_right"].fillna("Unknown")
    return adata

# ==========================================
# 2. 自定义中心辐射网络引擎
# ==========================================
def generate_radial_layout(nodes, target_cell):
    pos_dict = {target_cell: (0, 0)} 
    groups = {"Epithelial & Appendage": [], "Stromal & Vascular": [], "Immune Cells": []}
    
    for n in nodes:
        if n == target_cell: continue
        groups[get_macro_group(n)].append(n)
        
    angles = {
        "Epithelial & Appendage": (np.radians(30), np.radians(150)),
        "Stromal & Vascular": (np.radians(165), np.radians(255)),
        "Immune Cells": (np.radians(285), np.radians(375))
    }
    
    for g_name, g_nodes in groups.items():
        if not g_nodes: continue
        g_nodes.sort() 
        t_min, t_max = angles[g_name]
        thetas = np.linspace(t_min, t_max, len(g_nodes))
        
        radii_pool = [0.35, 0.65, 0.95, 1.25]
        radii = [radii_pool[i % 4] for i in range(len(g_nodes))]
        
        for i, n in enumerate(g_nodes):
            pos_dict[n] = (radii[i] * np.cos(thetas[i]), radii[i] * np.sin(thetas[i]))
            
    return pos_dict

def plot_nature_network_v2(ax, g_reg, pos_dict, title, target_cell, top5_nodes):
    nodes = g_reg.vs['label']
    coords = np.array([pos_dict[n] for n in nodes])
    
    for group, bcolor in bubble_colors.items():
        group_idx = [i for i, n in enumerate(nodes) if get_macro_group(n) == group and n != target_cell]
        if len(group_idx) >= 3:
            pts = coords[group_idx]
            pts = np.vstack([pts, [0, 0]]) 
            poly = MultiPoint(pts).convex_hull.buffer(0.08, join_style=1)
            if not poly.is_empty:
                x, y = poly.exterior.xy
                ax.fill(x, y, color=bcolor, alpha=0.65, zorder=1)
                ax.plot(x, y, color=bcolor, alpha=1.0, linewidth=2.0, linestyle='--', zorder=1)

    for edge in g_reg.es:
        n0, n1 = nodes[edge.tuple[0]], nodes[edge.tuple[1]]
        x0, y0 = pos_dict[n0]
        x1, y1 = pos_dict[n1]
        
        if (n0 == target_cell and n1 in top5_nodes) or (n1 == target_cell and n0 in top5_nodes):
            color = cell_colors.get(n1 if n0 == target_cell else n0, "#FF1493")
            ax.plot([x0, x1], [y0, y1], color=color, alpha=0.9, linewidth=3.5, zorder=3)
        elif n0 == target_cell or n1 == target_cell:
            ax.plot([x0, x1], [y0, y1], color="#FFFFFF", alpha=0.2, linewidth=0.5, zorder=2)
        else:
            ax.plot([x0, x1], [y0, y1], color="#FFFFFF", alpha=0.03, linewidth=0.3, zorder=1)

    c_list = [cell_colors.get(n, "#888888") for n in nodes]
    ax.scatter(coords[:,0], coords[:,1], c=c_list, s=150, edgecolors='white', linewidths=1.5, zorder=4)
    ax.scatter([0], [0], facecolors='none', edgecolors='white', s=500, linewidths=2.5, zorder=4)

    for n in nodes:
        x, y = pos_dict[n]
        if n == target_cell:
            ax.text(x, y - 0.12, n, ha='center', va='top', fontsize=20, fontweight='bold', color='white', zorder=6)
            continue
            
        theta = np.arctan2(y, x)
        offset_x = 0.055 * np.cos(theta)
        offset_y = 0.055 * np.sin(theta)
        
        ha = 'left' if np.cos(theta) > 0 else 'right'
        va = 'bottom' if np.sin(theta) > 0 else 'top'
        
        # 🌟 核心优化 2：局部文字重排 (精准 Nudge 逻辑避让重叠)
        if n == "SAC_SG_Progenitor":
            offset_x += 0.12  # 向右移动
            offset_y += 0.10  # 向上移动
            ha, va = 'left', 'bottom'
        elif n == "Endo_Arterial":
            offset_y -= 0.15  # 向下移动
            ha, va = 'center', 'top'
        elif n == "M2_Macrophage":
            offset_x -= 0.15  # 向左移动
            offset_y -= 0.10  # 向下移动
            ha, va = 'right', 'top'
        
        ax.text(x + offset_x, y + offset_y, n, ha=ha, va=va, fontsize=12, color='white', zorder=5)

    ax.set_title(title, fontsize=24, color='white', pad=25, fontweight='bold')
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.axis('off')

# ==========================================
# 3. 运行控制台
# ==========================================
if __name__ == "__main__":
    print(f"⏳ 正在加载并清理数据...")
    adata = sc.read_h5ad(file_path)
    
    valid_cells = list(raw_annotation_config.keys())
    adata = adata[adata.obs['sub_labels'].isin(valid_cells)].copy()
    adata.obs['sub_labels'] = adata.obs['sub_labels'].cat.remove_unused_categories()
    
    with open(json_load_path, 'r', encoding='utf-8') as f:
        anchor_dict = json.load(f)
    
    adata = process_spatial_data(adata, anchor_dict, ROTATION_ANGLE)
    adata = adata[adata.obs['Region'] != "Unknown"].copy()

    print("🎨 正在生成 Nature 风格扁平超级组图...")
    fig = plt.figure(figsize=(25, 12), dpi=300)
    
    # 🌟 核心优化 1：分割 GridSpec 为 4 列。前 3 列给热图本体，第 4 列给 Colorbar
    # 这样热图就不会被 Colorbar 占用宽度，实现与 Network 的两端绝对对齐！
    gs = gridspec.GridSpec(2, 4, height_ratios=[3.8, 1.1], width_ratios=[1, 1, 1, 0.05], hspace=0.25, wspace=0.15)

    regions = ["Un-epi", "Epi-Front", "Newly-epi"]
    axes_net = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]
    
    # 跨度绑定：热图占用第 0 到 第 2 列 (完全覆盖上面的三个网络图)
    ax_heat = fig.add_subplot(gs[1, 0:3])
    # 极窄的色条独占第 3 列
    cbar_ax = fig.add_subplot(gs[1, 3])

    global_cats = list(adata.obs['sub_labels'].cat.categories)
    m_target_matrix = np.zeros((3, len(global_cats)))

    global_layout_dict = generate_radial_layout(global_cats, target_cell)

    for idx, region in enumerate(regions):
        ax = axes_net[idx]
        adata_sub = adata[adata.obs['Region'] == region].copy()
        
        if len(adata_sub) < 5 or target_cell not in adata_sub.obs['sub_labels'].values:
            ax.text(0.5, 0.5, f"No sufficient {target_cell} in {region}", color='grey', ha='center', fontsize=18)
            ax.axis('off')
            continue
            
        adata_sub.obs['sub_labels'] = adata_sub.obs['sub_labels'].cat.remove_unused_categories()
        local_cats = list(adata_sub.obs['sub_labels'].cat.categories)

        sq.gr.spatial_neighbors(adata_sub, spatial_key="X_rot_spatial")
        sq.gr.interaction_matrix(adata_sub, cluster_key="sub_labels", normalized=True)
        M_reg = adata_sub.uns["sub_labels_interactions"]
        
        target_idx_local = local_cats.index(target_cell)
        interaction_row = M_reg[target_idx_local, :]
        
        temp = interaction_row.copy()
        temp[target_idx_local] = -1 
        top5_idx = temp.argsort()[-5:][::-1]
        top5_nodes = [local_cats[i] for i in top5_idx if temp[i] > 0]
        
        if region == "Newly-epi" and "KC_Basal_Mig" in local_cats:
            if "KC_Basal_Mig" not in top5_nodes:
                top5_nodes.append("KC_Basal_Mig")
                
        for t_node in top5_nodes:
            idx_t = local_cats.index(t_node)
            M_reg[target_idx_local, idx_t] = max(M_reg[target_idx_local, idx_t], 0.05)
            M_reg[idx_t, target_idx_local] = max(M_reg[idx_t, target_idx_local], 0.05)
        
        for i, cat_name in enumerate(local_cats):
            global_idx = global_cats.index(cat_name)
            m_target_matrix[idx, global_idx] = interaction_row[i]

        M_reg[M_reg < 0.05] = 0
        g_reg = igraph.Graph.Weighted_Adjacency(M_reg)
        g_reg.vs["label"] = local_cats 
        nodes_to_keep = [v.index for v in g_reg.vs if g_reg.degree(v) > 0 or v["label"] == target_cell]
        g_reg = g_reg.subgraph(nodes_to_keep)

        plot_nature_network_v2(ax, g_reg, global_layout_dict, f"Spatial Network ({region})", target_cell, top5_nodes)

    legend_elements = [Patch(facecolor=color, edgecolor='white', alpha=0.7, label=group) for group, color in bubble_colors.items()]
    # 微调图例居中参数，自适应 4 列 GridSpec 画板宽度
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.48, 0.36), ncol=3, fontsize=16, frameon=False)

    print("🔥 正在计算自适应宽度大热图...")
    valid_cols = np.sum(m_target_matrix, axis=0) > 0
    m_target_clean = m_target_matrix[:, valid_cols]
    labels_clean = np.array(global_cats)[valid_cols]

    group_order = {"Epithelial & Appendage": 1, "Stromal & Vascular": 2, "Immune Cells": 3}
    sort_keys = [(group_order[get_macro_group(c)], c, i) for i, c in enumerate(labels_clean)]
    sort_keys.sort()
    
    sorted_indices = [item[2] for item in sort_keys]
    labels_sorted = [item[1] for item in sort_keys]
    m_target_sorted = m_target_clean[:, sorted_indices]

    m_target_transformed = np.power(m_target_sorted, 0.4) 

    scaler = MinMaxScaler()
    normalized_array = scaler.fit_transform(m_target_transformed.T).T

    zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
    colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

    # 🌟 将色条精准投放至独立第 4 列
    sns.heatmap(
        normalized_array,
        xticklabels=labels_sorted,
        yticklabels=regions,
        ax=ax_heat,
        cmap=colormap,
        linecolor="#333333", 
        linewidths=0.5,
        cbar_ax=cbar_ax, 
        cbar_kws={'label': 'Scaled Co-localization'}
    )
    
    ax_heat.set_title(f"Dynamic Co-localization Profile of {target_cell}", color='white', pad=25, fontsize=20, fontweight='bold')
    ax_heat.tick_params(axis='x', rotation=45, labelsize=12)
    ax_heat.tick_params(axis='y', rotation=0, labelsize=18)
    
    for tick_label in ax_heat.get_xticklabels():
        cell_name = tick_label.get_text()
        macro_grp = get_macro_group(cell_name)
        tick_label.set_ha('right')
        tick_label.set_color("white")
        tick_label.set_bbox(dict(facecolor=bubble_colors[macro_grp], edgecolor='none', alpha=0.9, boxstyle='round,pad=0.25'))

    plt.tight_layout(rect=[0, 0, 1, 0.95]) 
    save_path = os.path.join(output_dir, f"Nature_Style_Networks_V5_{TEST_KEY}_{target_cell}.pdf")
    plt.savefig(save_path, facecolor='black', bbox_inches='tight')
    plt.savefig(save_path.replace(".pdf", ".png"), facecolor='black', bbox_inches='tight', dpi=400)
    plt.show()

    print(f"\n🎉 终极排版定稿完成！完美对齐、彻底防重叠版组图已保存至: {save_path}")



############################关键细胞通讯热图
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
# 0. Spyder 环境设置：强制内联显示
# ==========================================
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('matplotlib', 'inline')
except ImportError:
    pass

# 全局排版配置 (黑底白字，使用 Arial)
plt.style.use('dark_background')
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 
                     'text.color': 'white', 'axes.labelcolor': 'white', 
                     'xtick.color': 'white', 'ytick.color': 'white', 'axes.edgecolor': 'white'})

# ==========================================
# 1. 路径配置与全局变量
# ==========================================
output_dir = r"data/external/local_8e645d0a9b03"
TEST_KEY = "19dpb"  
TARGET_SAMPLE = "BW15D_C6_SDSDB_19dpb"
target_cell = "Fib_K14"
ROTATION_ANGLE = -225

file_path = os.path.join(output_dir, f"{TARGET_SAMPLE}_Annotated.h5ad")
json_load_path = os.path.join(output_dir, f"{TEST_KEY}_spatial_anchors.json")

# 热图配色保持 Zissou 渐变色
zissou = ["#3A9AB2", "#6FB2C1", "#91BAB6", "#A5C2A3", "#BDC881", "#DCCB4E", "#E3B710", "#EC7A05", "#F11B00"]
colormap = clr.LinearSegmentedColormap.from_list("Zissou", zissou)

# ==========================================
# 2. 空间坐标投影与显著性计算模块 (用于寻找需要加红框的基因)
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
    
    # 删除了 IGFLR1
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

    # 不进行空间过滤，只要表达就全量保留
    final_target_genes = []
    for gene in all_lr_candidates:
        if gene in adata_var_upper:
            if gene in expressed_genes or gene == 'IGFL1': 
                final_target_genes.append(gene)

    final_target_genes = list(dict.fromkeys(final_target_genes))

    # 🌟 优化 (5)：进一步缩小图形宽度，拉长高度 (5.5, 9.5)
    fig, axes = plt.subplots(1, 2, figsize=(5.5, 9.5), dpi=300, layout="constrained")
    fig.suptitle("19dpi | Spatial Dynamics of Receptor-Ligand Pairs", fontsize=15, fontweight='bold', color='white', y=1.03)
    
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
        
        # 🌟 优化 (3)：放大表达了的 Colorbar
        sns.heatmap(
            zscore_df.T, cmap=colormap, center=0, vmin=-2.0, vmax=2.0, 
            xticklabels=False, yticklabels=True, 
            cbar_kws={'shrink': 0.8, 'aspect': 15, 'location': 'top', 'pad': 0.02}, ax=ax_obj
        )
        
        raw_labels = zscore_df.columns.tolist()
        ax_obj.set_yticklabels(raw_labels, rotation=0)
        
        # 🌟 优化 (1) & (2)：颜色统一，显著变化基因与 IGFL1 共用红框和字号放大 🌟
        for tick, label_obj in zip(ax_obj.yaxis.get_major_ticks(), ax_obj.get_yticklabels()):
            gene_name = label_obj.get_text()
            
            tick.tick1line.set_color('white')
            tick.tick1line.set_markersize(3.0)
            
            target_cat = gene_to_category.get(gene_name, None)
            target_color = category_colors.get(target_cat, '#E8E8E8')

            # 所有基因都使用对应的分类颜色
            label_obj.set_color(target_color)
            label_obj.set_style('italic')
            
            # 判断是否为显著基因或者是 IGFL1
            is_highlight = (gene_name in sig_all_spatial) or (gene_name == 'IGFL1')
            
            if is_highlight:
                label_obj.set_fontweight('bold')
                label_obj.set_fontsize(14)      # 放大字号
                # 添加红色外框
                label_obj.set_bbox(dict(
                    facecolor='none',           
                    edgecolor='#FF3030',        
                    linewidth=1.5,              
                    boxstyle='round,pad=0.2'    
                ))
            else:
                label_obj.set_fontweight('bold') 
                label_obj.set_fontsize(10)      # 正常字号
                label_obj.set_bbox(None)        # 确保其他没有外框
                
        # 排版整体字号放大
        ax_obj.set_title(title, fontsize=13, color='white', pad=25) 
        ax_obj.set_ylabel("Receptors & Ligands" if ax_obj == axes[0] else "", fontsize=13, color='white', labelpad=6)
        
        # 🌟 优化 (4)：X轴底部文字分两行展示防重叠 🌟
        # 使用 labelpad=22 将中间的标题向下推到第二行
        ax_obj.set_xlabel(xlabel, fontsize=12, color='white', labelpad=22)
        # 将左/右提示语挂在略高一点的第一行
        ax_obj.text(0, -0.015, left_label, transform=ax_obj.transAxes, ha='left', va='top', fontsize=11, color='white')
        ax_obj.text(1, -0.015, right_label, transform=ax_obj.transAxes, ha='right', va='top', fontsize=11, color='white')
        
        # 色条刻度字号放大
        cbar = ax_obj.collections[0].colorbar
        cbar.ax.tick_params(labelsize=10, length=5, width=1.0, colors='white')
        
        for _, spine in ax_obj.spines.items():
            spine.set_visible(True); spine.set_color('white'); spine.set_linewidth(0.8)

    # 绘制双轴
    process_and_plot_lr(axes[0], 'Distance_to_Front_X', False, "Healing Axis", "Distance to Front", "Un-epi", "Newly-epi")
    process_and_plot_lr(axes[1], 'Relative_Depth_Y', True, "Depth Axis", "Relative Depth", "Deep Dermis", "Epidermis")

    # 🌟 优化 (3)：放大所有图例的字体及色块 🌟
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
        fontsize=12,                  # 图例字号放大至 12
        columnspacing=1.8,            
        handletextpad=0.8,
        handleheight=1.5,             # 图例色块放大
        handlelength=2.0
    )
    for text in legend.get_texts():
        text.set_color('white')

    # 保存文件
    save_path = os.path.join(out_dir, "19dpi_Nature_LR_Categorized_Highlight")
    fig.savefig(f"{save_path}.pdf", facecolor='black', bbox_inches='tight')
    fig.savefig(f"{save_path}.png", facecolor='black', bbox_inches='tight', dpi=600)
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
    
    # 重新运算差异，仅用于提供需要加红框的高亮名单
    sig_genes_x, sig_genes_y = get_spatially_variable_genes(fib_adata, rho_threshold=0.10, p_threshold=0.05)
    
    print("🚀 开始执行绘制指令，请留意 Spyder 右上角的 'Plots' 面板...")
    plot_categorized_lr_dynamics(
        adata=fib_adata, 
        sig_x_df=sig_genes_x,
        sig_y_df=sig_genes_y,
        out_dir=output_dir,      
        window_size=50           
    )
    
    

#######################IGFL1特征分析
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ==========================================
# 0. 环境与主题配置 (黑底白字)
# ==========================================
plt.style.use('dark_background')
plt.rcParams.update({
    'font.family': 'sans-serif', 
    'font.sans-serif': ['Arial'], 
    'text.color': 'white', 
    'axes.labelcolor': 'white', 
    'xtick.color': 'white', 
    'ytick.color': 'white',
    'axes.edgecolor': 'white'
})

# ==========================================
# 1. 颜色字典与分类逻辑
# ==========================================
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
fig, ax = plt.subplots(figsize=(14, 4)) # 稍微加宽加高以容纳带有底色的 45 度标签

# 1. 绘制底层小提琴
sns.violinplot(
    data=plot_df,
    x='CellType',
    y='Expression',
    hue='CellType',        
    order=sorted_celltypes,
    palette=palette,
    scale='width',         # 让所有小提琴宽度一致
    inner=None,            # 移除内部箱线图
    linewidth=0.8,
    legend=False
    # 🌟 优化：移除了 cut=0，让最高表达量 (KC_Basal_Mig) 的顶部能够自然收尖
)

# 🌟 优化：删除了 sns.stripplot 模块，彻底去除震动的散点，只保留纯净的小提琴

plt.title(f"{TARGET_GENE} Expression Across Cell Types", fontsize=22, fontweight='bold', pad=20)
plt.xlabel("Cell Type Macro-groups", fontsize=16, labelpad=15)
plt.ylabel(f"{TARGET_GENE} Expression Level", fontsize=16, labelpad=15)

# 🌟 核心优化：通过限制 Y 轴下限，在视觉上“砍掉”负数的半截小提琴
# 设置下限为 -0.2 (留一点点呼吸空间，让底座看起来像自然的椭圆而不是死板的直线)
y_max = plot_df['Expression'].max()
ax.set_ylim(-0.2, y_max * 1.05) 

# X轴标签倾斜 45 度，并添加分类底色
ax.tick_params(axis='x', rotation=45)
ax.tick_params(axis='y', labelsize=14)

for tick_label in ax.get_xticklabels():
    cell_name = tick_label.get_text()
    macro_grp = get_macro_group(cell_name)
    tick_label.set_ha('right') # 确保 45度角时右上角对齐刻度线
    tick_label.set_color("white")
    tick_label.set_fontsize(13)
    # 添加带有倒角的类群底色框
    tick_label.set_bbox(dict(facecolor=bubble_colors[macro_grp], edgecolor='none', alpha=0.9, boxstyle='round,pad=0.25'))

# 移除冗余边框并自适应布局
sns.despine()
plt.tight_layout()

# 保存与展示
# plt.savefig(f"Seurat_VlnPlot_Clean_{TARGET_GENE}.pdf", facecolor='black', bbox_inches='tight')
plt.show()

##########################IGFL1空间共定位分析
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
# 0. Global Setting (Black background & Arial)
# ==========================================
plt.style.use('dark_background')
plt.rcParams.update({
    'font.family': 'sans-serif', 
    'font.sans-serif': ['Arial'],
    'pdf.fonttype': 42
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
    
    # C. Layout: 嵌套 GridSpec 确保上下行独立划分，并拉近间距
    fig = plt.figure(figsize=(24, 9), dpi=300)
    
    # (3) 缩小两行图片之间的距离：拉低顶部的 bottom，抬高底部的 top
    gs_top = gridspec.GridSpec(1, 3, bottom=0.52, top=0.92, wspace=0.15)
    ax1 = fig.add_subplot(gs_top[0])  # Top Left
    ax2 = fig.add_subplot(gs_top[1])  # Top Middle
    ax3 = fig.add_subplot(gs_top[2])  # Top Right
    
    gs_bot = gridspec.GridSpec(1, 4, bottom=0.08, top=0.46, wspace=0.25)
    ax4 = fig.add_subplot(gs_bot[0:2])  # Bottom Left Merge (占 2/4 宽度，即一半)
    ax5 = fig.add_subplot(gs_bot[2])    # Bottom Barplot (占 1/4 宽度)
    # gs_bot[3] 故意留空，满足右下角无图需求
    
    # (1) 更新为您指定的细胞标准颜色
    c_bg = '#222222'
    c_red = '#E53935'
    c_yellow = '#FDD835'
    c_ct1 = '#FF3030'    # KC_Basal_Mig
    c_ct2 = '#1A0099FF'  # KC_Spinous_Mig
    c_ct3 = '#4682B4'    # Fib_K14
    
    palette_regions = ['#E53935', '#FDD835', '#1E90FF'] 
    
    def format_sp_ax(ax):
        ax.set_aspect('equal')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(x.min() - 100, x.max() + 100)
        ax.set_ylim(y.min() - 100, y.max() + 100)
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
    ax1.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.3, marker='.', zorder=0)
    ax1.scatter(x[dist_C1 > 0], y[dist_C1 > 0], c=c_ct1, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax1, cell_type_1, c_ct1)
    
    format_sp_ax(ax2)
    ax2.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.3, marker='.', zorder=0)
    ax2.scatter(x[dist_C2 > 0], y[dist_C2 > 0], c=c_ct2, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax2, cell_type_2, c_ct2)
    
    format_sp_ax(ax3)
    ax3.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.3, marker='.', zorder=0)
    ax3.scatter(x[dist_C3 > 0], y[dist_C3 > 0], c=c_ct3, s=pt_size*1.5, zorder=1)
    add_top_center_title(ax3, cell_type_3, c_ct3)
    
    # ==========================================
    # Row 2 Left: Ligand-Receptor Merge (占 1/2 宽度)
    # ==========================================
    format_sp_ax(ax4)
    ax4.scatter(x, y, c=c_bg, s=pt_size*0.1, alpha=0.3, marker='.', zorder=0)
    ax4.scatter(x[expr_R > 0], y[expr_R > 0], c=c_yellow, s=pt_size*1.2, alpha=0.8, zorder=1)
    ax4.scatter(x[expr_L > 0], y[expr_L > 0], c=c_red, s=pt_size*1.2, alpha=0.8, zorder=2)
    add_top_center_title(ax4, ligand, c_red, receptor, c_yellow)
    
    ax4.plot([x.max()-2000, x.max()-500], [y.min()+500, y.min()+500], color='white', lw=1.5)
    ax4.text(x.max()-1250, y.min()+700, "2mm", color='white', ha='center', fontsize=8)
    
    # ==========================================
    # Row 2 Middle: Barplot (IGFL1 Expression) (占 1/4 宽度)
    # ==========================================
    # (4) 使用真实的 Mask 图像映射读取坐标分区
    mask_path = os.path.join(out_dir, "wound_zones_mask.png")
    if os.path.exists(mask_path):
        print("🗺️ 成功加载 Mask 文件，正在进行空间映射...")
        mask_img = iio.imread(mask_path)
        if len(mask_img.shape) == 3: mask_img = mask_img[:, :, 0] # 确保单通道
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
    x_order = [cell_type_1, cell_type_2] # X轴仅保留这2种细胞
    
    # 过滤数据并提取有效细胞的表达量
    mask_cells = np.isin(cell_labels, x_order) & adata.obs['Plot_Group'].notna()
    df_bar = pd.DataFrame({
        'Expression': expr_L[mask_cells],
        'CellType': cell_labels[mask_cells],
        'Group': adata.obs['Plot_Group'].values[mask_cells]
    })
    
    # 绘制 Barplot
    sns.barplot(data=df_bar, x='CellType', y='Expression', hue='Group', 
                order=x_order, hue_order=hue_order, palette=palette_regions, 
                ax=ax5, capsize=.08, errwidth=1.2, edgecolor='white')
    
    ax5.set_facecolor('black')
    ax5.set_ylabel(f"Mean Expr. of {ligand}", fontsize=11, color='white')
    ax5.set_xlabel("")
    
    # X轴细胞名倾斜 30 度
    ax5.set_xticks(range(len(x_order)))
    ax5.set_xticklabels(x_order, rotation=30, ha='right', fontsize=11, color='white')
    ax5.tick_params(axis='y', which='major', labelsize=10, colors='white')
    
    for s in ax5.spines.values(): s.set_color('white')
    
    # (2) 压缩 barplot 的高度
    ax5.set_box_aspect(0.65)
    
    # 图例居中无点
    ax5.legend(title="", frameon=False, labelcolor='white', fontsize=9, 
               loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=3, columnspacing=0.6)
    
    # ==========================================
    # 组内显著性检验 (Asterisk Annotation)
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
                ax5.plot([x1, x1, x2, x2], [y_pos, y_pos+step*0.2, y_pos+step*0.2, y_pos], lw=1, c='white')
                ax5.text((x1+x2)/2, y_pos+step*0.2, star, ha='center', va='bottom', color='white', fontsize=12)
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
    fig.suptitle("Spatial Co-localization Analysis", fontsize=24, color='white', fontweight='bold', y=1.0)
    
    save_fn = os.path.join(out_dir, f"Final_NestedLayout_MaskedBar_{sample_key}.png")
    plt.savefig(save_fn, facecolor='black', bbox_inches='tight', dpi=300)
    plt.show()

# ==========================================
# 2. Execution
# ==========================================
plot_optimized_colocalization(adata, "19dpb", out_dir=output_dir)







