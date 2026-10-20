# ===================== visium =====================
import os
import gc
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
import spatialdata as sd
import spatialdata_io
import harmonypy as hm
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

# ==========================================
# 1. 路径与环境配置
# ==========================================
plt.style.use("dark_background")
base_dir = r"data/external/local_137b74f92ae2"
samples = ["BW9", "X002"]
adata_list = []


# ==========================================
# 2. 循环读取与预处理 (手动组装版)
# ==========================================
import cv2
import numpy as np

for s_name in samples:
    sample_path = os.path.join(base_dir, s_name)
    matrix_path = os.path.join(sample_path, "filtered_feature_bc_matrix")
    spatial_path = os.path.join(sample_path, "spatial")
    
    print(f"📖 正在手动读取 Visium 样本及背景图: {s_name}")
    
    try:
        # 1. 读取表达矩阵
        temp_ad = sc.read_10x_mtx(matrix_path, var_names='gene_symbols', cache=False)
        temp_ad.var_names_make_unique()
        
        # 2. 读取坐标
        pos_file = os.path.join(spatial_path, "tissue_positions.csv")
        if not os.path.exists(pos_file):
            pos_file = os.path.join(spatial_path, "tissue_positions_list.csv")
        positions = pd.read_csv(pos_file, header=None, index_col=0)
        temp_ad.obsm['spatial'] = positions.loc[temp_ad.obs_names, [4, 5]].values.astype(float)
        
        # 3. 增强型图片读取 (解决中文/长路径报错)
        img_path = os.path.join(spatial_path, "tissue_hires_image.png")
        if os.path.exists(img_path):
            # 关键改动：先用 numpy 读取二进制流，再用 cv2 解码
            img_array = np.fromfile(img_path, dtype=np.uint8)
            bgr_img = cv2.imdecode(img_array, cv2.IMREAD_COLOR)
            
            if bgr_img is not None:
                rgb_img = cv2.cvtColor(bgr_img, cv2.COLOR_BGR2RGB)
                temp_ad.uns['spatial'] = {
                    s_name: {
                        'images': {'lowres': rgb_img / 255.0},
                        'scalefactors': {
                            'tissue_lowres_scalefactor': 1.0, 
                            'spot_diameter_fullres': 30.0
                        }
                    }
                }
                print(f"🖼️  {s_name} 背景图成功加载并校色")
            else:
                print(f"⚠️  {s_name} 图片解码失败，请确认文件是否损坏")
        else:
            print(f"⚠️  {s_name} 路径下未找到 tissue_lowres_image.png")
        
        # 4. 标记样本 ID
        temp_ad.obs["sample_batch_new"] = s_name
        temp_ad.obs_names = f"{s_name}_" + temp_ad.obs_names
        
        adata_list.append(temp_ad)
        print(f"✅ {s_name} 组装成功，点数: {temp_ad.n_obs}")
        
    except Exception as e:
        print(f"❌ 跳过样本 {s_name}，读取失败: {e}")
    
    gc.collect()

# ==========================================
# 3. 数据合并与初步保存
# ==========================================
if len(adata_list) > 0:
    print(f"🚀 正在合并 {len(adata_list)} 个样本...")
    # join='inner' 保证基因集对齐，uns_merge="unique" 确保空间信息不冲突
    adata = ad.concat(adata_list, join='inner', uns_merge="unique")
    del adata_list
    gc.collect()

    # 设定输出路径 (适配您的 F 盘目录)
    output_dir = r"data/external/local_605e2dec15c8"
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    raw_h5ad = os.path.join(output_dir, "combined_raw_spatial.h5ad")
    pdf_path = os.path.join(output_dir, "Spatial_Harmony_Report.pdf")
    
    print(f"💾 正在保存原始合并数据至: {raw_h5ad}")
    adata.write_h5ad(raw_h5ad)

    # ==========================================
    # 4. 全局视觉设定 (黑底 Arial 风格)
    # ==========================================
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'figure.facecolor': 'black', 'axes.facecolor': 'black', 'axes.edgecolor': 'white',
        'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white',
        'text.color': 'white', 'legend.facecolor': 'black', 'legend.edgecolor': 'white'
    })
    sc.settings.set_figure_params(dpi=300, facecolor='black', fontsize=14, color_map='magma')

    # ==========================================
    # 5. QC 与 标准化
    # ==========================================
    pdf_path = r"data/external/local_b136a98947b5"
    pdf = PdfPages(pdf_path)
    
    print("🧹 正在进行 QC 计算...")
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], log1p=False, inplace=True)
    
    # 绘制 QC 小提琴图
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], 
                 groupby="sample_batch_new", rotation=90, multi_panel=True, show=False)
    plt.suptitle("Quality Control Metrics", color='white', y=1.05)
    pdf.savefig(plt.gcf(), bbox_inches='tight', facecolor='black')
    plt.close()

    print("🧬 正在进行标准化与高变基因筛选...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="sample_batch_new")
    sc.tl.pca(adata, svd_solver='arpack')

    # ==========================================
    # 6. Harmony 批次校正 (核心环节)
    # ==========================================
    print("🤝 正在运行 Harmony 整合...")
    X_pca = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    ho = hm.run_harmony(X_pca, adata.obs, 'sample_batch_new')
    
    # 处理 80 万细胞可能的非标准数组结构
    z_corr_raw = np.array(ho.Z_corr)
    if z_corr_raw.shape == (50, adata.n_obs):
        adata.obsm["X_pca_harmony"] = z_corr_raw.T
    elif z_corr_raw.shape == (adata.n_obs, 50):
        adata.obsm["X_pca_harmony"] = z_corr_raw
    else:
        adata.obsm["X_pca_harmony"] = np.stack(ho.Z_corr).T
    
    print(f"✅ Harmony 校正完毕，形状: {adata.obsm['X_pca_harmony'].shape}")
    del X_pca, ho
    gc.collect()

    # ==========================================
    # 7. 聚类与分样本可视化
    # ==========================================
    print("🗺️ 正在构建 UMAP 并聚类 (Leiden)...")
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony")
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=1.5, key_added="clusters")

    # 绘制分样本 UMAP
    samples_list = adata.obs["sample_batch_new"].unique().tolist()
    fig, axes = plt.subplots(1, len(samples_list), figsize=(6*len(samples_list), 5), facecolor='black')
    for i, s in enumerate(samples_list):
        sc.pl.umap(adata[adata.obs["sample_batch_new"] == s], color="clusters", ax=axes[i], 
                   show=False, title=f"UMAP: {s}", frameon=False, size=2)
    pdf.savefig(fig, bbox_inches='tight', facecolor='black')
    plt.close()

    # 绘制空间坐标分布 (Spatial Plot)
    fig, axes = plt.subplots(1, len(samples_list), figsize=(6*len(samples_list), 5), facecolor='black')
    for i, s in enumerate(samples_list):
        # 使用之前手动读入的 obsm['spatial']
        sc.pl.embedding(adata[adata.obs["sample_batch_new"] == s], basis="spatial", 
                        color="clusters", ax=axes[i], show=False, title=f"Spatial: {s}", frameon=False, size=1.5)
        axes[i].set_aspect('equal')
    pdf.savefig(fig, bbox_inches='tight', facecolor='black')
    plt.close()

    # 保存最终聚类数据
    final_h5ad = os.path.join(output_dir, "combined_clustered_800k.h5ad")
    adata.write_h5ad(final_h5ad)
    pdf.close()
    print(f"✨ 报告已生成: {pdf_path}")
    print(f"💾 聚类数据已保存: {final_h5ad}")

############################
# 1. 自动计算或手动设定缩放因子
# 对指定的样本进行中心对称处理
for s in ["BW9", "X002"]:
    if s in adata.uns['spatial']:
        img = adata.uns['spatial'][s]['images']['lowres']      
        # --- 方案 A：如果是“上下反了” (垂直翻转) ---
        # 使用 flipCode=0
        img_corrected = cv2.flip(img, 1) 
        #img_corrected = cv2.rotate(img, cv2.ROTATE_180)
        #img_corrected = cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)
        #img_corrected = cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
        # --- 方案 B：如果是“左右反了” (水平翻转) ---
        # 如果方案 A 不对，请注释掉上面一行，取消下面这行的注释
        # img_corrected = cv2.flip(img, 1) 
        # 将修正后的图片写回
        adata.uns['spatial'][s]['images']['lowres'] = img_corrected
        print(f"🔄 样本 {s} 已完成镜像翻转修正")

sample_to_fix = "BW9"
if sample_to_fix in adata.uns['spatial']:
    # 设定平移量（根据 image_98dd66 观察，点位需要向右下移动）
    # delta_x > 0 向右移；delta_y > 0 向下移
    # 建议先试：右移 150，下移 150 (单位为 fullres 像素)
    delta_x = -50 
    #delta_y = 150 
    
    # 仅针对 BW9 样本的坐标进行平移
    mask = adata.obs['sample_batch_new'] == sample_to_fix
    adata.obsm['spatial'][mask, 0] += delta_x  # 修改 X 坐标
    #adata.obsm['spatial'][mask, 1] += delta_y  # 修改 Y 坐标

# ------------------------------------------
# 重新绘图验证
# ------------------------------------------
sc.pl.spatial(
    adata[adata.obs['sample_batch_new'] == "BW9"], 
    library_id="BW9", 
    color="clusters", 
    spot_size=180, 
    alpha_img=0.7, 
    img_key="lowres",
    title="BW9: Epidermis-Dermis Realignment",
    show=True
)

# 验证 KRT14 是否精准覆盖在 HE 的表皮层上
sc.pl.spatial(
    adata[adata.obs['sample_batch_new'] == "BW9"], 
    library_id="BW9", 
    color="KRT14",  # 换成 KRT14 观察
    spot_size=180, 
    img_key="lowres",
    color_map="magma"
)


#################BW9
# sc.pl.spatial(adata, library_id="BW9", color="clusters",  spot_size=150, img_key="lowres", alpha_img=0.7)
# 同理绘制 X002
for s in ["X002"]:
    if s in adata.uns['spatial']:
        img = adata.uns['spatial'][s]['images']['lowres']      
        # --- 方案 A：如果是“上下反了” (垂直翻转) ---
        # 使用 flipCode=0
        img_corrected = cv2.flip(img, 1) 
        #img_corrected = cv2.rotate(img, cv2.ROTATE_180)
        #img_corrected = cv2.rotate(img, cv2.ROTATE_90_COUNTERCLOCKWISE)
        #img_corrected = cv2.rotate(img, cv2.ROTATE_90_CLOCKWISE)
        # --- 方案 B：如果是“左右反了” (水平翻转) ---
        # 如果方案 A 不对，请注释掉上面一行，取消下面这行的注释
        # img_corrected = cv2.flip(img, 1) 
        # 将修正后的图片写回
        adata.uns['spatial'][s]['images']['lowres'] = img_corrected
        print(f"🔄 样本 {s} 已完成镜像翻转修正")

sample = "X002"
if sample in adata.uns['spatial']:
    # 1. 调整缩放因子 (Scale Factor)
    # 增加该值会让点阵整体“变大”，请根据 image_9955a4 效果微调
    # 如果点阵还是比组织小，把 0.12 改成 0.13 或更高
    adata.uns['spatial'][sample]['scalefactors']['tissue_lowres_scalef'] = 0.10
    
    # 2. 坐标平移 (Translation)
    # 根据您的截图，点阵需要向右下移动
    # 建议值：X 向右移 120，Y 向下移 200 (单位为 fullres 像素)
    delta_x = 50 
    delta_y = 0 
    
    mask = adata.obs['sample_batch_new'] == sample
    adata.obsm['spatial'][mask, 0] += delta_x
    adata.obsm['spatial'][mask, 1] += delta_y


sc.pl.spatial(
    adata[adata.obs['sample_batch_new'] == "X002"], 
    library_id="X002", 
    color="clusters", 
    spot_size=180, 
    alpha_img=0.7, 
    img_key="lowres",
    title="X002: Clusters on HE",
    show=True
)


sample_to_fix = "X002"
if sample_to_fix in adata.uns['spatial']:
    # 设定平移量（根据 image_98dd66 观察，点位需要向右下移动）
    # delta_x > 0 向右移；delta_y > 0 向下移
    # 建议先试：右移 150，下移 150 (单位为 fullres 像素)
    delta_x = 150 
    delta_y = 150 
    
    # 仅针对 BW9 样本的坐标进行平移
    mask = adata.obs['sample_batch_new'] == sample_to_fix
    adata.obsm['spatial'][mask, 0] += delta_x  # 修改 X 坐标
    adata.obsm['spatial'][mask, 1] += delta_y  # 修改 Y 坐标

######################保存所有调整结果
import os
from datetime import datetime

# 1. 设置保存路径
output_dir = r"data/external/local_4e259f9e9282"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 2. 生成文件名（包含日期，防止覆盖旧版本）
current_time = datetime.now().strftime("%Y%m%d")
final_file = os.path.join(output_dir, f"Combined_Aligned_800k_{current_time}.h5ad")

# 3. 完整保存（包含 BW9 和 X002 的校准结果）
# 使用 gzip 压缩以节省空间，毕竟 800k 数据体积较大
adata.write_h5ad(final_file, compression='gzip')

print(f"✅ 成功保存合并后的全量校准数据：\n{final_file}")

#############################
import scanpy as sc

# 1. 重新读取
print("⏳ 正在重新加载数据，请稍候...")
adata = sc.read_h5ad(final_file)

# 2. 验证关键对齐参数是否存在
samples = ["BW9", "X002"]
for s in samples:
    if s in adata.uns['spatial']:
        scalef = adata.uns['spatial'][s]['scalefactors']['tissue_lowres_scalef']
        print(f"🔍 检查样本 {s}: 缩放因子为 {scalef:.4f}，坐标矩阵形状 {adata[adata.obs['sample_batch_new']==s].obsm['spatial'].shape}")

# 3. 快速可视化确认 (选一个样本看看 KRT14 是否还在原位)
sc.pl.spatial(
    adata[adata.obs['sample_batch_new'] == "X002"], 
    library_id="X002", 
    color="KRT1", 
    spot_size=180, 
    img_key="lowres",
    title="Final Verification: X002 Alignment",
    show=True
)

sc.pl.spatial(
    adata[adata.obs['sample_batch_new'] == "BW9"], 
    library_id="BW9", 
    color="KLF4", 
    spot_size=180, 
    img_key="lowres",
    title="Final Verification: X002 Alignment",
    show=True
)

##############################可视化
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import scanpy as sc
# ==========================================
# BW9 & X002 愈合 Marker 2x5 十 Panel 展示
# ==========================================

# 1. 配置参数
# 我们分别展示：BW9 (第一行), X002 (第二行)
sample_ids = ["BW9", "X002"] 
# 保持 5 个关键 Marker
genes_to_plot = ['KLF4', 'PDGFRA', 'IGFL1', 'KRT14', 'TP63']

# 使用您调试好的最佳参数
spot_size_base = 150  
alpha_img = 0.7

# 2. 创建 2x5 画布
# 调大 figsize 以容纳 10 张图
fig = plt.figure(figsize=(30, 12)) 
# 2 行 (samples), 5 列 (genes)
gs = gridspec.GridSpec(2, 5) 

# 3. 循环绘图
for row_idx, sample in enumerate(sample_ids):
    # 提取特定样本的数据
    adata_sample = adata[adata.obs['sample_batch_new'] == sample]
    
    for col_idx, gene in enumerate(genes_to_plot):
        # 计算在 gs 中的索引
        ax = fig.add_subplot(gs[row_idx, col_idx])
        
        # 调用 scanpy pl.spatial
        # show=False 必须设为 False，否则无法整合
        # basis 必须选择默认的 spatial，而不是 spatial_aligned (因为我们已经微调过 scalefactors)
        sc.pl.spatial(
            adata_sample,
            library_id=sample,
            color=gene,
            spot_size=spot_size_base,
            alpha_img=alpha_img,
            img_key="lowres",
            color_map="magma",  # 使用 magma 色阶，黄白色表达最高，清晰
            frameon=False,
            show=False,
            ax=ax
        )
        
        # 4. 手动设置标题 (样本 + 基因)
        ax.set_title(f"{sample}: {gene}", fontsize=18, fontweight='bold')

# 5. 调整布局，防止子图重叠
plt.tight_layout(pad=3.0)

# --- 关键：仅展示，不在此步保存 ---
# 6. 在 Spyder 的 'Plots' 窗口弹窗展示
# 确保 plt.show() 在最后
plt.show()

print("✅ 2x5 十 Panel 空间 Feature Plot 已展示在 IPython Plots 窗格中。")

#######################差异表达分析
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import os

# 1. 差异表达分析 (DEA)
print("⏳ 正在进行差异表达分析 (Wilcoxon)...")
# 假设您之前的聚类结果存在 'clusters' 中
sc.tl.rank_genes_groups(adata, groupby='clusters', method='wilcoxon')

# 2. 输出 DEA 结果为 CSV
output_dir = r"data/external/local_48efc8310eb6"
csv_path = os.path.join(output_dir, "All_Clusters_Markers_Ranked.csv")
result_df = sc.get.rank_genes_groups_df(adata, group=None)
result_df.to_csv(csv_path, index=False)
print(f"✅ 差异基因列表已保存至: {csv_path}")

# 3. 准备 Marker 基因列表 (从 txt 中提取并去重)
# 这里仅列出您提到的关键基因，如需全量可使用：markers = open('merged_cleaned_genes.txt').read().splitlines()
target_markers = open(r'data/external/local_828e7870eb7d').read().splitlines()
# 过滤掉 adata 中不存在的基因
valid_markers = [g for g in target_markers if g in adata.var_names]

# 4. 生成合并 PDF (Violin + 2x FeaturePlots)
pdf_path = os.path.join(output_dir, "Cell_Annotation_Atlas_Compressed.pdf")

print(f"⏳ 正在生成 PDF 报告 (共 {len(valid_markers)} 页)...")
with PdfPages(pdf_path) as pdf:
    for gene in valid_markers:
        # 创建画布：1行3列的比例排版
        fig = plt.figure(figsize=(22, 7))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1, 1])
        
        # --- Panel 1: 小提琴图 (全数据集) ---
        ax0 = fig.add_subplot(gs[0])
        sc.pl.violin(adata, keys=gene, groupby='clusters', rotation=90, show=False, ax=ax0)
        ax0.set_title(f"{gene} Expression by Cluster", fontsize=14)
        
        # --- Panel 2: BW9 Featureplot ---
        ax1 = fig.add_subplot(gs[1])
        sc.pl.spatial(
            adata[adata.obs['sample_batch_new'] == "BW9"],
            library_id="BW9", color=gene, img_key="lowres",
            spot_size=150, alpha_img=0.5, show=False, ax=ax1,
            title="BW9: Spatial Distribution", frameon=False
        )
        # 关键优化：光栅化斑点图层以极大减小 PDF 体积
        ax1.collections[0].set_rasterized(True)
        
        # --- Panel 3: X002 Featureplot ---
        ax2 = fig.add_subplot(gs[2])
        sc.pl.spatial(
            adata[adata.obs['sample_batch_new'] == "X002"],
            library_id="X002", color=gene, img_key="lowres",
            spot_size=180, alpha_img=0.5, show=False, ax=ax2,
            title="X002: Spatial Distribution", frameon=False
        )
        ax2.collections[0].set_rasterized(True)
        
        fig.suptitle(f"Marker Analysis: {gene}", fontsize=20, fontweight='bold', y=1.05)
        plt.tight_layout()
        
        # 保存当前页并关闭画布释放内存
        pdf.savefig(fig, bbox_inches='tight', dpi=150) # 降低 DPI 以压缩体积
        plt.close(fig)

print(f"🎉 报告已生成！请查看: {pdf_path}")

##############
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages

# 1. 设置保存路径
output_dir = r"data/external/local_4e259f9e9282"
pdf_path = os.path.join(output_dir, "Clusters_Spatial_Verification.pdf")

# 2. 准备绘图参数
samples = ["BW9", "X002"]
# 使用您之前调试完美的参数
spot_sizes = {"BW9": 180, "X002": 180} 

print(f"⏳ 正在生成 Clusters 空间图报告...")

with PdfPages(pdf_path) as pdf:
    for s in samples:
        # 创建画布
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # 提取子集并绘图
        sc.pl.spatial(
            adata[adata.obs['sample_batch_new'] == s],
            library_id=s,
            color="clusters",
            spot_size=spot_sizes[s],
            alpha_img=0.7,
            img_key="lowres",
            frameon=False,
            show=False,
            ax=ax,
            title=f"Sample {s}: Cluster Distribution"
        )
        
        # --- 关键：光栅化 Spot 图层 (极速减小 PDF 体积并防止卡顿) ---
        # sc.pl.spatial 生成的第一个 collection 是 spots
        ax.collections[0].set_rasterized(True)
        
        # 保存当前样本到 PDF 页
        pdf.savefig(fig, bbox_inches='tight', dpi=300)
        plt.close(fig)

print(f"✅ 聚类结果 PDF 已保存至：\n{pdf_path}")













# 1. 设置保存路径 (适配您的 F 盘工作目录)
deg_output_path = r"data/external/local_cb1dd8a39884"

# 2. 提取所有 Cluster 的统计数据并合并
result_list = []
for cluster_id in adata.obs['clusters'].cat.categories:
    # 提取特定 Cluster 的数据
    df = sc.get.rank_genes_groups_df(adata, group=str(cluster_id))
    df['cluster'] = cluster_id  # 添加一列标记所属 Cluster
    result_list.append(df)

# 3. 合并为一个大表并导出
full_deg_table = pd.concat(result_list)

# 4. 过滤掉不显著的基因 (可选，例如只保留 adj_p < 0.05)
# full_deg_table = full_deg_table[full_deg_table['pvals_adj'] < 0.05]

# 5. 写入 CSV
full_deg_table.to_csv(deg_output_path, index=False, encoding='utf_8_sig')

print(f"✅ 差异基因分析结果已保存至: {deg_output_path}")

# 根据您的 DEG 结果更新下面的映射关系
annotation_dict = {
    "0": "Basal Keratinocytes", 
    "1": "Fibroblasts", 
    "2": "Immune Cells",
    # ... 请根据您 PDF 中的 22 个 Cluster 编号补全
}

adata.obs['cell_type'] = adata.obs['clusters'].map(annotation_dict).astype('category')
sc.pl.embedding(adata, basis="spatial", color="cell_type", title="Annotated Tissues")


# 1. 空间分布图 (分样本展示)
sc.pl.embedding(adata, basis="spatial", color=["KRT14", "PDGFRA"], 
                batch_key="sample_batch_new", ncols=2, size=2, color_map="magma")

# 2. 共表达逻辑：利用散点图观察二者相关性
sc.pl.scatter(adata, x="KRT14", y="PDGFRA", color="clusters", title="Co-expression Correlation")

# 3. 针对“上皮化”观察：KRT14 (上皮) 与 PDGFRA (间质/成纤维) 的交界情况
# 我们可以创建一个简单的共表达评分存入 obs
adata.obs['KRT14_PDGFRA_score'] = adata[:, "KRT14"].X.toarray().flatten() * adata[:, "PDGFRA"].X.toarray().flatten()
sc.pl.embedding(adata, basis="spatial", color="KRT14_PDGFRA_score", title="KRT14 & PDGFRA Interaction Zones")



