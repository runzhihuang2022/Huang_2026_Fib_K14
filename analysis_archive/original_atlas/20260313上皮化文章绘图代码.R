cd external_asset_d552fae5454f

R

library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(SingleR)
library(reticulate)
library(cowplot)
library(stringr)
library(RColorBrewer)
library(ggplotify)
library(cowplot)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(tidyverse)
library(patchwork)
library(reshape2)
library(scales)
library(ggpubr)
library(pheatmap)
library(dplyr)
library(Seurat)
library(patchwork)
library(Nebulosa)
library(SCpubr)
library(ggplot2)
library(dplyr)
library(tidydr)
library(harmony)
library(dplyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(Matrix)
library(future)
library(R.utils)
library(pdftools)
library(SeuratDisk)  # 用于 .h5ad 文件读取
library(reticulate)  # 用于 Python 包支持

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

setwd("data/external/local_d552fae5454f")

pbmc = readRDS("data/external/local_d4358e4d1be9")
#saveRDS(pbmc, file = "pbmc_final.rds")

head(pbmc@meta.data)
unique(pbmc@meta.data$sub_labels)
unique(pbmc@meta.data$Grade)

# ==============================================================================
# 1. 定义排序后的 Label 列表 (已根据最新字典对齐)
# ==============================================================================

bulk_labels_sorted <- c(
  "Keratinocyte", "SAC_SG", "SAC_HF", "Melanocyte", "Fibroblast",
  "Endothelial", "Pericyte", "Neurocyte", "Mono_Macro", "Neutrophil",
  "Mast", "DC", "LGHS", "B", "CD4_T", "CD8_T", "NK"
)

sub_labels_sorted <- c(
  # --- 1. 上皮细胞 (Keratinocytes) ---
  "KC_Basal", "KC_Basal_Mig","KC_Basal_Prolif", "KC_Spinous",
  "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
  # --- 2. 汗腺 (Sweat Glands) ---
  "SAC_SG_Progenitor","SAC_SG_Clear","SAC_SG_Dark", "SAC_SG_Ductal",
  # --- 3. 毛囊与色素 (Hair Follicle & Melanocytes) ---
  "SAC_HF_HFSC", "SAC_HF_DP_DS", "SAC_HF_Matrix", "SAC_HF_ORS","SAC_HF_IRS", "Melanocyte",
  # --- 4. 间充质/成纤维细胞 (Fibroblasts) ---
  "Fib_Papi", "Fib_SFRP2", "Fib_Fasci", "Fib_EN1", "Fib_Inflama", "Fib_CD45", "Fib_Myo", "Fib_Prolif", "Fib_K14",
  # --- 5. 内皮与壁细胞 (Endothelial & Mural) ---
  "Endo_Arterial", "Endo_Venous","Endo_Lymphatic","Endo_Capillary", "Endo_Prolif", "Pericyte",
  # --- 6. 神经相关 (Neural) ---
  "Schwann", 
  # --- 7. 髓系免疫细胞 (Myeloid) ---
  "M1_Macrophage", "M2_Macrophage","LAM","Neutrophil",
  # --- 8. 树突状细胞与肥大细胞 (DC & Mast) ---
  "Mast","cDC1","cDC2","pDC","LGHS",
  # --- 9. B细胞谱系 (B Cells) ---
  "B_Mature","B_Plasma",
  # --- 10. T细胞与NK (T & NK Cells) ---
  "CD4_T","CD8_T","NK"
)

# ==============================================================================
# 2. 定义针对性颜色向量 (强制锁定字典配色)
# ==============================================================================

# --- Sub Labels 颜色映射 ---
celltype_colors_sub_labels <- c(
  # KC
  "KC_Basal"="#008941", "KC_Basal_Mig"="#FF3030", "KC_Basal_Prolif"="#FFD700", 
  "KC_Spinous"="#00FF00FF", "KC_Spinous_Mig"="#1A0099FF", "KC_Spinous_Mat"="#EE4C97FF", "KC_Granular"="#882D17",
  # SAC_SG
  "SAC_SG_Progenitor"="#0099CCFF", "SAC_SG_Clear"="#B7E4F9FF", "SAC_SG_Dark"="#193006FF", "SAC_SG_Ductal"="#8F7700FF",
  # SAC_HF & Melanocyte
  "SAC_HF_HFSC"="#00F5FF", "SAC_HF_DP_DS"="#CDBE70", "SAC_HF_Matrix"="#7FFFD4", "SAC_HF_ORS"="#A30059", "SAC_HF_IRS"="#CC9900FF", "Melanocyte"="#FF4500",
  # Fib
  "Fib_Papi"="#3DA873FF", "Fib_SFRP2"="#CCFF00FF", "Fib_Fasci"="#9370DB", "Fib_EN1"="#BDB76B", "Fib_Inflama"="#800000FF", "Fib_CD45"="#DA70D6", "Fib_Myo"="#B79762", "Fib_Prolif"="#EED2EE", "Fib_K14"="#4682B4",
  # Endo & Pericyte
  "Endo_Arterial"="#FF00FF", "Endo_Venous"="#4169E1", "Endo_Lymphatic"="#ADFF2F", "Endo_Capillary"="#CF4E9CFF", "Endo_Prolif"="#F08080", "Pericyte"="#FFFF00",
  # Neuro
  "Schwann"="#7CFC00",
  # Myeloid
  "M1_Macrophage"="#FF4040", "M2_Macrophage"="#54FF9F", "LAM"="#8B4513", "Neutrophil"="#FF7F50",
  # DC & Mast
  "Mast"="#FF1493", "cDC1"="#9B30FF", "cDC2"="#EE9A00", "pDC"="#D1EEEE", "LGHS"="#FFB5C5",
  # Lymphoid
  "B_Mature"="#0000CD", "B_Plasma"="#FAF0E6", "CD4_T"="#00BFFF", "CD8_T"="#1E90FF", "NK"="#00E5EE"
)

# --- Bulk Labels 颜色映射 (大类简化版) ---
celltype_colors_bulk_labels <- c(
  "Keratinocyte"="#008941", "SAC_SG"="#0099CCFF", "SAC_HF"="#00F5FF", "Melanocyte"="#FF4500",
  "Fibroblast"="#DCD300", "Endothelial"="#AD002AFF", "Pericyte"="#FFFF00", "Neurocyte"="#F38400",
  "Mono_Macro"="#5A0007", "Neutrophil"="#FF7F50", "Mast"="#FF1493", "DC"="#631879FF", "LGHS"="#E68FAC",
  "B"="#003C67FF", "CD4_T"="#C2B280", "CD8_T"="#882D17", "NK"="#A1C299"
)

# ==============================================================================
# 3. 执行元数据清洗与因子化
# ==============================================================================

# 1. 强制合并 Quies 亚群
pbmc@meta.data$sub_labels[pbmc@meta.data$sub_labels == "KC_Basal_Quies"] <- "KC_Basal"

# 2. 转换 Factor 并指定顺序
pbmc@meta.data$bulk_labels <- factor(pbmc@meta.data$bulk_labels, levels = bulk_labels_sorted)
pbmc@meta.data$sub_labels  <- factor(pbmc@meta.data$sub_labels,  levels = names(celltype_colors_sub_labels))

# 3. 分级信息
pbmc@meta.data$Grade <- factor(pbmc@meta.data$Grade, levels = c("DPTDI", "SPTDI", "Normal"))

# 4. 时间节点标准化为 dpi
pbmc@meta.data$Time_category <- gsub("dpt|dpb", "dpi", pbmc@meta.data$Time_category)
pbmc@meta.data$Time_category <- factor(pbmc@meta.data$Time_category, 
                                       levels = c("0-7dpi", "8-14dpi", "15-28dpi", "1-2mph", "Normal"))

# ==============================================================================
# 4. 绘图示例
# ==============================================================================
# DimPlot(pbmc, group.by = "sub_labels", cols = celltype_colors_sub_labels) + NoLegend()
# 建议在画图前，先对数据进行一次 Quies 合并：
#pbmc@meta.data$sub_labels[pbmc@meta.data$sub_labels == "KC_Basal_Quies"] <- "KC_Basal"
pbmc@meta.data$sub_labels <- factor(pbmc@meta.data$sub_labels, levels = names(celltype_colors_sub_labels))

##################细胞比例分析
# ==============================================================================
# 2. 修复后的比例计算 (按 Grade 统计)
# ==============================================================================

# 获取所有 Grade
tissues <- unique(pbmc@meta.data$Grade)

# 创建完整组合 (列名必须改为 Grade 以便 join)
bulk_combinations <- expand.grid(Grade = tissues, bulk_labels = bulk_labels_sorted)
sub_combinations <- expand.grid(Grade = tissues, sub_labels = sub_labels_sorted)

# --- 统计 Bulk Labels ---
bulk_prop_grade <- pbmc@meta.data %>%
  group_by(Grade, bulk_labels) %>%   # 修改：使用 Grade 分组
  summarise(n = n(), .groups = "drop") %>%
  right_join(bulk_combinations, by = c("Grade", "bulk_labels")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(Grade) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# --- 统计 Sub Labels ---
sub_prop_grade <- pbmc@meta.data %>%
  group_by(Grade, sub_labels) %>%    # 修改：使用 Grade 分组
  summarise(n = n(), .groups = "drop") %>%
  right_join(sub_combinations, by = c("Grade", "sub_labels")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(Grade) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# 保存
write.csv(bulk_prop_grade, file = "bulk_label_proportion_by_Grade.csv", row.names = FALSE)
write.csv(sub_prop_grade, file = "sub_label_proportion_by_Grade.csv", row.names = FALSE)

# ==============================================================================
# 3. 按样本 (orig.ident) 统计并绘制堆积图
# ==============================================================================

# 获取所有样本名
samples <- unique(pbmc@meta.data$orig.ident)

# --- 准备绘图数据 (Bulk Labels) ---
bulk_sample_comb <- expand.grid(orig.ident = samples, bulk_labels = bulk_labels_sorted)

bulk_prop_sample <- pbmc@meta.data %>%
  group_by(orig.ident, bulk_labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  right_join(bulk_sample_comb, by = c("orig.ident", "bulk_labels")) %>% # 补全缺失值为0
  mutate(n = replace_na(n, 0)) %>%
  group_by(orig.ident) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# --- 准备绘图数据 (Sub Labels) ---
sub_sample_comb <- expand.grid(orig.ident = samples, sub_labels = sub_labels_sorted)

sub_prop_sample <- pbmc@meta.data %>%
  group_by(orig.ident, sub_labels) %>%
  summarise(n = n(), .groups = "drop") %>%
  right_join(sub_sample_comb, by = c("orig.ident", "sub_labels")) %>%
  mutate(n = replace_na(n, 0)) %>%
  group_by(orig.ident) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

# 为了让样本在图中按组织类型排序（可选），我们可以关联 Grade 信息
# 如果不需要排序，可以跳过这一小段
sample_info <- unique(pbmc@meta.data[, c("orig.ident", "Grade")])
bulk_prop_sample <- left_join(bulk_prop_sample, sample_info, by = "orig.ident")
sub_prop_sample <- left_join(sub_prop_sample, sample_info, by = "orig.ident")

# ==============================================================================
# 4. 绘图代码
# ==============================================================================

# --- Plot 1: Bulk Labels Stacked Barplot ---
p1 <- ggplot(bulk_prop_sample, aes(x = orig.ident, y = freq, fill = bulk_labels)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) + # position="fill" 保证高度为1
  scale_fill_manual(values = celltype_colors_bulk_labels) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(title = "Cell Type Proportion by Sample (Bulk)", x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white") # 如果使用了分面
  ) 
# 如果想按 Grade 分面展示： + facet_grid(~Grade, scales = "free_x", space = "free")

# --- Plot 2: Sub Labels Stacked Barplot ---
p2 <- ggplot(sub_prop_sample, aes(x = orig.ident, y = freq, fill = sub_labels)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = celltype_colors_sub_labels) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(title = "Cell Type Proportion by Sample (Sub)", x = "Sample", y = "Proportion", fill = "Sub Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    panel.grid = element_blank()
  )

# ==============================================================================
# 5. 保存图片
# ==============================================================================
ggsave("Barplot_Proportion_Bulk_Sample.pdf", p1, width = 10, height = 6)
ggsave("Barplot_Proportion_Sub_Sample.pdf", p2, width = 14, height = 8)

# 同时保存样本级别的比例数据 CSV
write.csv(bulk_prop_sample, "bulk_label_proportion_by_Sample.csv", row.names = FALSE)
write.csv(sub_prop_sample, "sub_label_proportion_by_Sample.csv", row.names = FALSE)







##############Figure 1##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(showtext)
# ==============================================================================
# 2. 定义全局黑底主题函数 (优化：移动图例至下方)
# ==============================================================================
theme_black_umap <- function() {
  theme_void() + 
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b=15)),
      # --- 图例位置 ---
      legend.position = "bottom",          
      legend.box.margin = margin(t = 5),   
      # --- 【核心修改：缩小行间距】 ---
      # 1. 调整垂直间距：使用负值可以显著收紧行间距
      legend.spacing.y = unit(-0.1, "cm"), 
      # 2. 调整水平间距（如果需要）：控制标签之间的左右距离
      legend.spacing.x = unit(0.1, "cm"),
      # 3. 调整图例按行排列时的整体对齐间距
      legend.box.spacing = unit(0.2, "cm"),
      # --- 【核心修改：缩小图例图标高度】 ---
      # 减小 key 的高度也有助于行间距看起来更紧凑
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.text = element_text(size = 8, color = "white"), 
      legend.title = element_text(size = 10, face = "bold", color = "white")
      #plot.margin = margin(10, 10, 10, 10)
    )
}
# ==============================================================================
# 3. 绘制总体 UMAP (优化：图例分3行排布)
# ==============================================================================
p_global <- DimPlot(pbmc, 
                    group.by = "bulk_labels", 
                    reduction = "umap", 
                    cols = celltype_colors_bulk_labels, 
                    pt.size = 1.2,      
                    label = FALSE,      
                    raster = TRUE) +    
  theme_black_umap() +
  labs(title = "") +
  
  # --- 【修改：设置图例为 3 行】 ---
  guides(color = guide_legend(
    override.aes = list(size = 5), 
    nrow = 3,           # 强制设为 3 行
    byrow = TRUE        # 按行排列，视觉上更符合阅读习惯
  ))

# ==============================================================================
# 4. 保存 (建议：因为下方多了3行图例，高度可以适当增加一点)
# ==============================================================================
# 预览
print(p_global)

showtext_auto(FALSE)
# 保存为 PDF
# 提示：由于下方增加了图例空间，如果觉得图形被压缩，可以将 height 从 5 稍微调大一点，比如 6
ggsave("Global_UMAP_ALL_bulk_labels.pdf", 
       plot = p_global, 
       width = 6, 
       height = 7, 
       device = cairo_pdf,  # 必须使用 cairo 引擎以支持嵌入字体
       bg = "black",
       family = "Arial")    # 显式指定字体族

##############Figure 1 (针对 Fib_K14 连贯性与图例优化版)##############
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(showtext)

# ==============================================================================
# 2. 定义全局黑底主题函数 (优化：彻底移除图例相关设置)
# ==============================================================================
theme_black_umap_nolegend <- function() {
  theme_void() + 
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b=15)),
      # 👑 核心修改：确保没有任何图例空间残留
      legend.position = "none" 
    )
}

# ==========================================
# 💎 专家提示：如何让点看起来连起来？
# ==========================================
# 在 UMAP 图中，点是否连成线，取决于两个因素：
# 1. 细胞数量（我们无法改变）
# 2. pt.size (点的大小)。
# 对于 28 万细胞的大切片，提高 pt.size 是最直接的方法。

# ==============================================================================
# 3. 绘制总体 UMAP (优化：放大 pt.size，移除图例)
# ==============================================================================
p_global_fixed <- DimPlot(pbmc, 
                          group.by = "sub_labels", 
                          reduction = "umap", 
                          cols = celltype_colors_sub_labels, 
                          # 👑 核心修改：大幅提高点的大小（建议尝试 5-10），
                          # 让 Fib_K14 的深蓝色点连成清晰的线条
                          pt.size = 1.4,      
                          label = FALSE,      
                          raster = TRUE) +    # 保持栅格化防卡顿
  theme_black_umap_nolegend() +
  labs(title = "") +
  NoLegend() # 👑 核心修改：彻底强制删除图例

# ==============================================================================
# 4. 预览与高清保存
# ==============================================================================
# 关闭 showtext 防止字体干扰
showtext_auto(FALSE)

# 预览
print(p_global_fixed)

# 保存为 PDF
# 提示：因为没有了图例，画布可以恢复为接近正方形
ggsave("Global_UMAP_ALL_sub_labels_NoLegend.pdf", 
       plot = p_global_fixed, 
       width = 7, 
       height = 7, # 移除图例空间后，高度可缩小
       device = cairo_pdf, 
       bg = "black",
       family = "Arial")

print("✅ 修改 pt.size 并移除图例后的图谱已保存。")



# ==============================================================================
# Figure 1 专用：超紧凑型三大类横向图例 (解决堆叠与过宽间距)
# ==============================================================================
library(ggplot2)
library(showtext)
library(patchwork)

showtext_auto()

# 核心绘图函数：极大缩减了 plot.margin 并调整了内部间距
draw_compact_legend <- function(cell_list, title_text, col_count = 3) {
  sub_colors <- celltype_colors_sub_labels[cell_list]
  df <- data.frame(
    CellType = factor(names(sub_colors), levels = names(sub_colors)),
    y = 1
  )
  
  ggplot(df, aes(x = CellType, y = y, color = CellType)) +
    geom_point(size = 5) +
    scale_color_manual(values = sub_colors) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      legend.position = "bottom",
      legend.title = element_text(color = "white", face = "bold", size = 16),
      legend.text = element_text(color = "white", size = 12, family = "Arial"),
      # 💎 关键优化 1：消除绘图区域的所有边距，特别是左右边距
      plot.margin = margin(t = 5, r = 0, b = 5, l = 0), 
      # 💎 关键优化 2：缩紧图标和文字的距离
      legend.spacing.x = unit(0.1, "cm"),
      legend.key.width = unit(0.1, "cm")
    ) +
    guides(color = guide_legend(
      title = title_text,
      ncol = col_count,
      byrow = FALSE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(size = 6)
    ))
}

# ------------------------------------------------------------------------------
# 2. 定义三大类成员名单
# ------------------------------------------------------------------------------
list_epi <- c(
  "KC_Basal", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous", "KC_Spinous_Mig", 
  "KC_Spinous_Mat", "KC_Granular", "SAC_SG_Progenitor", "SAC_SG_Clear", 
  "SAC_SG_Dark", "SAC_SG_Ductal", "SAC_HF_HFSC", "SAC_HF_DP_DS", 
  "SAC_HF_Matrix", "SAC_HF_ORS", "SAC_HF_IRS", "Melanocyte"
)

list_stromal <- c(
  "Fib_Papi", "Fib_SFRP2", "Fib_Fasci", "Fib_EN1", "Fib_Inflama", "Fib_Myo", 
  "Fib_Prolif", "Fib_K14", "Fib_CD45", "Endo_Arterial", "Endo_Venous", 
  "Endo_Lymphatic", "Endo_Capillary", "Endo_Prolif", "Pericyte", "Schwann"
)

list_immune <- c(
  "M1_Macrophage", "M2_Macrophage", "LAM", "Neutrophil", "Mast", "cDC1", 
  "cDC2", "pDC", "LGHS", "B_Mature", "B_Plasma", "CD4_T", "CD8_T", "NK"
)

# 生成三个组件
p1 <- draw_compact_legend(list_epi, "Epithelial & Appendage")
p2 <- draw_compact_legend(list_stromal, "Stromal & Vascular")
p3 <- draw_compact_legend(list_immune, "Immune Cells")

# 💎 关键优化 3：使用负间距拼接 (Standard patchwork logic)
# widths 参数控制每个板块的比例，& 符号统一应用主题
# 通过 plot_layout 强制把三个板块拉近
combined_legend_final <- p1 + p2 + p3 + 
  plot_layout(widths = c(1, 1, 0.8)) & 
  theme(
    plot.margin = margin(1, 0, 1, 0),      # 剥夺所有不必要的上下左右外边距
    plot.spacing = unit(-1.5, "cm")        # 👑 杀手锏：板块之间的缝隙直接设为负数
  )

# 高清保存
# 💎 注意：由于间距缩减了，width 不再需要之前的 24，可以设为 16-18
ggsave("Figure1_Legend_UltraCompact.pdf", 
       plot = combined_legend_final, 
       width = 15,   
       height = 3,   
       device = cairo_pdf, 
       bg = "black")

showtext_auto(FALSE)


# ==============================================================================
# 1. 定义亚群分组（根据您的 bulk_labels 进行归类）
# ==============================================================================
target_groups <- list(
  "Keratinocyte" = "Keratinocyte",
  "SAC"          = c("SAC_SG", "SAC_HF"),
  "Fibroblast"   = "Fibroblast",
  "Endothelial"  = c("Endothelial", "Pericyte"),
  "Immune"       = c("Mono_Macro", "Neutrophil", "Mast", "DC", "LGHS", "B", "CD4_T", "CD8_T", "NK")
)

# 设定输出根目录
base_output_dir <- "./Subclustering_Results"
if (!dir.exists(base_output_dir)) dir.create(base_output_dir, recursive = TRUE)

# ==============================================================================
# 1. 调整后的黑底主题（图例置于右侧）
# ==============================================================================
theme_black_umap_right <- function() {
  theme_void() + 
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b=15)),
      
      # --- 修改：图例回归右侧 ---
      legend.position = "right", 
      legend.justification = "center",
      legend.text = element_text(size = 8, color = "white"), 
      legend.title = element_text(size = 10, face = "bold", color = "white"),
      legend.spacing.y = unit(0.1, "cm"),
      
      plot.margin = margin(20, 20, 20, 20)
    )
}

# ==============================================================================
# 2. 自动化循环处理逻辑
# ==============================================================================

# ==============================================================================
# 修改后的自动化循环处理逻辑 (强调高凝聚力与簇间拉开)
# ==============================================================================

showtext_auto(FALSE)

for (group_name in names(target_groups)) {
  
  message(sprintf("\n>>> 正在处理亚群 (高凝聚模式): %s <<<", group_name))
  
  sub_dir <- file.path(base_output_dir, group_name)
  if (!dir.exists(sub_dir)) dir.create(sub_dir)
  
  sub_obj <- subset(pbmc, subset = bulk_labels %in% target_groups[[group_name]])
  
  tryCatch({
    sub_obj <- NormalizeData(sub_obj) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
      ScaleData() %>%
      RunPCA(npcs = 30, verbose = FALSE) %>%
      # --- 【深度优化：模拟文章中的紧凑聚集感】 ---
      RunUMAP(dims = 1:20, 
              n.neighbors = 50,    # 增加邻居数：有助于平滑噪声，让簇内更圆润
              min.dist = 0.01,     # 极小距离：强制让属于同一簇的细胞紧贴在一起
              spread = 0.5,        # 减小散布：防止簇被无限拉大
              n.epochs = 500,      # 增加迭代次数：让流形布局更稳定
              repulsion.strength = 1.2, # 增加斥力：让不相关的簇弹得更远
              verbose = FALSE) %>%
      FindNeighbors(dims = 1:20, verbose = FALSE) %>%
      FindClusters(resolution = 0.8, verbose = FALSE)
    
    # --- 【绘图优化：增加点的大小并关闭 Raster 以保证清晰度】 ---
    p_sub <- DimPlot(sub_obj, 
                     group.by = "sub_labels", 
                     reduction = "umap", 
                     cols = celltype_colors_sub_labels, 
                     pt.size = 0.8,    # 适当增大点：点稍微重叠会产生更强的“色块”感
                     label = FALSE,
                     raster = FALSE) +  # 小样本亚群建议关闭 raster 以获得矢量质感
      theme_black_umap_right() + 
      labs(title = paste(group_name, "Sub-clusters")) +
      # 让图例圆点更大，方便辨认
      guides(color = guide_legend(
        override.aes = list(size = 4, alpha = 1), 
        ncol = 1
      ))
    
    # 2.5 保存
    file_base <- file.path(sub_dir, paste0(group_name, "_Compact_UMAP"))
    ggsave(paste0(file_base, ".pdf"), 
           plot = p_sub, width = 7, height = 6, 
           device = cairo_pdf, bg = "black")
    
    saveRDS(sub_obj, file = file.path(sub_dir, paste0(group_name, "_processed.rds")))
    
  }, error = function(e) {
    message(sprintf("X [错误]: %s 处理失败: %s", group_name, e$message))
  })
  
  rm(sub_obj, p_sub)
  gc()
}














########################时空组学分析



saveRDS(pbmc,"pbmc_final.rds")


##############Figure 2##############
########################拟时序分析
# ==============================================================================
# Figure 2 终极可视化：PAGA 连通图 + 拟时间演化图 (仅限 Fib_K14/SAC/KC & 原始拟时版)
# ==============================================================================
# ==============================================================================
# Figure 2 终极可视化：PAGA连通图 + 拟时间演化图 (完全体重构 + 底层手动映射版)
# ==============================================================================
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(ggraph)
library(igraph)
library(showtext)
library(patchwork)

# 如果报错找不到 as.cell_data_set，请确保加载了 SeuratWrappers
# library(SeuratWrappers) 

showtext_auto()
out_dir <- "Trajectory_Analysis_Output"
if(!dir.exists(out_dir)) dir.create(out_dir)

# ------------------------------------------------------------------
# 1. 定义出版级黑底主题
# ------------------------------------------------------------------
theme_black_publication <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "sans", size = 12),
      axis.text = element_blank(), 
      axis.title = element_text(color = "white", size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, color = "white", margin = margin(b = 15)),
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white", size = 10),
      legend.title = element_text(color = "white", size = 12, face = "bold"),
      panel.grid = element_blank()
    )
}

# ------------------------------------------------------------------
# 2. 颜色字典与全量背景数据准备
# ------------------------------------------------------------------
celltype_colors_sub_labels <- c(
  "KC_Basal"="#008941", "KC_Basal_Mig"="#FF3030", "KC_Basal_Prolif"="#FFD700", 
  "KC_Spinous"="#00FF00FF", "KC_Spinous_Mig"="#1A0099FF", "KC_Spinous_Mat"="#EE4C97FF", "KC_Granular"="#882D17",
  "SAC_SG_Progenitor"="#0099CCFF", "SAC_SG_Clear"="#B7E4F9FF", "SAC_SG_Dark"="#193006FF", "SAC_SG_Ductal"="#8F7700FF",
  "SAC_HF_HFSC"="#00F5FF", "SAC_HF_DP_DS"="#CDBE70", "SAC_HF_Matrix"="#7FFFD4", "SAC_HF_ORS"="#A30059", "SAC_HF_IRS"="#CC9900FF", "Melanocyte"="#FF4500",
  "Fib_Papi"="#3DA873FF", "Fib_SFRP2"="#CCFF00FF", "Fib_Fasci"="#9370DB", "Fib_EN1"="#BDB76B", "Fib_Inflama"="#800000FF", "Fib_CD45"="#DA70D6", "Fib_Myo"="#B79762", "Fib_Prolif"="#EED2EE", "Fib_K14"="#4682B4",
  "Endo_Arterial"="#FF00FF", "Endo_Venous"="#4169E1", "Endo_Lymphatic"="#ADFF2F", "Endo_Capillary"="#CF4E9CFF", "Endo_Prolif"="#F08080", "Pericyte"="#FFFF00",
  "Schwann"="#7CFC00",
  "M1_Macrophage"="#FF4040", "M2_Macrophage"="#54FF9F", "LAM"="#8B4513", "Neutrophil"="#FF7F50",
  "Mast"="#FF1493", "cDC1"="#9B30FF", "cDC2"="#EE9A00", "pDC"="#D1EEEE", "LGHS"="#FFB5C5",
  "B_Mature"="#0000CD", "B_Plasma"="#FAF0E6", "CD4_T"="#00BFFF", "CD8_T"="#1E90FF", "NK"="#00E5EE"
)
master_colors <- celltype_colors_sub_labels
master_colors["Other"] <- "grey15" 

sub_labels_char <- as.character(pbmc$sub_labels)
# 👑 锁定全量 Fib、KC、SAC 为目标谱系
target_mask <- grepl("^KC_|^SAC_|^Fib_", sub_labels_char)
target_celltypes <- unique(sub_labels_char[target_mask])

umap_coords <- Embeddings(pbmc, "umap")
bg_data <- data.frame(umap_1 = umap_coords[, 1], umap_2 = umap_coords[, 2], cluster = sub_labels_char)
bg_data$cluster[!bg_data$cluster %in% target_celltypes] <- "Other"

# ==============================================================================
# 3. 核心图 A: PAGA 连通图渲染
# ==============================================================================
message(">>> 正在渲染左图: PAGA 连通图...")
graph_name <- names(pbmc@graphs)[1]
adj_matrix <- pbmc@graphs[[graph_name]]
cell_cluster_matrix <- model.matrix(~ 0 + factor(pbmc$sub_labels))
colnames(cell_cluster_matrix) <- levels(as.factor(pbmc$sub_labels))

cluster_connectivity <- t(cell_cluster_matrix) %*% adj_matrix %*% cell_cluster_matrix
cluster_connectivity_norm <- log1p(as.matrix(cluster_connectivity))

valid_targets <- intersect(target_celltypes, colnames(cluster_connectivity_norm))
conn_target_only <- cluster_connectivity_norm[valid_targets, valid_targets]

net <- graph_from_adjacency_matrix(conn_target_only, mode = "undirected", weighted = TRUE, diag = FALSE)
net_filtered <- delete.edges(net, which(E(net)$weight < quantile(E(net)$weight, 0.80)))

centroids <- bg_data %>% filter(cluster %in% V(net_filtered)$name) %>% group_by(cluster) %>% summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))
layout_matrix <- as.matrix(centroids[match(V(net_filtered)$name, centroids$cluster), c("UMAP_1", "UMAP_2")])

p_paga <- ggraph(net_filtered, layout = layout_matrix) +
  geom_point(data = bg_data, aes(x = umap_1, y = umap_2, color = cluster), alpha = 0.5, size = 1, show.legend = FALSE) +
  geom_edge_link(aes(width = weight), color = "white", alpha = 0.4, show.legend = FALSE) +
  scale_edge_width(range = c(1, 5)) +
  geom_node_point(aes(color = name), size = 4, show.legend = FALSE) +
  geom_node_label(aes(label = name), repel = TRUE, size = 4, fontface = "bold", color = "white", fill = alpha("black", 0.7), label.size = 0, family = "sans") +
  scale_color_manual(values = master_colors) +
  theme_black_publication() + theme(legend.position = "none") +
  labs(title = "Epidermal-Mesenchymal Connectivity Graph", x = "UMAP_1", y = "UMAP_2")

# ==============================================================================
# 4. 核心图 B: Monocle3 底层重构与手动映射
# ==============================================================================
message(">>> 正在提取目标谱系重构 Monocle3 轨迹...")

# 提取子集
pbmc_sub <- subset(pbmc, cells = colnames(pbmc)[target_mask])

# 👑 核心修复：不依赖 SeuratWrappers，纯手工原生态构建 CDS 对象
# 1. 提取表达矩阵
counts_data <- GetAssayData(pbmc_sub, assay = "RNA", slot = "counts")
# 2. 提取细胞元数据
cell_meta <- pbmc_sub@meta.data
# 3. 构建基因元数据 (Monocle3 强制要求必须有一列叫 gene_short_name)
gene_meta <- data.frame(gene_short_name = rownames(counts_data), row.names = rownames(counts_data))

# 4. 生成 CDS 对象
cds_sub <- new_cell_data_set(counts_data,
                             cell_metadata = cell_meta,
                             gene_metadata = gene_meta)

# 5. 关键一步：把 Seurat 算好的 UMAP 坐标强行“继承”给 Monocle3，保证左右图形状绝对一致！
reducedDims(cds_sub)[["UMAP"]] <- Embeddings(pbmc_sub, reduction = "umap")

# 强制破除孤岛进行全域学习
cds_sub <- cluster_cells(cds_sub, reduction_method = "UMAP")
cds_sub <- learn_graph(cds_sub, use_partition = FALSE)

# 重新定根：Fib_K14 / SAC -> KC_Granular
start_cells <- colnames(cds_sub)[colData(cds_sub)$sub_labels %in% c("Fib_K14", "SAC_SG_Progenitor", "SAC_HF_HFSC")]
cds_sub <- order_cells(cds_sub, root_cells = start_cells)

message(">>> 正在将拟时间逆向投射至全量 UMAP...")
# 提取计算好的拟时间并合并到 bg_data 中
pt_values <- pseudotime(cds_sub)
bg_data$pseudotime <- NA
bg_data[names(pt_values), "pseudotime"] <- pt_values

# 👑 提取 Monocle3 的白色轨迹树突骨架 (完全脱离 plot_cells)
mst <- principal_graph(cds_sub)[["UMAP"]]
edge_list <- igraph::as_data_frame(mst, what = "edges")
dp_mst <- principal_graph_aux(cds_sub)[["UMAP"]]$dp_mst

edge_list$x <- dp_mst[1, edge_list$from]
edge_list$y <- dp_mst[2, edge_list$from]
edge_list$xend <- dp_mst[1, edge_list$to]
edge_list$yend <- dp_mst[2, edge_list$to]

message(">>> 正在渲染底层绘图 (Layer by Layer)...")
p_time <- ggplot() +
  # 第 1 层：绘制非目标细胞 (暗灰色背景)
  geom_point(data = filter(bg_data, is.na(pseudotime)), 
             aes(x = umap_1, y = umap_2), color = "grey15", size = 0.8, alpha = 0.5) +
  # 第 2 层：绘制有拟时间值的目标细胞 (Fib/KC/SAC 渐变色)
  geom_point(data = filter(bg_data, !is.na(pseudotime) & is.finite(pseudotime)), 
             aes(x = umap_1, y = umap_2, color = pseudotime), size = 1) +
  # 颜色标尺
  scale_color_viridis_c(option = "plasma", begin = 0.1, name = "Pseudotime") +
  # 第 3 层：手工绘制白色树突骨架
  geom_segment(data = edge_list, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               color = "white", linewidth = 1.2, lineend = "round") +
  theme_black_publication() + theme(legend.position = "right") +
  labs(title = "Pseudotime Progression", x = "UMAP_1", y = "UMAP_2")

# ==============================================================================
# 5. 拼图与保存
# ==============================================================================
message(">>> 正在合成最终高清 PDF...")
combined_plot <- p_paga + p_time + plot_layout(ncol = 2)

ggsave(file.path(out_dir, "Figure2_Manual_Mapped_Pseudotime.pdf"), 
       plot = combined_plot, width = 16, height = 8, device = cairo_pdf, bg = "black")

showtext_auto(FALSE)
message("✅ 完美！已跳过 Monocle3 绘图层，实现底层数据无损融合。")














# ==============================================================================
# 6. 拟时序热图 (黑底级联版)
# ==============================================================================
message(">>> 正在绘制拟时序级联热图 (Fib & KC Only)...")

# 提取目标细胞
pbmc_vis <- subset(pbmc, subset = Monocle3_Pseudotime != Inf & sub_labels %in% target_celltypes)
pbmc_vis$pseudotime_bin <- cut(pbmc_vis$Monocle3_Pseudotime, breaks = 100, labels = FALSE)

# 聚合表达量
expr_matrix <- GetAssayData(pbmc_vis, slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(as.matrix(expr_matrix)))
expr_df$bin <- pbmc_vis$pseudotime_bin
heatmap_data <- expr_df %>% group_by(bin) %>% summarise(across(everything(), mean)) %>% arrange(bin) %>% select(-bin)
heatmap_mat_scaled <- t(scale(as.matrix(heatmap_data)))

# 热图配色 (黑底建议使用冷暖对比强的色系)
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "black", "#FFCC00"))

# 绘图
ht <- Heatmap(heatmap_mat_scaled,
              name = "Z-score",
              col = col_fun,
              cluster_columns = FALSE, 
              cluster_rows = TRUE,     
              show_column_names = FALSE,
              row_title = "Lineage Driving Genes",
              column_title = "Lineage Differentiation (Pseudotime)",
              # 💎 核心修改：黑底标签
              row_names_gp = gpar(col = "white", fontsize = 12, fontfamily = "Arial"),
              column_title_gp = gpar(col = "white", fontsize = 14, fontfamily = "Arial", fontface = "bold"),
              row_title_gp = gpar(col = "white", fontsize = 14, fontfamily = "Arial"),
              heatmap_legend_param = list(
                title_gp = gpar(col = "white", fontsize = 12),
                labels_gp = gpar(col = "white", fontsize = 10)
              ),
              border = "grey30")

pdf(file.path(out_dir, "Figure2_Heatmap_Black.pdf"), width = 10, height = 8)
grid.rect(gp = gpar(fill = "black", col = NA)) # 强制铺黑底
draw(ht, newpage = FALSE)
dev.off()

message("✅ Figure 2 黑底高标准可视化完成！")





# ==============================================================================
# 0. 环境设置与准备
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
  library(scales)
})

cds = readRDS("data/external/local_dec344101783")

# ==============================================================================
# 1. 定义丢失的自定义函数与变量
# ==============================================================================

# 设定输出目录
out_dir <- "Trajectory_Analysis_Output"
if(!dir.exists(out_dir)) dir.create(out_dir)

get_valid_root_nodes <- function(cds, start_clusters, cluster_col = "sub_labels") {
  # 1.1 获取目标亚群的细胞
  cell_ids <- which(colData(cds)[, cluster_col] %in% start_clusters)
  if (length(cell_ids) == 0) stop("未在数据中找到指定的起始亚群名称。")
  
  # 1.2 获取离这些细胞最近的节点
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  candidate_nodes <- names(sort(table(closest_vertex[cell_ids, ]), decreasing = TRUE))
  
  # 1.3 核心修复：获取主图中实际存在的节点列表
  existing_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  
  # 1.4 仅保留存在的节点
  valid_nodes <- intersect(candidate_nodes, existing_nodes)
  
  if (length(valid_nodes) == 0) {
    message("警告：候选节点均不在主图中，尝试返回第一个主节点 Y_1")
    return("Y_1")
  }
  
  # 返回前 3 个最频繁且有效的节点
  return(head(valid_nodes, 3))
}

# 丢失函数 2: 黑底主题
theme_black_umap <- function() {
  theme_void() + 
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "Arial"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b=10)),
      legend.position = "bottom",
      legend.text = element_text(size = 8, color = "white"),
      legend.title = element_text(size = 10, face = "bold", color = "white")
    )
}

# ==============================================================================
# 1. 对象过滤与提取：严格限定 KC, SAC, Fibroblast 大类
# ==============================================================================
message(">>> 正在提取目标谱系 (KC, SAC, Fibroblast)...")

cds <- as.cell_data_set(pbmc)
fData(cds)$gene_short_name <- rownames(fData(cds))
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 利用正则表达式筛选 sub_labels 以 KC_, SAC_, Fib_ 开头的细胞
all_labels <- colData(cds)$sub_labels
target_mask <- grepl("^KC_|^SAC_|^Fib_", all_labels)

target_cells <- colnames(cds)[target_mask]
non_target_cells <- colnames(cds)[!target_mask] # 保留非目标细胞名单用于画灰色底图

# 切割出纯净的子集用于计算轨迹
cds_sub <- cds[, target_cells]
# 重新学习骨架（此时骨架只会在这三类细胞中穿梭，不会被免疫细胞干扰）
cds_sub <- learn_graph(cds_sub, use_partition = FALSE) 

# ==============================================================================
# 修复：防弹版节点寻找函数 (强制交叉验证节点合法性)
# ==============================================================================
get_multiple_principal_nodes_safe <- function(cds, start_clusters) {
  root_list <- c()
  
  # 1. 提取当前 cds 对象中【真实存在】的主图节点
  valid_graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
  
  # 2. 提取细胞到节点的投影矩阵，并确保用细胞 Barcode 作为名称匹配
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  # 强制转化为向量以便于匹配
  closest_vertex_vec <- as.character(closest_vertex[, 1])
  names(closest_vertex_vec) <- rownames(closest_vertex)
  
  for (sc in start_clusters) {
    # 3. 获取目标亚群的所有细胞 Barcode
    target_cells <- colnames(cds)[colData(cds)$sub_labels == sc]
    
    if(length(target_cells) > 0) {
      # 获取这些细胞对应的节点
      nodes_for_cluster <- closest_vertex_vec[target_cells]
      
      # 统计节点出现的频次，按从高到低排序
      node_freq <- sort(table(nodes_for_cluster), decreasing = TRUE)
      candidate_nodes <- names(node_freq)
      
      # 4. 💎 核心安全防线：取交集，确保候选节点必须在主图中！
      safe_nodes <- intersect(candidate_nodes, valid_graph_nodes)
      
      if (length(safe_nodes) > 0) {
        best_node <- safe_nodes[1] # 取合法节点中频次最高的那一个
        root_list <- c(root_list, best_node)
        message(sprintf("✅ 成功锁定 [%s] 的起点节点: %s", sc, best_node))
      } else {
        warning(sprintf("⚠️ [%s] 的细胞未能映射到主图的任何有效节点上！", sc))
      }
    } else {
      warning(sprintf("⚠️ 当前数据集中不存在亚群: [%s]", sc))
    }
  }
  
  # 去重返回，防止两个极其相近的亚群共用同一个起点节点
  return(unique(root_list))
}

# ==============================================================================
# 继续执行拟时序计算
# ==============================================================================
# 严格设定您的双起点
start_pops <- c("Fib_K14", "SAC_SG_Progenitor")

# 使用安全函数获取节点
root_nodes_safe <- get_multiple_principal_nodes_safe(cds_sub, start_pops)

# 打印最终找到的节点进行确认
print("选定的起点节点为：")
print(root_nodes_safe)

# 计算拟时序 (Pseudotime)
cds_sub <- order_cells(cds_sub, root_pr_nodes = root_nodes_safe)

message("🎉 拟时序计算完成，可以开始下一步的画图了！")

# ==============================================================================
# 3. 绘图：灰色背景兜底的精美展示
# ==============================================================================
showtext_auto()

# 提取灰色背景细胞的 UMAP 坐标
umap_bg <- as.data.frame(reducedDims(cds)[non_target_cells, ])
colnames(umap_bg) <- c("data_dim_1", "data_dim_2")

# --- 图 A: 亚群分布 (灰色底图 + 焦点染色) ---
p1_sub <- plot_cells(cds_sub, color_cells_by = "sub_labels",
                     label_cell_groups = FALSE, cell_size = 0.8) +
  scale_color_manual(values = celltype_colors) + 
  theme_black_umap() +
  labs(title = "Lineage Tracing: Fib & Appendage to KC")

# 👑 核心黑科技：将灰色背景细胞强行插入 ggplot 图层的最底层 (after = 0)
p1_sub$layers <- append(
  p1_sub$layers, 
  geom_point(data = umap_bg, aes(x = data_dim_1, y = data_dim_2), color = "#333333", size = 0.5, alpha = 0.4), 
  after = 0
)


# --- 图 B: 拟时序图 (同样注入灰色底图) ---
p2_time <- plot_cells(cds_sub, color_cells_by = "pseudotime",
                      label_cell_groups = FALSE, cell_size = 0.8) +
  scale_color_gradientn(colors = c("#2A6EBBFF", "#FAFD7CFF", "#D60C00FF"), name = "Pseudotime") +
  theme_black_umap() +
  theme(legend.position = "bottom")

p2_time$layers <- append(
  p2_time$layers, 
  geom_point(data = umap_bg, aes(x = data_dim_1, y = data_dim_2), color = "#333333", size = 0.5, alpha = 0.4), 
  after = 0
)


# --- 图 C: 轨迹骨架 ---
p3_graph <- plot_cells(cds_sub, label_leaves = TRUE, label_branch_points = TRUE,
                       cell_size = 0) +
  theme_black_umap() +
  labs(title = "Trajectory Backbone")

p3_graph$layers <- append(
  p3_graph$layers, 
  geom_point(data = umap_bg, aes(x = data_dim_1, y = data_dim_2), color = "#333333", size = 0.5, alpha = 0.4), 
  after = 0
)

# 保存函数保持不变...
save_trajectory_plots(p1_sub, "Focused_SubLabels_WithBG", out_dir)
save_trajectory_plots(p2_time, "Focused_Pseudotime_DualStart", out_dir)
save_trajectory_plots(p3_graph, "Focused_Backbone_WithBG", out_dir)

showtext_auto(FALSE)

# ==============================================================================
# 5. 保存 CDS 文件
# ==============================================================================
saveRDS(cds, file = file.path(out_dir, "Monocle3_Global_Trajectory_CDS.rds"))
message("√ 任务完成：PDF 图表与 CDS 对象已保存至 ", out_dir)










# ==============================================================================
# 1. Monocle3 准备（保留全图结构）
# ==============================================================================
message(">>> 正在基于全局 UMAP 结构构建 SAC, KC, Fibroblast 轨迹...")

# 直接使用全量 pbmc 对象，确保坐标不变
cds <- as.cell_data_set(pbmc)
fData(cds)$gene_short_name <- rownames(fData(cds))

# 关键：手动将 Seurat 的 UMAP 坐标赋值给 cds，确保结构完全一致
# 如果您已经运行过前面的全局 UMAP，Monocle3 会默认读取
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE) 

# 设定多个起点
start_clusters <- c("Fib_K14")
root_nodes <- get_valid_root_nodes(cds, start_clusters = start_clusters)

# 此时运行 order_cells 不再报错
cds <- order_cells(cds, root_pr_nodes = root_nodes)

# ==============================================================================
# 1. 对象过滤：仅保留 Keratinocyte 和 Fibroblast
# ==============================================================================
# 这样不仅能让图片更清爽，还能极大减小文件体积和计算压力
target_cells <- colnames(cds)[colData(cds)$bulk_labels %in% c("Keratinocyte", "Fibroblast")]
cds_sub <- cds[, target_cells]

# 过滤后建议重新运行轨迹推断的最后两步，以确保骨架适应新的细胞子集
cds_sub <- learn_graph(cds_sub, use_partition = FALSE) 

# ==============================================================================
# 2. 定向起点：以 Fib_K14 为起始
# ==============================================================================
# 自动寻找距离 Fib_K14 质心最近的轨迹节点作为起点
get_earliest_principal_node <- function(cds, time_bin="Fib_K14") {
  cell_ids <- which(colData(cds)$sub_labels == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[cell_ids, , drop=FALSE])
  root_nodes <- names(which.max(table(closest_vertex)))
  return(root_nodes)
}

root_node <- get_earliest_principal_node(cds_sub, "Fib_K14")
cds_sub <- order_cells(cds_sub, root_pr_nodes = root_node)

# ==============================================================================
# 3. 绘图与多格式保存 (PDF + TIFF)
# ==============================================================================
showtext_auto()

# 定义保存函数以减少重复代码
save_trajectory_plots <- function(plot_obj, name_base, out_dir) {
  # 保存 PDF (cairo设备，矢量但体积优化)
  ggsave(file.path(out_dir, paste0(name_base, ".pdf")), 
         plot = plot_obj, width = 6, height = 7, device = cairo_pdf, bg = "black")
  
  # 保存 TIFF (压缩位图，加载极快，投稿常用)
  # compression = "lzw" 是减小体积的关键
  tiff(file.path(out_dir, paste0(name_base, ".tiff")), 
       width = 6, height = 7, units = "in", res = 300, compression = "lzw")
  print(plot_obj)
  dev.off()
}

# --- 图 A: 亚群分布 (Focus) ---
p1_sub <- plot_cells(cds_sub, color_cells_by = "sub_labels",
                     label_cell_groups = FALSE, cell_size = 0.5) +
  scale_color_manual(values = celltype_colors) + # 直接使用您的49色表
  theme_black_umap() +
  labs(title = "Fib & KC Lineage")

# --- 图 B: 以 Fib_K14 为起点的拟时序 ---
p2_time <- plot_cells(cds_sub, color_cells_by = "pseudotime",
                      label_cell_groups = FALSE, cell_size = 0.5) +
  scale_color_gradientn(colors = c("#2A6EBBFF", "#FAFD7CFF", "#D60C00FF"), name = "Pseudotime") +
  theme_black_umap() +
  theme(legend.position = "bottom")

# --- 图 C: 轨迹骨架 (只看路径) ---
p3_graph <- plot_cells(cds_sub, label_leaves = TRUE, label_branch_points = TRUE,
                       cell_size = 0) +
  theme_black_umap() +
  labs(title = "Trajectory Backbone")

# 执行保存
save_trajectory_plots(p1_sub, "Focused_SubLabels", out_dir)
save_trajectory_plots(p2_time, "Focused_Pseudotime_FibStart", out_dir)
save_trajectory_plots(p3_graph, "Focused_Backbone", out_dir)

showtext_auto(FALSE)



# ==============================================================================
# 0. slingshot分析，加载必要的包
# ==============================================================================
library(slingshot)
library(ggplot2)
library(SingleCellExperiment)
library(RColorBrewer)

message(">>> 正在启动 Slingshot 轨迹推断...")

# ==============================================================================
# 1. 数据准备与子集提取
# ==============================================================================
# 假设您已经有了之前处理好的 cds 对象，我们从中提取 UMAP 和细胞标签
# (或者从 Seurat 对象 pbmc 中提取: umap_coords <- Embeddings(pbmc, "umap"))
umap_coords <- reducedDims(cds)[["UMAP"]]
cell_labels <- colData(cds)$sub_labels

# 定义目标群体掩码 (KC, SAC, Fibroblast)
target_mask <- grepl("^KC_|^SAC_|^Fib_", cell_labels)

# 拆分坐标：用于算轨迹的目标点 vs 用于做底图的灰色背景点
umap_target <- umap_coords[target_mask, ]
labels_target <- cell_labels[target_mask]

umap_bg <- as.data.frame(umap_coords[!target_mask, ])
colnames(umap_bg) <- c("UMAP_1", "UMAP_2")

# ==============================================================================
# 2. 运行 Slingshot (👑 核心：强行指定起点与终点)
# ==============================================================================
message(">>> 正在计算 Slingshot 主曲线，强行锚定起点与终点...")

# 运行 Slingshot 算法
# 允许多个起点 (start.clus) 和指定的终点 (end.clus)
sds <- slingshot(data = umap_target, 
                 clusterLabels = labels_target, 
                 start.clus = c("Fib_K14", "SAC_SG_Progenitor"), 
                 end.clus = "KC_Granular")

# 提取拟时间 (Pseudotime)
# Slingshot 会为每条分支计算一个伪时间矩阵，我们取平均值作为全局伪时间
pt_matrix <- slingPseudotime(sds)
global_pt <- rowMeans(pt_matrix, na.rm = TRUE)

# 将目标细胞的数据整理成数据框用于 ggplot
df_target <- data.frame(
  UMAP_1 = umap_target[, 1],
  UMAP_2 = umap_target[, 2],
  sub_labels = labels_target,
  pseudotime = global_pt
)

# 提取 Slingshot 拟合的平滑曲线 (Curves) 坐标，用于画轨迹线
curves <- slingCurves(sds)
curve_df <- do.call(rbind, lapply(names(curves), function(c_name) {
  c_obj <- curves[[c_name]]
  df <- as.data.frame(c_obj$s[c_obj$ord, ])
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$curve_id <- c_name
  return(df)
}))

# ==============================================================================
# 3. 黑底白字高定绘图 (ggplot2)
# ==============================================================================
message("🎨 正在渲染黑底白字轨迹图谱...")

# 定义纯黑主题函数
theme_pure_black <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white", family = "sans"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white", face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white", face = "bold"),
      plot.title = element_text(color = "white", face = "bold", size = 16, hjust = 0.5)
    )
}

# --- 图 A: 亚群分布 (加 Slingshot 轨迹) ---
p1_slingshot_clusters <- ggplot() +
  # 1. 铺设灰色兜底背景
  geom_point(data = umap_bg, aes(x = UMAP_1, y = UMAP_2), color = "#333333", size = 0.5, alpha = 0.3) +
  # 2. 绘制目标细胞 (按细胞类型着色)
  geom_point(data = df_target, aes(x = UMAP_1, y = UMAP_2, color = sub_labels), size = 0.8) +
  # 3. 绘制 Slingshot 主干轨迹线 (亮白色粗线)
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2, group = curve_id), 
            color = "white", size = 1.5, arrow = arrow(length = unit(0.2, "inches"))) +
  scale_color_manual(values = celltype_colors) + # 请确保您环境中依然有 celltype_colors 变量
  theme_pure_black() +
  labs(title = "Slingshot Lineage: Forced Endpoint to KC_Granular", color = "Cell Type") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# --- 图 B: 伪时间分布 (加 Slingshot 轨迹) ---
p2_slingshot_pt <- ggplot() +
  # 1. 铺设灰色兜底背景
  geom_point(data = umap_bg, aes(x = UMAP_1, y = UMAP_2), color = "#333333", size = 0.5, alpha = 0.3) +
  # 2. 绘制目标细胞 (按 Pseudotime 着色，使用岩浆色带)
  geom_point(data = df_target, aes(x = UMAP_1, y = UMAP_2, color = pseudotime), size = 0.8) +
  # 3. 绘制 Slingshot 主干轨迹线
  geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2, group = curve_id), 
            color = "white", size = 1.5, arrow = arrow(length = unit(0.2, "inches"))) +
  scale_color_viridis_c(option = "magma", name = "Pseudotime") +
  theme_pure_black() +
  labs(title = "Slingshot Pseudotime Progression")

# ==============================================================================
# 4. 保存为高分辨率 PDF (黑底)
# ==============================================================================
# 保存函数
save_trajectory_plots <- function(plot_obj, name_base, out_dir) {
  # 必须使用 cairo_pdf 否则可能背景会变成白色
  ggsave(file.path(out_dir, paste0(name_base, ".pdf")), 
         plot = plot_obj, width = 8, height = 7, device = grDevices::cairo_pdf, bg = "black")
}

save_trajectory_plots(p1_slingshot_clusters, "Slingshot_Clusters_ForcedEndpoint", out_dir)
save_trajectory_plots(p2_slingshot_pt, "Slingshot_Pseudotime_ForcedEndpoint", out_dir)

message("✅ Slingshot 轨迹计算与黑底可视化全部完成！")




########################################细胞比例折线图
# ==============================================================================
# 1. 数据准备 (保持您的逻辑)
# ==============================================================================
# ==============================================================================
# 全量标签比例趋势图绘制：黑底白字版
# ==============================================================================
library(tidyverse)

# 1. 自动绘图函数：针对黑底进行视觉优化
# 💡 这里的逻辑是计算每个 label 在对应 Time & Grade 组内的相对比例
plot_trends_black <- function(df, label_column, color_vector) {
  # 计算比例
  df %>%
    group_by(!!sym(label_column), Time_category, Grade) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(Time_category, Grade) %>%
    mutate(total_in_group = sum(cell_count),
           percentage = (cell_count / total_in_group) * 100) %>%
    ungroup() -> summary_df
  
  # 按照标签拆分为列表进行循环绘图
  label_list <- unique(summary_df[[label_column]])
  plot_list <- list()
  
  for (lb in label_list) {
    sub_df <- summary_df %>% filter(!!sym(label_column) == lb)
    
    p <- ggplot(sub_df, aes(x = Time_category, y = percentage, 
                            color = Grade, group = Grade)) +
      geom_line(size = 0.8) +
      geom_point(size = 3) +
      scale_color_manual(values = color_vector) +
      # 黑底主题
      theme_minimal() + 
      labs(title =lb,
           x = "Time Point", y = "Percentage (%)") +
      theme(
        plot.background = element_rect(fill = "black", color = NA),
        panel.background = element_rect(fill = "black", color = NA),
        panel.grid.major = element_line(color = "#333333"), 
        panel.grid.minor = element_blank(),
        text = element_text(color = "white", family = "sans"),
        axis.line = element_line(color = "white", size = 0.8),
        axis.text = element_text(color = "white", size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "white"),
        legend.position = "top",
        legend.text = element_text(color = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    plot_list[[lb]] <- p
  }
  return(plot_list)
}

# 2. 执行计算与绘图逻辑
# 确保数据已同步
df <- pbmc@meta.data
# 合并 Quies 到 Basal 并处理 Normal 时间节点标签
#df$sub_labels[df$sub_labels == "KC_Basal_Quies"] <- "KC_Basal"
df$Time_category <- as.character(df$Time_category)
df$Time_category[df$Grade == "Normal"] <- "Normal"

# 因子化
time_levels <- c("Normal", "0-7dpi", "8-14dpi", "15-28dpi", "1-2mph")
grade_levels <- c("DPTDI", "SPTDI", "Normal")
df$Time_category <- factor(df$Time_category, levels = time_levels)
df$Grade <- factor(df$Grade, levels = grade_levels)

# 执行绘图
all_sub_plots <- plot_trends_black(df, "sub_labels", Grade_colors)
all_bulk_plots <- plot_trends_black(df, "bulk_labels", Grade_colors)

# 3. 批量保存为多页 PDF
# 建议将 sub_labels 和 bulk_labels 分开保存，方便查阅
pdf("Bulk_Labels_Trends_BlackBG.pdf", width = 4.5, height = 3)
for (p in all_bulk_plots) { print(p) }
dev.off()

pdf("Sub_Labels_Trends_BlackBG.pdf", width = 4.5, height = 3)
for (p in all_sub_plots) { print(p) }
dev.off()

print("✅ 所有亚型比例趋势图（黑底版）已保存完成。")




















#############Fib_K14标志物
# ==============================================================================
# 1. 提取 Fibroblast 子集并进行差异分析
# ==============================================================================
message(">>> 正在提取 Fibroblast 子集并寻找亚群 Marker...")

# 提取成纤维细胞
fib_obj <- subset(pbmc, subset = bulk_labels == "Fibroblast")
# 确保 Identity 设置为 sub_labels
Idents(fib_obj) <- "sub_labels"

# 寻找所有亚群的 Marker 基因 (仅保留上调基因)
fib_markers <- FindAllMarkers(fib_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              verbose = FALSE)
# 定义文件保存路径
fib_csv_path <- file.path(out_dir, "Fibroblast_SubCluster_Markers_Full.csv")

# 执行保存
write.csv(fib_markers, file = fib_csv_path, row.names = FALSE)

message("√ 成纤维细胞 Marker 列表已保存至: ", fib_csv_path)

# 每个亚群筛选前 2-3 个基因用于展示
top_markers <- fib_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>% # 设为 2 可保持 DotPlot 简洁，若空间充裕可改为 3
  pull(gene) %>%
  unique()

message(">>> 筛选出的关键 Marker 基因: ", paste(top_markers, collapse = ", "))

# 根据细胞功能逻辑手工整理的 Marker 列表
# 建议每个亚群选取 2-3 个最具代表性且在其他群低表达的基因
fib_top_markers_organized <- list(
  "Fib_K14" = c("KRT14", "KRT5", "TACSTD2", "EPCAM"), 
  "Fib_Prolif" = c("MKI67", "TOP2A"),
  "Fib_Myo" = c("ACTA2", "TAGLN"),
  "Fib_CD45" = c("PTPRC", "CD74", "HLA-DRA"),
  "Fib_Inflama" = c("IL6", "CCL2", "CXCL5"),
  "Fib_Fasci" = c("PROCR", "DPT"),
  "Fib_EN1" = c( "MFAP4", "DLK1"),
  "Fib_SFRP2" = c("SFRP2", "CTHRC1"),
  "Fib_Papi" = c("COL6A1", "LAMA4"))


# 提取纯基因名向量（用于画图）
genes_to_plot <- unique(unlist(fib_top_markers_organized))

message(">>> 手工整理完成，共选取 ", length(genes_to_plot), " 个核心 Marker。")

# ==============================================================================
# 2. 绘制 DotPlot (黑底风格 + 渐变色)
# ==============================================================================

# 定义渐变色 (参考 image_2282de.png：紫-绿-黄)
dotplot_colors <- c("#440154FF", "#35B779FF", "#FDE725FF")

p_dot <- DotPlot(fib_obj, 
                 features = genes_to_plot, 
                 col.min = -1.5, 
                 col.max = 2.5,
                 dot.scale = 8) + 
  # 使用 gradientn 确保三色平滑过渡
  scale_color_gradientn(colors = dotplot_colors, name = "Average\nExpression") +
  # --- 【关键修复：点亮图例点】 ---
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      override.aes = list(color = "white", fill = "white") # 强制图例里的圆圈变白
    ),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  # 统一黑底主题优化
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "gray20", size = 0.2), 
    panel.border = element_rect(color = "white", fill = NA, size = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "white", family = "Arial", size = 10, face = "italic"),
    axis.text.y = element_text(color = "white", family = "Arial", size = 10),
    axis.title = element_blank(),
    legend.background = element_rect(fill = "black"),
    legend.key = element_rect(fill = "black"), # 确保图例背景也是黑的
    legend.text = element_text(color = "white", size = 8),
    legend.title = element_text(color = "white", size = 9, face = "bold"),
    plot.title = element_text(color = "white", hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title = "Fibroblast Sub-cluster Markers")

# ==============================================================================
# 3. 输出 Adobe 可编辑 PDF
# ==============================================================================
showtext_auto(FALSE)

# 建议增加宽度以适应基因名称横向排布
ggsave(file.path(out_dir, "Fibroblast_Markers_DotPlot.pdf"), 
       plot = p_dot, width = 7, height = 3, 
       device = cairo_pdf, bg = "black", family = "Arial")

showtext_auto(TRUE)

message("√ DotPlot 已生成并保存至: ", out_dir)


# ==============================================================================
# 1. 提取 Fib_K14 核心 Marker 基因
# ==============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# 从上传的 CSV 中筛选 Fib_K14 的显著 Marker (p_val_adj < 0.05 & avg_log2FC > 1)
markers_df <- read.csv("Fibroblast_SubCluster_Markers_Full.csv")
fib_k14_genes <- markers_df %>%
  filter(cluster == "Fib_K14", p_val_adj < 0.05, avg_log2FC > 1) %>%
  pull(gene)

# ==============================================================================
# 2. GO 功能富集分析 (重点关注 BP: Biological Process)
# ==============================================================================
ego_k14 <- enrichGO(gene          = fib_k14_genes,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)

# 筛选并突出显示：上皮化 (Epithelialization) 与 成纤维功能 (ECM/Wound)
# 我们手动挑选最能体现您科学假设的条目进行展示
terms_to_highlight <- c(
  "skin development", "epidermis development", "epithelial cell differentiation", # 上皮化
  "keratinocyte differentiation", "intermediate filament organization",          # 结构特征
  "extracellular matrix organization", "collagen fibril organization",           # 成纤维功能
  "wound healing", "cell-matrix adhesion"                                        # 创面修复
)

# 过滤富集结果
ego_res <- as.data.frame(ego_k14)
plot_df <- ego_res %>%
  filter(Description %in% terms_to_highlight | p.adjust < 1e-10) %>%
  slice_max(order_by = Count, n = 6) # 取前15个核心条目

enrich_out_dir <- file.path(sub_dir, "Enrichment_Results")
if (!dir.exists(enrich_out_dir)) dir.create(enrich_out_dir, recursive = TRUE)

# 4.1 保存完整原始数据表 (方便后续筛选和制作表格)
# ------------------------------------------------------------------------------
write.csv(as.data.frame(ego_k14), 
          file = file.path(enrich_out_dir, "Fib_K14_GO_BP_Full_Results.csv"), 
          row.names = FALSE)
# ==============================================================================
# 3. 黑底风格可视化 (DotPlot 风格)
# ==============================================================================
# 重新排序，让图表更有层级感
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))

p_enrich <- ggplot(plot_df, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradientn(colors = c("#FF0000", "#FFCC00", "#008941"), # 红-黄-绿 渐变
                        guide = guide_colorbar(reverse = TRUE)) +
  scale_size_continuous(range = c(5, 12)) +
  labs(title = "Fib_K14 Functional Enrichment (BP)",
       subtitle = "Highlighting Epithelialization & ECM Remodeling") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "grey20", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    text = element_text(color = "white", family = "Arial"),
    axis.text = element_text(color = "white", size = 12),
    axis.title = element_text(color = "white", size = 14),
    title = element_text(face = "bold", size = 16),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(color = "white")
  )

# 保存高质量 PDF
ggsave("Fib_K14_Enrichment_Black.pdf", plot = p_enrich, width = 10, height = 6)

# ==============================================================================
# 1. 提取并过滤指定的 GO 条目
# ==============================================================================
library(ggplot2)
library(dplyr)
library(showtext)

# 您指定的 6 条核心 GO ID
target_go_ids <- c(
  "GO:0008544", "GO:0043588", "GO:0030216", 
  "GO:0061436", "GO:0042060", "GO:0085029"
)

# 从富集结果对象中提取数据
# 假设 ego_k14 是您之前运行 enrichGO 生成的对象
plot_df <- as.data.frame(ego_k14) %>%
  filter(ID %in% target_go_ids) %>%
  # 计算 -log10(p.adjust) 用于绘图
  mutate(logP = -log10(p.adjust)) %>%
  # 按照 logP 排序，确保绘图时的层级感
  arrange(logP) %>%
  mutate(Description = factor(Description, levels = Description))

# ==============================================================================
# 2. 绘制黑底风格 Barplot
# ==============================================================================
p_bar <- ggplot(plot_df, aes(x = Description, y = logP, fill = logP)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  # 使用经典的红黄绿渐变或单色调，这里推荐使用您喜欢的渐变色
  scale_fill_gradientn(colors = c("#008941", "#FFCC00", "#FF0000")) +
  coord_flip() + # 横向柱状图
  labs(
    title = "Fib_K14 Enriched Functions",
    x = NULL,
    y = "-log10(Adjusted P-value)",
    fill = "-log10(p.adj)"
  ) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major.x = element_line(color = "grey30", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "white", family = "sans"),
    axis.text = element_text(color = "white", size = 11),
    axis.title = element_text(color = "white", size = 12),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    legend.background = element_rect(fill = "black"),
    legend.text = element_text(color = "white"),
    legend.position = "right"
  )

# ==============================================================================
# 3. 稳健保存 (解决 Font Type 报错)
# ==============================================================================
# 导出 PDF (矢量图，投稿首选)
ggsave("Fib_K14_Selected_GO_Barplot.pdf", 
       plot = p_bar, 
       width = 5, 
       height = 3, 
       device = cairo_pdf, 
       bg = "black")

# 导出 PNG (高分辨位图)
ggsave("Fib_K14_Selected_GO_Barplot.png", 
       plot = p_bar, 
       width = 5, 
       height = 3, 
       dpi = 300, 
       bg = "black")

message(">>> 指定 GO 条目的 Barplot 已保存成功。")














#setwd("data/external/local_58a7345f101a")
#pbmc = readRDS("pbmc_final.rds")

# 设置结果保存目录
out_dir <- "Trajectory_Analysis_Output"
if(!dir.exists(out_dir)) dir.create(out_dir)

message(">>> 步骤 1: 准备就绪，结果将保存至: ", out_dir)

# ==============================================================================
# 1. 数据提取与重新降维 (The Foundation)
# ==============================================================================

# ==============================================================================
# 2. 方法一: Monocle3 分析
# ==============================================================================
message(">>> 步骤 2: 运行 Monocle3 轨迹推断...")

# 2.1 转换对象
cds <- as.cell_data_set(pbmc)
# 修复基因名 (Monocle3 需要)
fData(cds)$gene_short_name <- rownames(fData(cds))

# 2.2 聚类 (Monocle3 内部需要分区)
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 2.3 学习轨迹图
# use_partition = FALSE 是连接跨胚层谱系的关键！
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# 2.4 自动寻找起点 (基于 Fib_EN1)
# 定义一个辅助函数来找起点
get_earliest_principal_node <- function(cds, start_cluster_name, cluster_col = "sub_labels"){
  cell_ids <- which(colData(cds)[, cluster_col] == start_cluster_name)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

# 设定您的起点亚群名称 (请确保这里与您的 meta.data 一致)
start_cluster <- "Fib_K14" 

# 检查是否存在，如果不存在则报错或使用 Fibroblast
if(sum(pbmc$sub_labels == start_cluster) > 0) {
  root_node <- get_earliest_principal_node(cds, start_cluster_name = start_cluster)
  cds <- order_cells(cds, root_pr_nodes = root_node)
  message("    已自动将起点设定为: ", start_cluster)
} else {
  message("    警告: 找不到 Fib_EN1，将使用交互模式或默认起点。尝试使用第一个主节点。")
  cds <- order_cells(cds, root_pr_nodes = "Y_1") # 备用方案
}

# 2.5 绘图
# 图 A: 轨迹线图
p_mono_traj <- plot_cells(cds, color_cells_by = "sub_labels", 
                          label_groups_by_cluster = FALSE, label_leaves = FALSE, 
                          label_branch_points = TRUE, graph_label_size = 3) +
  ggtitle("Monocle3: Trajectory Graph") + NoAxes()

# 图 B: 拟时序热图
p_mono_time <- plot_cells(cds, color_cells_by = "pseudotime", 
                          label_cell_groups = FALSE, label_leaves = FALSE, 
                          label_branch_points = FALSE) +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  ggtitle("Monocle3: Pseudotime") + NoAxes()

saveRDS(cds, file.path(out_dir, "Monocle3cds.rds"))

# 保存
ggsave(file.path(out_dir, "02_Monocle3_Trajectory.tiff"), p_mono_traj, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_dir, "03_Monocle3_Pseudotime.tiff"), p_mono_time, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 3. 方法二: Slingshot 分析
# ==============================================================================
message(">>> 步骤 3: 运行 Slingshot 平滑曲线分析...")

# 3.1 转换对象
sce <- as.SingleCellExperiment(pbmc)

# 3.2 运行 Slingshot
# start.clus 必须是 sub_labels 中的一个真实名称
sds <- slingshot(sce, clusterLabels = pbmc$sub_labels, 
                 reducedDim = "UMAP", 
                 start.clus = start_cluster) # 使用与 Monocle 相同的起点

# 3.3 提取拟时序并绘图
# Slingshot 的绘图通常用 base plot 比较方便，我们将其保存为 PDF
pdf(file.path(out_dir, "04_Slingshot_Curves.pdf"), width = 8, height = 8)

# 提取颜色
colors <- viridis::magma(100)[as.numeric(cut(slingPseudotime(sds)[,1], breaks = 100))]
# 处理 NA 值颜色 (有些细胞可能不在轨迹上)
colors[is.na(colors)] <- "grey90"

plot(reducedDims(sce)$UMAP, col = scales::alpha(colors, 0.6), pch = 16, cex = 0.5, 
     main = "Slingshot Trajectory", xlab = "UMAP_1", ylab = "UMAP_2", bty = "l")
lines(SlingshotDataSet(sds), lwd = 3, col = "black")
legend("topright", legend = paste0("Start: ", start_cluster), bty = "n")

dev.off()

# ==============================================================================
# 4. 方法三: 分化潜能验证 (Stemness Proxy)
# ==============================================================================
message(">>> 步骤 4: 运行分化潜能验证 (nFeature Proxy)...")

# 使用基因检出数 (nFeature_RNA) 作为干性/复杂度的代理指标
# 干细胞通常转录组更复杂(High nFeature)，分化细胞更特化(Low nFeature)
p_stem <- FeaturePlot(pbmc, features = "nFeature_RNA", pt.size = 0.5, order = TRUE) + 
  scale_color_viridis_c(option = "magma", direction = -1, name = "Complexity") + 
  ggtitle("Predicted Stemness (nFeature Proxy)") +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave(file.path(out_dir, "05_Stemness_Proxy.tiff"), p_stem, width = 8, height = 8, dpi = 300)

# ==============================================================================
# 5. 方法四: PAGA 风格图抽象分析 (Connectivity Graph)
# ==============================================================================
message(">>> 步骤 5: 运行 PAGA 风格图抽象分析...")

library(igraph)
library(ggraph)
library(reshape2)

# --- 5.1 计算亚群间的连接强度 (Connectivity Matrix) ---
# 获取 SNN 邻接矩阵 (基于之前 FindNeighbors 的结果)
# 注意：确保这里调用的图名称正确，通常是 "RNA_snn" 或 "integrated_snn"
graph_name <- names(pbmc@graphs)[1] 
message("    使用邻接图: ", graph_name)

adj_matrix <- pbmc@graphs[[graph_name]]
clusters <- pbmc$sub_labels
unique_clusters <- levels(as.factor(clusters))

# 创建分群矩阵 (One-hot encoding)
cell_cluster_matrix <- model.matrix(~ 0 + factor(clusters))
colnames(cell_cluster_matrix) <- gsub("factor\\(clusters\\)", "", colnames(cell_cluster_matrix))

# 核心计算：Cluster-to-Cluster Connectivity
# 公式：亚群A与亚群B之间的连接数 = (A的细胞) * SNN矩阵 * (B的细胞)
cluster_connectivity <- t(cell_cluster_matrix) %*% adj_matrix %*% cell_cluster_matrix

# 归一化 (可选，为了消除亚群大小差异的影响)
# 这里我们使用简单的 Jaccard-like 归一化或直接使用原始强度
# 为了视觉效果，我们取对数平滑一下强连接
cluster_connectivity_norm <- log1p(as.matrix(cluster_connectivity))

# --- 5.2 构建 PAGA 网络图 ---
# 转换为 igraph 对象
net <- graph_from_adjacency_matrix(cluster_connectivity_norm, mode = "undirected", weighted = TRUE, diag = FALSE)

# 过滤弱连接 (只保留最强的骨架，类似 PAGA 的阈值过滤)
# 阈值策略：保留权重在前 20% 的边，或者大于平均值的边
threshold <- quantile(E(net)$weight, 0.75) # 只显示最强的 25% 连接
net_filtered <- delete.edges(net, which(E(net)$weight < threshold))

# --- 5.3 计算节点坐标 (UMAP Centroids) ---
# 让 PAGA 图叠加在 UMAP 坐标系上，这是 PAGA 的精髓
umap_coords <- Embeddings(pbmc, "umap")
centroids <- data.frame(umap_coords) %>%
  mutate(cluster = clusters) %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(umap_1), UMAP_2 = median(umap_2))

# 确保 centroids 顺序与 net 节点顺序一致
ordered_centroids <- centroids[match(V(net_filtered)$name, centroids$cluster), ]
layout_matrix <- as.matrix(ordered_centroids[, c("UMAP_1", "UMAP_2")])

# --- 5.4 绘图 (Publication Quality) ---
# 提取颜色 (匹配之前的 sub_labels 颜色)
# 注意：需要确保 colors 向量包含所有出现的亚群
colors_use <- celltype_colors_sub_labels[names(celltype_colors_sub_labels) %in% unique_clusters]
head(umap_coords)

p_paga <- ggraph(net_filtered, layout = layout_matrix) +
  # 1. 绘制背景点 (淡淡的 UMAP 轮廓，作为参考)
  geom_point(data = data.frame(umap_coords, cluster=clusters), 
             aes(x = umap_1, y = umap_2, color = cluster), 
             alpha = 0.1, size = 0.5, show.legend = FALSE) +
  
  # 2. 绘制边 (连接线) - 粗细代表连接强度
  geom_edge_link(aes(width = weight), color = "grey30", alpha = 0.6, show.legend = FALSE) +
  scale_edge_width(range = c(0.5, 3)) + # 控制线条粗细范围
  
  # 3. 绘制节点 (亚群中心)
  geom_node_point(aes(color = name), size = 8, show.legend = FALSE) +
  
  # 4. 添加文本标签 (带白底)
  geom_node_label(aes(label = name), repel = TRUE, size = 4, 
                  fontface = "bold", fill = alpha("white", 0.7), label.size = 0) +
  
  scale_color_manual(values = celltype_colors_sub_labels) +
  theme_void() +
  ggtitle("PAGA-like Connectivity Map") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# --- 5.5 保存 ---
ggsave(file.path(out_dir, "07_PAGA_Connectivity.tiff"), p_paga, width = 12, height = 12, dpi = 300)
ggsave(file.path(out_dir, "07_PAGA_Connectivity.pdf"), p_paga, width = 12, height = 12)


# ==============================================================================
# 6. 拟时序基因表达可视化 (仅展示特定细胞类型)
# ==============================================================================
#cds <- readRDS("Monocle3cds.rds")
message(">>> 步骤 6: 绘制基因随拟时序变化趋势 (仅展示上皮/间充质)...")

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(Seurat)

# ------------------------------------------------------------------
# 6.0 定义目标细胞与颜色 (用于可视化筛选)
# ------------------------------------------------------------------
target_celltypes <- c(
  # --- 1. 上皮细胞 (Keratinocytes) ---
  "KC_Basal", "KC_Basal_Quies", "KC_Basal_Mig", "KC_Basal_Prolif", "KC_Spinous",
  "KC_Spinous_Mig", "KC_Spinous_Mat", "KC_Granular",
  # --- 4. 间充质/成纤维细胞 (Fibroblasts) ---
  "Fib_K14"
)

# 完整的颜色列表 (假设顺序与您的因子 levels 对应，或者它是命名向量)
celltype_colors_sub_labels <- c(
  "#008941","#3B4992FF","#D60C00FF","#FAFD7CFF","#00FF00FF","#1A0099FF","#EE4C97FF","#882D17","#2A6EBBFF")

# ------------------------------------------------------------------
# 6.1 准备数据：将 Monocle3 的拟时序整合回 Seurat (对所有细胞)
# ------------------------------------------------------------------
# 提取拟时序值 (基于所有细胞计算的结果)
pseudotime_vals <- pseudotime(cds, reduction_method = "UMAP")

# 将拟时序添加到 Seurat 元数据中
pbmc <- AddMetaData(pbmc, metadata = pseudotime_vals, col.name = "Monocle3_Pseudotime")

# 过滤掉没有拟时序值的细胞 (Infinite)
# pbmc_traj 此时仍然包含所有类型的细胞 (免疫、血管等)，只要它们在轨迹上
pbmc_traj <- subset(pbmc, subset = Monocle3_Pseudotime != Inf)

message("    已提取全局拟时序信息，当前总细胞数: ", ncol(pbmc_traj))

# ------------------------------------------------------------------
# 6.2 读取并清洗目标基因列表
# ------------------------------------------------------------------
gene_file <- "validate_gene.txt"
if(file.exists(gene_file)){
  gene_df <- read.table(gene_file, header = TRUE, stringsAsFactors = FALSE)
  target_genes <- unique(as.character(gene_df[,1]))
  
  # 基因名检查
  valid_genes <- c()
  for(g in target_genes){
    if(g %in% rownames(pbmc_traj)){
      valid_genes <- c(valid_genes, g)
    } else {
      g_title <- stringr::str_to_title(g)
      if(g_title %in% rownames(pbmc_traj)){
        valid_genes <- c(valid_genes, g_title)
      }
    }
  }
  target_genes <- unique(valid_genes)
  message("    成功加载并匹配目标基因: ", length(target_genes), " 个")
} else {
  stop("未找到 validate_gene.txt 文件")
}

# ==============================================================================
# 6.3 可视化 A: 基因趋势折线图 (双图组合版 - PNG输出)
# ==============================================================================
message(">>> 正在绘制双图组合趋势图 (PNG格式)...")

library(patchwork) # 用于拼图
library(ggplot2)
library(dplyr)
library(viridis)

# 1. 设置 PNG 保存目录 (避免文件过多杂乱)
png_dir <- file.path(out_dir, "Gene_Trends_PNGs")
if(!dir.exists(png_dir)) dir.create(png_dir)

# 2. 提取所有细胞的数据
# (这里沿用您之前的逻辑，先提取再过滤)
plot_data_all <- FetchData(pbmc_traj, vars = c("sub_labels", "Monocle3_Pseudotime", target_genes))

# 3. 过滤数据框，只保留目标细胞类型
plot_data_filtered <- plot_data_all %>% 
  filter(sub_labels %in% target_celltypes)

# 检查
if(nrow(plot_data_filtered) == 0) stop("过滤后没有剩余细胞，请检查 target_celltypes！")

# 4. 循环绘图
for(gene in target_genes){
  
  # --- 准备单基因数据 ---
  df_plot <- plot_data_filtered[, c("sub_labels", "Monocle3_Pseudotime", gene)]
  colnames(df_plot) <- c("CellType", "Pseudotime", "Expression")
  
  # --- 左图: 按细胞类型 (Cell Type) 着色 ---
  # 展示：是谁在表达？
  p1 <- ggplot(df_plot, aes(x = Pseudotime, y = Expression)) +
    # 散点：按细胞类型着色
    geom_point(aes(color = CellType), size = 1.0, alpha = 0.6) +
    # 拟合线
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 1.2, fill = "grey80") +
    # 颜色
    scale_color_manual(values = celltype_colors_sub_labels) +
    theme_classic() +
    labs(title = paste0(gene, " by CellType"),
         x = "Pseudotime", y = "Expression") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_blank() # 去掉图例标题更简洁
    )
  
  # --- 右图: 按表达量 (Expression Intensity) 着色 ---
  # 展示：表达强度的分布和总体趋势 (不被细胞类型干扰)
  p2 <- ggplot(df_plot, aes(x = Pseudotime, y = Expression)) +
    # 散点：按表达量本身着色 (模拟 FeaturePlot 风格)
    # 使用 order = order(df_plot$Expression) 确保高表达点画在上面
    geom_point(aes(color = Expression), size = 1.0, alpha = 0.8) +
    # 拟合线
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 1.2, fill = "grey80") +
    # 颜色：使用 Magma 色阶 (黑/紫 -> 亮黄) 代表强度
    scale_color_viridis_c(option = "magma", direction = 1, name = "Expr") +
    theme_classic() +
    labs(title = paste0(gene, " Expression Trend"),
         x = "Pseudotime", y = "Expression") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    )
  
  # --- 拼图 ---
  # 左右排列，p1 需要图例空间大一些，所以宽度给宽一点
  combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(1.2, 1))
  
  # --- 保存为 PNG ---
  # 文件名包含基因名
  file_name <- file.path(png_dir, paste0(gene, "_Trend_Combined.png"))
  
  # 保存 (宽12英寸，高5英寸，分辨率300)
  ggsave(filename = file_name, plot = combined_plot, width = 10, height = 4, dpi = 300)
  
  message("  -> Saved: ", basename(file_name))
}

message("✅ 所有 PNG 图片已保存至: ", png_dir)










# ------------------------------------------------------------------
# 6.4 可视化 B: 拟时序热图 (仅展示 Target Cells)
# ------------------------------------------------------------------
message("    正在绘制拟时序级联热图 (已过滤细胞类型)...")

# 1. 【关键步骤】创建一个只包含目标细胞的临时 Seurat 对象
pbmc_vis <- subset(pbmc_traj, sub_labels %in% target_celltypes)

# 2. 将拟时序切分为 100 个区间 (Bins)
# 注意：我们使用 subset 后的对象进行 binning，这样热图反映的是目标细胞群内的变化
pbmc_vis$pseudotime_bin <- cut(pbmc_vis$Monocle3_Pseudotime, breaks = 100, labels = FALSE)

# 3. 计算每个 Bin 中目标基因的平均表达量
expr_matrix <- GetAssayData(pbmc_vis, slot = "data")[target_genes, ]
expr_df <- as.data.frame(t(as.matrix(expr_matrix)))
expr_df$bin <- pbmc_vis$pseudotime_bin

# 按 bin 聚合求均值
heatmap_data <- expr_df %>%
  group_by(bin) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>% # 增加 na.rm 防止某些 bin 为空
  arrange(bin) %>%
  select(-bin)

# 转置为 (基因 x 时间)
heatmap_mat <- t(as.matrix(heatmap_data))

# 4. 基因标准化 (Z-score)
# 处理可能的 NA (如果有某些 bin 没有目标细胞)
heatmap_mat[is.na(heatmap_mat)] <- 0 
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# 5. 绘制热图
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

ha <- HeatmapAnnotation(Pseudotime = anno_empty(border = FALSE, height = unit(4, "mm"))) 

ht <- Heatmap(heatmap_mat_scaled,
              name = "Z-score",
              col = col_fun,
              cluster_columns = FALSE, 
              cluster_rows = TRUE,     
              show_column_names = FALSE,
              show_row_names = TRUE,
              column_title = "Pseudotime (Epi/Mes Only)",
              row_names_gp = gpar(fontsize = 8),
              top_annotation = ha)

pdf(file.path(out_dir, "07_Gene_Cascade_Heatmap_Filtered.pdf"), width = 8, height = max(6, length(target_genes)*0.15))
draw(ht)
decorate_annotation("Pseudotime", {
  grid.lines(c(0, 1), c(0.5, 0.5), gp = gpar(lwd = 2, col = "black"))
  grid.polygon(x = c(0.98, 1, 0.98), y = c(0.2, 0.5, 0.8), gp = gpar(fill = "black", col = NA))
})
dev.off()

message(">>> 可视化完成 (Filtered)！")

# ==============================================================================
# 6.5 组合可视化: FeaturePlot + TrendPlot (仅展示上皮/间充质)
# ==============================================================================
message(">>> 步骤 6.5: 生成 FeaturePlot + TrendPlot 组合图 (已过滤)...")

##################################
library(patchwork)
library(ggplot2)
library(Seurat)
library(viridis)

# 1. 【关键步骤】创建一个仅包含目标细胞的 Seurat 子集对象
#    这确保了 FeaturePlot 只画出目标细胞，同时保留了它们在全局 UMAP 中的原始坐标
pbmc_vis <- subset(pbmc_traj, sub_labels %in% target_celltypes)

# 2. 设置保存路径
save_dir <- file.path(out_dir, "Gene_Combined_Plots_Filtered")
if(!dir.exists(save_dir)) dir.create(save_dir)

# 3. 创建多页 PDF
pdf_file <- file.path(save_dir, "All_Genes_Feature_Trend_Combined.pdf")
pdf(pdf_file, width = 14, height = 6) # 宽一点，容纳两个图

for(gene in target_genes){
  
  message("    Plotting: ", gene)
  
  # ------------------------------------------------------------------
  # 左图: FeaturePlot (展示基因在 UMAP 轨迹上的表达)
  # ------------------------------------------------------------------
  # 注意：这里使用的是 pbmc_vis (子集)，所以免疫细胞不会显示
  p_feature <- FeaturePlot(pbmc_vis, features = gene, pt.size = 1.0, order = TRUE) +
    scale_color_viridis_c(option = "magma", direction = 1, name = "Expr") +
    ggtitle(paste0(gene, " Expression\n(Epi/Mes Only)")) +
    theme_void() + # 去除坐标轴，突出轨迹形态
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # ------------------------------------------------------------------
  # 右图: Trend Plot (展示基因随拟时序的变化)
  # ------------------------------------------------------------------
  # 提取数据 (只包含 pbmc_vis 中的目标细胞)
  df_trend <- FetchData(pbmc_vis, vars = c("sub_labels", "Monocle3_Pseudotime", gene))
  colnames(df_trend) <- c("CellType", "Pseudotime", "Expression")
  
  p_trend <- ggplot(df_trend, aes(x = Pseudotime, y = Expression)) +
    # 1. 散点 (按细胞类型着色)
    geom_point(aes(color = CellType), size = 1.5, alpha = 0.6) +
    # 2. 拟合曲线 (黑色实线)
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 1.2, fill = "grey85") +
    # 3. 样式调整
    scale_color_manual(values = celltype_colors_sub_labels) + # 使用全局颜色表
    labs(title = paste0(gene, " Temporal Dynamics"), 
         x = "Pseudotime (Start -> End)", 
         y = "Expression Level") +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(color = "black")
    ) +
    # 强制让 X 轴从 0 开始，更符合时间逻辑
    scale_x_continuous(expand = c(0, 0)) 
  
  # ------------------------------------------------------------------
  # 拼图: 左 Feature + 右 Trend
  # ------------------------------------------------------------------
  # 布局：左图占 1 份宽，右图占 1.5 份宽 (右图有复杂的图例)
  combined_plot <- p_feature + p_trend + plot_layout(ncol = 2, widths = c(1, 1.5))
  
  print(combined_plot)
}

dev.off()

message("✅ 组合图绘制完成！文件位于: ", pdf_file)

##################Figure 4 
#################SCENIC
########################SCENIC
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(SingleR)
library(reticulate)
library(cowplot)
library(stringr)

setwd("data/external/local_d552fae5454f")
#setwd("data/external/local_1d8bfb594ee9")

pbmc = readRDS("pbmc_final.rds")
#pbmc = readRDS("pbmc.rds")

dge <- GetAssayData(object = pbmc, slot = "counts")
#rm(pbmc,dge)

pbmc  <- FindVariableFeatures(pbmc,selection.method = "vst", nfeatures = 20000)

gene = as.factor(unique(pbmc@assays$RNA@meta.data$var.features))

dge <- dge[which(rownames(dge)%in%gene),]

dim(dge)

#metadata=read.csv("metadata_final.csv",sep=",",header=T,check.names=F,row.names = 1)

#pbmc@meta.data = metadata
#colnames(metadata)

cell.info = pbmc@meta.data[,c("sub_labels","Grade","Time_category")]

cell.info[1:5,1:3]

colnames(cell.info)[1] = "CellType"

colnames(cell.info)[2] = "Grade"

colnames(cell.info)[3] = "Time_category"

cell.info$CellState <- paste(cell.info$CellType, cell.info$Grade, sep = "_")

rm(pbmc)


as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

dge = as_matrix(dge)

# 查看维度确认
dim(dge)


data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(
  org="hgnc", # 物种名 or hgnc, or dmel or mgi
  dbDir="data/external/local_288f93454abb", # RcisTarget databases location
  dbs="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", # file name of motif database
  datasetTitle="SCENIC on Human Cell Atlas", # choose a name for your analysis
  nCores=1
)

genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)

length(genesKept)

dge <- dge[genesKept, ]
dim(dge)

source("data/external/local_b39572d9d069")
avg20.rep1 <- AvgN(dge, cell.info, seed=1)
avg20.rep2 <- AvgN(dge, cell.info, seed=2)
avg20.rep3 <- AvgN(dge, cell.info, seed=3)

#########Save results
if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")

##########Avg20 DGE (grouped UMI)
source("data/external/local_ab96fac0c1d3")

saveLoom <- function(exprMat, output){
  ## prepare cell cluster info
  cellInfo <- data.frame(
    row.names = colnames(exprMat),
    cellType = sapply(strsplit(colnames(exprMat), split = "\\."), function(x) x[1])
  )
  ## 存储loom文件
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo)
  close_loom(loom)
}

saveLoom(avg20.rep1, "output/s1_avg20_rep1.loom")
saveLoom(avg20.rep2, "output/s1_avg20_rep2.loom")
saveLoom(avg20.rep3, "output/s1_avg20_rep3.loom")

###Filtered DGE (UMI)
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)

########
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)


###############################退出R语言
####################
cd external_asset_d552fae5454f

conda activate py38_env

export PATH=/root/.local/bin:$PATH

bash external_asset_40c5852e01ea avg20_rep1
bash external_asset_40c5852e01ea avg20_rep2
bash external_asset_40c5852e01ea avg20_rep3

############################
/root/anaconda3/envs/scenicpy_protocol/bin/pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple cffi
/root/anaconda3/envs/scenicpy_protocol/bin/pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple cmake==3.18.4
/root/anaconda3/envs/scenicpy_protocol/bin/pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple MulticoreTSNE
/root/anaconda3/envs/scenicpy_protocol/bin/pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple --ignore-installed llvmlite

#############################
source ~/miniconda3/etc/profile.d/conda.sh

#conda create -n py38_env python=3.8
conda activate py38_env

pip install umap-learn -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple cffi
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple cmake==3.18.4
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple MulticoreTSNE
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple --ignore-installed llvmlite
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple pandas==1.3.5
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple loompy
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple pyscenic

pip uninstall numba llvmlite numpy -y

pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple numpy==1.19.5
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple numba==0.55.0
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple llvmlite==0.38.1
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple "dask[array]" "dask[dataframe]"
pip install --user -i https://pypi.tuna.tsinghua.edu.cn/simple dask==2021.6.2 distributed==2021.6.2

f_loom_path_scenic=output/s2_avg20_rep1.pyscenic.loom
ctx_output=output/s2_avg20_rep1.reg.tsv
sample_name=output/s3_avg20_rep1
threads=10
min_regulon_size=10
/root/miniconda3/envs/py38_env/bin/python external_asset_0d1c015c08f4 $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

f_loom_path_scenic=output/s2_avg20_rep2.pyscenic.loom
ctx_output=output/s2_avg20_rep2.reg.tsv
sample_name=output/s3_avg20_rep2
threads=10
min_regulon_size=10
/root/miniconda3/envs/py38_env/bin/python external_asset_0d1c015c08f4 $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

f_loom_path_scenic=output/s2_avg20_rep3.pyscenic.loom
ctx_output=output/s2_avg20_rep3.reg.tsv
sample_name=output/s3_avg20_rep3
threads=10
min_regulon_size=10
/root/miniconda3/envs/py38_env/bin/python external_asset_0d1c015c08f4 $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size


#######pySCENIC及调控网络构建
##############################作图
cd external_asset_0458a1039a3f

R

library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(SingleR)
library(reticulate)
library(cowplot)
library(stringr)

setwd("data/external/local_0458a1039a3f")

pbmc = readRDS("pbmc_final.rds")

#metadata=read.csv("metadata.csv",sep=",",header=T,check.names=F,row.names = 1)

#pbmc@meta.data = metadata

###########Figure 1
library(ggplot2)
library(plyr)
library(ggsci)

#Load data
emb.tsne <- read.table("output/s3_avg20_rep1.tsne.txt", sep = "\t", row.names = 1, header = T)
colnames(emb.tsne) <- paste0("tSNE_", 1:2)

library(tidyverse)
emb.umap = pbmc@reductions$umap@cell.embeddings 

cell.info <- readRDS("output/s1_cell.info.rds")

cell.info <- cbind(cell.info, emb.tsne[rownames(cell.info), ])
cell.info <- cbind(cell.info, emb.umap[rownames(cell.info), ])

head(cell.info)

saveRDS(cell.info, "output/s4_cell.info.rds")
cell.info <- readRDS("output/s4_cell.info.rds")

## calculating the position of cluster labels
get_label_pos <- function(data, emb = "tSNE", group.by="ClusterID") {
  new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by)]
  colnames(new.data) <- c("x","y","cluster")
  clusters <- names(table(new.data$cluster))
  new.pos <- lapply(clusters, function(i) {
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = i)
  })
  do.call(rbind, new.pos)
}

#tSNE
pdf(file="tSNE.pdf",height=4,width=6.5)
ggplot(cell.info, aes(tSNE_1, tSNE_2, color=as.character(CellType))) + 
  geom_point(size=.3) + 
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "tSNE"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color="black")
  )
dev.off()

#UMAP
pdf(file="UMAP.pdf",height=4,width=6.5)
ggplot(cell.info, aes(umap_1, umap_2, color=as.character(CellType))) + 
  geom_point(size=.3) + 
  geom_text(inherit.aes = F, data = get_label_pos(cell.info, emb = "umap"), aes(x,y,label=label), size=3) + 
  theme_bw(base_size = 15) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color="black")
  )
dev.off()

##########tSNEPlot.tissue
tSNEPlot.tissue <- function(data, tissue.use, topn) {
  clusters.use <- table(subset(data, Tissue == tissue.use)$ClusterID)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$Tissue == tissue.use & data$ClusterID %in% clusters.use, data$CellType, "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  
  ggplot(data, aes(tSNE_1, tSNE_2, color=new.cluster)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1, 
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(tissue.use) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}

##########tSNEPlot.tissue
UAMPPlot.tissue <- function(data, tissue.use, topn) {
  clusters.use <- table(subset(data, Tissue == tissue.use)$ClusterID)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$Tissue == tissue.use & data$ClusterID %in% clusters.use, data$CellType, "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.05, 0.2)
  
  ggplot(data, aes(UMAP_1, UMAP_2, color=new.cluster)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1, 
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(tissue.use) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}


#tSNEPlot.cellType
tSNEPlot.cellType <- function(data, cell.type, topn) {
  tissues.use <- table(subset(data, CellType == cell.type)$Tissue)
  if(topn > length(tissues.use)) {
    topn <- length(tissues.use)
  }
  tissues.use <- names(head(sort(tissues.use, decreasing = T), n=topn))
  data$new.tissue <- ifelse(data$CellType == cell.type & data$Tissue %in% tissues.use, data$Tissue, "Others")
  data$new.tissue <- factor(data$new.tissue, levels = c(tissues.use, "Others"))
  data$pt.size <- ifelse(data$new.tissue == "Others", 0.05, 0.2)
  
  ggplot(data, aes(tSNE_1, tSNE_2, color=new.tissue)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1,
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(cell.type) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}

UAMPPlot.cellType <- function(data, cell.type, topn) {
  tissues.use <- table(subset(data, CellType == cell.type)$Tissue)
  if(topn > length(tissues.use)) {
    topn <- length(tissues.use)
  }
  tissues.use <- names(head(sort(tissues.use, decreasing = T), n=topn))
  data$new.tissue <- ifelse(data$CellType == cell.type & data$Tissue %in% tissues.use, data$Tissue, "Others")
  data$new.tissue <- factor(data$new.tissue, levels = c(tissues.use, "Others"))
  data$pt.size <- ifelse(data$new.tissue == "Others", 0.05, 0.2)
  
  ggplot(data, aes(UMAP_1, UMAP_1, color=new.tissue)) + 
    geom_point(size=data$pt.size) + 
    scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1,
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(cell.type) + 
    theme(legend.justification = c(0,0), 
          legend.position = c(0,0),
          legend.title = element_blank(),
          legend.key = element_rect(fill = alpha("white", 0)),
          legend.background = element_rect(fill=alpha('white', 0)),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}


#####Figure 2
library(data.table)
library(pbapply)
library(plyr)
# philentropy包用于计算JS散度，关于这个package的详细介绍在此：
# http://blog.fens.me/r-entropy/
#BiocManager::install("philentropy")
#install.packages("latex2exp")
#install.packages("dendextend")

library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

#读入RAS矩阵
#rasMat: rows=cellID, cols=regulon, values=Regulon Acticity Score

rasMat <- fread("output/s3_avg20_rep1.AUCell.txt", sep = "\t", header = T, data.table = F) # data.frame
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
saveRDS(rasMat, "output/s5_avg20_rep1.rasMat.rds")
#读入细胞类型矩阵
cell.info <- readRDS("output/s4_cell.info.rds")
cell.types <- names(table(cell.info$CellType))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$CellType == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)

#计算RSS矩阵(Regulon Specificity Score)
rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
saveRDS(rssMat, "output/s5_avg20_rep1.rssMat.rds")
rssMat <- readRDS("output/s5_avg20_rep1.rssMat.rds")
binMat <- read.table("output/s3_avg20_rep1.binary_mtx.txt", sep = "\t", header = T, row.names = 1, check.names = FALSE)
colnames(binMat) <- sub("(+)", "", colnames(binMat), fixed = T)

##Regulon Rank Plot
source("data/external/local_396f4bcc3a07")
options(stringsAsFactors = FALSE)

write.table(rssMat,file="rssMat.txt",quote=F,sep="\t",col.names=T)

#DimPlot
cell.info <- readRDS("output/s4_cell.info.rds")
cell.info <- cbind(cell.info, binMat[rownames(cell.info), ])

cell.info[1:5,1:10]

source("data/external/local_9b049f44409e")

#########combine
# 加载必要的包
library(data.table)
library(cowplot)

# 读取 rssMat.txt 文件
rssMat <- read.table("rssMat.txt", sep = "\t", row.names = 1, header = T)

# 获取所有细胞亚型列名
cell.types <- colnames(rssMat)

# 找出每个细胞亚型中打分最高的转录因子
top.regulons <- apply(rssMat, 2, function(x) names(which.max(x)))

# 加载你自己的绘图函数
PlotRegulonRank <- function(rssMat, cell.type, topn = 3) {
  # 拿到每个转录因子的 RSS 得分
  scores <- rssMat[, cell.type]
  
  # 生成 data frame，保留行名
  data <- data.frame(
    Regulons = rownames(rssMat),
    RSS = scores,
    stringsAsFactors = FALSE
  )
  
  # 按 RSS 降序排列
  data <- data[order(-data$RSS), ]
  
  # 提取前 topn 个高分 regulon 用于标注
  data.label <- head(data, n = topn)
  
  # 添加颜色列：前 topn 为高亮，其余为背景色
  data$pt.col <- ifelse(data$Regulons %in% data.label$Regulons, "#007D9B", "#BECEE3")
  
  # 只展示前 200 个点
  data <- head(data, 200)
  
  # 绘图
  ggplot(data, aes(x = reorder(Regulons, -RSS), y = RSS)) + 
    geom_point(size = 2, aes(color = pt.col)) +
    scale_color_identity() +
    ggrepel::geom_text_repel(
      data = data.label, aes(label = Regulons),
      size = 6, color = "#007D9B"
    ) + 
    ggtitle(cell.type) + 
    ylab("Specificity score") + xlab("Regulons") +
    theme_bw(base_size = 12) + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )
}

# 定义绘图组合函数
fig2Plot <- function(cell.type, regulon) {
  p.list <- list(
    PlotRegulonRank(rssMat, cell.type),
    DimPlot(cell.info, cell.type = cell.type),
    DimPlot(cell.info, regulon = regulon)
  )
  cowplot::plot_grid(plotlist = p.list, ncol = 3, rel_widths = c(3, 5, 5))
}

# 生成 PDF：每页一个图像
pdf("SCENIC_each_celltype_page.pdf", width = 12, height = 5)
for (cell.type in cell.types) {
  regulon <- top.regulons[cell.type]
  print(fig2Plot(cell.type, regulon))
}
dev.off()


###########利用SEEK对regulon进行验证
Lmo2.seek <- read.csv("../seek_results/Rara_SEEK.txt", sep = "\t", header = T)
Irf4.seek <- read.csv("../seek_results/Hey2_SEEK.txt", sep = "\t", header = T)
Sox10.seek <- read.csv("../seek_results/Usf1_SEEK.txt", sep = "\t", header = T)
Nkx2.1.seek <- read.csv("../seek_results/Hoxd9_SEEK.txt", sep = "\t", header = T)
SeekPlot <- function(seek.res, key.words){
  seek.res$Related <- grepl(key.words, seek.res$Description, ignore.case = T)
  
  data.use <- data.frame(
    Dataset = seek.res$Rank,
    P = -log10(seek.res$Coexpression.PValue),
    Related = seek.res$Related & seek.res$Coexpression.PValue < 0.01
  )
  m = sum(seek.res$Related)
  n = sum(data.use$Related)
  M = nrow(seek.res)
  N = sum(seek.res$Coexpression.PValue < 0.01)
  fisher.res <- fisher.test(matrix(c(n,m-n,N-n,M-N-m+n), ncol = 2), alternative = "greater")
  max.p <- max(data.use$P[is.finite(data.use$P)])
  data.use$P <- ifelse(is.finite(data.use$P), data.use$P, max.p)
  
  ggplot(data.use, aes(Dataset, P)) + 
    geom_point(color="#4590CE", size=3) + 
    geom_point(inherit.aes = F, data = subset(data.use, Related), 
               aes(Dataset, P), color="#E2AE2D", size=3) + 
    geom_hline(yintercept=2, color="grey", size=1) + 
    annotate("text", x=Inf, y=Inf, hjust=1.3, vjust=8.5, label=paste0("(", n, " out of ", m, ")")) + 
    annotate("text", x=Inf, y=Inf, hjust=1.4, vjust=10, label=paste0("p=", signif(fisher.res$p.value, 3))) + 
    ylab(TeX("-log_{10}(p-value)")) + 
    ggtitle("") +
    theme_bw(base_size = 15) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}
p1.1 <- SeekPlot(Lmo2.seek, "(hematopoietic)|(erythroid)")
p2.1 <- SeekPlot(Irf4.seek, "(B cell)|(lympho)")
p3.1 <- SeekPlot(Sox10.seek, "(CNS)|(Oligodend)|(brain)|spinal cord")
p4.1 <- SeekPlot(Nkx2.1.seek, "(lung)|(Alveolar)")
p1.2 <- cowplot::plot_grid(p1, p1.1, rel_widths = c(13,3))
p2.2 <- cowplot::plot_grid(p2, p2.1, rel_widths = c(13,3))
p3.2 <- cowplot::plot_grid(p3, p3.1, rel_widths = c(13,3))
p4.2 <- cowplot::plot_grid(p4, p4.1, rel_widths = c(13,3))

cowplot::plot_grid(p1.2, p2.2, p3.2, p4.2, ncol = 1)

#####
if(!dir.exists("../figure")) {
  dir.create("../figure")
}
pdf("../figures/Figure2.pdf", width = 16, height = 20, units = "in")
cowplot::plot_grid(p1.2, p2.2, p3.2, p4.2, ncol = 1)
dev.off()

########Figure 3
library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)

##Load data
rasMat <- readRDS("output/s5_avg20_rep1.rasMat.rds")
regulons <- read.table("output/s3_avg20_rep1.regulons.txt", sep = "\t")
regulon.names <- sapply(strsplit(regulons$V1, split = "\\("), function(x) x[1])
regulon.sizes <- sapply(strsplit(regulons$V3, split = ","), function(x) length(x))
regulon.names <- regulon.names[regulon.sizes>=10]

rasMat <- rasMat[, regulon.names]
dim(rasMat)

pccMat <- cor(rasMat)
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}
csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)
round(pccMat[1:10,1:10], 2)

csiMat.binary <- matrix(as.numeric(csiMat >= 0.7), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
csiMat.binary[1:10,1:10]

saveRDS(csiMat, "output/s6_avg20_rep1.csiMat.rds")
write.table(csiMat.binary, "output/s6_avg20_rep1.csiMat.binary.txt", sep = "\t")

#Figure 3A
mat = readRDS("output/s6_avg20_rep1.csiMat.rds")

#Dendrogram
h = 8
row_dend = as.dendrogram(hclust(dist(mat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))

pdf("Dendrogram.pdf", width = 10, height = 6)
plot(row_dend)
dev.off()

#Heatmap
pdf("regulation_heatmap.pdf", width = 8, height = 8)
col_range = c(0.7, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

ht <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

lgd <- Legend(
  col_fun = col_fun, 
  title = "", 
  at = col_range, 
  labels = c("low", "high"), 
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
decorate_heatmap_body("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2), 
            hjust = 0, vjust = 0, default.units = "npc", 
            gp = gpar(fill=NA, col="#FCB800", lwd=3))
  grid.text(label = paste0("M",clusters),
            x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
            default.units = "npc",
            hjust = 1, vjust = 0.5,
            gp = gpar(fontsize=12, fontface="bold"))
})
decorate_column_dend("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
            default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
})
dev.off()

#Figure S3C
tree = column_dend(ht)
ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
clusters <- names(table(ind))
regulon.clusters <- data.frame(regulon=names(ind), cluster=paste0("M",ind))

write.table(regulon.clusters, "output/s6_avg20_rep1.regulon_clusters.txt", sep = "\t", quote = F, row.names = F)
k = length(clusters)
cell.info <- readRDS("output/s4_cell.info.rds")
moduleRasMat <- lapply(paste0("M",1:k), function(x){
  regulon.use <- subset(regulon.clusters, cluster == x)$regulon
  rowMeans(rasMat[, regulon.use])
})
names(moduleRasMat) <- paste0("M",1:k)
moduleRasMat <- do.call(cbind, moduleRasMat)
cell.info <- cbind(cell.info, moduleRasMat[rownames(cell.info), ])
p.list <- lapply(paste0("M",1:k), function(module){
  data.use <- cell.info
  expression.color <- c("darkblue", "lightblue", "green", "yellow", "red")
  max.val <- quantile(data.use[, module], 0.99)
  low.val <- quantile(data.use[, module], 0.1)
  data.use[, module] <- ifelse(data.use[, module] > max.val, max.val, data.use[, module])
  ggplot(data.use, aes(umap_1, umap_2, color=get(module))) + 
    geom_point(size=0.05) + 
    theme_bw(base_size = 15) + 
    ggtitle(module) + 
    scale_color_gradientn(name = NULL, colors = expression.color) + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black")
    )
})

pdf("regulation_module_featureplot.pdf", width = 9, height = 6)
cowplot::plot_grid(plotlist = p.list, ncol = 3)
dev.off()

###Cell Type
cell.info1 = cell.info[,c("ClusterID","CellType")]

colnames(cell.info1)[1] = "Cluster"

colnames(cell.info1)[2] = "Cell_Type"

cellType.info <- cell.info1
## group module score by clusterID
clusters <- names(table(cell.info$ClusterID))

moduleScoreInCellType <- lapply(clusters, function(x) {
  cells.use <- rownames(subset(cell.info, ClusterID == x))
  colMeans(cell.info[cells.use, paste0("M",1:k)])
})
names(moduleScoreInCellType) <- clusters
moduleScoreInCellType <- do.call(rbind, moduleScoreInCellType)
moduleScoreInCellType <- as.data.frame(moduleScoreInCellType)

moduleScoreInCellType$CellType <- mapvalues(
  x = rownames(moduleScoreInCellType),
  from = cellType.info$Cluster,
  to = cellType.info$Cell_Type
)

plotCellTypeRank <- function(data, module, topn=5){
  data.use <- data
  data.use <- data.use[order(data.use[, module], decreasing = TRUE), ]
  data.use$Rank <- 1:nrow(data.use)
  data.use$pt.col <- ifelse(data.use$Rank <= topn, "#007D9B", "#BECEE3")
  data.label <- head(data.use, n = topn)
  data.label$delta <- c(Inf, abs(diff(data.label[, module])))
  
  ggplot(data.use, aes(Rank, get(module))) + 
    geom_point(size=3, color=data.use$pt.col) + 
    geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Rank, get(module), label=CellType), size=4, max.iter = 2e4) + 
    ggtitle(module) + 
    ylab("Regulon activity score") + 
    xlab("Cell type") + 
    theme_bw(base_size = 12) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5, face = "bold")
    )
}
p.list <- lapply(paste0("M",1:k), function(x) plotCellTypeRank(moduleScoreInCellType, module = x, topn = 10))

pdf("regulation_module_plotlist.pdf", width = 9, height = 6)
cowplot::plot_grid(plotlist = p.list, ncol = 3)
dev.off()

































##################Figure 5

########################cellchat
# 确保已经加载必要的包
####################CellChat
####cellchat全流程####
# install.packages("magrittr")
# devtools::install_github('satijalab/seurat-data')

#devtools::install_github("sqjin/CellChat") # 安装CellChat

# http://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.1.0.tar.gz
#install.packages('data/external/local_ba5b961c5282', repos=NULL, type="source")
library(Seurat)
library(CellChat)
library(tidyverse)
library(magrittr)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(RColorBrewer)
library(dplyr)

# 读取 pbmc 对象
#pbmc <- readRDS("data/external/local_d82bdbd11022")

# 可选：设置工作目录
# workDir = "data/external/local_28c5260b4ac4"
# setwd(workDir)

# 提取相关数据
# 这里假设pbmc是已经创建好的Seurat对象，直接进行后续分析
# 若需查看对象结构，使用 `pbmc` 来检查。
#pbmc <- subset(
#  pbmc,
#  subset = tissue == "Skin" 
  # Grade %in% c("SPTDI_Non_Scar", "DPTDI_Scar", "Normal") &
  #bulk_labels == "Endothelial"
#)
# 数据预处理：过滤细胞，去除低质量细胞
# 删除Seurat对象中的图形数据（邻接矩阵和SNN图等）
pbmc@graphs <- list()

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 500 & mt_percent < 15)    # 过滤低质量细胞

# 按照 pbmc 的细胞顺序对 metadata 排序
cell_names_pbmc <- colnames(pbmc)
metadata_sorted <- pbmc@meta.data[cell_names_pbmc, ]
write.csv(metadata_sorted, "metadata_sorted.csv")

# 获取所有 tissue 的唯一类型
tissue <- unique(pbmc@meta.data$Grade)

# 定义生成 CellChat 对象的函数
generate_cellchat <- function(tissue, pbmc, metadata_sorted) {
  message("Processing: ", tissue)
  
  spe_cells <- rownames(metadata_sorted[metadata_sorted$Grade == tissue, ])
  if (length(spe_cells) < 50) {
    message("Too few cells for: ", tissue, " — Skipping.")
    return(NULL)
  }
  
  data_input <- GetAssayData(pbmc, slot = "data")[, spe_cells]
  meta <- metadata_sorted[spe_cells, c("sub_labels", "Grade")]
  meta$sub_labels <- as.vector(meta$sub_labels)
  
  cellchat <- createCellChat(object = data_input)
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "sub_labels")
  cellchat@DB <- CellChatDB.human
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # cellchat <- projectData(cellchat, PPI.human) # 可选
  
  cellchat <- computeCommunProb(cellchat)
  # cellchat <- filterCommunication(cellchat, min.cells = 10) # 可选
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat, paste0(tissue, "_ALL_cellchat.rds"))
  
  return(cellchat)
}

# 批量执行
cellchat_list <- lapply(tissue, function(tg) generate_cellchat(tg, pbmc, metadata_sorted))
names(cellchat_list) <- tissue


###########################
library(Seurat)
library(CellChat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(pheatmap)

# 所有需要比较的组对
group_pairs <- list(
  c("DPTDI", "Normal"),
  c("SPTDI", "Normal"),
  c("DPTDI", "SPTDI"))

# 函数：处理单个 CellChat 对象
process_cellchat <- function(cellchat) {
  cellchat@images <- list()
  cellchat@options$datatype <- "RNA"
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- aggregateNet(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  
  return(cellchat)
}

# 函数：处理一对对象并保存
compare_and_save <- function(group1, group2) {
  file1 <- paste0(group1, "_ALL_cellchat.rds")
  file2 <- paste0(group2, "_ALL_cellchat.rds")
  
  if (!file.exists(file1) || !file.exists(file2)) {
    message("Skipping pair: ", group1, " vs ", group2, " — one or both files missing.")
    return(NULL)
  }
  
  message("Processing pair: ", group1, " vs ", group2)
  
  obj1 <- readRDS(file1)
  obj2 <- readRDS(file2)
  
  obj1 <- process_cellchat(obj1)
  obj2 <- process_cellchat(obj2)
  
  object.list <- list()
  object.list[[group1]] <- obj1
  object.list[[group2]] <- obj2
  
  merged <- mergeCellChat(object.list, add.names = names(object.list))
  
  out_name <- paste0(group1, "_vs_", group2, "_merged_cellchat.rds")
  saveRDS(merged, out_name)
  message("Saved merged file: ", out_name)
}

# 批量执行
walk(group_pairs, ~compare_and_save(.x[1], .x[2]))

# ==============================================================================
# 1. 环境准备
# ==============================================================================
library(CellChat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)

# 确保输出主目录存在
output_base_dir <- "CellChat_Single_Sample_Visualizations"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)

# 获取所有 RDS 文件
rds_files <- list.files(pattern = "_cellchat.rds$")

# ==============================================================================
# 2. 辅助函数定义
# ==============================================================================

# 黑底主题
theme_cellchat_black <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      panel.grid = element_line(color = "grey20"),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5, face = "bold", size = 16)
    )
}

# Base Plot 类型的黑底 Circle 图保存函数
save_circle_black <- function(file_path, mat, title_name) {
  # 使用 tiff 或 pdf 均可，这里以 tiff 为例
  tiff(file_path, width = 2500, height = 2500, res = 300, bg = "black")
  par(bg = "black", mar = c(1,1,1,1))
  netVisual_circle(mat, weight.scale = TRUE, title.name = title_name, 
                   label.edge = FALSE, edge.curved = 0.2, shape = "circle")
  dev.off()
}

# ==============================================================================
# 3. 循环执行可视化
# ==============================================================================
for (f in rds_files) {
  sample_name <- gsub(".rds", "", f)
  message("\n>>> 正在处理样本: ", sample_name)
  
  cellchat <- readRDS(f)
  
  # 自动获取通路
  sample_pathways <- cellchat@netP$pathways
  outdir <- file.path(output_base_dir, sample_name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 判断是否为合并对象 (Merged object)
  is_merged <- "datasets" %in% slotNames(cellchat)
  
  tryCatch({
    
    # --- (1) Bubble Plot ---
    pdf(file.path(outdir, "01_Bubble_Visual.pdf"), width = 12, height = 10)
    if (is_merged) {
      # 合并对象：展示组间对比气泡图 (默认对比 dataset 1 和 2)
      p <- netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45) + 
        theme_cellchat_black() + ggtitle(paste(sample_name, "Comparison Bubble"))
      print(p)
    } else {
      # 单样本：按 Source 循环
      for (cell in levels(cellchat@idents)) {
        p <- try(netVisual_bubble(cellchat, sources.use = cell, title = paste("Source:", cell)), silent = TRUE)
        if (!inherits(p, "try-error")) { # 修复为 inherits
          print(p + theme_cellchat_black())
        }
      }
    }
    dev.off()
    
    # --- (2) Circle 图 (仅限单样本或取合并样本的总强度) ---
    if (!is_merged) {
      save_circle_black(file.path(outdir, "02_Circle_Interaction_Count.tiff"), cellchat@net$count, "Count")
      save_circle_black(file.path(outdir, "03_Circle_Interaction_Weight.tiff"), cellchat@net$weight, "Weight")
    }
    
    # --- (3) Pathway 分析 ---
    if (length(sample_pathways) > 0) {
      
      # Chord Plots
      pdf(file.path(outdir, "04_Pathway_Chord_Network.pdf"), width = 10, height = 10, bg = "black")
      par(bg = "black", fg = "white")
      for (pw in sample_pathways) {
        try({
          # 如果是合并对象，netVisual_aggregate 需要指定 dataset
          if (is_merged) {
            # 合并对象默认画第一组的，或者您可以循环数据集
            netVisual_aggregate(cellchat, signaling = pw, layout = "chord", dataset = 1)
          } else {
            netVisual_aggregate(cellchat, signaling = pw, layout = "chord")
          }
          title(main = pw, col.main = "white")
        }, silent = TRUE)
      }
      dev.off()
      
      # Contribution & Expression
      pdf(file.path(outdir, "05_Pathway_Contribution_Expression.pdf"), width = 12, height = 9)
      for (pw in sample_pathways) {
        # 贡献度在合并对象中通常展示对比图
        p_con <- try(netAnalysis_contribution(cellchat, signaling = pw), silent = TRUE)
        if (!inherits(p_con, "try-error")) {
          print(p_con + theme_cellchat_black() + ggtitle(paste(pw, "Contribution")))
        }
        
        p_exp <- plotGeneExpression(cellchat, signaling = pw) + theme_cellchat_black()
        print(p_exp)
      }
      dev.off()
    }
    
    message("√ 样本 ", sample_name, " 处理成功。")
    
  }, error = function(e) {
    message("X 样本 ", sample_name, " 失败: ", e$message)
  })
  
  rm(cellchat); gc()
}

message("\n========== 所有样本常规可视化任务已结束 ==========")


# ==============================================================================
# 1. 环境准备与数据读取 (DPTDI vs SPTDI)
# ==============================================================================
library(CellChat)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)

# 强制加载 cairo 设备以保证 PDF 矢量图的字体兼容性
library(Cairo)

# 定义输出主目录
output_base_dir <- "CellChat_Single_Sample_Visualizations"
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)

# --- 请黄医生精确修改以下路径 ---
# 指向您 DPTDI 组的 .rds 文件
DPTDI_path <- "DPTDI_ALL_cellchat.rds" 
# 指向您 SPTDI 组的 .rds 文件
SPTDI_path <- "SPTDI_ALL_cellchat.rds" 

message("▶ 正在读取并预处理 DPTDI 与 SPTDI 数据...")

# 读取并确保中心度已计算 (这是进行贡献度分析的前提)
cellchat.DPTDI <- readRDS(DPTDI_path)
cellchat.DPTDI <- netAnalysis_computeCentrality(cellchat.DPTDI, slot.name = "netP")

cellchat.SPTDI <- readRDS(SPTDI_path)
cellchat.SPTDI <- netAnalysis_computeCentrality(cellchat.SPTDI, slot.name = "netP")

# 构建对比列表，并赋予明确的组名
object.list <- list(DPTDI = cellchat.DPTDI, SPTDI = cellchat.SPTDI)

# 合并对象
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# 清理内存
rm(cellchat.DPTDI, cellchat.SPTDI); gc()

message("√ DPTDI_vs_SPTDI 合并对象构建成功。")


# ==============================================================================
# 2. IGFL1 专项差异可视化 (DPTDI vs SPTDI)
# ==============================================================================

# 定义专项输出目录
igfl1_out_dir <- file.path(output_base_dir, "Comparison_DPTDI_vs_SPTDI_IGFL1")
dir.create(igfl1_out_dir, showWarnings = FALSE, recursive = TRUE)

message("\n▶ 正在生成 IGFL1 组间差异可视化图集 (DPTDI vs SPTDI)...")

# --- 定义黑底主题函数 ---
theme_cellchat_black <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      panel.grid = element_line(color = "grey20"),
      legend.background = element_rect(fill = "black"),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      plot.title = element_text(color = "white", hjust = 0.5, face = "bold")
    )
}

# ----------------------------------------------------------------------------
# (1) 总体强度对比：IGFL1 在深浅烧伤中的总通讯概率
# ----------------------------------------------------------------------------
# 此图用于 SCI 文章展示组间总体的信号升降
pdf(file.path(igfl1_out_dir, "01_IGFL1_Overall_Strength_Comparison.pdf"), width = 7, height = 6)
p1 <- rankNet(cellchat_merged, mode = "comparison", stacked = T, signaling = c("IGFL1")) +
  theme_cellchat_black() +
  # 设定标志性对比色：DPTDI 深红，SPTDI 橙黄 (代表不同深度火)
  scale_fill_manual(values = c("#A10000FF", "#FAFD7CFF")) + 
  ggtitle("IGFL1 Pathway Strength")
print(p1)
dev.off()

# ----------------------------------------------------------------------------
# (2) 气泡图：展示所有相关细胞亚群的通讯差异
# ----------------------------------------------------------------------------
# 此图是最详尽的通讯网络细节，重点观察您关注的 Fib_CD45 或 Fib_K14
pdf(file.path(igfl1_out_dir, "02_IGFL1_Bubble_Interactions_DPTDI_vs_SPTDI.pdf"), width = 11, height = 9)
p2 <- netVisual_bubble(cellchat_merged, signaling = c("IGFL1"), comparison = c(1, 2), angle.x = 45) +
  theme_cellchat_black() +
  ggtitle("DPTDI vs SPTDI: IGFL1 Interactions")
print(p2)
dev.off()

# ----------------------------------------------------------------------------
# (3) 独立网络图：对比 IGFL1 的信号来源分布
# ----------------------------------------------------------------------------
# 并列展示两组，观察 DPTDI 中是否出现了“信号中断”或“位置漂移”
pdf(file.path(igfl1_out_dir, "03_IGFL1_Network_Source_Comparison.pdf"), width = 14, height = 7, bg = "black")
par(mfrow = c(1, 2), bg = "black", fg = "white") # 关键设置

for (i in 1:length(object.list)) {
  pathway_present <- "IGFL1" %in% object.list[[i]]@netP$pathways
  if (pathway_present) {
    # 使用 layout = "circle" 展示亚群排布
    netVisual_aggregate(object.list[[i]], signaling = "IGFL1", layout = "circle",
                        title.name = paste("IGFL1 -", names(object.list)[i]))
  } else {
    plot.new()
    text(0.5, 0.5, paste("IGFL1 inactive in", names(object.list)[i]), col = "white", cex = 1.3)
  }
}
dev.off()

# ----------------------------------------------------------------------------
# (4) 贡献度分析：核心配受体对的组间差异
# ----------------------------------------------------------------------------
# 用于确认是否是特定的配体 (如 IGFL1) 的下调导致了 DPTDI 通讯减弱
pdf(file.path(igfl1_out_dir, "04_IGFL1_L-R_Pair_Contribution.pdf"), width = 9, height = 7)
p4 <- netAnalysis_contribution(cellchat_merged, signaling = "IGFL1") +
  theme_cellchat_black() +
  ggtitle("L-R Contribution: IGFL1")
print(p4)
dev.off()

# ----------------------------------------------------------------------------
# (5) 靶向基因表达验证：点图展示
# ----------------------------------------------------------------------------
# 将通讯概率差异回推到基因表达差异，增加文章说服力
pdf(file.path(igfl1_out_dir, "05_IGFL1_Gene_Expression_Validation.pdf"), width = 11, height = 6)
# 颜色渐变使用灰到红，代表表达量
p5 <- plotGeneExpression(cellchat_merged, signaling = "IGFL1", type = "dot", 
                         color.use = c("grey15", "#D60C00FF")) +
  theme_cellchat_black()
print(p5)
dev.off()

# 保存合并后的对象备查
saveRDS(cellchat_merged, file.path(igfl1_out_dir, "CellChat_Merged_DPTDI_vs_SPTDI.rds"))

message("✅ IGFL1 专项差异分析图集 (DPTDI vs SPTDI) 已成功保存至: ", igfl1_out_dir)










# -------- GSVA --------
library(Seurat)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)

# 导入 GMT 文件
geneSets <- getGmt("data/external/local_d7185c3af581")

as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# 获取总的细胞数量
total_cells <- ncol(pbmc)
batch_size <- 20000  # 每次处理6万个细胞
num_batches <- ceiling(total_cells / batch_size)  # 计算批次数量

# 存储所有 GSVA 结果的列表
all_gsva_results <- list()

# 逐批进行 GSVA 分析
for (batch in 1:num_batches) {
  start_cell <- (batch - 1) * batch_size + 1
  end_cell <- min(batch * batch_size, total_cells)
  
  # 获取当前批次的细胞名称
  batch_cells <- colnames(pbmc)[start_cell:end_cell]
  
  # 提取当前批次细胞的表达矩阵（直接通过列索引）
  df.data <- GetAssayData(pbmc, slot = "data")[, batch_cells]
  
  # 将稀疏矩阵转换为普通矩阵
  df.data1 = as_matrix(df.data)
  
  # 运行 GSVA 分析
  gsvaResult <- gsva(df.data1, 
                     geneSets,
                     method = "ssgsea",
                     ssgsea.norm = TRUE,
                     min.sz = 1, 
                     max.sz = 5000, 
                     verbose = TRUE,
                     parallel.sz = 4)
  
  # 将当前批次的 GSVA 结果存入列表
  all_gsva_results[[batch]] <- gsvaResult
}

# 合并所有 GSVA 结果
gsva_combined <- do.call(cbind, all_gsva_results)

# 因为是normalized后的数据，所以使用默认参数就可以了。根据自己的电脑cpu来调整parallel.sz
gsvascore = gsva_combined
dim(gsva_combined)
gsvascore[1:5, 1:5]

df.group <- data.frame(umi = names(Idents(pbmc)), 
                       cluster = as.character(pbmc@meta.data$bulk_labels), 
                       stringsAsFactors = F)


#????GSVA????
#????GSVA????
meta_data=pbmc@meta.data
gsvaOut0 = t(gsva_combined)
gsvaOut1 = cbind(meta_data,gsvaOut0)

dim(meta_data)
dim(gsvaOut0)

head(gsvascore)

gsvaOut2=rbind(id=colnames(gsvaOut1), gsvaOut1)
write.table(gsvaOut2, file="GSVA.result.txt", sep="\t", quote=F, col.names=T)

gsvascore=rbind(id=colnames(gsvascore), gsvascore)
write.table(gsvascore, file="gsvascore.txt", sep="\t", quote=F, col.names=F)

#############GSVA可视化
#library(Seurat)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)
library(patchwork)
library(grid) # 确保加载 grid

# ==============================================================================
# 1. 核心配置：应用您的颜色与排序 (完全保持一致)
# ==============================================================================

# --- 更新 Metadata 因子水平 ---
# 这一步至关重要，决定了所有绘图的顺序
pbmc@meta.data$bulk_labels <- factor(pbmc@meta.data$bulk_labels, levels = bulk_labels_sorted)
pbmc@meta.data$sub_labels <- factor(pbmc@meta.data$sub_labels, levels = sub_labels_sorted)
pbmc@meta.data$Grade       <- factor(pbmc@meta.data$Grade,       levels = c("DPTDI", "SPTDI", "Normal"))

# 设置默认 Identity 为 sub_labels
Idents(pbmc) <- "sub_labels"

# ==============================================================================
# 2. 准备 GSVA 数据
# ==============================================================================
# 假设 gsva_combined 已经计算好
library(Seurat)
library(dplyr)

# ==============================================================================
# 1. 读取 GSVA 结果文件
# ==============================================================================
message(">>> 正在读取 gsvascore.txt ...")

# check.names = FALSE 非常重要！防止 R 把细胞名的 "-" 自动改成 "."
gsva_data <- read.table("gsvascore.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# 转换为矩阵
gsva_matrix <- as.matrix(gsva_data)

# ==============================================================================
# 2. 数据清洗与对齐 (解决 subscript out of bounds 报错的核心)
# ==============================================================================
message(">>> 正在对齐细胞 ID ...")

# 1. 获取 pbmc 和 gsva 数据的共有细胞
common_cells <- intersect(colnames(pbmc), colnames(gsva_matrix))

if(length(common_cells) == 0) {
  stop("❌ 错误：pbmc 中的细胞名与 gsvascore.txt 中的细胞名完全不匹配！\n请检查 txt 文件中的列名格式 (例如是否把 '-' 变成了 '.')")
}

message(paste("  - pbmc 细胞数:", ncol(pbmc)))
message(paste("  - GSVA 文件细胞数:", ncol(gsva_matrix)))
message(paste("  - 共有细胞数:", length(common_cells)))

# 2. 筛选并重排矩阵，使其与 pbmc 完全一致
# 这一步会自动丢弃那些在 pbmc 中已经被过滤掉的细胞
gsva_matrix_aligned <- gsva_matrix[, common_cells]

# 如果 pbmc 里有些细胞在 GSVA 文件里没有（比较少见，但为了代码健壮性），我们需要只保留共有的
# 注意：这会临时把 pbmc 缩小到共有细胞的范围，通常我们希望 pbmc 保持不变
# 更好的做法是：只取交集部分的数据放入 Assay，Seurat 会自动处理对齐

# ==============================================================================
# 3. 创建 Assay 并存入 Seurat
# ==============================================================================
message(">>> 正在创建 GSVA Assay ...")

# 创建新的 Assay
# 注意：Seurat v5+ 建议使用 CreateAssay5Object，但在 v4 或通用情况下 CreateAssayObject 即可
# 警告 Feature names cannot have underscores 是正常的，Seurat 会自动处理
gsva_assay <- CreateAssayObject(counts = gsva_matrix_aligned)

# 将 Assay 添加到 pbmc (Seurat 会自动根据细胞名对齐，缺失的细胞会填 NA 或报错，但我们上面已经筛选了交集)
# 为了安全，我们先取子集确保 pbmc 只包含有 GSVA 数据的细胞（如果必须一一对应）
# 或者直接赋值：
pbmc[["GSVA"]] <- gsva_assay

# 如果 pbmc 中有细胞在 GSVA 里没有，赋值后这些细胞在 GSVA Assay 里就没有数据。
# 为了后续画图不报错，我们通常只保留有数据的细胞：
pbmc <- subset(pbmc, cells = common_cells)

# ==============================================================================
# 4. 标准化与准备
# ==============================================================================
# 切换默认 Assay
DefaultAssay(pbmc) <- "GSVA"

# ScaleData (热图可视化必须步骤)
message(">>> 正在进行 ScaleData ...")
pbmc <- ScaleData(pbmc, features = rownames(gsva_matrix_aligned))

# 输出目录
out_dir <- "GSVA_Visualization_Consistent_Colors"
if(!dir.exists(out_dir)) dir.create(out_dir)

pathways <- rownames(gsva_matrix)

# ==============================================================================
# 3. 绘制汇总图：Heatmap & DotPlot (修复版)
# ==============================================================================
message(">>> 正在绘制汇总 Heatmap 和 Dotplot...")

# --- 关键修复 1: 确保 data slot 有数据 ---
# 防止 data slot 为空导致 AverageExpression 报错
if(nrow(GetAssayData(pbmc, assay = "GSVA", slot = "data")) == 0){
  pbmc <- SetAssayData(pbmc, assay = "GSVA", slot = "data", 
                       new.data = GetAssayData(pbmc, assay = "GSVA", slot = "counts"))
}

# --- 关键修复 2: 改用 data slot 计算平均值，然后手动 Scale ---
# 1. 计算平均值 (使用原始分数的平均)
avg_gsva <- AverageExpression(pbmc, assays = "GSVA", group.by = "sub_labels", slot = "data")$GSVA

# 2. 确保列顺序与 sub_labels_sorted 一致 (并移除没有细胞的亚群)
present_levels <- intersect(sub_labels_sorted, colnames(avg_gsva))
avg_gsva <- avg_gsva[, present_levels, drop = FALSE]

# 3. 手动进行行标准化 (Row Z-score)
# t(scale(t(x))) 是标准的行归一化方法
avg_gsva_scaled <- t(scale(t(avg_gsva)))

# 处理可能产生的 NA (如果某一行标准差为0)
avg_gsva_scaled[is.na(avg_gsva_scaled)] <- 0

# --- 3.1 绘制 Heatmap ---

# 创建顶端颜色注释
# 确保颜色向量只包含存在的亚群，防止报错
used_colors <- celltype_colors_sub_labels[present_levels]

ha_top <- HeatmapAnnotation(
  CellType = present_levels,
  col = list(CellType = used_colors),
  show_legend = FALSE, # 图例太多可隐藏，或改为 TRUE
  show_annotation_name = FALSE
)

# 动态计算尺寸
h_calc <- max(8, nrow(avg_gsva_scaled) * 0.4)
w_calc <- max(10, ncol(avg_gsva_scaled) * 0.6)

# 绘制
p_heatmap <- grid.grabExpr(draw(
  Heatmap(avg_gsva_scaled,
          name = "Z-score",
          top_annotation = ha_top, 
          col = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")), # 蓝白红配色
          cluster_rows = TRUE,      # 通路聚类
          cluster_columns = FALSE,  # 亚群不聚类，保持自定义顺序
          show_row_names = TRUE,
          show_column_names = TRUE,
          column_names_rot = 45,
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          column_title = "Average GSVA Score (Z-score)"
  ),
  padding = unit(c(1, 1, 1, 1), "cm")
))

# ==============================================================================
# 3.2 绘制 DotPlot (GSVA 专用优化版)
# ==============================================================================
message(">>> 正在绘制 GSVA DotPlot ...")

# 1. 确保 DefaultAssay 正确
DefaultAssay(pbmc) <- "GSVA"

# 2. 检查 features 是否存在
valid_features <- intersect(rownames(gsva_data), rownames(pbmc))
if(length(valid_features) == 0) stop("❌ 错误：Seurat 对象中找不到对应的 GSVA 通路名。")

# 3. 绘图
# scale = TRUE: 对每一行(通路)进行 Z-score 标准化，这样可以清晰看到该通路在哪个亚群高
# col.min/max: 限制 Z-score 范围，防止极端值影响颜色分布
p_dot <- DotPlot(pbmc, 
                 features = valid_features, 
                 group.by = "sub_labels", 
                 scale = TRUE,  # 关键：按行标准化，突出差异
                 col.min = -2, 
                 col.max = 2) +
  coord_flip() + # 翻转坐标轴，让通路在 Y 轴，亚群在 X 轴
  
  # 使用红蓝配色 (Blue=低, White=0, Red=高)
  scale_color_gradientn(colours = c("#377EB8", "white", "#E41A1C")) +
  
  # 美化主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8), # 防止通路名太挤
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 1)
  ) +
  labs(title = "GSVA Scores DotPlot (Z-score)", y = NULL, x = NULL)

# 4. 打印预览
print(p_dot)

message("✅ GSVA DotPlot 绘制完成。注意：由于 GSVA 数据的稠密性，点的大小(Pct)可能差异不明显，请主要关注颜色差异。")

# ==============================================================================
# 4. 循环绘制单通路图 (Feature, Violin-CellType, Violin-Grade)
# ==============================================================================
message(">>> 正在绘制单通路详细图表...")

plot_list_feature <- list()
plot_list_vln <- list()
plot_list_vln_split <- list()

# 定义 Grade 比较组
my_comparisons <- list(c("DPTDI", "Normal"), c("SPTDI", "Normal"), c("DPTDI", "SPTDI"))

for (pathway in pathways) {
  safe_name <- gsub("[^a-zA-Z0-9]", "_", pathway)
  
  # --- FeaturePlot ---
  p_feat <- FeaturePlot(pbmc, features = pathway, order = T, min.cutoff = "q1", max.cutoff = "q99") +
    scale_color_viridis_c(option = "plasma") +
    ggtitle(paste0(pathway)) +
    theme(plot.title = element_text(size = 11))
  plot_list_feature[[safe_name]] <- p_feat
  
  # --- ViolinPlot (按 CellType，应用您的配色) ---
  # 使用 scale_fill_manual 应用自定义颜色
  p_vln <- VlnPlot(pbmc, features = pathway, group.by = "sub_labels", pt.size = 0, raster = FALSE) +
    scale_fill_manual(values = celltype_colors_sub_labels) + 
    geom_boxplot(width=0.2, fill="white", outlier.shape = NA, alpha = 0.5) +
    ggtitle("Score by Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", axis.title.x = element_blank())
  plot_list_vln[[safe_name]] <- p_vln
  
  # --- ViolinPlot (按 Grade 分组) ---
  # 提取数据绘图
  plot_data <- FetchData(pbmc, vars = c(pathway, "sub_labels", "Grade"))
  colnames(plot_data)[1] <- "Score"
  
  # Grade 配色 (建议使用简单区分色，以免和细胞类型混淆)
  disease_cols <- c("DPTDI" = "#E41A1C", "SPTDI" ="#4DAF4A", "Normal" =  "#377EB8")
  
  p_vln_split <- ggplot(plot_data, aes(x = sub_labels, y = Score, fill = Grade)) +
    geom_violin(scale = "width", trim = TRUE, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, fill="white", alpha=0.5) +
    stat_compare_means(aes(group = Grade), label = "p.signif", method = "kruskal.test", 
                       label.y.npc = "top", hide.ns = TRUE) + 
    scale_fill_manual(values = disease_cols) +
    theme_bw() +
    labs(y = "GSVA Score", x = NULL, title = "Score by Grade Status") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_list_vln_split[[safe_name]] <- p_vln_split
}

# ==============================================================================
# 5. 保存输出
# ==============================================================================
message(">>> 正在保存 PDF 和 PNG...")

# --- PDF ---
combined_pdf <- file.path(out_dir, "GSVA_Analysis_Consistent_Colors.pdf")
pdf(combined_pdf, width = 14, height = 9)

# Page 1: Heatmap
grid.newpage()
grid.draw(p_heatmap)

# Page 2: Dotplot
print(p_dot + theme(plot.margin = margin(1, 1, 1, 1, "cm")))

# Page 3+: Loop
for (name in names(plot_list_feature)) {
  # 布局设计：左上 Feature，右上 CellType Vln，下方 Grade Vln
  layout_design <- "
  AAABBB
  CCCCCC
  "
  combined_plot <- plot_list_feature[[name]] + plot_list_vln[[name]] + plot_list_vln_split[[name]] +
    plot_layout(design = layout_design, heights = c(1, 1.2)) +
    plot_annotation(title = paste0("Pathway: ", name), theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  print(combined_plot)
}
dev.off()

# --- PNG (创建子文件夹) ---
sub_dirs <- c("FeaturePlot", "Violin_CellType", "Violin_Grade")
for(d in sub_dirs) dir.create(file.path(out_dir, d), showWarnings = F)

# 保存汇总图
ggsave(file.path(out_dir, "0_Summary_Heatmap.png"), p_heatmap, width = w_calc, height = h_calc, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "0_Summary_Dotplot.png"), p_dot, width = 12, height = max(8, length(pathways)*0.4), bg = "white", limitsize = FALSE)

# 保存单图
for (name in names(plot_list_feature)) {
  ggsave(file.path(out_dir, "FeaturePlot", paste0(name, ".png")), plot_list_feature[[name]], width = 6, height = 6, bg = "white")
  ggsave(file.path(out_dir, "Violin_CellType", paste0(name, ".png")), plot_list_vln[[name]], width = 10, height = 6, bg = "white")
  ggsave(file.path(out_dir, "Violin_Grade", paste0(name, ".png")), plot_list_vln_split[[name]], width = 12, height = 6, bg = "white")
}

# 恢复默认 Assay
#DefaultAssay(pbmc) <- "RNA"
message(paste0("✅ 分析完成！颜色已校准。结果保存在: ", out_dir))



############################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggVennDiagram)
  library(ggplot2) 
})

############################
## 1. 读取差异基因结果
############################

deg_files <- list.files(
  "DEG_txt",
  pattern = "_DEGs.txt$",
  full.names = TRUE
)

for (deg_file in deg_files) {
  
  message("=== Processing DEG file: ", basename(deg_file))
  
  ## ---------- 解析文件名 ----------
  fname <- basename(deg_file)
  fname_core <- str_remove(fname, "_DEGs.txt$")
  
  parts <- str_split(fname_core, "_", simplify = TRUE)
  
  tissue     <- parts[1]
  label_type <- parts[length(parts)]
  comparison <- paste(parts[2:(length(parts) - 1)], collapse = "_")
  
  ## ---------- 创建输出目录 ----------
  out_base <- file.path("Output", tissue, comparison, label_type)
  
  dir.create(file.path(out_base, "TF_overlap"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_base, "Ligand_overlap"),
             recursive = TRUE, showWarnings = FALSE)
  
  ## ---------- 读取 DEG ----------
  deg <- read.table(
    deg_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  ## 仅分析 rss 中存在的细胞
  cell_types <- intersect(unique(deg$cluster), colnames(rss))
  
  ## ---------- 细胞级循环 ----------
  for (cell in cell_types) {
    
    message("  -> Cell: ", cell)
    
    deg_genes <- deg %>%
      filter(cluster == cell) %>%
      pull(gene) %>%
      unique()
    
    ## ===== TF =====
    tf_top <- get_top_tf(cell)
    tf_intersect <- intersect(deg_genes, tf_top)
    
    if (length(tf_intersect) > 0) {
      p_tf <- ggVennDiagram(
        list(DEG = deg_genes, TF = tf_top),
        label_alpha = 0
      ) + ggtitle(paste(cell, "DEG ∩ TF"))
      
      ggsave(
        file.path(out_base, "TF_overlap", paste0(cell, "_TF_venn.pdf")),
        p_tf, width = 4, height = 4
      )
      
      deg %>%
        filter(cluster == cell, gene %in% tf_intersect) %>%
        write.table(
          file.path(out_base, "TF_overlap", paste0(cell, "_TF_DEG.tsv")),
          sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
    
    ## ===== Ligand =====
    ligands <- get_cellchat_ligands(cell)
    ligand_intersect <- intersect(deg_genes, ligands)
    
    if (length(ligand_intersect) > 0) {
      p_lig <- ggVennDiagram(
        list(DEG = deg_genes, Ligand = ligands),
        label_alpha = 0
      ) + ggtitle(paste(cell, "DEG ∩ Ligand"))
      
      ggsave(
        file.path(out_base, "Ligand_overlap", paste0(cell, "_Ligand_venn.pdf")),
        p_lig, width = 4, height = 4
      )
      
      deg %>%
        filter(cluster == cell, gene %in% ligand_intersect) %>%
        write.table(
          file.path(out_base, "Ligand_overlap", paste0(cell, "_Ligand_DEG.tsv")),
          sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
  }
}

message("✅ All DEG files processed with new directory structure.")

#######################################关键受体配体对可视化
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggVennDiagram)
  library(ggplot2)
  library(Seurat)
  library(readr)
  library(gridExtra)
})

# ==============================================================================
# 0. 参数与路径设置
# ==============================================================================

# 基础路径
base_dir <- "data/external/local_d552fae5454f"
deg_dir  <- file.path(base_dir, "DEG_txt")
cellchat_dir <- file.path(base_dir, "CellChat_Single_Sample_Visualizations")

# 目标文件
target_deg_file <- file.path(deg_dir, "Skin_SPTDI_VS_DPTDI_sub_labels_DEGs.txt")
tf_file <- "data/external/local_93307c7e46a2"
cre_file <- "data/external/local_984f38af71fc"

# 颜色定义
celltype_colors_sub_labels <- c("#008941","#3B4992FF","#D60C00FF","#FAFD7CFF","#00FF00FF","#1A0099FF","#EE4C97FF","#882D17",
                                "#0099CCFF","#B7E4F9FF","#193006FF","#8F7700FF",
                                "#2B3D26","#99991EFF","#69C8ECFF","#A30059","#CC9900FF","#171717FF",
                                "#3DA873FF","#CCFF00FF","#E7C76FFF","#DCD300","#800000FF","#C50084FF","#1B1919FF","#8787FFFF","#2A6EBBFF",
                                "#AD002AFF","#2600D1FF","#76EE00","#CF4E9CFF","#E64B35FF","#FFDC91FF",
                                "#F38400","#84BD00FF",
                                "#CD534CFF","#3DA873FF","#5F559BFF","#DC0000FF",
                                "#e9148f","#631879FF","#FFC20AFF","#A1CAF1","#E68FAC",
                                "#003C67FF","#3D3B25FF",
                                "#C2B280", "#882D17","#A1C299")

grade_colors <- c("#e5192c","#E1C62F","#3a77b7")

# 输出目录
out_dir_root <- file.path("Output", "Skin", "SPTDI_VS_DPTDI", "sub_labels", "Intersection_Analysis")
dir.create(out_dir_root, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. 读取参考基因集 (Reference Gene Lists)
# ==============================================================================
message(">>> Loading reference gene lists...")

# 1.1 Load TFs
# 假设是单列或第一列是基因名
tfs_df <- read.table(tf_file, header = FALSE, stringsAsFactors = FALSE) # 根据实际文件调整 header
ref_tfs <- unique(tfs_df[,1]) 
message("   Loaded ", length(ref_tfs), " TFs.")

# 1.2 Load Cre Genes
cre_df <- read.table(cre_file, header = FALSE, stringsAsFactors = FALSE)
ref_cre <- unique(cre_df[,1])
message("   Loaded ", length(ref_cre), " Cre-related genes.")

# 1.3 Load CellChat Interactions (Recursive Search)
message("   Scanning CellChat interactions...")
cc_files <- list.files(cellchat_dir, pattern = "Interactions.csv", recursive = TRUE, full.names = TRUE)
ref_cellchat <- c()

if(length(cc_files) > 0){
  for(f in cc_files){
    tmp <- read.csv(f, stringsAsFactors = FALSE)
    # 假设 Interations.csv 中包含 ligand 和 receptor 列
    # 根据实际 CellChat 输出调整列名，通常是 'ligand' 和 'receptor' 或 'interaction_name_2'
    if(all(c("ligand", "receptor") %in% colnames(tmp))){
      genes <- c(tmp$ligand, tmp$receptor)
    } else {
      # 备选方案：如果是 CellChat 对象导出的 data.frame，可能列名不同，这里做简单处理
      # 尝试提取 gene symbol 格式的字符串
      genes <- unlist(tmp) 
    }
    ref_cellchat <- c(ref_cellchat, genes)
  }
  ref_cellchat <- unique(ref_cellchat)
  # 清洗一下，去除非基因名的杂质
  ref_cellchat <- ref_cellchat[grep("^[A-Z0-9]+$", ref_cellchat)] 
  message("   Loaded ", length(ref_cellchat), " unique CellChat interaction genes.")
} else {
  warning("   No Interactions.csv files found!")
}

# ==============================================================================
# 2. 读取差异基因 (Target DEGs)
# ==============================================================================
message(">>> Loading DEGs...")
if(!file.exists(target_deg_file)) stop("DEG file not found: ", target_deg_file)

deg_data <- read.table(target_deg_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
# 确保列名正确，通常 Seurat FindMarkers 输出包含 'gene' 和 'cluster'
if(!"gene" %in% colnames(deg_data)) deg_data$gene <- rownames(deg_data)

unique_clusters <- unique(deg_data$cluster)

# ==============================================================================
# 3. 定义绘图与交集函数
# ==============================================================================

# 通用函数：处理交集并生成 Venn PDF 和 基因列表
process_intersection <- function(deg_df, ref_list, ref_name, output_subdir) {
  
  save_dir <- file.path(out_dir_root, output_subdir)
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  venn_plots <- list()
  intersect_genes_list <- list() # 用于后续可视化
  
  for(cell in unique_clusters) {
    # 提取该细胞类型的 DEGs (通常筛选显著的，例如 p_val_adj < 0.05)
    cell_degs <- deg_df %>% filter(cluster == cell, p_val_adj < 0.05) %>% pull(gene)
    
    if(length(cell_degs) == 0) next
    
    # 计算交集
    overlap <- intersect(cell_degs, ref_list)
    intersect_genes_list[[cell]] <- overlap
    
    # 只有当有交集时才画图
    if(length(overlap) > 0) {
      # 1. 绘制 Venn
      p <- ggVennDiagram(
        list(DEG = cell_degs, Ref = ref_list),
        category.names = c(paste0(cell, "\nDEGs"), ref_name),
        label_alpha = 0
      ) + 
        scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
        ggtitle(paste0(cell, ": DEG vs ", ref_name)) +
        theme(plot.title = element_text(hjust = 0.5, size = 10))
      
      venn_plots[[cell]] <- p
      
      # 2. 保存交集基因列表
      out_table <- deg_df %>% filter(cluster == cell, gene %in% overlap)
      write.table(out_table, 
                  file.path(save_dir, paste0("Genelist_", cell, "_intersect.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  
  # 3. 输出 Venn PDF (多页)
  if(length(venn_plots) > 0) {
    pdf(file.path(save_dir, paste0("All_Cells_", ref_name, "_Venn.pdf")), width = 5, height = 4)
    for(p in venn_plots) print(p)
    dev.off()
    message("   Generated Venn PDF for ", ref_name)
  }
  
  return(intersect_genes_list)
}

# ==============================================================================
# 4. 执行交集分析
# ==============================================================================

message(">>> Processing TF Intersections...")
genes_tf <- process_intersection(deg_data, ref_tfs, "TF", "TF_Overlap")

message(">>> Processing Cre Intersections...")
genes_cre <- process_intersection(deg_data, ref_cre, "Cre", "Cre_Overlap")

message(">>> Processing CellChat Intersections...")
genes_cc <- process_intersection(deg_data, ref_cellchat, "CellChat", "CellChat_Overlap")


# ==============================================================================
# 5. 表达量可视化 (Violin & Line Plots)
# ==============================================================================
# 前提：需要 pbmc 对象存在
if(exists("pbmc")) {
  message(">>> Generating expression plots using Seurat object...")
  
  # --- 准备时间轴数据 ---
  # 如果 pbmc@meta.data 中没有 Time 列，尝试解析 Grade
  # 假设 Grade 格式类似 "Skin_D1", "Skin_D3" 等
  if(!"Time" %in% colnames(pbmc@meta.data)) {
    message("   'Time' column not found. Attempting to parse from 'Grade'...")
    # 这是一个示例解析逻辑，请根据实际 Grade 命名修改
    # 提取数字部分作为时间
    pbmc$Time <- str_extract(pbmc$Grade, "\\d+") 
    # 如果提取失败，或者 Time 是无序的，这里需要手动指定因子水平
    pbmc$Time <- factor(pbmc$Time, levels = stringr::str_sort(unique(pbmc$Time), numeric = TRUE))
  }
  
  # 确保 sub_labels 是因子以匹配颜色
  pbmc$sub_labels <- as.factor(pbmc$sub_labels)
  
  # ---------------------------------------------------------
  # 可视化函数
  # ---------------------------------------------------------
  generate_visualization <- function(target_genes_list, category_name) {
    
    save_dir <- file.path(out_dir_root, paste0(category_name, "_Overlap"), "Visualization")
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    
    pdf_file <- file.path(save_dir, paste0("Integrated_Plots_", category_name, ".pdf"))
    pdf(pdf_file, width = 12, height = 10) # 大画布
    
    processed_genes <- c() # 避免重复画同一个基因
    
    # 遍历每个细胞类型及其交集基因
    for(cell in names(target_genes_list)) {
      genes <- target_genes_list[[cell]]
      # 限制基因数量，避免 PDF 过大，每个细胞最多画 Top 10 (按 logFC 排序)
      # 需重新从 deg_data 获取 logFC 排序
      top_genes <- deg_data %>% 
        filter(cluster == cell, gene %in% genes) %>% 
        arrange(desc(abs(avg_log2FC))) %>% 
        head(5) %>% # 【注意】这里只取前5个最显著的基因进行作图，以免太慢
        pull(gene)
      
      if(length(top_genes) == 0) next
      
      message("   Plotting ", category_name, " for ", cell)
      
      # 子集化 Seurat 对象，只保留该细胞类型（用于 Plot 2 和 3）
      # 注意：Plot 1 需要用全数据集
      
      for(gene in top_genes) {
        
        # --- 图 1: 不同 sub_labels 水平的 ViolinPlot ---
        # 展示该基因在所有细胞类型中的表达
        p1 <- VlnPlot(pbmc, features = gene, group.by = "sub_labels", pt.size = 0, cols = celltype_colors_sub_labels) +
          ggtitle(paste0(gene, " (Global Expression)")) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), legend.position = "none") +
          xlab("")
        
        # --- 图 2: 该细胞中不同 Grade 的 ViolinPlot ---
        # 需要先 subset 到该细胞类型
        sub_obj <- subset(pbmc, sub_labels == cell)
        
        p2 <- VlnPlot(sub_obj, features = gene, group.by = "Grade", pt.size = 0.1, cols = grade_colors) +
          ggtitle(paste0(gene, " in ", cell, " by Grade")) +
          theme(legend.position = "none") +
          xlab("")
        
        # --- 图 3: 折线图 (Time_category vs Expression by Grade) ---
        # 提取数据
        expr_data <- FetchData(sub_obj, vars = c(gene, "Grade", "Time_category"))
        colnames(expr_data) <- c("Expression", "Grade", "Time_category")
        
        # 聚合数据：计算每个 Time点、每个 Grade 的均值
        # 注意：这里假设 Grade 和 Time_category 是关联的。如果每个 Grade 只有一个 Time_category 点，则只有一条线。
        # 如果 Grade 是分组（如 WT vs KO），Time_category 是时间点，则会有多条线。
        # 这里按用户要求：不同的折线颜色表示不同的 Grade，x轴为时间
        
        plot_dat <- expr_data %>%
          group_by(Grade, Time_category) %>%
          summarise(Mean_Exp = mean(Expression), .groups = "drop")
        
        p3 <- ggplot(plot_dat, aes(x = Time_category, y = Mean_Exp, group = Grade, color = Grade)) +
          geom_line(linewidth = 1) +
          geom_point(size = 3) +
          scale_color_manual(values = grade_colors) +
          theme_classic() +
          labs(title = paste0(gene, " Dynamics in ", cell), y = "Avg Expression") +
          theme(legend.position = "right")
        
        # --- 组合三张图到一页 ---
        # 布局：上层放 global vln，下层放 local vln 和 line plot
        combined_plot <- gridExtra::arrangeGrob(
          p1, 
          gridExtra::arrangeGrob(p2, p3, ncol = 2), 
          nrow = 2, heights = c(1, 1),
          top = paste0("Gene: ", gene, " | Cell: ", cell, " | Cat: ", category_name)
        )
        
        grid::grid.newpage()
        grid::grid.draw(combined_plot)
      }
    }
    dev.off()
    message("   Saved visualization PDF to: ", pdf_file)
  }
  
  # 执行可视化
  message(">>> Visualizing TF overlaps...")
  generate_visualization(genes_tf, "TF")
  
  message(">>> Visualizing Cre overlaps...")
  generate_visualization(genes_cre, "Cre")
  
  message(">>> Visualizing CellChat overlaps...")
  generate_visualization(genes_cc, "CellChat")
  
} else {
  warning("Seurat object 'pbmc' not found. Skipping expression visualization steps.")
}

message("✅ All tasks completed successfully.")




