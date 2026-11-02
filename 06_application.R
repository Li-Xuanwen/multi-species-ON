source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

setwd(r_work_dir)
lapply(c('Seurat', 'dplyr', 'patchwork', 'ggplot2', 'scplotter', 'SCNT'), require, character.only = TRUE)
my_colors <- c(
  "OLG"    = "#A3C9DD",
  "AS"     = "#2272A9",
  "FB"     = "#ADD487",
  "EC"     = "#369939",
  "OPC"    = "#F5BA6E",
  "Mural"  = "#EF7C1B",
  "MG"     = "#F19695",
  "MP"     = "#0FBC82"
)

# 02 读取外部模型预测对象并整理小鼠 sample ----

# pnas
sce_pnas <- GetSeurat(file.path(scarches_dir, "pnas/pnas_peri.h5ad"))
sce_pnas_merge <- GetSeurat(file.path(scarches_dir, "pnas_merge/4sp_pnas.h5ad"))

# onc 数据
sce_onc <- GetSeurat(file.path(scarches_dir, "onc/onc.h5ad"))
sce_onc_merge <- GetSeurat(file.path(scarches_dir, "onc_merge/onc_merge.h5ad"))

# mouse
sce_mouse <- GetSeurat(file.path(scarches_dir, "mouse/mouse.h5ad"))

sce_5sp <- GetSeurat(file.path(scarches_dir, "5sp/5sp_cp.h5ad"))

sce_5sp@meta.data <- sce_5sp@meta.data %>%
  mutate(
    # 我们要修改的是 sample 列
    sample = case_when(
      sp == "mouse" & orig.ident == "MOS" ~ "Mouse3",
      sp == "mouse" & orig.ident == "MOD" ~ "Mouse4",
      orig.ident == "Ms_ON_sheath_Control1" ~ "Mouse1",
      orig.ident == "Ms_ON_sheath_Control2" ~ "Mouse2",
      # 关键：对于不满足以上所有条件的情况，保持 sample 列的原值
      TRUE ~ sample
    )
  )

# 03 外部辅助函数移除离群细胞 ----

library(dbscan)
sce_pnas_outlier_res <- remove_outliers_dbscan(sce_pnas, "orig_cell_type", "umap",minPts = 50)

sce_mouse_outlier_res <- remove_outliers_dbscan(sce_mouse, "cell_type", "umap",minPts = 50)

sce_onc_outlier_res <- remove_outliers_dbscan(sce_onc, "cell_type", "umap",minPts = 50)

sce_onc_outlier_res <- remove_outliers_dbscan(sce_onc, "cell_type", "umap",minPts = 50)

# 04 三组 UMAP 对照与九宫格导出 ----

## 3*3 拼图 ##########################################
p1 <- CellDimPlot(sce_pnas_merge,
                  group_by = "orig.ident",  reduction = "umap",legend.position = "none",
                  palcolor = c('grey','grey','grey','grey','red','red','red','red','red','grey','grey','red','red',
                  'grey','grey','grey','grey','grey','grey')
) + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                axis.text = element_blank(),
                                                panel.grid=element_blank(),
                                                panel.border = element_blank())

p2 <- CellDimPlot(sce_pnas_outlier_res$sce_clean,
                  group_by = 'orig_cell_type',legend.position = 'none',
                  palcolor = my_colors[sce_pnas$orig_cell_type   %>% table() %>% names()]) +
  labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                axis.text = element_blank(),
                                                panel.grid=element_blank(),
                                                panel.border = element_blank())

# 3. 创建第二个图 (p2)，同样存入变量
p3 <- CellDimPlot(sce_pnas,
                  group_by = 'predictions_scanvi',legend.position = 'none',
                  highlight = 'predictions_scanvi == "Mural"',
                  highlight_stroke = 1.5,
                  highlight_size = 3,
                  palcolor = my_colors[sce_pnas$predictions_scanvi %>% table() %>% names()]) +
  labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                axis.text = element_blank(),
                                                panel.grid=element_blank(),
                                                panel.border = element_blank())

p4 <- CellDimPlot(sce_5sp,
                  group_by = 'sp',reduction = 'umap',legend.position = 'none',
                  palcolor = c('grey','grey','grey','red','grey')
) + xlab('') + ylab('') + theme(axis.ticks = element_blank(), axis.text = element_blank()) +
  labs(subtitle = NULL)  + labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                                         axis.text = element_blank(),
                                                                         panel.grid=element_blank(),
                                                                         panel.border = element_blank())
p5 <- CellDimPlot(sce_mouse_outlier_res$sce_clean,
                        group_by = "cell_type", reduction = "umap", alpha = 0.1,
                  legend.position = "none",
                        palcolor = my_colors[levels(sce_mouse$cell_type)],
                        label_insitu = F
) + labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                axis.text = element_blank(),
                                                panel.grid=element_blank(),
                                                panel.border = element_blank())

p6 <- CellDimPlot(sce_mouse,
                  group_by = "predictions_scanvi", reduction = "umap", alpha = 0.1,
                  legend.position = "none",
                  palcolor = my_colors[levels(sce_mouse$predictions_scanvi)],
                  label_insitu = F
) + labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                  axis.text = element_blank(),
                                                  panel.grid=element_blank(),
                                                  panel.border = element_blank())

p7 <- CellDimPlot(sce_onc_merge,
            group_by = 'orig.ident',reduction = 'umap',legend.position = 'none',
            palcolor = c('grey','grey','grey','grey','grey','grey','grey','grey','red','grey',
          'grey','grey','grey')
) + xlab('') + ylab('')+
  labs(subtitle = NULL) +theme(axis.ticks = element_blank(),
                                                  axis.text = element_blank(),
                                                  panel.grid=element_blank(),
                                                  panel.border = element_blank())

p8 <- CellDimPlot(sce_onc_outlier_res$sce_clean,
                  group_by = "cell_type", reduction = "umap", alpha = 0.1,
                  legend.position = 'none',
                  palcolor = my_colors[levels(sce_onc$cell_type)],
                  label_insitu = F
) + labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                                         axis.text = element_blank(),
                                                                         panel.grid=element_blank(),
                                                                         panel.border = element_blank())

p9 <- CellDimPlot(sce_onc_outlier_res$sce_clean,
                  group_by = "predictions_scanvi", reduction = "umap", alpha = 0.1,
                  legend.position = 'none',
                  highlight = 'predictions_scanvi == "MG"',
                  highlight_stroke = 1.5,
                  highlight_size = 3,
                  palcolor = my_colors[levels(sce_onc$predictions_scanvi)],
                  label_insitu = F
) + labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                                         axis.text = element_blank(),
                                                                         panel.grid=element_blank(),
                                                                         panel.border = element_blank())
library(ggplot2)
library(patchwork)

# 1. 创建并修改好您的最终图形对象
umap9_plot <- (p1 + p2 + p3) /
              (p4 + p5 + p6) /
              (p7 + p8 + p9)

final_plot <- umap9_plot + theme(plot.margin = margin(1, 1, 1, 1, unit = "mm"))

# 2. 使用 ggsave() 一行代码完成保存
ggsave(
  filename = file.path(figure_dir, "fig3/9umap_new.pdf"),
  plot = final_plot,
  width = 5,
  height = 5,
  units = "in" # 明确指定单位是英寸，更清晰
)

# 05 继承注释对象的 marker 检查 ----

DotPlot(object =scesinglet_list_filter[[1]],
        features = human_markers ,
        cols = c("lightgrey", "blue"),
        group.by = 'cell_type',
        #assay = "RNA",
        dot.scale = 4, 
        cluster.idents = FALSE, 
        scale.by = "radius",  
        scale = TRUE,   
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of human marker")

#save.image(file.path(r_work_dir, "3appli.RData"))

# 06 按原筛选口径计算预测准确率 ----

mouse_normal <- data.frame(cell_type = sce_mouse$cell_type,
                           cell_type_peri = sce_mouse$predictions_scanvi)
mouse_onc <- data.frame(cell_type = subset(sce_onc_merge, orig.ident == 'crush2.1wpi')$cell_type,
                        cell_type_peri = subset(sce_onc_merge, orig.ident == 'crush2.1wpi')$predictions_scanvi) %>%
  filter(cell_type_peri != 'MG')
pnas <- data.frame(cell_type = sce_pnas$orig_cell_type ,
                       cell_type_peri = sce_pnas$predictions_scanvi) %>%
  filter(cell_type_peri != 'Mural')

# --- 3. 核心步骤：计算按细胞类型分类的准确率和总准确率 ---
calculate_accuracy <- function(df, sample_name) {
  # 预处理：将因子转为字符
  df_char <- df %>%
    mutate(
      cell_type = as.character(cell_type),
      cell_type_peri = as.character(cell_type_peri)
    )

  # 计算按细胞类型分组的准确率
  grouped_acc <- df_char %>%
    group_by(cell_type) %>%
    summarise(accuracy = sum(cell_type == cell_type_peri) / n(), .groups = 'drop')

  # 计算总准确率
  total_acc <- df_char %>%
    summarise(accuracy = sum(cell_type == cell_type_peri) / n()) %>%
    mutate(cell_type = "Total")

  # 合并分组结果和总计结果
  bind_rows(total_acc, grouped_acc) %>%
    mutate(sample = sample_name)
}

# 应用函数到每个数据集
accuracy_normal <- calculate_accuracy(mouse_normal, "Mouse Normal")
accuracy_onc <- calculate_accuracy(mouse_onc, "Mouse Onc")
accuracy_pnas <- calculate_accuracy(pnas, "PNAS")

# 07 汇总及绘制准确率柱图 ----

combined_accuracy <- bind_rows(accuracy_normal, accuracy_onc, accuracy_pnas)

all_cell_types <- unique(combined_accuracy$cell_type)
ordered_levels <- c("Total", sort(all_cell_types[all_cell_types != "Total"]))
combined_accuracy$cell_type <- factor(combined_accuracy$cell_type, levels = ordered_levels)
sample_ordered_levels <- c('PNAS','Mouse Normal','Mouse Onc')
combined_accuracy$sample <- factor(combined_accuracy$sample, levels = sample_ordered_levels)

# 将 sample 列也转换为因子，以便控制绘图时的顺序和颜色
combined_accuracy$sample <- factor(combined_accuracy$sample, levels = c( "PNAS","Mouse Normal", "Mouse Onc"))

combined_accuracy <- combined_accuracy %>%
  mutate(
    label_y_position = case_when(
      sample == "Mouse Normal" ~ 1 + 0.04,
      sample == "Mouse Onc"    ~ 1 + 0.11,
      sample == "PNAS"         ~ 1 + 0.18,
      TRUE                     ~ 1      
    )
  )

library(scales)
pdf(file = file.path(figure_dir, "fig3/Accuray_bar_new.pdf"),width = 5.7,height = 3.7)
# --- 5. 使用 ggplot2 绘制分组柱状图 ---
ggplot(combined_accuracy, aes(x = cell_type, y = accuracy, fill = sample)) +
  # 使用 geom_bar 并设置 stat = "identity" 来绘制已有y值（而不是计数）的柱状图
  # position = "dodge" 使不同样本的柱子并排显示
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +

  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey40") +
  geom_text(
    aes(y = label_y_position, label = percent(accuracy, accuracy = 0.01)),
    position = position_dodge(width = 0.9),
    size = 3.5
  ) +

  scale_fill_brewer(palette = "Pastel1") +

  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.20), expand = c(0, 0)) +

  labs(
    title = "Prediction Accuracy by Cell Type and Sample",
    #subtitle = "Comparison of predicted cell types vs. ground truth, including total accuracy",
    x = "Cell Type",
    y = "Accuracy",
    fill = "Sample" 
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30", linewidth = 0.5) +

  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text( hjust = 1, size = 11, color = 'black'), # 旋转x轴标签以防重叠
    axis.text.y = element_text(size = 11,color = 'black'),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "top" # 将图例放在顶部
)
dev.off()

# 08 保存工作区并查看五物种预测 ----
save.image(file.path(r_work_dir, "3appli.RData"))

CellDimPlot(sce_5sp,
                  group_by = 'predictions_scanvi',
                  palcolor = my_colors[sce_5sp$predictions_scanvi %>% table() %>% names()]) +
  labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                axis.text = element_blank(),
                                                panel.grid=element_blank(),
                                                panel.border = element_blank())
