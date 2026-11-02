source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 环境与依赖 ----
setwd(r_work_dir)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(ggplot2)

# 02 同一物种内合并样本并记录 sample ----
cat("步骤 1: 正在合并各物种内的样本并创建 'sample' 列...\n")
sce_merged_list <- list()
for (species_name in names(scRNAlist)) {
  sample_list <- scRNAlist[[species_name]]

  merged_obj <- NULL # 初始化一个空对象

  # 如果该物种有多个样本，则进行合并
  if (length(sample_list) > 1) {
    cat("  -> 正在合并", species_name, "的", length(sample_list), "个样本。\n")

    # 核心改动：在合并前，为每个独立的样本对象预先添加 'sample' 列
    for (sample_name in names(sample_list)) {
      sample_list[[sample_name]]$sample <- sample_name
    }

    # 现在进行合并。由于每个对象都已包含 'sample' 列，Seurat会自动保留此信息。
    merged_obj <- merge(x = sample_list[[1]], y = sample_list[-1])

  # 如果该物种只有一个样本，直接使用
  } else if (length(sample_list) == 1) {
    cat("  ->", species_name, "只有一个样本，无需合并。\n")
    merged_obj <- sample_list[[1]]

    # 对于单个样本，其 orig.ident 可能是原始长名称，
    # 因此我们直接使用其在列表中的简洁名称来创建 sample 列
    merged_obj$sample <- names(sample_list)
  }

  # 将处理好的对象存入新列表
  if (!is.null(merged_obj)) {
    sce_merged_list[[species_name]] <- merged_obj
  }
}
cat("所有物种样本合并完成, 'sample' 列已成功创建。\n\n")

# 03 计算 QC 指标 ----

# 目标：为每个细胞计算线粒体、血红蛋白、核糖体基因表达比例等关键QC指标。
# 优化：使用一个配置列表来管理不同物种的基因前缀，避免在循环中使用复杂的if/else结构。

cat("步骤 2: 正在计算QC指标...\n")

# --- QC基因模式配置 ---
# 将物种特异性的信息放在配置中，使代码更具可读性和可扩展性
qc_patterns <- list(
  Human  = list(mt = "^MT-", hb = "^HB[^(P)]", rb = "^RP[SL]"),
  Monkey = list(mt = "^MT-", hb = "^HB[^(P)]", rb = "^RP[SL]"),
  # 大鼠(Rat)的线粒体基因通常是 "mt-" (小写)
  Rat    = list(mt = "^mt-", hb = "^Hb[^(p)]", rb = "^Rp[sl]"),
  Goat   = list(mt = "^MT-", hb = "^HB[^(P)]", rb = "^RP[SL]") # 假设Goat与Human类似
)

for (species_name in names(sce_merged_list)) {
  sce <- sce_merged_list[[species_name]]
  patterns <- qc_patterns[[species_name]]

  # 计算各种基因集的百分比并存入元数据
  sce[["mt_percent"]] <- PercentageFeatureSet(sce, pattern = patterns$mt)
  sce[["hb_percent"]] <- PercentageFeatureSet(sce, pattern = patterns$hb)
  sce[["rb_percent"]] <- PercentageFeatureSet(sce, pattern = patterns$rb)

  # 新增质控指标：细胞复杂度 (log10(基因数) / log10(UMI数))
  # 复杂度低的细胞可能是低质量的，复杂度高的可能暗示双细胞
  sce[["log10_genes_per_umi"]] <- log10(sce$nFeature_RNA) / log10(sce$nCount_RNA)

  sce_merged_list[[species_name]] <- sce
}
cat("QC指标计算完成。\n\n")

# 04 过滤前 QC 可视化 ----

# 目标：通过可视化检查原始数据的质量分布，为设定过滤阈值提供依据。

cat("步骤 3: 正在生成质控前的可视化图表...\n")

# --- 小提琴图 (Violin Plots) ---
# 展示主要QC指标的分布
qc_features <- c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi")
violin_plots_before <- lapply(names(sce_merged_list), function(species) {
  VlnPlot(sce_merged_list[[species]], features = qc_features, pt.size = 0, ncol = 6, group.by = "sample") +
    plot_annotation(title = paste(species, "- 质控前分布"))
})

# --- 散点图 (Feature Scatter Plots) ---
# 查看基因数和UMI数的相关性，以及线粒体比例与其他指标的关系
scatter_plots_before <- lapply(names(sce_merged_list), function(species) {
  p1 <- FeatureScatter(sce_merged_list[[species]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample") + NoLegend()
  p2 <- FeatureScatter(sce_merged_list[[species]], feature1 = "nCount_RNA", feature2 = "mt_percent", group.by = "sample") + NoLegend()
  (p1 + p2) + plot_annotation(title = paste(species, "- 质控前相关性"))
})

# --- 新增：密度图 (Density Plots) ---
# 更清晰地查看 nFeature_RNA 和 log10_genes_per_umi 的分布，并标示计划的过滤阈值
density_plots_before <- lapply(names(sce_merged_list), function(species) {
  metadata <- sce_merged_list[[species]]@meta.data

  # 绘制 nFeature_RNA (基因数) 的密度分布
  p1 <- metadata %>%
      ggplot(aes(color=sample, x=nFeature_RNA, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      scale_x_log10() +
      geom_vline(xintercept = 300, linetype="dashed") +
      labs(title="基因数分布 (log10 scale)")

  # 绘制 log10_genes_per_umi (细胞复杂度) 的密度分布
  p2 <- metadata %>%
      ggplot(aes(x=log10_genes_per_umi, color=sample, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8, linetype="dashed") +
      labs(title="细胞复杂度分布")

  (p1 / p2) + plot_annotation(title = paste(species, "- 质控前密度图"))
})

# 打印第一个物种的图表作为示例
print(violin_plots_before[[1]])
print(scatter_plots_before[[1]])
print(density_plots_before[[1]])
cat("质控前图表已生成（已显示第一个物种作为示例）。\n\n")

# 05 按原阈值过滤并汇总细胞数 ----

# 目标：根据设定的阈值，移除低质量细胞。
# 优化：将过滤阈值定义为变量，方便调整。并创建一个质控总结表。

cat("步骤 4: 正在根据阈值过滤细胞...\n")

# --- 定义过滤参数 ---
# 将阈值放在一起，便于修改和理解
qc_thresholds <- list(
  min_features = 300,      # 每个细胞最少检测到的基因数
  max_features_quantile = 0.95, # 基因数上限，使用分位数可自适应不同样本
  max_mt_percent = 2.5,      # 线粒体基因比例上限
  min_log10_genes_per_umi = 0.85 # 细胞复杂度下限
)

# 记录质控前后的细胞数
qc_summary <- data.frame(
  Species = names(sce_merged_list),
  Cells_Before = sapply(sce_merged_list, ncol)
)

# --- 执行过滤 ---
sce_filtered_list <- list()
for (species_name in names(sce_merged_list)) {
  sce <- sce_merged_list[[species_name]]

  # 动态计算基因数上限
  max_features <- quantile(sce$nFeature_RNA, probs = qc_thresholds$max_features_quantile)

  cat("  -> 过滤", species_name, "使用阈值: nFeature_RNA >", qc_thresholds$min_features,
      "& nFeature_RNA <", round(max_features),
      "& mt_percent <", qc_thresholds$max_mt_percent, "%\n")

  # 使用 subset 函数进行过滤
  sce_filtered_list[[species_name]] <- subset(sce,
    subset = nFeature_RNA > qc_thresholds$min_features &
             nFeature_RNA < max_features &
             mt_percent < qc_thresholds$max_mt_percent &
             log10_genes_per_umi > qc_thresholds$min_log10_genes_per_umi &
                hb_percent < 0.25 &
                rb_percent < 2.5 &
                nCount_RNA < quantile(nCount_RNA,0.95) &
                nCount_RNA > 500
  )
}

# 更新总结表
qc_summary$Cells_After <- sapply(sce_filtered_list, ncol)
qc_summary$Cells_Removed <- qc_summary$Cells_Before - qc_summary$Cells_After
qc_summary$Percent_Removed <- scales::percent(qc_summary$Cells_Removed / qc_summary$Cells_Before, accuracy = 0.1)

# 打印质控总结表
cat("\n--- 质控总结 ---\n")
print(qc_summary)
cat("\n")

# 06 过滤后 QC 可视化 ----

# 目标：生成与质控前相同的图表，以检验过滤效果。

cat("步骤 5: 正在生成质控后的可视化图表...\n")

# --- 质控后小提琴图 ---
violin_plots_after <- lapply(names(sce_filtered_list), function(species) {
  VlnPlot(sce_filtered_list[[species]], features = qc_features, pt.size = 0, ncol = 6, group.by = "sample") +
    plot_annotation(title = paste(species, "- 质控后分布"))
})

# --- 质控后散点图 ---
scatter_plots_after <- lapply(names(sce_filtered_list), function(species) {
  p1 <- FeatureScatter(sce_filtered_list[[species]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample") + NoLegend()
  p2 <- FeatureScatter(sce_filtered_list[[species]], feature1 = "nCount_RNA", feature2 = "mt_percent", group.by = "sample") + NoLegend()
  (p1 + p2) + plot_annotation(title = paste(species, "- 质控后相关性"))
})

# --- 新增：质控后密度图 ---
density_plots_after <- lapply(names(sce_filtered_list), function(species) {
  metadata <- sce_filtered_list[[species]]@meta.data

  # 绘制 nFeature_RNA (基因数) 的密度分布
  p1 <- metadata %>%
      ggplot(aes(color=sample, x=nFeature_RNA, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      scale_x_log10() +
      geom_vline(xintercept = 300, linetype="dashed") +
      labs(title="基因数分布 (log10 scale)")

  # 绘制 log10_genes_per_umi (细胞复杂度) 的密度分布
  p2 <- metadata %>%
      ggplot(aes(x=log10_genes_per_umi, color=sample, fill=sample)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8, linetype="dashed") +
      labs(title="细胞复杂度分布")

  (p1 / p2) + plot_annotation(title = paste(species, "- 质控后密度图"))
})

# 打印第一个物种的质控后图表作为示例
print(violin_plots_after[[1]])
print(scatter_plots_after[[1]])
print(density_plots_after[[1]])
cat("质控后图表已生成（已显示第一个物种作为示例）。\n\n")

cat("所有质控流程已完成！\n")

# 07 导出 QC PDF ----

# 目标: 将所有生成的图表按 "质控前" 和 "质控后" 分别保存到两个PDF文件中。
# 方法: 使用 pdf() 函数打开设备，循环打印所有图表，然后用 dev.off() 关闭设备。

cat("步骤 6: 正在将所有可视化结果保存到PDF文件...\n")

# --- 保存质控前的图表 ---
pdf(file.path(figure_dir, "qc/QC_Plots_Before_Filtering.pdf"), width = 16, height = 12)
# 循环遍历每个物种的索引 (1, 2, 3, 4)
for (i in seq_along(sce_merged_list)) {
  cat("  -> 正在保存", names(sce_merged_list)[i], "的质控前图表...\n")
  # 在PDF中为每个物种打印三种图
  print(violin_plots_before[[i]])
  print(scatter_plots_before[[i]])
  print(density_plots_before[[i]])
}
dev.off() # 关闭并保存PDF文件
cat("质控前图表已保存到 'QC_Plots_Before_Filtering.pdf'\n")

# --- 保存质控后的图表 ---
pdf(file.path(figure_dir, "qc/QC_Plots_After_Filtering.pdf"), width = 16, height = 12)
# 再次循环遍历每个物种的索引
for (i in seq_along(sce_filtered_list)) {
  cat("  -> 正在保存", names(sce_filtered_list)[i], "的质控后图表...\n")
  # 在PDF中为每个物种打印三种图
  print(violin_plots_after[[i]])
  print(scatter_plots_after[[i]])
  print(density_plots_after[[i]])
}
dev.off() # 关闭并保存PDF文件
cat("质控后图表已保存到 'QC_Plots_After_Filtering.pdf'\n\n")

cat("所有质控流程已完成！\n")

# 08 保存当前会话快照 ----
save.image(file.path(r_work_dir, "qc.RData"))
