source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 载入原有依赖 ----
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(clustree)
library(DoubletFinder)
library(scDblFinder)
library(SingleCellExperiment)

# 02 各物种独立预处理与 PCA ----

# 这是所有下游分析的基础。对每个物种独立进行处理。

cat("步骤 1: 正在对每个物种进行标准化、降维...\n")
sce_processed_list <- list()
for (species_name in names(sce_filtered_list)) {
  cat("  -> 正在处理:", species_name, "\n")
  sce <- sce_filtered_list[[species_name]]

  # 1. 标准化数据: LogNormalize 方法
  sce <- NormalizeData(sce)

  # 2. 寻找高变基因: vst 方法，寻找 2000 个高变基因
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

  # 3. 缩放数据: 对所有基因进行缩放
  sce <- ScaleData(sce)

  # 4. 运行PCA
  sce <- RunPCA(sce)

  sce_processed_list[[species_name]] <- sce
}
cat("所有物种的标准化和PCA已完成。\n\n")

# 03 多样本对象执行 CCA 整合 ----

# 目标: 对于那些由多个样本合并而成的对象，使用CCA方法进行整合。

cat("步骤 1.5: 正在整合多样本数据以去除批次效应...\n")
for (species_name in names(sce_processed_list)) {
    sce <- sce_processed_list[[species_name]]

    # 检查对象中是否包含多个样本 (通过 'sample' 列的唯一值数量判断)
    if (length(unique(sce$sample)) > 1) {
        cat("  -> 发现", species_name, "有多个样本，正在进行CCA整合...\n")

        # 使用 IntegrateLayers V5 方法进行整合
        sce <- IntegrateLayers(
            object = sce,
            method = CCAIntegration,
            orig.reduction = "pca",
            new.reduction = "integrated.cca",
            verbose = FALSE
            # 按 'sample' 列进行批次校正
            #group.by = "sample"
        )
        sce_processed_list[[species_name]] <- sce
    } else {
        cat("  ->", species_name, "只有一个样本，跳过整合。\n")
    }
}
cat("数据整合步骤完成。\n\n")

# 04 固定维度、分辨率探索与初始聚类 ----

# PC的选择对下游聚类和UMAP结果有决定性影响。

cat("步骤 2: 正在确定显著PC数并进行初步聚类...\n")
elbow_plots <- list()
jackstraw_plots <- list()
clustree_plots <- list()

# 确定用于下游分析的PC维度
# 通常选择15-20个PC，但这取决于数据的复杂性。这里我们统一使用前20个PC。
DIMS_TO_USE <- 1:20

for (species_name in names(sce_processed_list)) {
  cat("  -> 正在分析:", species_name, "\n")
  sce <- sce_processed_list[[species_name]]

  # 动态选择使用的降维方法：如果进行了整合，则用 'integrated.cca'，否则用 'pca'
  reduction_method <- ifelse("integrated.cca" %in% Reductions(sce), "integrated.cca", "pca")
  cat("     将使用 '", reduction_method, "' 进行下游分析。\n")

  # 1. 可视化PC的贡献度 (Elbow Plot) - 始终基于原始PCA
  elbow_plots[[species_name]] <- ElbowPlot(sce, ndims = 30, reduction = "pca")

  # 2. JackStraw 分析 (计算量大，但更精确) - 始终基于原始PCA
  # sce <- JackStraw(sce, num.replicate = 100, dims = 30)
  # sce <- ScoreJackStraw(sce, dims = 1:30)
  # jackstraw_plots[[species_name]] <- JackStrawPlot(sce, dims = 1:30)

  # 3. 基于选定的PC进行细胞聚类 (使用正确的reduction)
  sce <- FindNeighbors(sce, dims = DIMS_TO_USE, reduction = reduction_method)
  # 测试多个分辨率，以找到最合适的分群粒度  -  直接使用 0.5
  sce_tree<- FindClusters(sce, resolution = seq(0.1, 1.0, by = 0.1))

  # 4. 使用 clustree 可视化不同分辨率下的聚类关系
  clustree_plots[[species_name]] <- clustree(sce_tree)

  # 5. 选择一个中间分辨率 (e.g., 0.5) 作为默认聚类结果
  sce <- FindClusters(sce, resolution = 0.5)

  # 6. 运行 UMAP 和 t-SNE (使用正确的reduction)
  sce <- RunUMAP(sce, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "umap.uncorrected")
  sce <- RunTSNE(sce, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "tsne.uncorrected")

  sce <- JoinLayers(sce)
  sce_processed_list[[species_name]] <- sce
  rm(sce,sce_tree)
}
cat("初步聚类和可视化已完成。\n\n")

# 05 双细胞检测：两种算法分别计算 ----

# 定义一个函数来批量处理双细胞检测，结合 scDblFinder 和 DoubletFinder 的结果。
# 注意：双细胞检测通常在整合前进行，因此仍使用 'pca' 结果。

cat("步骤 3: 正在进行双细胞检测...\n")

batch_doublet_filter <- function(sce_list, dims_use) {
  library(DoubletFinder)
  library(scDblFinder)
  doublet_plots <- list()
  sce_with_doublets <- list()

  for (species_name in names(sce_list)) {
    message(paste("  -> 正在检测", species_name, "的双细胞..."))
    sce_tmp <- sce_list[[species_name]]

    # --- Part A: 使用 scDblFinder ---
    sce_sce <- as.SingleCellExperiment(sce_tmp)
    sce_sce <- scDblFinder(sce_sce)
    sce_tmp$scDblFinder.class <- sce_sce$scDblFinder.class

    # --- Part B: 使用 DoubletFinder (基于未整合的PCA)---
    # 1. 寻找最佳 pK 值 (无 ground-truth)
    sweep.res <- paramSweep(sce_tmp, PCs = dims_use, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    best_pk_info <- find.pK(sweep.stats)
    best_pk <- as.numeric(as.character(best_pk_info$pK[which.max(best_pk_info$BCmetric)]))

    # 2. 估计双细胞比例 (基于10X Genomics的公式: 约 5% / 10,000 个细胞)
    nExp_poi <- round(0.05 * (ncol(sce_tmp) / 10000) * ncol(sce_tmp))

    # 3. 运行 DoubletFinder
    sce_tmp <- doubletFinder(sce_tmp, PCs = dims_use, pN = 0.25, pK = best_pk, nExp = nExp_poi, sct = FALSE)

    # 4. 清理元数据列名
    df_col_name <- colnames(sce_tmp@meta.data)[grepl("^DF.classifications", colnames(sce_tmp@meta.data))]
    sce_tmp$DF.classifications <- sce_tmp[[df_col_name]]
    sce_tmp[[df_col_name]] <- NULL

    # --- Part C: 生成可视化图 ---
    p1 <- DimPlot(sce_tmp, group.by = "scDblFinder.class", reduction = "umap.uncorrected") + labs(title="scDblFinder Predictions")
    p2 <- DimPlot(sce_tmp, group.by = "DF.classifications", reduction = "umap.uncorrected") + labs(title="DoubletFinder Predictions")
    doublet_plots[[species_name]] <- p1 + p2

    sce_with_doublets[[species_name]] <- sce_tmp
    rm(sce_tmp,sce_sce)
  }

  return(list(doublet_plots = doublet_plots, sce_with_doublets = sce_with_doublets))
}

# 执行双细胞检测
doublet_results <- batch_doublet_filter(sce_processed_list, dims_use = DIMS_TO_USE)
cat("双细胞检测完成。\n\n")

# 06 按 scDblFinder 结果过滤并重新聚类 ----

cat("步骤 4: 正在过滤双细胞并重新运行UMAP/t-SNE...\n")
scesinglet_list_filter <- list()
for (species_name in names(doublet_results$sce_with_doublets)) {
    cat("  -> 正在过滤:", species_name, "\n")
    sce <- doublet_results$sce_with_doublets[[species_name]]

    # 根据 scDblFinder 的判定保留 singlet。
    singlets <- subset(sce, subset = scDblFinder.class  == "singlet")

    #  基于选定的PC进行细胞聚类 (使用正确的reduction)
    singlets <- FindNeighbors(singlets, dims = DIMS_TO_USE, reduction = reduction_method)

    # 重新聚类
    singlets <- FindClusters(singlets, resolution = 0.5)

    # 再次动态选择正确的降维方法
    reduction_method <- ifelse("integrated.cca" %in% Reductions(singlets), "integrated.cca", "pca")

    # 移除双细胞后，细胞间的相对关系可能改变，因此需要重新运行UMAP/t-SNE
    singlets <- RunUMAP(singlets, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "umap")
    singlets <- RunTSNE(singlets, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "tsne")

    scesinglet_list_filter[[species_name]] <- singlets
}
cat("双细胞过滤和最终可视化已完成。\n\n")

# 07 检视结果与导出 UMAP ----
cat("步骤 5: 生成最终的可视化图表...\n")
# 示例: 打印第一个物种的双细胞检测图和过滤后的UMAP图
if (length(doublet_results$doublet_plots) > 0) {
    print(doublet_results$doublet_plots[[1]])
}
if (length(scesinglet_list_filter) > 0) {
    print(DimPlot(scesinglet_list_filter[[1]], reduction = "umap", group.by = "seurat_clusters", label = TRUE))
}

cat("所有流程已完成！\n")

p1 <- DimPlot(scesinglet_list_filter[[1]], reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  DimPlot(scesinglet_list_filter[[2]], reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  DimPlot(scesinglet_list_filter[[3]], reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  DimPlot(scesinglet_list_filter[[4]], reduction = "umap", group.by = "seurat_clusters", label = TRUE)

pdf(file = file.path(figure_dir, "tmp/4sp_sampe_cluster.pdf"),width = 8,height = 8)
print(p1)
dev.off()

p2 <- DimPlot(scesinglet_list_filter[[1]], reduction = "umap", group.by = "sample", label = TRUE) +
  DimPlot(scesinglet_list_filter[[2]], reduction = "umap", group.by = "sample", label = TRUE) +
  DimPlot(scesinglet_list_filter[[3]], reduction = "umap", group.by = "sample", label = TRUE) +
  DimPlot(scesinglet_list_filter[[4]], reduction = "umap", group.by = "sample", label = TRUE)
pdf(file = file.path(figure_dir, "tmp/4sp_sampe_umap.pdf"),width = 7,height = 7)
print(p2)
dev.off()

# 08 保存当前会话快照 ----
save.image(file.path(r_work_dir, "cluster.RData"))
