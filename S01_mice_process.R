rm(list =ls())
source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

setwd(r_work_dir)

# 加载所需包 ####
lapply(c("Seurat","dplyr","ggplot2","scplotter","tidyr","ggpubr","stringr","patchwork","SeuratObject","scDblFinder","DoubletFinder","dbscan"), require, character.only = TRUE)

# 01 读取小鼠原始矩阵 ----

mouse_path = file.path(matrix_dir, "mouse")

matrix_paths <- list.files(
    path = mouse_path,
    pattern = "filtered_feature_bc_matrix",
    recursive = TRUE,    # 递归查找，深入所有子目录
    full.names = TRUE,   # 返回完整路径
    include.dirs = TRUE  # 确保返回的是目录
  )

 # 初始化一个临时的、属于当前物种的列表
sce_mouse_list <- list()

  # 遍历找到的所有样本路径
for (data_path in matrix_paths) {

    # 提取样本名 (通常是 filtered_feature_bc_matrix 的上一级目录名)
    sample_name <- basename(dirname(data_path))
    cat("  -> 正在读取样本:", sample_name, "\n")

    # 读取 10x 数据
    counts <- Read10X(data.dir = data_path)

    # 创建 Seurat 对象
    seurat_obj <- CreateSeuratObject(counts, project = sample_name, min.cells = 3, min.features = 300)

    # 将创建好的对象添加到物种列表中
    sce_mouse_list[[sample_name]] <- seurat_obj
  names(sce_mouse_list) <- paste0('mouse', seq_along(sce_mouse_list))
}

names(sce_mouse_list) <- paste0('mouse', seq_along(sce_mouse_list))

# 02 合并样本与计算质控指标 ----

sce_mice <- merge(x = sce_mouse_list[[1]], y = sce_mouse_list[-1])

# 计算质控指标
# --- QC基因模式配置 ---
# 将物种特异性的信息放在配置中，使代码更具可读性和可扩展性
qc_patterns <- list(list(mt = "^Mt-", hb = "^Hb[^(p)]", rb = "^Rp[sl]"))

# 计算各种基因集的百分比并存入元数据
sce_mice[["mt_percent"]] <- PercentageFeatureSet(sce_mice, pattern = qc_patterns)
sce_mice[["hb_percent"]] <- PercentageFeatureSet(sce_mice, pattern = qc_patterns)
sce_mice[["rb_percent"]] <- PercentageFeatureSet(sce_mice, pattern = qc_patterns)

  # 新增质控指标：细胞复杂度 (log10(基因数) / log10(UMI数))
  # 复杂度低的细胞可能是低质量的，复杂度高的可能暗示双细胞
sce_mice[["log10_genes_per_umi"]] <- log10(sce_mice$nFeature_RNA) / log10(sce_mice$nCount_RNA)

##### 质控前可视化 #####

# --- 小提琴图 (Violin Plots) ---
# 展示主要QC指标的分布
qc_features <- c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi")
violin_plots_before <- VlnPlot(sce_mice, features = qc_features, pt.size = 0, ncol = 3) +
    labs(title = paste('mouse', "- 质控前分布"))

# --- 散点图 (Feature Scatter Plots) ---
# 查看基因数和UMI数的相关性，以及线粒体比例与其他指标的关系

p1 <- FeatureScatter(sce_mice, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
p2 <- FeatureScatter(sce_mice, feature1 = "nCount_RNA", feature2 = "mt_percent") + NoLegend()
scatter_plots_before <- (p1 + p2) + plot_annotation(title = paste("mouse", "- 质控前相关性"))

# --- 新增：密度图 (Density Plots) ---
# 更清晰地查看 nFeature_RNA 和 log10_genes_per_umi 的分布，并标示计划的过滤阈值

metadata <- sce_mice@meta.data

# 绘制 nFeature_RNA (基因数) 的密度分布
p1 <- metadata %>%
      ggplot(aes(color=orig.ident, x=nFeature_RNA, fill=orig.ident)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      scale_x_log10() +
      geom_vline(xintercept = 300, linetype="dashed") +
      labs(title="基因数分布 (log10 scale)")

# 绘制 log10_genes_per_umi (细胞复杂度) 的密度分布
p2 <- metadata %>%
      ggplot(aes(x=log10_genes_per_umi, color=orig.ident, fill=orig.ident)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8, linetype="dashed") +
      labs(title="细胞复杂度分布")

  (p1 / p2) + plot_annotation(title = paste("mouse", "- 质控前密度图"))

# 打印第一个物种的图表作为示例
print(violin_plots_before[[1]])
print(scatter_plots_before[[1]])
print(density_plots_before[[1]])
cat("质控前图表已生成（已显示第一个物种作为示例）。\n\n")

# 03 按阈值筛选细胞 ----

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
  Species = "mouse",
  Cells_Before = ncol(sce_mice)
)

# --- 执行过滤 ---
sce_mice_filtered

  # 动态计算基因数上限
  max_features <- quantile(sce_mice$nFeature_RNA, probs = qc_thresholds$max_features_quantile)

  cat("  -> 过滤", species_name, "使用阈值: nFeature_RNA >", qc_thresholds$min_features,
      "& nFeature_RNA <", round(max_features),
      "& mt_percent <", qc_thresholds$max_mt_percent, "%\n")

  # 使用 subset 函数进行过滤
  sce_mice_filtered <- subset(sce_mice,
    subset = nFeature_RNA > qc_thresholds$min_features &
             nFeature_RNA < max_features &
             mt_percent < qc_thresholds$max_mt_percent &
             log10_genes_per_umi > qc_thresholds$min_log10_genes_per_umi &
                hb_percent < 0.25 &
                rb_percent < 2.5 &
                nCount_RNA < quantile(nCount_RNA,0.95) &
                nCount_RNA > 500
  )

# 更新总结表
qc_summary$Cells_After <- ncol(sce_mice_filtered)
qc_summary$Cells_Removed <- qc_summary$Cells_Before - qc_summary$Cells_After
qc_summary$Percent_Removed <- scales::percent(qc_summary$Cells_Removed / qc_summary$Cells_Before, accuracy = 0.1)

# 打印质控总结表
cat("\n--- 质控总结 ---\n")
print(qc_summary)
cat("\n")

#### 5. 质控后可视化 ####
# 目标：生成与质控前相同的图表，以检验过滤效果。

cat("步骤 5: 正在生成质控后的可视化图表...\n")

# --- 质控后小提琴图 ---
violin_plots_after <- VlnPlot(sce_mice_filtered, features = qc_features, pt.size = 0, ncol = 6) +
    plot_annotation(title = paste(species, "- 质控后分布"))

# --- 质控后散点图 ---
p1 <- FeatureScatter(sce_mice_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
  p2 <- FeatureScatter(sce_mice_filtered, feature1 = "nCount_RNA", feature2 = "mt_percent") + NoLegend()
scatter_plots_after <- (p1 + p2) + plot_annotation(title = paste("mouse", "- 质控后相关性"))

# --- 新增：质控后密度图 ---

  metadata <- sce_mice_filtered@meta.data

  # 绘制 nFeature_RNA (基因数) 的密度分布
  p1 <- metadata %>%
      ggplot(aes(color=orig.ident, x=nFeature_RNA, fill=orig.ident)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      scale_x_log10() +
      geom_vline(xintercept = 300, linetype="dashed") +
      labs(title="基因数分布 (log10 scale)")

  # 绘制 log10_genes_per_umi (细胞复杂度) 的密度分布
  p2 <- metadata %>%
      ggplot(aes(x=log10_genes_per_umi, color=orig.ident, fill=orig.ident)) +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = 0.8, linetype="dashed") +
      labs(title="细胞复杂度分布")

density_plots_after <- (p1 / p2) + plot_annotation(title = paste("mouse", "- 质控后密度图"))

# 04 标准化、整合与聚类 ----

sce_mice_filtered_processed <-sce_mice_filtered %>%
 NormalizeData() %>%
 FindVariableFeatures(selection.method = "vst", nfeatures = 2000)%>%
  ScaleData() %>%
  RunPCA()

##### cca 合并 #####
cat("步骤 1.5: 正在整合多样本数据以去除批次效应...\n")

sce_mice_filtered_processed_cca <- IntegrateLayers(
        object = sce_mice_filtered_processed,
        method = CCAIntegration,
        orig.reduction = "pca",
        new.reduction = "integrated.cca",
        verbose = FALSE
        # 按 'sample' 列进行批次校正
         #group.by = "sample"
)

elbow_plots
jackstraw_plots <- list()
clustree_plots <- list()

# 通常选择15-20个PC，但这取决于数据的复杂性。这里我们统一使用前20个PC。
DIMS_TO_USE <- 1:20

# 动态选择使用的降维方法：如果进行了整合，则用 'integrated.cca'，否则用 'pca'
reduction_method <- ifelse("integrated.cca" %in% Reductions(sce_mice_filtered_processed_cca), "integrated.cca", "pca")

# 1. 可视化PC的贡献度 (Elbow Plot) - 始终基于原始PCA
elbow_plots <- ElbowPlot(sce_mice_filtered_processed_cca, ndims = 30, reduction = "pca")

sce_mice_filtered_processed_cca <- FindNeighbors(sce_mice_filtered_processed_cca, dims = DIMS_TO_USE, reduction = reduction_method)
# 测试多个分辨率，以找到最合适的分群粒度  -  直接使用 0.5
sce_tree <- FindClusters(sce_mice_filtered_processed_cca, resolution = seq(0.1, 1.0, by = 0.1))

# 4. 使用 clustree 可视化不同分辨率下的聚类关系
clustree_plots <- clustree(sce_tree)

# 5. 选择一个中间分辨率 (e.g., 0.5) 作为默认聚类结果
sce_mice_filtered_processed_cca <- FindClusters(sce_mice_filtered_processed_cca, resolution = 0.5)

# 6. 运行 UMAP 和 t-SNE (使用正确的reduction)
sce_mice_filtered_processed_cca <- RunUMAP(sce_mice_filtered_processed_cca, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "umap.uncorrected")
sce_mice_filtered_processed_cca <- RunTSNE(sce_mice_filtered_processed_cca, dims = DIMS_TO_USE, reduction = reduction_method, reduction.name = "tsne.uncorrected")
sce_mice_filtered_processed_cca <- JoinLayers(sce_mice_filtered_processed_cca)

# 05 双细胞检测与过滤 ----

sce_list <- list(
  mouse = sce_mice_filtered_processed_cca
)
doublet_results <- batch_doublet_filter(sce_list, dims_use = DIMS_TO_USE)

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

# 06 人工细胞类型注释 ----

DimPlot(scesinglet_list_filter[[1]], reduction = "umap", group.by = "orig.ident", label = TRUE) +
  DimPlot(scesinglet_list_filter[[1]], reduction = "umap",  label = TRUE)
mouse_markers <- list(Mural = c("Acta2", "Myh11" ),
                      OPC = c("Pcdh15","Tnr",'Pdgfra'), # noncanonical
                      EC = c("Flt1", "Erg"),
                      MG = c('P2ry12','Hexb','Siglech'),
                      #Skeletal_Muscle_Cells = c('Myh4','Tnnt3'),
                      #Schwann_Cells = c('Mpz'),
                      #Bcell = c('Cd19','Ms4a1','Cd22','Bank1','Cd79a','Cd20'),
                      #Tcell = c('Cd8','Cd3',"Cd4",'Trac'),
                      MP = c('Mrc1','Cd163','F13a1'),
                      #LECs = c('Lyve1','Prox1'),
                      #Epithelial_Cells = c('Cdh1','Pkp2','Cxadr'),
                      FB = c('Dcn','Col1a1','Col3a1','Bicc1','Eya2'),
                      #VLMC = c('Slc47a1', 'Sema3g','Slit3','Slit2'),
                      AS = c("Gfap",'Aqp4'),
                      #Adipocyte = c("Pnpla2", "Plin1"),
                      OLG = c("Mobp","Mog"))

DotPlot(object = scesinglet_list_filter[[1]],
        features = mouse_markers ,
        cols = c("lightgrey", "blue"),
        #assay = "RNA",
        dot.scale = 4,   # 控制点的大小比例
        cluster.idents = FALSE,   # 是否将细胞群按照树状结构聚类
        scale.by = "radius",          # 控制点大小的标准化方式
        scale = TRUE,                 # 是否对数据进行标准化
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of mouse marker")

VlnPlot(scesinglet_list_filter[[1]], features =  c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi"), pt.size = 0, group.by = "RNA_snn_res.0.5")

mouse_cell_type_map <- c(
  "0" = "OLG",
  "1" = "OLG",
  "2" = "AS",
  "3" = "FB",
  "4" = "FB",
  "5" = "FB", 
  "6" = "FB",
  "7" = "EC",
  "8" = "EC",
  "9" = "MG",
  "10" = "OPC",
  "11" = "MP",
  "12" = "EC",
  "13" = "None", # 删除SMC
  "14" = "None", # 删除 Epithelial_Cells
  "15" = "Mural", #
  "16" = "None",  # 删除 LECs
  "17" = "None", # 删除 mt群
  "18" = "None", # Neuron
  "19" = "None", # Schwann_Cells
  "20" = "None",
  "21" = "FB", #
  "22" = "None", # Adipocyte
  "23" = "None" # 太少
)

##### findmarkers 确定存疑的亚群 #####
markers_list <- lapply(scesinglet_list_filter, function(df) {
   FindAllMarkers(df, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
})

##### 确定并查看结果 #####
scesinglet_list_filter[[1]]$cell_type <- setNames(mouse_cell_type_map[as.character(scesinglet_list_filter[[1]]$RNA_snn_res.0.5)],
 Cells(scesinglet_list_filter[[1]]))

# 07 异常群与离群细胞过滤 ----

remove_cluster_list <- list(
  mouse = c(13,14,16,17,18,19,20,22,23)
)

sce_mice_clean <- subset(scesinglet_list_filter[[1]], subset = seurat_clusters %in% remove_cluster_list$mouse, invert = TRUE)

DotPlot(object = sce_mice_clean,features = mouse_markers ,group.by = 'cell_type',
 dot.scale = 4, cluster.idents = FALSE, scale.by = "radius", scale = TRUE, col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of mouse marker") +
  DimPlot(sce_mice_clean,group.by = 'cell_type')

#### 4. 去除离群值 ####
sce_rmoutlies_list <- remove_outliers_dbscan(sce_mice_clean, "cell_type", "umap",minPts = 50)

# 08 汇总并保存小鼠对象 ----

qc_summary <- data.frame(
  Species = 'mouse',
  Cells_Before = ncol(sce_mice)
)

# 更新总结表
qc_summary$Cells_After <- ncol(sce_rmoutlies_list$sce_clean)
qc_summary$Cells_Removed <- qc_summary$Cells_Before - qc_summary$Cells_After
qc_summary$Percent_Removed <- scales::percent(qc_summary$Cells_Removed / qc_summary$Cells_Before, accuracy = 0.1)

# 打印质控总结表
cat("\n--- 质控总结 ---\n")
print(qc_summary)
cat("\n")

saveRDS(sce_rmoutlies_list$sce_clean,file = file.path(r_work_dir, "mice_seu.RDS"))
save.image(file.path(r_work_dir, "mice.RData"))
