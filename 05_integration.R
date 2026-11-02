source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 环境及外部小鼠对象 ----
setwd(file.path(r_work_dir, "revise"))
#load(file.path(r_work_dir, "revise/integration/integrate.RData"))
lapply(c("biomaRt","Seurat","dplyr","ggplot2","scplotter","tidyr","ggpubr","stringr","patchwork",
"SeuratObject","scDblFinder","DoubletFinder","dbscan"), require, character.only = TRUE)

sce_rmoutlie_list$Mouse <- readRDS(file.path(r_work_dir, "mice_seu.RDS"))

# 02 从归档数据库获取同源关系 ----

##### 设置物种的数据库 #####
human_ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host = 'https://dec2021.archive.ensembl.org/') # https://jan2024.archive.ensembl.org/
mouse_ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl",host = 'https://dec2021.archive.ensembl.org/')
goat_ensembl <- useMart("ensembl",dataset="chircus_gene_ensembl",host = 'https://dec2021.archive.ensembl.org/')
Macaque_ensembl <- useMart("ensembl",dataset="mmulatta_gene_ensembl",host = 'https://dec2021.archive.ensembl.org/')
Rat_ensembl <- useMart("ensembl",dataset="rnorvegicus_gene_ensembl",host = 'https://dec2021.archive.ensembl.org/')

# 大鼠基因符号（RGD）映射到人类基因符号（HGNC）。
RtoH <- getLDS(attributes = "rgd_symbol",  
               filters = "rgd_symbol", 
               mart = Rat_ensembl,
               values = rownames(sce_rmoutlie_list[[3]]),
               attributesL = "hgnc_symbol",  
               martL = human_ensembl,
               uniqueRows = TRUE) %>%
  dplyr::filter(!grepl("-", HGNC.symbol)) %>%
  dplyr::group_by(RGD.symbol) %>%
  dplyr::filter(n() == 1) %>% 
  dplyr::ungroup() %>% dplyr::group_by(HGNC.symbol) %>%
  dplyr::slice(1) %>% 
  dplyr::ungroup()%>% as.data.frame()

GtoH <- getLDS(attributes = "hgnc_symbol", 
               filters = "hgnc_symbol", 
               mart = goat_ensembl, 
               values = rownames(sce_rmoutlie_list[[4]]), 
               attributesL = "hgnc_symbol", 
               martL = human_ensembl,
               uniqueRows = TRUE)

listAttributes(Macaque_ensembl) %>% dplyr::filter(grepl("symbol", name)) 
MatoH <- getLDS(attributes = "hgnc_symbol", 
               filters = "hgnc_symbol",
               mart = Macaque_ensembl, 
               values = rownames(sce_rmoutlie_list[[2]]), 
               attributesL = "hgnc_symbol", 
               martL = human_ensembl, 
               uniqueRows = TRUE)

mtoH <- getLDS(attributes = "mgi_symbol", 
               filters = "mgi_symbol", 

               mart = mouse_ensembl, 
               values = rownames(sce_rmoutlie_list[[5]]),
               attributesL = "hgnc_symbol", 
               martL = human_ensembl, 
               uniqueRows = TRUE) %>%
  dplyr::filter(!grepl("-", HGNC.symbol)) %>%
  dplyr::group_by(MGI.symbol) %>%
  dplyr::filter(n() == 1) %>% 
  dplyr::ungroup() %>% dplyr::group_by(HGNC.symbol) %>%
  dplyr::slice(1) %>% 
  dplyr::ungroup()%>% as.data.frame()

# 03 同源映射应用与预处理 ----

rownames(GtoH) <- GtoH$HGNC.symbol
GoattoH_sce <- subset(sce_rmoutlie_list[[4]], features = GtoH$HGNC.symbol)
rownames(GoattoH_sce) <- GtoH[rownames(GoattoH_sce),] %>% pull('HGNC.symbol.1')

rownames(RtoH) <- RtoH$RGD.symbol
RattoH_sce <- subset(sce_rmoutlie_list[[3]], features = RtoH$RGD.symbol)
rownames(RattoH_sce) <- RtoH[rownames(RattoH_sce),] %>% pull('HGNC.symbol')

rownames(mtoH) <- mtoH$MGI.symbol
MousetoH_sce <- subset(sce_rmoutlie_list[[5]], features = mtoH$MGI.symbol)
rownames(MousetoH_sce) <- mtoH[rownames(MousetoH_sce),] %>% pull('HGNC.symbol')

GoattoH_sce <- GoattoH_sce %>% NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = T)

RattoH_sce <- RattoH_sce %>% NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%
  ScaleData() %>% 
  RunPCA(npcs = 50, verbose = T)

# 04 建立五物种列表并取基因交集 ----

sce_list <- list(
  Mouse = MousetoH_sce,
  Rat = RattoH_sce,
  Goat = GoattoH_sce,
  Macaque = sce_rmoutlie_list$Monkey,
  Human = sce_rmoutlie_list$Human
)
cat("Your data has been loaded and structured into a list of 5 objects.\n\n")

# 2.5. FIND COMMON GENES AND SUBSET OBJECTS
# --------------------------------------------------------------------------
# To ensure compatibility, we will find the intersection of genes across all
# datasets and subset each object to only contain these shared genes.
cat("Step 2: Finding common genes across all 5 datasets...\n")

# Get a list of all gene names from each object
gene_lists <- lapply(sce_list, rownames)

# Find the intersection of all gene lists
common_genes <- Reduce(intersect, gene_lists)
cat(paste("Found", length(common_genes), "common genes.\n"))

# Subset each Seurat object in the list to keep only the common genes
sce_list_filtered <- lapply(sce_list, function(obj) {
  obj <- subset(obj, features = common_genes)
  return(obj)
})
cat("All Seurat objects have been subset to the common gene set.\n\n")

# 05 山羊 H5 导出及中间合并 ----
sce_goat_cp <- GoattoH_sce
sce_goat_cp[["RNA"]] <- as(object = sce_goat_cp[["RNA"]], Class = "Assay")
dior::write_h5(sce_goat_cp,file=file.path(sccaf_dir, "scRNA_goat_trans.h5"))
rm(sce_goat_cp)

sce_merge <- merge(x=RattoH_sce,y=c(GoattoH_sce,sce_rmoutlie_list$Monkey,sce_rmoutlie_list$Human),
             add.cell.ids = c("Rat", "Goat","Macaque","Human"))

# 06 整合方法依赖与公共预处理 ----

# ---
# Title: "Single-Cell Data Integration and Benchmarking"
# Author: Gemini
# Date: 2024-09-18
# Description: This script benchmarks various Seurat integration methods
#              (CCA, RPCA, Harmony, FastMNN, LIGER) by timing their execution.
# ---

# 1. INSTALL AND LOAD LIBRARIES
# --------------------------------------------------------------------------
# Ensure all necessary packages are installed
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(Seurat, SeuratWrappers, harmony, rliger, dplyr, patchwork)

# 2. LOAD AND PREPARE YOUR DATA
# --------------------------------------------------------------------------
# This section has been modified to use your actual data objects.
# IMPORTANT: Ensure your Seurat objects (RattoH_sce, GoattoH_sce, etc.)
# are already loaded into your R environment before running this script.
# For example, you might load them from .Rds files like this:
# RattoH_sce <- readRDS("path/to/your/RattoH_sce.rds")
# GoattoH_sce <- readRDS("path/to/your/GoattoH_sce.rds")
# sce_rmoutlie_list <- readRDS("path/to/your/sce_rmoutlie_list.rds")

# 3. STANDARD PRE-PROCESSING (COMMON FOR MOST METHODS)
# --------------------------------------------------------------------------
# This pre-processing is required for methods like Harmony, RPCA, and CCA (v5).
cat("Step 2: Performing standard pre-processing...\n")
sce_merge <- merge(x = sce_list_filtered[[5]], y = c(sce_list_filtered[[4]], sce_list_filtered[[3]], sce_list_filtered[[2]]))
sce_merge <- NormalizeData(sce_merge) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE)
cat("Standard pre-processing complete.\n\n")

# 07 各整合方法运行及计时 ----

# --------------------------------------------------------------------------
# Initialize a list to store the timing results for each method.
timing_results <- list()

# --- Method 1: Seurat v3 CCA ---
cat("Running Integration Method 1: Seurat v3 CCA...\n")
timing_results$CCA_v3 <- system.time({
  # Pre-process individual objects for this method
  sce_list_cca <- lapply(sce_list_filtered[2:5], function(obj) {
      NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=F)
  })

  anchors <- FindIntegrationAnchors(object.list = sce_list_cca,
                                    dims = 1:30, anchor.features = 2000, reduction = 'cca')
  integrated_cca_v3 <- IntegrateData(anchorset = anchors)

  # Post-integration analysis
  DefaultAssay(integrated_cca_v3) <- "integrated"
  integrated_cca_v3 <- ScaleData(integrated_cca_v3, verbose = FALSE)
  integrated_cca_v3 <- RunPCA(integrated_cca_v3, npcs = 30, verbose = FALSE)
  integrated_cca_v3 <- RunUMAP(integrated_cca_v3, reduction = "pca", dims = 1:30)
})
cat("Seurat v3 CCA complete.\n\n")

# --- Method 2: Seurat v5 CCAIntegration ---
cat("Running Integration Method 2: Seurat v5 CCA...\n")
timing_results$CCA_v5 <- system.time({
  sce_merge_cca5 <- IntegrateLayers(object = sce_merge, method = CCAIntegration,
                                    orig.reduction = "pca", new.reduction = "integrated.cca",
                                    verbose = FALSE)
  sce_merge_cca5 <- FindNeighbors(sce_merge_cca5, reduction = "integrated.cca", dims = 1:20)
  sce_merge_cca5 <- RunUMAP(sce_merge_cca5, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
})
cat("Seurat v5 CCA complete.\n\n")

# --- Method 3: Harmony ---
cat("Running Integration Method 3: Harmony...\n")
timing_results$Harmony <- system.time({
  sce_merge_harmony <- IntegrateLayers(object = sce_merge, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.harmony",
                                       verbose = FALSE)
  sce_merge_harmony <- RunUMAP(sce_merge_harmony, reduction = "integrated.harmony", dims = 1:20, reduction.name = "umap.harmony")
})
cat("Harmony complete.\n\n")

# --- Method 4: FastMNN ---
cat("Running Integration Method 4: FastMNN...\n")

sce_fastmnn <- merge(sce_list_filtered[[5]], y = c(sce_list_filtered[[4]], sce_list_filtered[[3]], sce_list_filtered[[2]]),
             add.cell.ids = c("Rat", "Goat","Macaque","Human"))

timing_results$FastMNN <- system.time({

  sce_fastmnn <- sce_fastmnn %>% NormalizeData() %>% FindVariableFeatures()
  sce_fastmnn_res <- RunFastMNN(object.list = SplitObject(sce_fastmnn,split.by = 'orig.ident'))
  sce_fastmnn_res <- RunUMAP(sce_fastmnn_res, reduction = "mnn", dims = 1:20)

})

cat("FastMNN complete.\n\n")

# --- Method 5: RPCA ---
cat("Running Integration Method 5: RPCA...\n")
# Increase memory limit if needed
options(future.globals.maxSize = 16 * 1024^3) # 16GB
timing_results$RPCA <- system.time({
  sce_merge_rpca <- IntegrateLayers(object = sce_merge, method = RPCAIntegration,
                                    orig.reduction = "pca", new.reduction = "integrated.rpca",
                                    verbose = FALSE)
  sce_merge_rpca <- RunUMAP(sce_merge_rpca, reduction = "integrated.rpca", dims = 1:20, reduction.name = "umap.rpca")
})
cat("RPCA complete.\n\n")

# --- Method 6: LIGER ---
cat("Running Integration Method 6: LIGER...\n")
timing_results$LIGER <- system.time({
  # LIGER requires its own preprocessing steps
  sce_liger <- sce_merge %>%
    normalize() %>%
    selectGenes() %>%
    scaleNotCenter() %>%
    runINMF(k = 20) %>%
    quantileNorm()
    # To visualize with UMAP in Seurat, you would need to add the liger reduction back
    # sce_merge[['inmf']] <- create.dim.reduc(sce_liger@H.norm, key = "inmf_")
    # sce_merge <- RunUMAP(sce_merge, reduction = "inmf", dims = 1:20)
  sce_liger <- RunUMAP(sce_liger, reduction = "inmfNorm", dims = 1:20)
})
cat("LIGER complete.\n\n")

# 08 汇总实测运行时间与检查嵌入 ----

# --------------------------------------------------------------------------
cat("=====================================================\n")
cat("          Integration Benchmarking Results           \n")
cat("=====================================================\n\n")

# Convert the list of timings to a more readable data frame
timing_df <- do.call(rbind, lapply(timing_results, function(x) {
  data.frame(
    User_Time_s = x['user.self'],
    System_Time_s = x['sys.self'],
    Elapsed_Time_s = x['elapsed']
  )
}))

# scvi 582.49
# scanvi 855.19

rownames(timing_df) <- names(timing_results)
print(timing_df)

cat("\n'Elapsed_Time_s' is the total wall-clock time taken for the computation.\n")

# You can also generate plots to visualize the results
# For example, a DimPlot for one of the methods:
# DimPlot(sce_merge_harmony, reduction = 'umap.harmony', group.by = 'orig.ident')

DimPlot(integrated_cca_v3,,group.by = 'sample') +
  DimPlot(sce_merge_cca5,group.by = 'sample') +
  DimPlot(sce_merge_harmony,group.by = 'sample') +
  DimPlot(sce_fastmnn_res,group.by = 'sample') +
  DimPlot(sce_merge_rpca,group.by = 'sample') +
  DimPlot(sce_liger,group.by = 'sample')

# 09 导出五物种 H5 文件 ----

colnames(sce_list_filtered[[5]]) <- paste(sce_list_filtered[[5]]$sample,colnames(sce_list_filtered[[5]]),sep = '_')
colnames(sce_list_filtered[[4]]) <- paste(sce_list_filtered[[4]]$sample,colnames(sce_list_filtered[[4]]),sep = '_')
colnames(sce_list_filtered[[3]]) <- paste(sce_list_filtered[[3]]$sample,colnames(sce_list_filtered[[3]]),sep = '_')
colnames(sce_list_filtered[[2]]) <- paste(sce_list_filtered[[2]]$sample,colnames(sce_list_filtered[[2]]),sep = '_')

sce_list_filtered_cp <- sce_list_filtered
sce_list_filtered_cp <- lapply(sce_list_filtered_cp, function(obj) {

  obj[['RNA']] <- as(object = obj[['RNA']], Class = "Assay")

  return(obj)
})

dior::write_h5(sce_list_filtered_cp[[5]],file=file.path(h5_dir, "scRNA_human.h5"),)
dior::write_h5(sce_list_filtered_cp[[4]],file=file.path(h5_dir, "scRNA_monkey.h5"))
dior::write_h5(sce_list_filtered_cp[[3]],file=file.path(h5_dir, "scRNA_goat_trans.h5"))
dior::write_h5(sce_list_filtered_cp[[2]],file=file.path(h5_dir, "scRNA_rat_trans.h5"))
dior::write_h5(sce_list_filtered_cp[[1]],file=file.path(h5_dir, "scRNA_mice_trans.h5"))

# 10 定义并调用 counts H5AD 导出函数 ----

# 从 RNA assay 提取原始 counts，构建用于导出的对象。
write_h5ad_counts_only <- function(seu_obj, file_path) {


  raw_counts <- seu_obj@assays$RNA@counts

  if (is.null(raw_counts) || ncol(raw_counts) == 0) {
    stop("错误：提取到的 counts 矩阵为空，请检查原始对象是否被过度瘦身或破坏！")
  }

  seu_clean <- Seurat::CreateSeuratObject(counts = raw_counts, meta.data = seu_obj@meta.data)

  #seu_clean[["RNA"]] <- as(object = seu_clean[["RNA"]], Class = "Assay")

  GetH5ad(
  seu_clean,
  output_path = file_path,
  mode = "sc",
  assay = "RNA")
}

# 导出各物种的 counts H5AD。

write_h5ad_counts_only(sce_list_filtered_cp[[5]],file=file.path(h5_dir, "scRNA_human.h5ad"))
write_h5ad_counts_only(sce_list_filtered_cp[[4]], file = file.path(h5_dir, "scRNA_monkey.h5ad"))
write_h5ad_counts_only(sce_list_filtered_cp[[3]], file = file.path(h5_dir, "scRNA_goat_trans.h5ad"))
write_h5ad_counts_only(sce_list_filtered_cp[[2]], file = file.path(h5_dir, "scRNA_rat_trans.h5ad"))
write_h5ad_counts_only(sce_list_filtered_cp[[1]], file = file.path(h5_dir, "scRNA_mice_trans.h5ad"))

# 11 导出整合 UMAP 坐标 ----

write.csv(data.frame(integrated_cca_v3@reductions$umap@cell.embeddings,
row.names = c(colnames(sce_list_filtered_cp[[2]]),colnames(sce_list_filtered_cp[[3]]),
colnames(sce_list_filtered_cp[[4]]),colnames(sce_list_filtered_cp[[5]]))), file.path(benchmark_dir, "files/umap_cca3.csv"), row.names = T)

write.csv(data.frame(sce_merge_cca5@reductions$umap@cell.embeddings,
row.names = c(colnames(sce_list_filtered_cp[[5]]),colnames(sce_list_filtered_cp[[4]]),
colnames(sce_list_filtered_cp[[3]]),colnames(sce_list_filtered_cp[[2]]))), file.path(benchmark_dir, "files/umap_cca5.csv"), row.names = T)

write.csv(data.frame(sce_fastmnn_res@reductions$umap@cell.embeddings,
row.names = c(colnames(sce_list_filtered_cp[[5]]),colnames(sce_list_filtered_cp[[4]]),
colnames(sce_list_filtered_cp[[3]]),colnames(sce_list_filtered_cp[[2]]))), file.path(benchmark_dir, "files/umap_fasmnn.csv"), row.names = T)

write.csv(data.frame(sce_liger@reductions$umap@cell.embeddings,
row.names = c(colnames(sce_list_filtered_cp[[5]]),colnames(sce_list_filtered_cp[[4]]),
colnames(sce_list_filtered_cp[[3]]),colnames(sce_list_filtered_cp[[2]]))), file.path(benchmark_dir, "files/umap_liger.csv"), row.names = T)

# 12 绘制固定的运行时间比较 ----

# 1. 准备原始数据
# 创建数据框
data_tm <- data.frame(
  Method = c("CCA_v3", "CCA_v5", "Harmony", "RPCA", "LIGER", "FastMNN"),
  Elapsed_Time_s = c(7717.019, 767.200, 783.695, 439.124, 847.139, 601.151)
)

# 2. 添加 SCVI 和 SCANVI 数据
# 假设您提供的数字是 Elapsed_Time_s
new_data <- data.frame(
  Method = c("SCVI", "SCANVI"),
  Elapsed_Time_s = c(577.68, 832.26)
)

# 合并数据
final_df <- rbind(data_tm, new_data)

# 3. 绘图
ggplot(final_df, aes(x = Elapsed_Time_s,
                     # 关键点：使用 reorder 进行排序
                     # -Elapsed_Time_s 表示取负值排序，这样数值小的（负值大）排在坐标轴上方
                     y = reorder(Method, -Elapsed_Time_s))) +

  # 绘制柱状图，stat="identity"表示使用实际数值，fill设置灰度颜色
  geom_bar(stat = "identity", fill = "grey40", width = 0.7) +

  # 添加数值标签（可选，增加可读性）
  geom_text(aes(label = round(Elapsed_Time_s, 1)),
            hjust = -0.1,  # 标签向右偏移一点
            size = 3.5) +

  # 设置坐标轴和标题
  labs(
    title = "Execution Time Comparison (Elapsed Time)",
    x = "Elapsed Time (seconds)",
    y = "Method"
  ) +

  # 扩大 X 轴范围，防止最长的柱子标签被切掉
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +

  # 使用简洁主题
  theme_classic() +

  # 美化细节
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(color = "black")
  )
