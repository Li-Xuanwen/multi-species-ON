# 01 环境清理与工作目录 ----
rm(list =ls())
source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

setwd(project_dir)

# 02 载入原有依赖 ----

lapply(c("Seurat","dplyr","ggplot2","tidyr","ggpubr","stringr","patchwork","SeuratObject"), require, character.only = TRUE)

# 数据预处理 ####
library(Seurat)
library(dplyr)

# 03 原始数据目录与物种定义 ----

# 使用一个列表来清晰地管理每个物种的名称和其对应的数据源路径。
# 这样做便于未来添加或修改物种，而无需更改核心处理逻辑。
species_data_sources <- list(
  Human  = list(path = file.path(shared_matrix_dir, "Huma")),
  Monkey = list(path = file.path(shared_matrix_dir, "Monkey")),
  Rat    = list(path = file.path(matrix_dir, "rat")),
  Goat   = list(path = file.path(matrix_dir, "goat"))
)

# 04 逐物种、逐样本创建 Seurat 对象 ----

# 初始化一个最终的空列表，用于存储所有物种的数据
scRNAlist <- list()

# 遍历在上面定义的每个物种
for (species_name in names(species_data_sources)) {

  # 获取当前物种的数据路径
  species_path <- species_data_sources[[species_name]]$path

  cat("========================================\n")
  cat("开始处理物种:", species_name, "\n")
  cat("数据路径:", species_path, "\n")

  # 检查物种的基础目录是否存在
  if (!dir.exists(species_path)) {
    warning(paste("指定的路径不存在，跳过物种:", species_path))
    next # 'next'会跳过本次循环，继续处理下一个物种
  }

  # 使用 list.files 递归查找所有 'filtered_feature_bc_matrix' 文件夹
  # 这是最关键的一步，它能自动处理不同深度的目录结构 (如Human vs Monkey)
  matrix_paths <- list.files(
    path = species_path,
    pattern = "filtered_feature_bc_matrix",
    recursive = TRUE,    # 递归查找，深入所有子目录
    full.names = TRUE,   # 返回完整路径
    include.dirs = TRUE  # 确保返回的是目录
  )

  # 如果没有找到任何样本，则打印警告并跳到下一个物种
  if (length(matrix_paths) == 0) {
    warning(paste("在", species_name, "的路径下未找到任何 'filtered_feature_bc_matrix' 文件夹。"))
    next
  }

  # 初始化一个临时的、属于当前物种的列表
  species_seurat_list <- list()

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
    species_seurat_list[[sample_name]] <- seurat_obj
  }

  # --- 按要求重命名小列表中的样本 ---
  # 将样本名重命名为 "物种名 + 数字" (例如: Human1, Human2)
  if (length(species_seurat_list) > 0) {
    names(species_seurat_list) <- paste0(species_name, seq_along(species_seurat_list))
  }

  # 将处理完的当前物种的列表，作为一个整体，添加到最终的大列表中
  scRNAlist[[species_name]] <- species_seurat_list

  cat("物种", species_name, "处理完成，共找到", length(species_seurat_list), "个样本。\n\n")
}

# 05 结束提示与临时变量清理 ----

cat("所有指定的物种均已处理完毕。\n")
# 清理不再需要的临时变量
rm(species_data_sources, species_name, species_path, matrix_paths, data_path, sample_name, counts, seurat_obj, species_seurat_list)
