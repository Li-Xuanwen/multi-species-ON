source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

setwd(r_work_dir)
lapply(c('Seurat', 'dplyr', 'patchwork', 'ggplot2', 'scplotter', 'stringr','tidyr','SingleCellExperiment','ggsignif',
'ewceData','EWCE','MAGMA.Celltyping'), require, character.only = TRUE)
#save.image(file.path(rdata_dir, "subtype.RData"))

# 02 读取并整理青光眼相关基因表 ----

#  355 gene  from NG - Prof. Han  Table Sheet5
glu_gene <- readxl::read_xlsx(file.path(disease_dir, "41588_2023_1428_MOESM5_ESM.xlsx"), sheet = 5)[3:357,10] %>%
  # 重命名列名（假设列名为"Gene"）
  setNames("Gene") %>%
  # 拆分分号（支持前后空格、全角分号等）
  separate_rows(Gene, sep = "\\s*[;；]\\s*") %>%
  # 清理基因名：去空格、转大写、过滤空值
  mutate(Gene = str_trim(Gene),  # 去除首尾空格
         Gene = toupper(Gene)) %>%  # 统一大写
  filter(!is.na(Gene) & Gene != "") %>%  # 移除空值
  # 严格去重
  distinct(Gene)

# 03 构建主分析 CTD ----
counts_matrix <- GetAssayData(sce_glia_final, assay = "RNA", slot = "counts")
# 1. 提取细胞类型的注释
cell_annotations <- sce_glia_final@meta.data$cell_subtype
names(cell_annotations) <- colnames(sce_glia_final) # 确保细胞名和矩阵列名对应

# 2. 生成CellTypeData (CTD) 对象
ctd <- generate_celltype_data(
  exp = counts_matrix,
  annotLevels = list(level1 = cell_annotations),
  groupName = "MyDiseaseAnalysis",
  savePath = file.path(r_work_dir, "revise/glu_score") # 保存路径，"."代表当前工作目录
)

# 04 按表达比例筛选疾病基因 ----

genes_to_check <- intersect(glu_genes, rownames(sce_glia_final))
# 3. 检查是否有基因可供检查
if (length(genes_to_check) == 0) {

  warning("您提供的 'glu_genes' 均不在 'sce_glia_final' 对象中。")
  filtered_glu_genes <- c() # 返回一个空向量

} else {

  # 4. (关键) 计算 DotPlot 数据
  #    我们使用 'group.by' 参数来指定按 'cell_subtype' 分组
  #    我们只提取 '$data' 部分，这里面包含了所有计算好的百分比

  message("正在计算各分群的表达百分比...")
  dp_data <- DotPlot(sce_glia_final, features = genes_to_check, group.by = "cell_subtype")$data

  # 5. 查看 dp_data 结构 (可选, 供您理解)
  #    它是一个数据框 (data frame)，包含:
  #    - dp_data$pct.exp      : 表达百分比 (0-100 的值)
  #    - dp_data$features.plot: 基因名称
  #    - dp_data$id           : cell_subtype 的分群名称

  # 6. 筛选出在 *任何* 分群中表达 > 5% 的行
  #    注意: DotPlot 的 'pct.exp' 是 0-100 的值，所以 5% 就是 5
  passed_rows <- dp_data[dp_data$pct.exp > 5, ]

  # 7. 从这些 "通过" 的行中，提取出所有不重复的基因名
  filtered_glu_genes <- unique(passed_rows$features.plot)

  # 8. (可选) 将结果从 factor 转换为 character 向量，更易于使用
  filtered_glu_genes <- as.character(filtered_glu_genes)

}

# 05 背景基因及主 EWCE 分析 ----

# 检查一下你的基因有多少个能匹配到单细胞数据中
# 这是一个很好的习惯，可以避免后续出错
background_genes <- rownames(sce_glia_final)
good_genes <- intersect(filtered_glu_genes , background_genes)
print(sprintf("你的基因集中有 %d 个基因, 其中 %d 个可以在单细胞数据中找到。",
              length(glu_genes), length(good_genes)))
# 184
write.table(good_genes,file = file.path(r_work_dir, "revise/good_genes.txt"))
write.table(background_genes,file = file.path(r_work_dir, "revise/bk_genes.txt"))

print(sprintf("你的基因集中有 %d 个基因, 其中 %d 个可以在单细胞数据中找到。",
              length(glu_genes), length(good_genes_2)))
load(file.path(r_work_dir, "revise/glu_score/ctd_MyDiseaseAnalysis.rda"))
full_results = EWCE::bootstrap_enrichment_test(sct_data=ctd,
                                               hits=good_genes,
                                               bg=background_genes,
                                               genelistSpecies = 'human',
                                               sctSpecies = 'human',
                                               sctSpecies_origin = 'human',
                                         reps=10000,annotLevel=1
                                         )

print(full_results$results[order(full_results$results$p),3:5][1:6,])
print(ewce.plot(full_results$results,mtc_method="BH"))

save_dir = file.path(r_work_dir, "revise/glu_score")
EWCE::generate_bootstrap_plots(sct_data=ctd,hits=good_genes,bg=background_genes,
                               reps=10000,annotLevel=1,
                               full_results=full_results,
                               listFileName="VignetteGraphs",
                               genelistSpecies = 'human',
                               sctSpecies = 'human',
                               sctSpecies_origin = 'human',
                               save_dir )

# 06 主 EWCE 可视化 ----
library(ggplot2)
# full_results 是你的EWCE运行结果
results_df <- full_results$results

# 增加一列用于标记显著性
results_df$significant <- ifelse(results_df$q < 0.05, "Yes", "No")

ggplot(results_df, aes(x = reorder(CellType, -sd_from_mean), y = sd_from_mean, fill = significant)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
  labs(
    title = "EWCE Enrichment Score",
    x = "Cell Type",
    y = "Enrichment Score (sd_from_mean)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 07 AS / OPC 相关基因筛选 ----

#### 筛选差异的 glu - related gene ####
avg_expr_hsa_sub <- AverageExpression(sce_glia_final, assays = "RNA", slot = "data", group.by = "cell_subtype")$RNA
mat_z_opc_glu <- avg_expr_hsa_sub[filtered_glu_genes,c(6:8)] %>% t() %>% scale() %>% t() %>% as.data.frame()
opc2_glu <- which(avg_expr_hsa_sub[mat_z_opc_glu %>% filter(OPC2 >0.7) %>% rownames(),2] > 2)

mat_z_AS_glu <- avg_expr_hsa_sub[filtered_glu_genes,c(1:2)] %>% t() %>% scale() %>% t() %>% as.data.frame()
AS1_glu <- which(avg_expr_hsa_sub[mat_z_AS_glu %>% filter(AS1 > 0.7) %>% rownames(),1] > 4)

DotPlot(subset(sce_glia_final,cell_subtype %in% c('OPC1','OPC2','OPC3')),
                                  features = names(opc2_glu),
                        group.by = "cell_subtype",col.min = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DotPlot(subset(sce_glia_final,cell_subtype %in% c('AS1','AS2')),
                                  features = names(AS1_glu[c(2,3,4,5,7,8,9,18)]),
                        group.by = "cell_subtype",col.min = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# 08 读取外部眼组织参照并合并 counts ----

load('eye_glaucoma.RData')
seurat_TM <- subset(sce_anterior,cell_type == 'TM_Fibro')
seurat_TM$cell_subtype <- seurat_TM$cell_type
sce_retina_AS$cell_subtype <- sce_retina_AS$cell_type

### 矩阵
common_genes <- Reduce(
  intersect,
  list(
    rownames(sce_glia_final),
    rownames(seurat_TM),
    rownames(sce_retina_AS)
  )
)

mat1 <- GetAssayData(sce_glia_final, assay = "RNA", layer = "counts")[common_genes, ]
mat2 <- GetAssayData(seurat_TM, assay = "RNA", layer = "counts")[common_genes, ]
mat3 <- GetAssayData(sce_retina_AS, assay = "RNA", layer = "counts")[common_genes, ]

counts_matrix_tmp <- cbind(mat1, mat2, mat3)

## 注释
cell_annotations_tmp <- c(
  as.character(sce_glia_final$cell_subtype),
  as.character(seurat_TM$cell_subtype),
  as.character(sce_retina_AS$cell_subtype)
)

cell_annotations_tmp <- trimws(cell_annotations_tmp)
names(cell_annotations_tmp) <- colnames(counts_matrix_tmp)

stopifnot(ncol(counts_matrix_tmp) == length(cell_annotations_tmp))
stopifnot(identical(colnames(counts_matrix_tmp), names(cell_annotations_tmp)))
stopifnot(!any(is.na(cell_annotations_tmp)))

# 09 构建外部比较 CTD 及候选基因 ----

ctd_tmp <- generate_celltype_data(
  exp = counts_matrix_tmp,
  annotLevels = list(level1 = cell_annotations_tmp),
  groupName = "MyDiseaseAnalysis",
  savePath = file.path(r_work_dir, "revise/glu_score_tmp") # 保存路径，"."代表当前工作目录
)

#  1. 筛选出 'glu_genes' 中真实存在于对象中的基因
genes_to_check_tmp <- intersect(glu_genes, common_genes)
# 3. 检查是否有基因可供检查
if (length(genes_to_check_tmp) == 0) {

  warning("没有共有基因哦")
  filtered_glu_genes <- c() # 返回一个空向量

} else {

  # 4. (关键) 计算 DotPlot 数据
  #    我们使用 'group.by' 参数来指定按 'cell_subtype' 分组
  #    我们只提取 '$data' 部分，这里面包含了所有计算好的百分比

  message("正在计算各分群的表达百分比...")
  dp_data <- DotPlot(sce_glia_final, features = genes_to_check_tmp, group.by = "cell_subtype")$data

  # 5. 查看 dp_data 结构 (可选, 供您理解)
  #    它是一个数据框 (data frame)，包含:
  #    - dp_data$pct.exp      : 表达百分比 (0-100 的值)
  #    - dp_data$features.plot: 基因名称
  #    - dp_data$id           : cell_subtype 的分群名称

  # 6. 筛选出在 *任何* 分群中表达 > 5% 的行
  #    注意: DotPlot 的 'pct.exp' 是 0-100 的值，所以 5% 就是 5
  passed_rows <- dp_data[dp_data$pct.exp > 5, ]

  # 7. 从这些 "通过" 的行中，提取出所有不重复的基因名
  filtered_glu_genes <- unique(passed_rows$features.plot)

  # 8. (可选) 将结果从 factor 转换为 character 向量，更易于使用
  filtered_glu_genes <- as.character(filtered_glu_genes)

}

# 10 外部比较 EWCE 与结果图 ----

# 检查一下你的基因有多少个能匹配到单细胞数据中
# 这是一个很好的习惯，可以避免后续出错
background_genes <- rownames(sce_glia_final)
good_genes <- intersect(filtered_glu_genes , background_genes)
print(sprintf("你的基因集中有 %d 个基因, 其中 %d 个可以在单细胞数据中找到。",
              length(glu_genes), length(good_genes)))
# 184
write.table(good_genes,file = file.path(r_work_dir, "revise/good_genes.txt"))
write.table(background_genes,file = file.path(r_work_dir, "revise/bk_genes.txt"))

print(sprintf("你的基因集中有 %d 个基因, 其中 %d 个可以在单细胞数据中找到。",
              length(glu_genes), length(good_genes_2)))
load(file.path(r_work_dir, "revise/glu_score_tmp/ctd_MyDiseaseAnalysis.rda"))
full_results_tmp = EWCE::bootstrap_enrichment_test(sct_data=ctd,
                                               hits=good_genes,
                                               bg=common_genes,
                                               genelistSpecies = 'human',
                                               sctSpecies = 'human',
                                               sctSpecies_origin = 'human',
                                         reps=10000,annotLevel=1
                                         )

print(full_results$results[order(full_results$results$p),3:5][1:6,])
print(ewce.plot(full_results$results,mtc_method="BH"))

# full_results 是你的EWCE运行结果
results_df_tmp <- full_results_tmp$results

# 增加一列用于标记显著性
results_df_tmp$significant <- ifelse(results_df_tmp$q < 0.05, "Yes", "No")

ggplot(results_df_tmp, aes(x = reorder(CellType, -sd_from_mean), y = sd_from_mean, fill = significant)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Yes" = "red", "No" = "grey")) +
  labs(
    title = "EWCE Enrichment Score",
    x = "Cell Type",
    y = "Enrichment Score (sd_from_mean)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 11 保存指定结果 ----

fig6_objects <- c(
  "sce_glia_final",
  "sce_retina_AS",
  # gene list
  "glu_gene",
  "glu_genes",
  "genes_to_check",
  "filtered_glu_genes",
  "background_genes",
  "good_genes",

  # main EWCE
  "counts_matrix",
  "cell_annotations",
  "ctd",
  "full_results",
  "results_df",

  # subtype differential glaucoma-related genes
  "avg_expr_hsa_sub",
  "mat_z_opc_glu",
  "opc2_glu",
  "mat_z_AS_glu",
  "AS1_glu",

  # external comparison: retina-AS / TM_Fibro
  "seurat_TM",
  "common_genes",
  "counts_matrix_tmp",
  "cell_annotations_tmp",
  "ctd_tmp",
  "genes_to_check_tmp",
  "full_results_tmp",
  "results_df_tmp"
)

# 只保存当前环境中真实存在的对象，避免因为某些对象不存在而报错
fig6_objects_exist <- fig6_objects[fig6_objects %in% ls()]

save(
  list =fig6_objects_exist,
  file = file.path(rdata_dir, "fig6.RData")
)

# 12 额外对象读取 ----
seurat_5sp <- readRDS(file.path(rdata_dir, "seurat_5sp.rds"))
