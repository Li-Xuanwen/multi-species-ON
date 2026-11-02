source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 载入原有依赖 ----
lapply(c("biomaRt","Seurat","dplyr","ggplot2","scplotter","tidyr","ggpubr","stringr","patchwork",
"SeuratObject"), require, character.only = TRUE)

# 02 原始小鼠分组映射与物种样本构成 ----
mouse_orig_map <- c(
  "MOD" = "OLG",
  "MOS" = "OLG",
  "Ms_ON_sheath_Control1" = "FB",
  "Ms_ON_sheath_Control2" = "FB"
)

sce_list_filtered$Mouse$sample <- setNames(mouse_orig_map[as.character(sce_list_filtered$Mouse$orig.ident)],
 Cells(sce_list_filtered$Mouse))

pdf(file.path(figure_dir, "sup2/stat_sample.pdf"),height = 5,width = 20)
CellStatPlot(sce_list_filtered$Human,ident = 'sample', group_by = "cell_type", frac = "group",
             swap = TRUE, position = "stack") + RotatedAxis() + ggtitle('human') + xlab(NULL) |
  CellStatPlot(sce_list_filtered$Macaque,ident = 'sample', group_by = "cell_type", frac = "group",
               swap = TRUE, position = "stack") + RotatedAxis() + ggtitle('macaque')+ xlab(NULL) |
  CellStatPlot(sce_list_filtered$Rat,ident = 'sample', group_by = "cell_type", frac = "group",
                swap = TRUE, position = "stack") + RotatedAxis() + ggtitle('rat') + xlab(NULL) |
  CellStatPlot(sce_list_filtered$Mouse,ident = 'sample', group_by = "cell_type", frac = "group",
                  swap = TRUE, position = "stack") + RotatedAxis()+ ggtitle('mouse') + xlab(NULL) |
  CellStatPlot(sce_list_filtered$Goat,ident = 'sample', group_by = "cell_type", frac = "group",
               swap = TRUE, position = "stack") + RotatedAxis()+ ggtitle('goat') + xlab(NULL)
dev.off()

# 03 修改现有五物种对象的样本字段 ----

is_mouse <- sce_5sp$sp == "mouse"

# 2. 使用这个逻辑向量来对指定的列进行赋值
sce_5sp$MOS[is_mouse] <- "Mouse1"
sce_5sp$MOD[is_mouse] <- "Mouse2"

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

sce_5sp$sample <- factor(sce_5sp$sample,levels = c('Human1',  'Human2', 'Monkey1', 'Monkey2','Goat1','Goat2','Goat3','Goat4',
 'Rat1','Rat2','Rat3','Rat4','Mouse1','Mouse2','Mouse3','Mouse4'))

FeatureStatPlot(sce_5sp, features = c("nFeature_RNA", "nCount_RNA"),
                ident = "sample", facet_scales = "free_y")

# 04 重新载入另一五物种对象并处理细胞 ID ----

library(SCNT)
sce_5sp <- GetSeurat(file.path(scarches_dir, "5sp/5sp.h5ad"))

new_cell_ids <- sub("-\\d+$", "", colnames(sce_5sp))
sce_5sp <- RenameCells(object = sce_5sp, new.names = new_cell_ids)

# 05 读取外部邻接矩阵和 latent 表 ----

connectivities <- readMM(file.path(scarches_dir, "5sp/connectivities.mtx"))
distances <- readMM(file.path(scarches_dir, "5sp/distances.mtx"))
latent_df <- read.csv(file.path(scarches_dir, "5sp/latent_representation_20.csv"),row.names = 1)

rownames(latent_df) <- sub("-\\d-\\d+$", "",rownames(latent_df))

rownames(connectivities) <- rownames(latent_df)
colnames(connectivities) <- rownames(latent_df)
rownames(distances) <- rownames(latent_df)
colnames(distances) <- rownames(latent_df)

# 06 写入图结构与整理潜在空间列名 ----

sce_5sp@graphs$RNA_snn_ann <- as.Graph(connectivities[colnames(sce_5sp),])
sce_5sp@graphs$RNA_dist_ann <- as.Graph(distances[colnames(sce_5sp),])

##### latent 对齐 #####
# 转换为矩阵
latent_matrix <- as.matrix(latent_df)

# 假设 latent_matrix 当前列名为 "latent_0", "latent_1", ...
colnames(latent_matrix) <- as.character(1:ncol(latent_matrix))  # 改为 "1", "2", ..., "10"

# 07 写入外部 UMAP 并读取标签表 ----

umap_loci <- read.csv(file.path(scarches_dir, "5sp/5sp_umap_scarches.csv"))[,2:3]
rownames(umap_loci) <- rownames(latent_df)
# 将新的 UMAP 坐标添加到 Seurat 对象的 reductions 中
sce_5sp[["umap_latent"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_loci[colnames(sce_5sp),]),
  key = "umap_",
  assay = DefaultAssay(sce_5sp)
)

peri_celltype <- read.csv(file.path(scarches_dir, "5sp/peri_cell_type.csv"))
rownames(peri_celltype) <- rownames(latent_df)

# 08 计算标签交叉表及按行匹配率 ----

cross_table <- table(sce_5sp$cell_type, sce_5sp$cell_type_peri)
print(cross_table)

# 计算每行（即每种 cell_type）的匹配百分比
# 对角线上的数值表示匹配的细胞数
matching_by_type <- diag(cross_table)

# 计算每行（即每种 cell_type）的总数
total_by_type <- rowSums(cross_table)

# 计算每种 cell_type 的匹配百分比
matching_percentage_by_type <- (matching_by_type / total_by_type) * 100

# 将结果转换为数据框以便更好地展示
matching_df <- data.frame(
  cell_type = names(matching_percentage_by_type),
  matching_percentage = matching_percentage_by_type
)

print(matching_df)
