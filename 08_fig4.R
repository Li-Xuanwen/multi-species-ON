source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

setwd(file.path(revised_figure_dir, "fig4"))
lapply(c('Seurat', 'dplyr', 'patchwork', 'ggplot2', 'tibble','scplotter', 'ROGUE','SCNT','ComplexUpset','tidyr','MetaNeighbor','SingleCellExperiment'), require, character.only = TRUE)
sce_glia <- subset(sce_5sp,predictions_scanvi %in% c('AS','OLG','OPC','MG'))

# 02 胶质细胞 UMAP 查看 ----

p_5umap <- CellDimPlot(sce_glia, group_by = 'predictions_scanvi',facet_by = 'species', facet_ncol = 5,
palcolor = my_colors[sce_pnas$orig_cell_type   %>% table() %>% names()] )+
  labs(x = "", y = '', subtitle = NULL) + theme(axis.ticks = element_blank(),
                                                  axis.text = element_blank(),
                                                  panel.grid=element_blank(),
                                                  panel.border = element_blank())

p_umap <- CellDimPlot(sce_glia, group_by = 'predictions_scanvi',
palcolor = my_colors[sce_pnas$orig_cell_type   %>% table() %>% names()])+
  labs(x = "", y = '', subtitle = NULL)

# 03 各物种 marker 与集合交集 ----

## 差异基因 top 热图 + pathway ####
markers_list_glia_celltype <- lapply(sce_list_filtered, function(df) {
  Idents(df) <- 'cell_type'
  df_glia <- subset(df,cell_type %in% c('AS','OLG','OPC','MG'))
  FindAllMarkers(df_glia, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
})

markers_list_glia_celltype_filter <- lapply(markers_list_glia_celltype, function(df) {
  df %>%
    mutate(Difference = pct.1 - pct.2) %>%
    dplyr::filter(avg_log2FC > 2.5 & Difference > 0.5)
})

# 合并数据（确保markers_list_celltype_filter中包含对应的数据框）
merged_df <- bind_rows(
  markers_list_glia_celltype$Mouse %>% mutate(species = "Mouse"),
  markers_list_glia_celltype$Rat %>% mutate(species = "Rat"),
  markers_list_glia_celltype$Goat %>% mutate(species = "Goat"),
  markers_list_glia_celltype$Macaque %>% mutate(species = "Macaque"),
  markers_list_glia_celltype$Human %>% mutate(species = "Human")
)
#merged_df$gene <- toupper(merged_df$gene)
#write.csv(merged_df,file = "independent_markers_with_pvalue_logfc.csv")
# 去重，并准备presence列
merged_df1 <- merged_df %>%
  mutate(gene = as.character(gene)) %>%
  distinct(species, cluster, gene) %>%
  mutate(presence = 1)

# 宽数据变换：以species为列，cluster+gene为行，presence为值
upset_data <- merged_df1 %>%
  pivot_wider(
    names_from = species,
    values_from = presence,
    values_fill = list(presence = 0),
    values_fn = list(presence = sum)
  )

# 假设 upset_data 结构如下：
# upset_data: data.frame with columns: cluster, gene, chx, mcc, rno, hsa
# chx, mcc, rno, hsa 列为0/1表示该 (cluster, gene) 是否在对应物种有表达

# 1. 对cluster进行汇总，检查该cluster是否在所有物种中都有至少一个基因出现
cluster_presence <- upset_data %>%
  group_by(cluster) %>%
  summarise(
    chx_present = any(Goat == 1),
    mcc_present = any(Macaque == 1),
    rno_present = any(Rat == 1),
    hsa_present = any(Human == 1),
    mmu_present = any(Mouse == 1)
  )

# 2. 筛选出在所有物种中都有基因出现的cluster
clusters_in_all <- cluster_presence %>%
  filter(chx_present & mcc_present & rno_present & hsa_present & mmu_present) %>%
  pull(cluster)

# 若 clusters_in_all 为空，说明无满足条件的cluster
if (length(clusters_in_all) == 0) {
  stop("No clusters are present in all species. Check your data or relax the criteria.")
}

# 3. 过滤 upset_data，只保留上述符合条件的cluster
upset_data_filtered <- upset_data %>%
  filter(cluster %in% clusters_in_all)

# 4. 准备绘图列
set_cols <- colnames(upset_data_filtered)[!(colnames(upset_data_filtered) %in% c("cluster", "gene"))]

# 5. 绘制UpSet图
library(ComplexUpset)
upset(
  upset_data_filtered %>%
      filter(cluster %in% c('OLG','AS','MG','OPC')),
  intersect = set_cols,
  base_annotations = list(
    'Gene number' = intersection_size(
      counts = FALSE,
      mapping = aes(fill = cluster)
    ) +
      scale_fill_manual(
        name="cell type",
        values=my_colors
      )
  ),
  width_ratio = 0.1,
  set_sizes = upset_set_size()
) +
theme(
  axis.text.x = element_text(angle=90),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# 04 第二组 UpSet 筛选 ----

upset_data_after_filtering <- upset_data_filtered %>%
  # 使用 `!` 来选择不满足括号内条件的行
  filter(
    !(
      # 条件1：排除所有单个物种独有的基因 (行的和为 1)
      (rowSums(across(all_of(set_cols))) == 1) |

        # 条件2：排除 mcc-hsa 共有且仅属于这两个物种的基因
        # (当数据是 0/1 时，需要用 == 1 和 == 0 进行显式比较)
        (mcc == 1 & hsa == 1 & chx == 0 & rno == 0 & mmu == 0) |

        # 条件3：排除 mmu-rno 共有且仅属于这两个物种的基因
        (mmu == 1 & rno == 1 & chx == 0 & hsa == 0 & mcc == 0)
    )
  )
upset(
  upset_data_after_filtering %>% filter(cluster %in% c('OLG','AS','MG','OPC')),
  intersect = set_cols,
  base_annotations = list(
    'Gene numer' = intersection_size(
      counts = FALSE,
      mapping = aes(fill = cluster)
    ) +
      scale_fill_manual(
        name = "cell type",  # 设置图例标题
        values = c("FB"     = "#E8BFA2",
                   "EC"     = "#6CC4D2",
                   "Mural"  = "#EF7C1B",
                   "AS" = "#D32528",
                   "OLG" = "#3778AE",
                   "OPC" = "#4DA649",
                   "MG" = "#8E4C96")
      )
  ),
  width_ratio = 0.1,
  set_sizes = (
    upset_set_size() +
      theme(axis.text.x = element_text(angle = 90),
            panel.background = element_blank(),    # 去除背景
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  )
)+
  labs(title = "UpSet Plot of Clusters Present in All Five Species")

# 05 保守 marker、Venn 与表格导出 ----
library(VennDiagram)

# 筛选在所有物种中均表达的基因
filtered_genes <- upset_data_filtered %>%
  filter(Goat == 1, Macaque == 1,Rat == 1, Human == 1, Mouse == 1)

# 按 cluster 分组，生成类似 mouse_markers 的列表
con_markers <- split(filtered_genes$gene, filtered_genes$cluster)

# 定义你想要的顺序
desired_order <- c("OPC", "OLG", "EC", "MG", "Mural", "AS", "FB")

# 使用 match() 函数对 con_markers 进行排序
# 只保留 desired_order 中存在的元素
con_markers_sorted <- con_markers[desired_order]

# 调整边距 绘制基因的 dotplot  ####
sce_5sp$cell_type <- factor(sce_5sp$cell_type,
                                                levels = c("OLG", "OPC", "AS", "FB", "EC", "Mural","MG"))
sce_5sp <- sce_5sp[, order(sce_5sp$cell_type)]
par(mar = c(5, 4, 4, 2) + 0.1)  # 设置边距
cairo_pdf("co-celltypes_markers.pdf", family = "Arial",width = 12, height = 5)
markerPlot(sce_5sp,'cell_type')
dev.off()

# 绘制 Venn 图，基于四个物种中的所有基因集合
gene_sets <- list(
  Goat = upset_data_filtered$gene[upset_data_filtered$Goat == 1] %>% unique(),
  Macaque = upset_data_filtered$gene[upset_data_filtered$Macaque == 1] %>% unique(),
  Rat = upset_data_filtered$gene[upset_data_filtered$Rat == 1] %>% unique(),
  Human = upset_data_filtered$gene[upset_data_filtered$Human == 1] %>% unique(),
  Mouse = upset_data_filtered$gene[upset_data_filtered$Mouse == 1] %>% unique()
)

write.csv(upset_data_filtered,file = "marker_of_5_species_independent.csv")
# 绘制 Venn 图
venn.plot <- venn.diagram(
  x = gene_sets,            # 基因集合
  filename = NULL,          # 不保存到文件，直接绘制
  col = "transparent",
  fill = c("#F08F89", "#C9B31B", "#C29EC9", "#89B4E0","#CE90BE"), # 各物种颜色
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#F08F89", "#C9B31B", "#C29EC9", "#89B4E0","#CE90BE"),
  margin = 0.2
)

# 显示 Venn 图
dev.new()

grid.draw(venn.plot)

# 06 MetaNeighborUS 物种和细胞类型 AUROC ----

#### MettaNeighbor ####
# 高变基因
highly_variable_genes = rownames(sce_pnas_merge)
sce_5sp_sce <- as.SingleCellExperiment(subset(sce_5sp,cell_type %in% c('AS','MG','OLG','OPC')))

# 使用 MetaNeighbor 计算细胞类型之间的相关性
Aurocs_matrix <- MetaNeighborUS(var_genes = highly_variable_genes,
                                    dat = sce_5sp_sce,
                                    study_id = sce_5sp_sce$species,
                                    cell_type = sce_5sp_sce$predictions_scanvi,
                                    fast_version = TRUE)

# 热图展示 MetaNeighbor 结果
par(mar = c(5, 4, 4, 2) + 0.1)  # 设置边距
cairo_pdf("MetaNeighbor_sp.pdf", family = "Arial",width = 7, height = 7)
Heatmap(Aurocs_matrix,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRampPalette(c("#3288BD", "white", "#D53E4F"))(50),
        rect_gp = gpar(col = "black"),
       )  #cell_fun = function(j, i, x, y, width, height, fill) {
         #grid.text(format(Aurocs_matrix[i, j], digits = 2), x, y,
        #           gp = gpar(col = "black", cex = 0.8))
       # }
dev.off()

# 重新排列矩阵
row_names <- rownames(Aurocs_matrix)
split_names <- strsplit(row_names, "\\|")
species <- sapply(split_names, `[`, 1)
cell_type <- sapply(split_names, `[`, 2)

# 定义排序规则
species_order <- c("Human", "Macaque", "Goat", "Rat", "Mouse")
cell_type_order <- unique(cell_type)  # 保持原始细胞类型顺序

# 转换为因子
species_factor <- factor(species, levels = species_order)
cell_type_factor <- factor(cell_type, levels = cell_type_order)

# 生成排序索引
new_order <- order(cell_type_factor, species_factor)
sorted_matrix <- Aurocs_matrix[new_order,
   new_order]

# 热图展示 MetaNeighbor 结果
par(mar = c(5, 4, 4, 2) + 0.1)  # 设置边距
cairo_pdf("MetaNeighbor_cell_type.pdf", family = "Arial",width = 7, height = 7)
Heatmap(sorted_matrix,
        cluster_rows = F,
        cluster_columns = F,
        col = colorRampPalette(c("#3288BD", "white", "#D53E4F"))(50),
        rect_gp = gpar(col = "black"),
        #cell_fun = function(j, i, x, y, width, height, fill) {
        #  grid.text(format(sorted_matrix[i, j], digits = 2), x, y,
        #            gp = gpar(col = "black", cex = 0.8))
        #}
        )
dev.off()

# 07 基于 GO 基因集的 MetaNeighbor ----

# 根据整理的人的GO通路查看 某个细胞的保守程度
# 构建合适的数据格式
sce_5sp <- NormalizeData(sce_5sp)
Meta_df <- SingleCellExperiment(
  #Matrix(sce_all_final_ann@assays$RNA@layers$data),
  #as.SingleCellExperiment(sce_all_final_ann),
  assays = list(logcounts = GetAssayData(sce_5sp, assay = "RNA", slot = "data")),
  colData = data.frame(study_id = sce_5sp$species,
                       cell_type = sce_5sp$predictions_scanvi)
)
celltype_matrix = model.matrix(~as.SingleCellExperiment(sce_5sp)$predictions_scanvi - 1)
AUROC_scores  = MetaNeighbor(dat = Meta_df,
                             experiment_labels  =as.numeric(factor(Meta_df$study_id)),
                             celltype_labels = celltype_matrix,
                             genesets = GOhuman,
                             bplot = TRUE,
                             fast_version = TRUE)

# 转换为数据框并添加行名（细胞类型）
# 1. 去除列名中的 'Meta_df$cell_type' 部分
colnames(AUROC_scores) <- gsub("as\\.SingleCellExperiment\\(sce_5sp\\)\\$predictions_scanvi", "", colnames(AUROC_scores))
# 2. 转换成长格式
AUROC_df <-  AUROC_scores %>% as.data.frame() %>%
  rownames_to_column(var = "GO_term") %>%
  pivot_longer(cols = -GO_term, names_to = "Cell_Type", values_to = "AUROC")

# 3. 按中位数排序 cell_type    ! 逆序
AUROC_df$Cell_Type <- factor(AUROC_df$Cell_Type, levels = c('OLG','MG','EC','AS','MP','FB','Mural','OPC'))

# 4 . 设置颜色对应的细胞类型
cell_type_colors <- my_colors

# 绘制排序后的箱线图
levels(AUROC_df$Cell_Type) <- c('OLG','MG','AS','OPC','EC','FB','Mural','MP')
ggplot(AUROC_df %>% filter(Cell_Type %in% c('OLG','MG','AS','OPC')), aes(x = Cell_Type, y = AUROC,fill = Cell_Type)) +
  geom_boxplot( color = "black",outlier.shape = NA) +
  labs(
    title = "AUROC Scores by Cell Type (Sorted by Median)",
    x = "Cell Type",
    y = "AUROC"
  ) +
  scale_y_continuous(limits = c(0.5, 1.05), breaks = seq(0.5, 1.05, by = 0.1)) + # 设置Y轴范围和单位
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8),
        axis.text.y = element_text(angle = 45, hjust = 1,size = 8),
        axis.title = element_text(,size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = cell_type_colors) + # 设置颜色
  labs(x = "Cell Type", y = "AUROC", title = "AUROC Scores by Cell Type")

# 08 ROGUE 计算及AUROC 再绘图片段 ----

# 没用
expr <-  GetAssayData(sce_5sp, layer = "counts") %>% as.matrix()

# 获取分类信息
meta <- sce_5sp@meta.data[,c('species','predictions_scanvi','orig.ident')]

# 表达熵模型
ent.res <- SE_fun(expr)
SEplot(ent.res)

# marker 检查与固定基因表达矩阵 ----

sce_glia <- sce_glia %>% NormalizeData()
markers_glia <- FindAllMarkers(sce_glia,group.by = 'cell_type',min.pct = 0.5, only.pos =T)

top5_markers <- markers_glia %>%
  group_by(cluster) %>%              # 1. 根据 'cluster' 列分组
  #arrange(desc(avg_log2FC)) %>%    # 2. 在每个组内，按 'avg_log2FC' 降序排列
  slice_head(n = 5)

features <- list(
  AS = c( "DCLK1","SORBS1","NEBL","PLCE1","GFAP"),
  MG = c( "INPP5D","DOCK8","APBB1IP","P2RY12", "PTPRC"),
  OPC = c("TNR","PCDH15","DSCAM","VCAN",'CSPG4'),
  OLG = c("SLC24A2","MOBP","MOG","ENPP2","ST18")
)

# 定义Y轴从下到上的顺序
desired_order <- c("AS", "MG", "OPC", "OLG")
sce_glia$cell_type <- factor(x = sce_glia$cell_type, levels = desired_order)

#cairo_pdf(file.path(project_dir, "figure3/5markerplot.pdf"), family = "Arial",width = 12, height = 5)
# markerPlot(sce_clean_glia,features, 'cell_type')
#dev.off()

features_in <- c("PBX1","SORBS1","NEBL","PLCE1","GFAP",'IKZF1',"INPP5D","DOCK8","P2RY12", "PTPRC","SOX6","TNR","PCDH15","VCAN",'CSPG4',
  "ST18","SLC24A2","MOBP","MOG","ENPP2")

species_list <- c('Goat','Human','Macaque','Rat','Mouse')
classes <- c('OLG','OPC','AS','MG')

# 1. 提前生成所有可能的列名
all_colnames <- unlist(lapply(species_list, function(species) {
  paste0(species, "_", classes)
}))

# 2. 使用所有基因和所有可能的列名来创建空的、预分配好内存的矩阵
#    行是所有基因，列是所有物种-细胞类型的组合
#    用 NA (Not Available) 来填充
full_log <- matrix(NA,
                   nrow = length(features_in),
                   ncol = length(all_colnames),
                   dimnames = list(features_in, all_colnames))

full_classnorm <- matrix(NA,
                         nrow = length(features_in),
                         ncol = length(all_colnames),
                         dimnames = list(features_in, all_colnames))

for(i in 1:length(species_list)){
  print(i)
  species <- species_list[i]
  obj <- subset(sce_glia, subset = species == species_list[i])
  Idents(obj) <- "cell_type"
  obj <- subset(obj, idents = classes)
  obj@meta.data$cell_class0 = factor(obj@meta.data$cell_type, levels = classes)
  levels(obj@meta.data$cell_class0) =c( 'OLG','OPC','AS','MG')
  Idents(obj) <- "cell_class0"
  # Calculate average expression for all genes by class
  #features_in <- features[features %in% rownames(obj)]
  avg_exp <- AverageExpression(obj, features = features_in, assays = "RNA", slot = "count")$RNA
  colnames(avg_exp) <- paste0(species, "_", colnames(avg_exp))
  # Normalize
  log_exp <- log1p(avg_exp)
  class_norm <- avg_exp / apply(avg_exp, 1, max)
  # save data
  full_log[features_in, colnames(log_exp)] <-  as.matrix(log_exp[features_in,])
  full_classnorm[features_in, colnames(class_norm)] <-  as.matrix(class_norm[features_in,])
}

# 10 整理热图展示次序并导出 ----
order <- c("Human_AS","Macaque_AS","Rat_AS","Mouse_AS","Goat_AS",
           "Human_MG","Macaque_MG","Rat_MG","Mouse_MG","Goat_MG",
           "Human_OPC","Macaque_OPC","Rat_OPC","Mouse_OPC","Goat_OPC",
           "Human_OLG","Macaque_OLG","Rat_OLG","Mouse_OLG","Goat_OLG")
library(reshape2)
file_table <- full_classnorm[,order]
melt_table <- melt(file_table, id.vars = "gene")

par(mar = c(5, 4, 4, 2) + 0.1)  # 设置边距
cairo_pdf(file.path(revised_figure_dir, "fig4/5marker_new.pdf"), family = "Arial",width = 7, height = 5)
ggplot(data = melt_table, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(colour="white", size=0.2) +
  guides(fill=guide_legend(title="Norm. Exp.")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank()) +
  scale_fill_gradient(low = "white",  high = "#447241", na.value = "gray70")
dev.off()

# 11 保存指定分析对象 ----
save(
  # 输入依赖
  sce_5sp,
  sce_pnas_merge,
  my_colors,
  GOhuman,
 # cell_type_colors_new,

  # glia UMAP
  sce_glia,
  p_5umap,
  p_umap,

  # conserved marker
  markers_list_glia_celltype,
  markers_list_glia_celltype_filter,
  merged_df,
  merged_df1,
  upset_data,
  cluster_presence,
  clusters_in_all,
  upset_data_filtered,
  set_cols,
  #upset_data_after_filtering,
  filtered_genes,
  con_markers,
  desired_order,
  #con_markers_sorted,
  gene_sets,
  venn.plot,

  # MetaNeighbor
  highly_variable_genes,
  sce_5sp_sce,
  Aurocs_matrix,
  row_names,
  split_names,
  species,
  cell_type,
  species_order,
  cell_type_order,
  species_factor,
  cell_type_factor,
  new_order,
  sorted_matrix,

  # GO AUROC
  Meta_df,
  celltype_matrix,
  AUROC_scores,
  AUROC_df,
  cell_type_colors,

  # marker heatmap
  markers_glia,
  top5_markers,
  features,
  features_in,
  species_list,
  classes,
  all_colnames,
  full_log,
  full_classnorm,
  order,
  file_table,
  melt_table,

  file = file.path(rdata_dir, "fig4.RData")
)
