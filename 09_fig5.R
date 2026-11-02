source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 环境与依赖 ----
setwd(file.path(revised_figure_dir, "fig5"))
#load(file.path(rdata_dir, "fig4.RData"))
lapply(c('Seurat', 'dplyr', 'patchwork', 'ggplot2', 'scplotter','tidyr','SCNT'), require, character.only = TRUE)
#save.image(file.path(rdata_dir, "subtype.RData"))

# 02 亚型 marker 检查 ----
Idents(sce_glia_final) <- 'cell_subtype'
sce_glia_final <- NormalizeData(sce_glia_final)

markers_ASsub <- FindAllMarkers(
      subset(sce_glia_final,cell_subtype %in% c('AS1','AS2','AS3')),
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.5,
      logfc.threshold = 1
    )

markers_ASsub$diff <- markers_ASsub$pct.1 -markers_ASsub$pct.2
DotPlot(subset(sce_glia_final,cell_subtype %in% c('AS1','AS2','AS3')),
features = filter(markers_ASsub,cluster == 'AS3',diff >0.2)$gene, scale =F,col.min =0)

markers_OPCsub <- FindAllMarkers(
      subset(sce_glia_final,cell_subtype %in% c('OPC1','OPC2','OPC3')),
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.5,
      logfc.threshold = 1
    )

markers_OPC3<- FindMarkers(
      subset(sce_glia_final,cell_subtype %in% c('OPC1','OPC2','OPC3')),
      assay = "RNA",
      ident.1 = 'OPC3',
      only.pos = TRUE,
      min.pct = 0.5,
      logfc.threshold = 1
    )
markers_OPCsub$diff <- markers_OPCsub$pct.1 -markers_OPCsub$pct.2
markers_OLGsub <- FindAllMarkers(
      subset(sce_glia_final,cell_subtype %in% c('OLG1','OLG2')),
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.5,
      logfc.threshold = 1
    )
markers_OLGsub$diff <- markers_OLGsub$pct.1 -markers_OLGsub$pct.2

# 03 亚型合并与标签交换 ----
sce_glia_final$cell_subtype <- factor(
  x,
  levels = c(
    "AS1", "AS2",
    "MG",
    "OLG1", "OLG2",
    "OPC1", "OPC2", "OPC3"
  )
)

markers_ASsub_new <- FindAllMarkers(
      subset(sce_glia_final,cell_subtype %in% c('AS1','AS2')),
      assay = "RNA",
      only.pos = TRUE,
      min.pct = 0.5,
      logfc.threshold = 1
    )

# 04 亚型 marker 图与表达检查 ----
pdf(file = file.path(revised_figure_dir, "fig5/gene_dotplot_sub.pdf"),width = 10,height = 3.7)
markerPlot(subset(sce_glia_final_all ,cell_subtype %in% c('AS1','AS2','OLG1','OLG2','OPC1','OPC2','OPC3')),
  features = list(AS = c("AQP4","GFAP",'HIF3A','PAX8','NLGN1','SLIT2','TNC','CAPS2'),
               OLG= c('MOBP','MOG','SLC5A11','SAMD12'),
              OPC = c("PDGFRA","TNR",'ADAM33','CRB1','TCF7L2','SIRT2')))
dev.off()

FeatureStatPlot(subset(sce_glia_final_all,cell_subtype %in% c('OPC1','OPC2','OPC3')), features = c(
       "PDGFRA", "OLIG2", "CRB1",
       "SOX6"
   ), ident = "cell_subtype", add_bg = TRUE, stack = TRUE,
   legend.position = "top", legend.direction = "horizontal")

# 05 补齐指定样本的 species ----

human_origins <- c(
  "Hu220OSON", "Hu221OSON", "Hu222OSON",
  "Hu235ODON", "HuPt14ODON", "Pt6ON", "Pt14OSON"
)

# 仅将上述样本中 species 为 NA 的细胞标记为 Human
sce_glia_final$species <- as.character(sce_glia_final$species)

idx <- is.na(sce_glia_final$species) &
  sce_glia_final$orig.ident %in% human_origins

table(sce_glia_final$orig.ident[idx])  # 修改前核对

sce_glia_final$species[idx] <- "Human"

# 06 物种构成柱图 ----

plot_data <- sce_glia_final@meta.data %>%
  # --- 第一步：计算背景（分母） ---
  # 计算每个【主群】(cell_type) 中，各【物种】(sp) 的总细胞数
  # 这一步是为了消除不同物种在不同主群采样深度的偏差
  group_by(predictions_scanvi, species) %>%
  mutate(total_sp_in_major = n()) %>%
  ungroup() %>%

  # --- 第二步：计算观测值（分子） ---
  # 统计每个【亚群】(cell_subtype) 中各物种的数量
  # 保留 cell_type 和 total_sp_in_major 用于后续计算
  group_by(predictions_scanvi, cell_subtype, species, total_sp_in_major) %>%
  summarise(count = n(), .groups = 'drop') %>%

  # --- 第三步：初次矫正 (Observed / Expected) ---
  # 这里的逻辑是：如果该物种在这个主群里很多，我们要除以它的总数进行“惩罚”
  mutate(norm_weight = count / total_sp_in_major) %>%

  # --- 第四步：二次归一化 (Re-scale to 100%) ---
  # 为了画堆叠图，我们需要让每个亚群内的比例加起来等于 1
  group_by(cell_subtype) %>%
  mutate(final_prop = norm_weight / sum(norm_weight)) %>%
  ungroup()

# ============================================================
# 3. 绘图 (All in One Stacked Barplot)
# ============================================================

# 为了让图更好看，可以按照主群对亚群进行排序（可选）
# 这一步让属于同一个 cell_type 的亚群在 X 轴上挨在一起
plot_data$cell_subtype <- factor(plot_data$cell_subtype,
                                 levels = unique(plot_data$cell_subtype[order(plot_data$predictions_scanvi)]))

p <- ggplot(plot_data, aes(x = cell_subtype, y = final_prop, fill = species)) +
  geom_col(position = "stack", width = 0.8) +

  # 坐标轴格式：百分比
  scale_y_continuous(labels = scales::percent_format(), expand = c(0,0)) +

  # 如果你想在图上方标注这个亚群属于哪个主群，可以用 facet_grid
  # scales = "free_x" 和 space = "free" 很关键，能让宽度自适应
  facet_grid(~ predictions_scanvi, scales = "free_x", space = "free") +

  labs(
    title = "Corrected Species Composition across All Subtypes",
    subtitle = "Normalized by species abundance within each Major Group",
    y = "Balanced Proportion",
    x = "Cell Subtypes",
    fill = "Species"
  ) +

  theme_classic() +
  theme(
    # X轴文字旋转45度，防止重叠
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    # 去掉分面背景，让图看起来更干净
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    # 调整图例位置
    legend.position = "right"
  )

print(p)

# 07 样本排除与 Ro/e 分布 ----

library(Startrac)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(sscVis)
library(circlize)

sce_glia_5sp <- subset(x =sce_glia_final,subset = orig.ident %in% c("8W-2","8W-3","MOD","ROD","ROS","goat-normal-ON1",
                       "goat-normal-ON2","goat_on_control2","goat_on_control3","human1", "human2" ))
keep_cells <- colnames(sce_glia_final)[
  !(sce_glia_final$sample %in% c("Human1", "Human2"))
]

sce_glia_final_del <- subset(
  sce_glia_final,
  cells = keep_cells
)

sce_glia_final_del$sample %>% table()
R_oe <- calTissueDist(sce_glia_final_del@meta.data,
                      byPatient = F,
                      colname.cluster = "cell_subtype",
                      #colname.patient = "patient",
                      colname.tissue = "species",
                      method = "chisq",
                      min.rowSum = 0)
#R_oe <- R_oe[c(1,2,3,7),]

# 定义颜色函数（这里将数值从 0 映射到蓝色，从 2 映射到红色）
col_fun <- colorRamp2(c(0, 2), c("#FFFDE7", "#E65100"))

ht <- Heatmap(
  mat,
  name = "R_o/e",
  show_heatmap_legend = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "right",
  show_column_names = TRUE,
  show_row_names = TRUE,
  col = col_fun,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "R_o/e",
    at = seq(0, 2, by = 0.5),
    labels = seq(0, 2, by = 0.5),
    legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", mat[i, j]), x, y,
              gp = gpar(fontsize = 8))
  }
)

pdf(file.path(revised_figure_dir, "fig5/R_oe_heatmap.pdf"), width = 10, height = 8)
draw(ht, newpage = TRUE)
dev.off()

# 08 OPC 外部数据读入及 counts 导出 ----

seurat_NSD_OPC <- readRDS(file.path(reference_dir, "GSE198951_COPs/GSE198951_NSD1_OPC_Seurat.RDS"))
write_h5ad_counts_only(seurat_NSD_OPC, file.path(subtype_export_dir, "human_COPs.h5ad"))
write_h5ad_counts_only(sce_glia_final, file.path(subtype_export_dir, "human_OPC_count_final.h5ad"))

y_opc <- UpdateSeuratObject(seurat_NSD_OPC)
y_opc$cell_subtype <- y_opc$newannot
if (!inherits(y_opc[["RNA"]], "Assay5")) {
  y_opc[["RNA"]] <- as(y_opc[["RNA"]], Class = "Assay5")
}
DefaultAssay(y_opc) <- 'RNA'
y_opc[["SCT"]] <- NULL
y_opc <- subset(
  y_opc,
  features = rownames(sce_human_glia_all_gene_new)
)

# 9 构建 OPC counts 对象 ----

y_counts <- LayerData(
    y_opc,
    assay = "RNA",
    layer = "counts"
)

y_opc_clean <- CreateSeuratObject(
    counts = y_counts,
    meta.data = y_opc@meta.data,
    project = "opc"
)

x_opc <- subset(sce_human_glia_all_gene_new,cell_subtype %in% c('OPC1','OPC2','OPC3'))
counts <- LayerData(
    x_opc,
    assay = "RNA",
    layer = "counts"
)

x_opc_clean <- CreateSeuratObject(
    counts = counts,
    meta.data = x_opc_clean@meta.data
)
x_opc_clean <- subset(
  x_opc_clean,
  features = rownames(y_opc)
)

# 10 OPC Harmony / CCA 整合 ----

opc_merged <- merge(
  x = x_opc_clean,
  y = y_opc_clean,
  add.cell.ids = c("glia", "NSD"),
  merge.data = FALSE,
  project = "opc_integration"
) %>% JoinLayers()

opc_merged[["RNA"]] <- split(
  opc_merged[["RNA"]],
  f = opc_merged$orig.ident
)

opc_preprocessed <- NormalizeData(
  opc_merged,
  verbose = FALSE
)

opc_preprocessed <- FindVariableFeatures(
  opc_preprocessed,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

opc_preprocessed <- ScaleData(
  opc_preprocessed,
  verbose = FALSE
)

opc_preprocessed <- RunPCA(
  opc_preprocessed,
  npcs = 30,
  verbose = FALSE
)

opc_preprocessed <- IntegrateLayers(
    object = opc_preprocessed,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.harmony",
    verbose = FALSE
  )

opc_preprocessed <- RunUMAP(
    opc_preprocessed,
    reduction = "integrated.harmony",
    dims = 1:20,
    reduction.name = "umap.harmony"
)
opc_preprocessed <- RunTSNE(opc_preprocessed, reduction = "integrated.harmony", dims = 1:20, reduction.name = "tsne.harmony")

opc_preprocessed <- IntegrateLayers(object = opc_preprocessed, method = CCAIntegration,
                                    orig.reduction = "pca", new.reduction = "integrated.cca",
                                        k.weight = 50,
                                    verbose = FALSE)
opc_preprocessed <- FindNeighbors(opc_preprocessed, reduction = "integrated.cca", dims = 1:20)
opc_preprocessed <- RunUMAP(opc_preprocessed, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")
opc_preprocessed <- RunTSNE(opc_preprocessed, reduction = "integrated.cca", dims = 1:20, reduction.name = "tsne.cca")

CellDimPlot(opc_preprocessed,
  group_by = "cell_subtype", reduction = "umap.harmony",
  highlight = 'cell_subtype == "COP"'
)

CellDimPlot(opc_preprocessed,
  group_by = "cell_subtype", reduction = "umap.cca",
  highlight = 'cell_subtype == "OPC3"'
) +
CellDimPlot(opc_preprocessed,
  group_by = "cell_subtype", reduction = "umap.cca",
  highlight = 'cell_subtype == "COP"'
)

CellDimPlot(opc_preprocessed,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "OPC3"'
) +
CellDimPlot(opc_preprocessed,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "COP"'
)

# 11 离群过滤及候选类群展示 ----

# 查看结果
opc_hdb <- dbscan_outlier(
  sobj = opc_preprocessed,
  reduction = "umap.cca",
  ident = "cell_subtype",
  minPts = 5,
  eps = 0.5
)

####  可视化候选细胞
plot_order <- c(
  "OPC1",
  "OPC",
  "OPC2",
  "OPC3",
  "COP"

)

opc_hdb <- opc_hdb[, order(match(opc_hdb$cell_subtype, plot_order))]
CellDimPlot(opc_hdb,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "OPC3"',
  show_stat = FALSE,      # 不显示每类细胞数量 / legend 中的 n
  add_mark = FALSE       # 不画分群外边框 / hull / ellipse
) +
CellDimPlot(opc_hdb,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "COP"',
  show_stat = FALSE,      # 不显示每类细胞数量 / legend 中的 n
  add_mark = FALSE       # 不画分群外边框 / hull / ellipse
)

CellDimPlot(opc_hdb,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "OPC1"',
  show_stat = FALSE,      # 不显示每类细胞数量 / legend 中的 n
  add_mark = FALSE       # 不画分群外边框 / hull / ellipse
)
CellDimPlot(opc_hdb,
  group_by = "cell_subtype", reduction = "tsne.cca",
  highlight = 'cell_subtype == "OPC2"',
  show_stat = FALSE,      # 不显示每类细胞数量 / legend 中的 n
  add_mark = FALSE       # 不画分群外边框 / hull / ellipse
)

# 12 AS1 嵌入与特征图 ----
seurat_AS <- subset(sce_glia_final,predictions_scanvi == 'AS')
seurat_AS[["RNA"]] <- split(
  seurat_AS[["RNA"]],
  f = seurat_AS$orig.ident
)
seurat_AS <- RunUMAP(seurat_AS, reduction = "latent", dims = 1:10)

seurat_AS <- NormalizeData(seurat_AS) %>% 
    FindVariableFeatures(selection.method = "vst",nfeatures = 1000) %>% 
    ScaleData() %>%
    RunPCA(npcs = 50, verbose = T)

seurat_AS_harmony <- IntegrateLayers(seurat_AS,method = HarmonyIntegration,
                                 orig.reduction = "pca",new.reduction = "harmony",
                                 theta = 5,
                                 max.iter.harmony = 20,
                                 tau  =T,
                                 verbose =F)
seurat_AS_harmony <- RunUMAP(seurat_AS_harmony, reduction = "harmony", dims = 1:20)
seurat_AS_harmony <- RunTSNE(seurat_AS_harmony, reduction = "harmony", dims = 1:20)
DimPlot(seurat_AS_harmony,group.by = "cell_subtype")
AS1_tsne <- DimPlot(seurat_AS_harmony,group.by = "cell_subtype",reduction = 'tsne')
ggsave(
  filename =  file.path(revised_figure_dir, "fig5/AS1_tsne.pdf"),
  plot = AS1_tsne,
  width = 6,
  height = 5
)

AS_feature <- FeaturePlot(seurat_AS_harmony,reduction = 'tsne',features = c('CAPS2','TNC'))
ggsave(
  filename =  file.path(revised_figure_dir, "fig5/AS1_feature.pdf"),
  plot = AS_feature,
  width = 10,
  height = 5
)

seurat_AS_cca <- IntegrateLayers(
    object = seurat_AS,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )
seurat_AS_cca <- RunUMAP(seurat_AS_cca, reduction = "integrated.cca", dims = 1:20)
seurat_AS_cca <- RunTSNE(seurat_AS_cca, reduction = "integrated.cca", dims = 1:20)
DimPlot(seurat_AS_cca,group.by = "cell_subtype")
DimPlot(seurat_AS_cca,group.by = "cell_subtype",reduction = 'tsne')

# 13 OPC 表达检查与保存现有对象 ----

seurat_OPC[["RNA"]] <- split(
  seurat_OPC[["RNA"]],
  f = seurat_OPC$orig.ident
)
OPC_feature <- FeaturePlot(seurat_OPC,features = c('CRB1','GPR17'))

save_AS1_tsne_objects <- function(
    file = file.path(rdata_dir, "AS1_tsne.RData")
){

  # 需要保存的对象
  obj_list <- c(

    # 原始对象
    "seurat_glia_final",

    # AS subtype
    "seurat_AS",
    "seurat_AS_harmony",
    "seurat_AS_cca",
    "AS1_tsne",
    "AS_feature",

    # OPC subtype
    "seurat_OPC",
    "seurat_OPC_harmony"
  )

  # 只保存当前环境存在的对象
  obj_exist <- obj_list[obj_list %in% ls(envir = .GlobalEnv)]

  save(
    list = obj_exist,
    file = file,
    compress = "xz",
    envir = .GlobalEnv
  )

  message(
    "Saved objects:\n",
    paste(obj_exist, collapse = "\n")
  )

}

# 执行保存
save_AS1_tsne_objects()