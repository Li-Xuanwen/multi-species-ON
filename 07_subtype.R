source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

lapply(c('Seurat', 'dplyr', 'patchwork', 'ggplot2', 'scplotter','tidyr'), require, character.only = TRUE)
save.image('subtype_tmp.RData')
sce_all <- GetSeurat(file.path(scarches_dir, "5sp_pnas/5sp_pnas.h5ad"))

# 定义所有细胞类型
cell_types <- unique(sce_all$predictions_scanvi)

# 02 基于潜在表示逐类群聚类 ----

# 对每个细胞类型进行亚群分析
sce_subtype_list_20 <- list()
markers_subtype_list_20 <- list()
sce_all <- NormalizeData(sce_all)
for (type in cell_types) {
  # 提取子集（自动保留 scanvi_latent）
  sce_sub <- subset(sce_all, predictions_scanvi == type)

  # 可选：补充局部降维（若潜变量维度不足）
  #sce_sub <- FindVariableFeatures(sce_sub, selection.method = "vst", nfeatures = 2000)
  sce_sub <- ScaleData(sce_sub)
  #sce_sub <- RunPCA(sce_sub)

  # 构建邻域图（优先使用 scanvi_latent）
  sce_sub <- FindNeighbors(
    sce_sub,
    reduction = "latent",  # 通过scanvi
    dims = 1:20                  # 使用所有潜变量维度
    #graph.name = paste0('new', "_snn")
  )

  # 聚类分析
  sce_sub <- FindClusters(
    sce_sub,
    #graph.name = paste0('new', "_snn"),
    resolution = 1
  )

  # 在当前细胞类群内寻找亚群标记基因。
  markers <- FindAllMarkers(
    sce_sub,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.5,
    logfc.threshold = 1
  )
  sce_sub <- RunUMAP(sce_sub, reduction = "latent", dims = 1:20,reduction.key = 'UMAPsub_')

  # 存储结果
  sce_subtype_list_20[[type]] <- sce_sub
  markers_subtype_list_20[[type]] <- markers
}

markers_subtype_list_20_filtered <- lapply(markers_subtype_list_20,FUN = function(x) {
  x %>%
    mutate(Difference = pct.1 - pct.2) %>%
    filter(pct.1 >0.6 & Difference > 0.4)
})

# 03 OPC 亚群注释 ----

# OPC
markers_subtype_list_20_filtered[['OPC']]$cluster %>% table()

sce_subtype_list_20[['OPC']]$RNA_snn_res.1 %>% table()
#   0    1    2    3    4    5    6    7    8    9   10   11   12
# 1098  935  864  849  846  718  655  495  442  145   97   60   29
DotPlot(sce_subtype_list_20[['OPC']],features =c('TNR','PDGFRA',unique(markers_subtype_list_20_filtered[['OPC']]$gene)),col.min = 0)
DotPlot(sce_subtype_list_20[['OPC']],features = unique(c('TNR','PDGFRA','PCDH15')),
  group.by = 'RNA_snn_res.1',col.min = 0,scale = F)
  DimPlot(sce_subtype_list_20[['OPC']],label = T)

DotPlot(sce_subtype_list_20[['OPC']],features = c('BNC2','CEMIP'))

CellStatPlot(sce_subtype_list_20[['OPC']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')
markers_subtype_list_20_filtered[['OPC']]

# 删除12   4+8  7  其他
OPC_subtype <- c(
  "0" = "OPC1",
  "1" = "OPC1",
  "2" = "OPC1", # WEAK - SLC47A1
  "3" = "OPC1",
  "4" = "OPC2", # WEAK
  "5" = "OPC1", # 低表细胞 删除
  "6" = "OPC1",
  "7" = "OPC3", #
  "8" = "OPC2",
  "9" = "OPC1",
  "10" = "OPC3",
  "11" = "OPC1"
)

# 删除低质量群  并重命名
sce_subtype_list_20[['OPC']] <- subset(sce_subtype_list_20[['OPC']], subset = RNA_snn_res.1 == '12', invert = TRUE)
sce_subtype_list_20[['OPC']]$RNA_snn_res.1 <- droplevels(sce_subtype_list_20[['OPC']]$RNA_snn_res.1)
sce_subtype_list_20[['OPC']]$cell_subtype <- setNames(OPC_subtype[as.character(sce_subtype_list_20[['OPC']]$RNA_snn_res.1)], Cells(sce_subtype_list_20[['OPC']]))

DotPlot(sce_subtype_list_20[['OPC']],features = c('TNR','PDGFRA','BNC2','CEMIP','UTRN','ITGA6'),group.by = 'cell_subtype',
  col.min=0,scale = F)

CellStatPlot(sce_subtype_list_20[['OPC']], ident= 'cell_subtype',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')+
  DimPlot(sce_subtype_list_20[['OPC']], group.by = 'cell_subtype')

sce_tree_OPC <- FindClusters(sce_subtype_list_20[['OPC']], resolution = seq(0.1, 1.0, by = 0.1))
clustree_plots_OPC <- clustree(sce_tree_OPC)

# 04 AS 亚群注释 ----

markers_subtype_list_20_filtered[['AS']]$cluster %>% table()

DimPlot(sce_subtype_list_20[['AS']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
  DotPlot(sce_subtype_list_20[['AS']],features = c('TNC','CAPS2','GFAP','AQP4'),col.min = 0,scale =F)
CellStatPlot(sce_subtype_list_20[['AS']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')
sce_subtype_list_20[['AS']]$RNA_snn_res.1 %>% table()

DotPlot(sce_subtype_list_20[['AS']],features = c('TNC','CAPS2','GFAP','AQP4'),col.min = 0,scale =F)

sce_tree_AS <- FindClusters(sce_subtype_list_20[['AS']], resolution = seq(0.1, 1.0, by = 0.1))
clustree_plots <- clustree(sce_tree_AS)

# 定义群
AS_subtype <- c(
  "0" = "AS1",
  "1" = "AS2",
  "2" = "AS1", # WEAK - SLC47A1
  "3" = "AS1",
  "4" = "AS1", # WEAK
  "5" = "AS1", # 低表细胞 删除
  "6" = "AS2",
  "7" = "AS2", #
  "8" = "AS1",
  "9" = "AS2",
  "10" = "AS1",
  "11" = "AS1"
)

# 重命名
sce_subtype_list_20[['AS']]$cell_subtype <- setNames(AS_subtype[as.character(sce_subtype_list_20[['AS']]$RNA_snn_res.1)], Cells(sce_subtype_list_20[['AS']]))
CellStatPlot(sce_subtype_list_20[['AS']], ident= 'cell_subtype',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')+
  DimPlot(sce_subtype_list_20[['AS']], group.by = 'cell_subtype')

DotPlot(sce_subtype_list_20[['AS']],features = c('GFAP','AQP4','TNC','CAPS2'),group.by = 'cell_subtype',
  col.min=0,scale = F)

# 05 EC 与 Mural 表达检查 ----

sce_subtype_list_20[['EC']]$RNA_snn_res.1 %>% table()
#   0    1    2    3    4    5    6    7    8    9   10   11   12
# 1232 1079  913  902  891  723  689  682  631  347  301  124  100

markers_subtype_list_20_filtered[['EC']]
DimPlot(sce_subtype_list_20[['EC']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
  DotPlot(sce_subtype_list_20[['EC']],features = c("FLT1", "ERG"),group.by = 'RNA_snn_res.1',col.min = 0,scale = F)

DotPlot(sce_subtype_list_20[['EC']],features = unique(markers_subtype_list_20_filtered[['EC']]$gene),group.by = 'RNA_snn_res.1',col.min = 0)
CellStatPlot(sce_subtype_list_20[['EC']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')

# Mural - 没有亚群
markers_subtype_list_20_filtered[['Mural']]$cluster %>% table()
# 0  1  2  3  4  5  6  7  8
# 0 10  0  0  2  1  3 12  0
DimPlot(sce_subtype_list_20[['Mural']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
  DotPlot(sce_subtype_list_20[['Mural']],features = c("ACTA2","MYH11"),group.by = 'RNA_snn_res.1',col.min = 0,scale = F)

DotPlot(sce_subtype_list_20[['Mural']],features = unique(markers_subtype_list_20_filtered[['Mural']]$gene),group.by = 'RNA_snn_res.1',col.min = 0)
llStatPlot(sce_subtype_list_20[['Mural']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')

# 06 MG 亚群注释 ----

markers_subtype_list_20_filtered[['MG']]$cluster %>% table()
# 0  1  2  3  4  5  6  7  8  9 10 11
# 0  1 10  0 44  0  0  4  3 22  0  5

DimPlot(sce_subtype_list_20[['MG']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
    DotPlot(sce_subtype_list_20[['MG']],features =c('P2RY12','C3','CD163','MRC1','F13A1'),group.by = 'RNA_snn_res.1',col.min = 0,scale = F)

DotPlot(sce_subtype_list_20[['MG']],features =unique(markers_subtype_list_20_filtered[['MG']]$gene),col.min = 0)

# 定义群
#
MG_subtype <- c(
  "0" = "MG1",
  "1" = "MG1",
  "2" = "MG1", # WEAK - SLC47A1
  "3" = "MG2",
  "4" = "MG2", # WEAK
  "5" = "MG1", # 低表细胞 删除
  "6" = "MG1",
  "7" = "MG1", #
  "8" = "MG1",
  "9" = "MG3",
  "10" = "MG1",
  "11" = "MG1"
)

# 重命名
sce_subtype_list_20[['MG']]$cell_subtype <- setNames(MG_subtype[as.character(sce_subtype_list_20[['MG']]$RNA_snn_res.1)], Cells(sce_subtype_list_20[['MG']]))
# 可视化
CellStatPlot(sce_subtype_list_20[['MG']], ident= 'cell_subtype',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')+
  DimPlot(sce_subtype_list_20[['MG']], group.by = 'cell_subtype')

# 07 OLG 亚群注释 ----

markers_subtype_list_20_filtered[['OLG']]

sce_subtype_list_20[['OLG']]$RNA_snn_res.1 %>% table()
#    0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16
#10694 10139  7242  7041  6045  5885  5700  5016  4418  3504  3016  2947  1965  1241   153   133    67

DimPlot(sce_subtype_list_20[['OLG']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
    DotPlot(sce_subtype_list_20[['OLG']],features =c('MOBP','MOG'),group.by = 'RNA_snn_res.1',col.min = 0,scale = F)

CellStatPlot(sce_subtype_list_20[['OLG']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')

DotPlot(sce_subtype_list_20[['OLG']],features =unique(markers_subtype_list_20_filtered[['OLG']]$gene),col.min = 0,scale = F)

# 定义群
OLG_subtype <- c(
  "0" = "OLG1",
  "1" = "OLG1",
  "2" = "OLG1",
  "3" = "OLG2",
  "4" = "OLG3",
  "5" = "OLG1",
  "6" = "OLG1",
  "7" = "OLG1",
  "8" = "OLG1",
  "9" = "OLG3",
  "10" = "OLG2",
  "11" = "OLG3",
  "12" = "OLG1",
  "13" = "OLG1",
  "14" = "OLG1",
  "15" = "OLG1"
)

# 删除低质量群  并重命名
sce_subtype_list_20[['OLG']] <- subset(sce_subtype_list_20[['OLG']], subset = RNA_snn_res.1 == '16', invert = TRUE)
sce_subtype_list_20[['OLG']]$RNA_snn_res.1 <- droplevels(sce_subtype_list_20[['OLG']]$RNA_snn_res.1)
sce_subtype_list_20[['OLG']]$cell_subtype <- setNames(OLG_subtype[as.character(sce_subtype_list_20[['OLG']]$RNA_snn_res.1)], Cells(sce_subtype_list_20[['OLG']]))

# 可视化
CellStatPlot(sce_subtype_list_20[['OLG']], ident= 'cell_subtype',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')+
  DimPlot(sce_subtype_list_20[['OLG']], group.by = 'cell_subtype')

# 08 FB 亚群检查 ----

markers_subtype_list_20_filtered[['FB']]$cluster %>% table()

sce_subtype_list_20[['FB']]$RNA_snn_res.1 %>% table()

DimPlot(sce_subtype_list_20[['FB']],reduction = 'umap',group.by = 'RNA_snn_res.1',label = T) +
    DotPlot(sce_subtype_list_20[['FB']],features =c('DCN','COL1A2','COL1A1'),group.by = 'RNA_snn_res.1',col.min = 0,scale = F)

CellStatPlot(sce_subtype_list_20[['FB']], ident= 'RNA_snn_res.1',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')

# 09 合并胶质细胞亚群标签 ----

sce_glia <- subset(sce_all, cells = c(Cells(sce_subtype_list_20[['OLG']]),
                                        Cells(sce_subtype_list_20[['OPC']]),
                                      Cells(sce_subtype_list_20[['MG']]),
                                    Cells(sce_subtype_list_20[['AS']])))

subtype_order <- c(sce_subtype_list_20[['OLG']]$cell_subtype,
                                        sce_subtype_list_20[['OPC']]$cell_subtype,
                                      sce_subtype_list_20[['MG']]$cell_subtype,
                                    sce_subtype_list_20[['AS']]$cell_subtype)
sce_glia$cell_subtype <- subtype_order[colnames(sce_glia)]

CellDimPlot(sce_glia,
              group_by = c("predictions_scanvi","cell_subtype"), reduction = "umap",label = T)
CellStatPlot(sce_glia, ident= 'cell_subtype',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident')

# 10 物种分布与 Ro/e ----

library(Startrac)
library(ggplot2)
#library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
library(circlize)
library(grid)

R_oe <- calTissueDist(sce_glia@meta.data,
                      byPatient = F,
                      colname.cluster = "cell_subtype",
                      #colname.patient = "patient",
                      colname.tissue = "sp",
                      method = "chisq",
                      min.rowSum = 0)
R_oe <- R_oe[c(1,2,3,7),]

rownames(R_oe) <- c("OLG", "OPC", "AS", "Immune")
colnames(R_oe) <- c("Goat", "Human", "Macaque", "Rat")

# 定义颜色函数（这里将数值从 0 映射到蓝色，从 2 映射到红色）
col_fun <- colorRamp2(c(0, 2), c("#FFFDE7", "#E65100"))

# 绘制热图
Heatmap(as.matrix(R_oe),
        name = "R_o/e",                       # 图例标题
        show_heatmap_legend = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = 'right',
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
          grid.text(sprintf("%.2f", R_oe[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)

# 11 外部注释与亚群 marker 展示 ----

env <- new.env()
load(file.path(pnas_workspace_dir, "pnas_work.RData"),envir =  env)

# 按照手动亚群绘制比例
CellStatPlot(sce_subtype_list_20[['OPC']], ident= 'cell_type_sub_manul',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident') /
  CellDimPlot(sce_subtype_list_20[['OPC']],
              group_by = "cell_type_sub_manul", reduction = "umap",label = T)

CellStatPlot(sce_subtype_list_20[['FB']], ident= 'cell_type_sub_manul',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident') /
  CellDimPlot(sce_subtype_list_20[['FB']],
              group_by = "cell_type_sub_manul", reduction = "umap",label = T)

CellStatPlot(sce_subtype_list_20[['OLG']], ident= 'cell_type_sub_manul',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident') /
  CellDimPlot(sce_subtype_list_20[['OLG']],
              group_by = "cell_type_sub_manul", reduction = "umap",label = T)

CellStatPlot(sce_subtype_list_20[['MG']], ident= 'cell_type_sub_manul',plot_type = "bar",
             group_by = "sp", position = "stack",frac = 'ident') /
  CellDimPlot(sce_subtype_list_20[['MG']],
              group_by = "cell_type_sub_manul", reduction = "umap",label = T)

# 手动亚群的 marker gene 展示
FeatureStatPlot(sce_subtype_list_20[['OPC']], features = split(makrer_sub_list[['OLG']]$gene,makrer_sub_list[['OLG']]$cluster), ident = "cell_type_sub_manul", plot_type = "dot")

tmp <- makrer_sub_list[['OPC']] %>% group_by(cluster) %>% top_n(15)
markerPlot(sce_subtype_list_20[['OPC']],tmp$gene,'cell_type_sub_manul') + coord_flip()

tmp <- makrer_sub_list[['FB']] %>% group_by(cluster) %>% top_n(5)
markerPlot(sce_subtype_list_20[['FB']],tmp$gene,'cell_type_sub_manul') + coord_flip()

tmp <- makrer_sub_list[['OLG']] %>% group_by(cluster) %>% top_n(5)
markerPlot(sce_subtype_list_20[['OLG']],tmp$gene,'cell_type_sub_manul') + coord_flip()

tmp <- makrer_sub_list[['MG']] %>% group_by(cluster) %>% top_n(5)
markerPlot(sce_subtype_list_20[['MG']],tmp[-20,]$gene,'cell_type_sub_manul') + coord_flip()

# 12 GO 与 KEGG 结果展示 ----

data(enrich_example, package = "plotthis")

df_sorted <- go_sub_list[['MG']] %>%
  group_by(group) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()

df_sorted <- df_sorted %>%
  mutate(Description = factor(Description, levels = unique(Description)))
ggplot() +
  geom_col(data = df_sorted %>% group_by(group) %>% top_n(3), aes(x = -log10(pvalue), y = Description, fill = group), width = 0.75) +
  theme_classic() +
  #scale_fill_manual(name = "Regulation", values = mycol, labels = c("Down", "Up")) +
  #scale_color_manual(name = "Regulation", values = mycol, labels = c("Down", "Up")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(#axis.text.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_line(size = 1.2),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 16)) +
  labs(
    title = 'GO Enriched pathways',
    x = 'log10(pvalue)'
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  geom_vline(xintercept = 0, color = "black", size = 1.2) +
  theme(legend.position = "right",
        axis.title.y = element_blank())

df_sorted <- kegg_sub_list[['MG']] %>%
  group_by(group) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()

df_sorted <- df_sorted %>%
  mutate(Description = factor(Description, levels = unique(Description)))
ggplot() +
  geom_col(data = df_sorted %>% group_by(group) %>% top_n(3), aes(x = -log10(pvalue), y = Description, fill = group), width = 0.75) +
  theme_classic() +
  #scale_fill_manual(name = "Regulation", values = mycol, labels = c("Down", "Up")) +
  #scale_color_manual(name = "Regulation", values = mycol, labels = c("Down", "Up")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(#axis.text.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line.x = element_line(size = 1.2),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16)) +
  labs(
    title = 'KEGG Enriched pathways',
    x = 'log10(pvalue)'
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) +
  geom_vline(xintercept = 0, color = "black", size = 1.2) +
  theme(legend.position = "right",
        axis.title.y = element_blank())

# 13 五物种对象的另一亚群分析分支 ----
sce_subtype_list2 <- list()
markers_subtype_list2 <- list()

for (type in cell_types) {
  # 提取子集（自动保留 scanvi_latent）
  sce_sub <- subset(sce_5sp, cell_type == type)

  # 可选：补充局部降维（若潜变量维度不足）
  #sce_sub <- RunPCA(sce_sub, npcs = 20)

  # 构建邻域图（优先使用 scanvi_latent）
  sce_sub <- FindNeighbors(
    sce_sub,
    reduction = "latent",  # 通过scanvi
    dims = 1:20,                  # 使用所有潜变量维度
    graph.name = paste0(type, "_snn")
  )

  # 聚类分析
  sce_sub <- FindClusters(
    sce_sub,
    graph.name = paste0(type, "_snn"),
    resolution = 1
  )

  # 在当前细胞类群内寻找亚群标记基因。

  sce_sub <- ScaleData(sce_sub)
  markers <- FindAllMarkers(
    sce_sub,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.5,
    logfc.threshold = 1
  )
  markers %>%
    mutate(Difference = pct.1 - pct.2) %>%
    filter(pct.1 >0.6 & Difference > 0.4)
  sce_sub <- RunUMAP(sce_sub, reduction = "latent", dims = 1:20,reduction.key = 'UMAPsub_')

  # 存储结果
  sce_subtype_list2[[type]] <- sce_sub
  markers_subtype_list2[[type]] <- markers
}
