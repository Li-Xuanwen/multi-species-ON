source(file.path(Sys.getenv("MULTI_ON_CODE_DIR", unset = getwd()), "config", "paths.R"))

# 01 载入原有依赖与工作目录 ----
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scplotter)
library(dbscan)
setwd(r_work_dir)

# 02 人类 marker 与人工注释 ----

human_markers <- list(
  OPC = c("PCDH15","TNR"), # noncanonical
  Mural = c("ACTA2","MYH11"), # ,'CALD1'
  MP = c('CD163','MRC1','F13A1'),
  MG = c('P2RY12','C3'),
  #VLMC = c('SLC47A1', 'SEMA3G','SLIT3'),
  EC = c("FLT1", "ERG"),
  Tcell = c('CD8A','CD3D'), # ,'CD8' 'CD3','CD8B',
  #Adipocyte = c("PNPLA2", "PLIN1"),
  FB = c('DCN','COL1A2'),
  #Epithelial_Cells = c('CDH1','PKP2'),
  AS = c("AQP4","GFAP"),  #
  OLG = c("MBP", "MOBP"))

#### 2. 确认分群 ####

##### human #####
DotPlot(object =scesinglet_list_filter[[1]],
        features = human_markers,
        cols = c("lightgrey", "blue"),
        #assay = "RNA",
        dot.scale = 4,  
        cluster.idents = FALSE,  
        scale.by = "radius",
        scale = TRUE, 
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of human marker")

human_cell_type_map <- c(
  "0" = "OLG",
  "1" = "AS",
  "2" = "FB",
  "3" = "OLG",
  "4" = "FB",
  "5" = "None",
  "6" = "MG",
  "7" = "FB", #
  "8" = "EC",
  "9" = "Epithelial_Cells",
  "10" = "MP",
  "11" = "Mural",
  "12" = "FB",
  "13" = "OPC",
  "14" = "None",
  "15" = "Adipocyte"
)

# 查看数据 确认5 是低质量
VlnPlot(scesinglet_list_filter[[1]], features =  c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi"), pt.size = 0, group.by = "RNA_snn_res.0.5")

# 确认群
scesinglet_list_filter[[1]]$cell_type <- setNames(human_cell_type_map[as.character(scesinglet_list_filter[[1]]$RNA_snn_res.0.5)], Cells(scesinglet_list_filter[[1]]))
# 确认结果
CellDimPlot(scesinglet_list_filter[[1]], reduction = 'umap', group_by = 'cell_type')
DotPlot(object =scesinglet_list_filter[[1]],
        features = human_markers ,
        cols = c("lightgrey", "blue"),
        group.by = 'cell_type',
        #assay = "RNA",
        dot.scale = 4,
        cluster.idents = FALSE,
        scale.by = "radius",
        scale = TRUE,
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of human marker")

# 03 猴 marker、人工注释与 MG 重标 ----

monkey_markers <- list(
  Mural = c("ACTA2","MYH11"), # ,'CALD1'
  #Adipocyte = c("PNPLA2", "PLIN1"),
  #Glial_Cells = c("RELN","S100B","SORCS1"),
  EC = c("FLT1", "ERG"),
  OPC = c("PCDH15","TNR",'PDFGRA'), #
  #VLMC = c('SLC47A1', 'SEMA3G','SLIT3'),
  # Epithelial_Cells = c('CDH1','PKP2'),
  #Pericyte = c('PDGFRB'),
  FB = c('DCN','COL1A2','COL1A1'),
  MG = c('TMEM119','C3','P2RY12'),
  #Bcell = c('CD19','MS4A1','CD22'),
  Tcell = c('CD8A','CD3D','TRAC'), # ,'CD8' 'CD3','CD8B',
  MP = c('CD163','MRC1','F13A1'),
  AS = c("AQP4","GFAP"),  #
  OLG = c("MBP", "MOBP"))

DotPlot(object =scesinglet_list_filter[[2]],
        features = monkey_markers ,
        cols = c("lightgrey", "blue"),
        #assay = "RNA",
        dot.scale = 4,
        cluster.idents = FALSE,
        scale.by = "radius",  
        scale = TRUE, 
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of macaque marker")

macaque_cell_type_map <- c(
  "0" = "OLG",
  "1" = "OLG",
  "2" = "MG",
  "3" = "OLG",
  "4" = "AS",
  "5" = "OPC",
  "6" = "OLG",
  "7" = "FB",
  "8" = "EC",
  "9" = "FB", # SLC47A1
  "10" = "FB",
  "11" = "MP",
  "12" = "Mural",
  "13" = "Epithelial_Cells",
  "14" = "Tcell",
  "15" = "OPC"
)

# 查看数据 确认是否存在低质量的群
VlnPlot(scesinglet_list_filter[[2]], features =  c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi"), pt.size = 0, group.by = "RNA_snn_res.0.5")

scesinglet_list_filter[[2]]$cell_type <- setNames(macaque_cell_type_map[as.character(scesinglet_list_filter[[2]]$RNA_snn_res.0.5)], Cells(scesinglet_list_filter[[2]]))
# 将11群中PYRY12 表达大于1的定义为MG的定义为MP
p2ry12_expression <- GetAssayData(scesinglet_list_filter[[2]], assay = "RNA", layer = "data")["P2RY12", ]
cells_to_relabel_as_MG <- (scesinglet_list_filter[[2]]$RNA_snn_res.0.5 == 11) & (p2ry12_expression > 2.5)
scesinglet_list_filter[[2]]$cell_type[which(cells_to_relabel_as_MG == T) %>% names()] <- "MG"

FeaturePlot(scesinglet_list_filter[[2]],features = c('CD163','MRC1'),reduction = 'umap')

# 结果确认
CellDimPlot(scesinglet_list_filter[[2]], reduction = 'umap', group_by = 'cell_type',
highlight = 'cell_type == "MP"')
DotPlot(object =scesinglet_list_filter[[2]],
        features = monkey_markers ,
        cols = c("lightgrey", "blue"),
        group.by = 'cell_type',
        #assay = "RNA",
        dot.scale = 4,   # 控制点的大小比例
        cluster.idents = FALSE,   # 是否将细胞群按照树状结构聚类
        scale.by = "radius",          # 控制点大小的标准化方式
        scale = TRUE,                 # 是否对数据进行标准化
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of macaque marker")

# 04 山羊 marker 与人工注释 ----

goat_markers <- list(Mural = c("ACTA2","MYH11"), # ,'CALD1'
  #Adipocyte = c("PNPLA2", "PLIN1"),
  #Glial_Cells = c("RELN","S100B","SORCS1"),
  EC = c("FLT1", "ERG"),
  #Bcell = c('BANK1', 'PRKCB','VAV3 '),#'CD19','MS4A1','CD22','CD20','CD79A','CD79B'),
  OPC = c("TNR","PDGFRA"), #
  #VLMC = c('SLC47A1', 'SEMA3G','SLIT3'),
  #Epithelial_Cells = c('CDH1','PKP2'),
  #Pericyte = c('PDGFRB'),
  FB = c('DCN','COL1A2','COL1A1'),
  MG = c('TMEM119','P2RY12','C3'),
  #Bcell = c('CD19','MS4A1','CD22'),
  Tcell = c('CD8A','CD3D','TRAC'), # ,'CD8' 'CD3','CD8B',
  #MP = c('CD14','MRC1','F13A1'),
  AS = c("AQP4","GFAP"),  #
  OLG = c("MBP", "MOBP"))

DotPlot(object =scesinglet_list_filter[[4]],
        features = goat_markers ,
        cols = c("lightgrey", "blue"),
        #assay = "RNA",
        dot.scale = 4,
        cluster.idents = FALSE, 
        scale.by = "radius",  
        scale = TRUE, 
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of goat marker")

goat_cell_type_map <- c(
  "0" = "OLG",
  "1" = "AS",
  "2" = "OLG",
  "3" = "OLG",
  "4" = "AS",
  "5" = "MG",
  "6" = "FB",
  "7" = "OPC",
  "8" = "FB",# SLC47A1
  "9" = "EC",
  "10" = "OLG",
  "11" = "None", 
  "12" = "OLG",
  "13" = "FB",
  "14" = "Mural",
  "15" = "OPC",
  "16" = "AS"
)

# 查看数据 确认是否存在低质量的群
VlnPlot(scesinglet_list_filter[[2]], features =  c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi"), pt.size = 0, group.by = "RNA_snn_res.0.5")

scesinglet_list_filter[[4]]$cell_type <- setNames(goat_cell_type_map[as.character(scesinglet_list_filter[[4]]$RNA_snn_res.0.5)], Cells(scesinglet_list_filter[[4]]))

# 确认MP 和 MG
FeaturePlot(scesinglet_list_filter[[4]],features = c('P2RY12','MRC1'),reduction = 'umap')

# 05 大鼠 marker 与人工注释 ----

mouse_markers <- list(Mural = c("Acta2", "Myh11" ),
                      OPC = c("Pcdh15","Tnr",'Pdgfra'),
                      EC = c("Flt1", "Erg"),
                      MG = c('P2ry12','Hexb','C3','Sall1','Siglech','Fcrls'),
                      #Skeletal_Muscle_Cells = c('Myh4','Tnnt3'),
                      #Schwann_Cells = c('Mpz'),
                      #Bcell = c('Cd19','Ms4a1','Cd22','Bank1','Cd79a','Cd20'),
                      #Tcell = c('Cd8','Cd3',"Cd4",'Trac'),
                      MP = c('Mrc1','Cd163','F13a1'),
                      #Epithelial_Cells = c('Cdh1','Pkp2','Cxadr'),
                      #VLMC = c('Slc47a1', 'Sema3g','Slit3','Atp1a2','Slit2'),
                      #LECs = c('Lyve1','Prox1'),
                      FB = c('Dcn','Col1a1','Col3a1'),
                      AS = c("Gfap",'Aqp4'),
                      OLG = c("Mobp","Mog"))

DotPlot(object = scesinglet_list_filter[[3]],
        features = mouse_markers ,
        cols = c("lightgrey", "blue"),
        #assay = "RNA",
        dot.scale = 4, 
        cluster.idents = FALSE, 
        scale.by = "radius",  
        scale = TRUE,  
        col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of rat marker")

rat_cell_type_map <- c(
  "0" = "OLG",
  "1" = "OLG",
  "2" = "FB",
  "3" = "FB",
  "4" = "None",  # Skeletal_Muscle_Cells
  "5" = "FB",
  "6" = "FB",
  "7" = "AS",
  "8" = "EC",
  "9" = "MP",
  "10" = "OPC",
  "11" = "Mural",
  "12" = "EC",
  "13" = "OLG", #
  "14" = "EC", 
  "15" = "MG", # ?
  "16" = "None", #Epithelial_Cells
  "17" = "None", #Schwann_Cells
  "18" = "None", # mt 太高
  "19" = "None", # Neuron
  "20" = "None", # 过少
  "21" = "None" # 过少
)

# 查看数据 确认是否存在低质量的群
VlnPlot(scesinglet_list_filter[[3]], features =  c("nFeature_RNA", "nCount_RNA", "mt_percent", "hb_percent", "rb_percent", "log10_genes_per_umi"), pt.size = 0, group.by = "RNA_snn_res.0.5")

scesinglet_list_filter[[3]]$cell_type <- setNames(rat_cell_type_map[as.character(scesinglet_list_filter[[3]]$RNA_snn_res.0.5)],
 Cells(scesinglet_list_filter[[3]]))

# 06 计算各群 marker ----

markers_list <- lapply(scesinglet_list_filter, function(df) {
   FindAllMarkers(df, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
})

# 07 按预设列表移除异常群 ----

sce_final_list <- list()

remove_cluster_list <- list(
  Human = c(5,9,15,14),
  Macaque = c(13),
  Goat = c(11),
  Rat = c(4,16:21)
)

for (species_name in names(scesinglet_list_filter)) {
    cat("  -> 正在处理:", species_name, "\n")
    sce <- scesinglet_list_filter[[species_name]]

    # 检查当前物种是否在移除列表中
    if (species_name %in% names(remove_cluster_list)) {
        clusters_to_filter <- remove_cluster_list[[species_name]]
        cat("     将移除以下聚类: ", paste(clusters_to_filter, collapse = ", "), "\n")

        # 使用 subset 函数高效移除细胞，invert = TRUE 表示保留不在列表中的细胞
        # 'seurat_clusters' 是 FindClusters(resolution=0.5) 默认保存聚类结果的列
        sce_clean <- subset(sce, subset = seurat_clusters %in% clusters_to_filter, invert = TRUE)

        sce_final_list[[species_name]] <- sce_clean

    } else {
        # 如果当前物种不在移除列表中，则直接将原对象存入最终列表
        cat("     无需移除聚类。\n")
        sce_final_list[[species_name]] <- sce
    }
}

# 08 结果检查 ----

DotPlot(object = sce_rmoutlie_list[[1]],features = human_markers ,group.by = 'cell_type',
 dot.scale = 4, cluster.idents = FALSE, scale.by = "radius", scale = TRUE, col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of human marker")
  DotPlot(object = sce_rmoutlie_list[[2]],features = monkey_markers ,group.by = 'cell_type',
 dot.scale = 4, cluster.idents = FALSE, scale.by = "radius", scale = TRUE, col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of macaque marker")
  DotPlot(object = sce_rmoutlie_list[[3]],features = mouse_markers ,group.by = 'cell_type',
 dot.scale = 4, cluster.idents = FALSE, scale.by = "radius", scale = TRUE, col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of rat marker")
  DotPlot(object = sce_rmoutlie_list[[4]],features = goat_markers ,group.by = 'cell_type',
 dot.scale = 4, cluster.idents = FALSE, scale.by = "radius", scale = TRUE, col.min = 0
)+ RotatedAxis() + ggtitle("High Expression Genes of goat marker")

# 09 按细胞类型移除 UMAP 离群值 ----

umap_rmoutlie_plot <- list()
sce_rmoutlie_list <- list()
for (species_name in names(sce_final_list)) {
  cat("  -> 正在处理:", species_name, "\n")
  sce <- sce_final_list[[species_name]]

  # 移除
  result <- remove_outliers_dbscan(sce, "cell_type", "umap",minPts = 50)
  umap_rmoutlie_plot[[species_name]] <- result$umap_plot
  sce_rmoutlie_list[[species_name]] <- result$sce_clean
  rm(sce,result)
  cat("  -> 处理完成:", species_name, "\n")
}

# 10 汇总保留细胞数并保存 ----

# 记录质控前后的细胞数
qc_summary <- data.frame(
  Species = names(sce_merged_list),
  Cells_Before = sapply(sce_merged_list, ncol)
)

# 更新总结表
qc_summary$Cells_After <- sapply(sce_rmoutlie_list, ncol)
qc_summary$Cells_Removed <- qc_summary$Cells_Before - qc_summary$Cells_After
qc_summary$Percent_Removed <- scales::percent(qc_summary$Cells_Removed / qc_summary$Cells_Before, accuracy = 0.1)

# 打印质控总结表
cat("\n--- 质控总结 ---\n")
print(qc_summary)
cat("\n")

save.image("annotation.RData")
