library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(SeuratData)
library(readxl)
library(presto)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(SeuratWrappers)
library(glmGamPoi)
library(harmony)
options(future.globals.maxSize = 1e9)


setwd("/.../.../")
GData <- readRDS("merged_QC_split_SCT_merged_harmony-azimuth_dims27_2nd-MADs_scDblFinder_NoStressed.rds")

DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = TRUE, pt.size = 0.5, label.size = 4, repel = TRUE, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Cluster_celltypes", label = TRUE, pt.size = 0.5, label.size = 4, repel = TRUE, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 0.5, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Stimulation", cols = c("grey", "darkred"), label = F, pt.size = 0.5, raster = FALSE)


Idents(GData) <- GData$harmony_clusters
Cluster_celltypes <- c("Cytotoxic T/NK cells", "Memory T cells", "Naive CD4 T cells", "B cells", "Naive CD8 T cells", "Monocytes",
                       "Dendritic cells")
names(Cluster_celltypes) <- levels(GData)
GData <- RenameIdents(GData, Cluster_celltypes)
GData$Cluster_celltypes <- as.character(Idents(GData))  # or factor(Idents(GData))
table(Idents(GData))



#### Monocytes ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData<-readRDS("Monocytes.rds")

GData <- subset(GData, subset = Stimulation == "Unstim")
GData <- subset(GData, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 4, raster = T)


#### Dendritic cells ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData<-readRDS("Dendritic.rds")

GData <- subset(GData, subset = Stimulation == "Unstim")
GData <- subset(GData, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 4, raster = T)


#### B cells ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData_B <-readRDS("Bcells.rds")

GData <- subset(GData_B, subset = Stimulation == "Unstim")
GData <- subset(GData_B, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 2, raster = T)



#### Cytotoxic T/NK ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData_Cyto <- readRDS("Cytotoxic.rds")

GData <- subset(GData_Cyto, subset = Stimulation == "Unstim")
GData <- subset(GData_Cyto, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 2, raster = T)


#### Memory T ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData_MemoryT <- readRDS("MemoryT.rds")

GData <- subset(GData_MemoryT, subset = Stimulation == "Unstim")
GData <- subset(GData_MemoryT, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 2, raster = T)


#### Naive CD4 ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData_NaiveCD4 <- readRDS("NaiveCD4T.rds")

GData <- subset(GData_NaiveCD4, subset = Stimulation == "Unstim")
GData <- subset(GData_NaiveCD4, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 2, raster = T)


#### Naive CD8 ####

setwd("/scratch/ts4594/RDS_files/Final_Final/New_UMAPs_celltype-level/")
GData_NaiveCD8 <- readRDS("NaiveCD8.rds")

#GData_NaiveCD8 <- subset(GData, subset = Cluster_celltypes == "Naive CD8 T cells")
#saveRDS(GData_NaiveCD8, file = "NaiveCD8.rds")

GData <- subset(GData_NaiveCD8, subset = Stimulation == "Unstim")
GData <- subset(GData_NaiveCD8, subset = Stimulation == "LPS")

GData <- SCTransform(GData, verbose = FALSE)
DefaultAssay(GData) <- "SCT"
GData <- FindVariableFeatures(GData, selection.method = "vst", nfeatures = 2000)
invisible(gc())
GData <- RunPCA(GData, npcs = 100, features = VariableFeatures(object = GData))
ElbowPlot(GData, ndims = 100)
dims_use <- 1:30
GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                             dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                             verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 80)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 80, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle", cols = c("#CD142B", "#355EA9"), label = F, pt.size = 2, raster = T)

