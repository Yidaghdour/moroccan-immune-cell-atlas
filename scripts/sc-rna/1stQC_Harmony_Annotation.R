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
GData<-readRDS("merged_pools.rds")

#### QC and Filtering ####

GData[["percent.mt"]] <- PercentageFeatureSet(GData, pattern = "^MT-")

# Generate violin and feature scatter plots per pool 
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
Idents(GData) <- "Pool"

plots <- lapply(features, function(feat) {
  VlnPlot(GData, features = feat, pt.size = 0) +            # no jittered points
    geom_boxplot(width = 0.08,                              # slim box so violin shows
                 outlier.shape = 16, outlier.size = 0.7,    # show only outliers
                 fill = NA, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23,    # mean (diamond-ish)
                 size = 2, stroke = 0.6) +
    stat_summary(fun = median, geom = "crossbar",           # median line
                 width = 0.4, fatten = 0, size = 0.5)
})
wrap_plots(plots, ncol = 3)

p1 <- FeatureScatter(GData, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(GData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

## Lean filtering per pool to remove cells with very low number of genes expressed and number of UMIs and very high mitochondrial reads percentage
pool_col <- "Pool" 
nmads <- 7

md <- GData@meta.data

# Work on transformed metrics
md$log10_nCount   <- log10(md$nCount_RNA + 1)
md$log10_nFeature <- log10(md$nFeature_RNA + 1)

# Compute per-group thresholds
mad_thresholds <- function(x, g, nmads=7) {
  # returns a data.frame with group, med, mad, lower, upper
  groups <- unique(g)
  out <- lapply(groups, function(gi){
    xi <- x[g == gi]
    m  <- median(xi, na.rm = TRUE)
    s  <- mad(xi, constant = 1, na.rm = TRUE)  # constant=1 => raw MAD; scale is absorbed by nmads
    data.frame(
      group = gi,
      med = m,
      mad = s,
      lower = m - nmads*s,
      upper = m + nmads*s
    )
  })
  do.call(rbind, out)
}

g <- md[[pool_col]]

thr_lib  <- mad_thresholds(md$log10_nCount,   g, nmads)
thr_feat <- mad_thresholds(md$log10_nFeature, g, nmads)
thr_mt   <- mad_thresholds(md$percent.mt,     g, nmads)

# map thresholds back to cells
idx <- match(g, thr_lib$group)
md$lib_lower  <- thr_lib$lower[idx];  md$lib_upper  <- thr_lib$upper[idx]
idx <- match(g, thr_feat$group)
md$feat_lower <- thr_feat$lower[idx]; md$feat_upper <- thr_feat$upper[idx]
idx <- match(g, thr_mt$group)
md$mt_upper   <- thr_mt$upper[idx]   

# flag outliers (low lib, low feat, high mito)
md$discard_mad <- (md$log10_nCount   < md$lib_lower)  |
  (md$log10_nFeature < md$feat_lower) |
  (md$percent.mt     > md$mt_upper)


GData@meta.data <- md
keep_cells <- rownames(md)[!md$discard_mad]

GData <- subset(GData, cells = keep_cells)


#### SCTransform per sample ####

GData$Pool <- as.character(GData$Pool)
GData$library <- paste(GData$Sample.ID, paste0("pool", GData$Pool), sep="_")
GData$library <- as.factor(GData$library)
invisible(gc())

GData <- SplitObject(GData, split.by = "library")
invisible(gc())
libs <- lapply(X = GData, FUN = SCTransform, method = "glmGamPoi", conserve.memory = TRUE, verbose = FALSE)
#### Merge and Set HVGs ####
GData <- merge(libs[[1]], libs[-1])
DefaultAssay(GData) <- "SCT"
features <- SelectIntegrationFeatures(object.list = libs, nfeatures = 3000)
VariableFeatures(GData) <- features
libs <- NULL
invisible(gc())


#### Run pca and harmony ####
set.seed(12)
options(uwot_parallel = FALSE)

GData <- RunPCA(GData, npcs = 100, verbose = TRUE)
dims_use <- 1:22
GData <- FindNeighbors(GData, dims = dims_use, reduction = "pca", k.param = 40, verbose = TRUE)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "unintegrated_clusters", random.seed  = 12345, verbose = TRUE)
GData <- RunUMAP(GData, dims = dims_use, reduction = "pca", reduction.name = "umap.unintegrated", n.neighbors = 40, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())

GData <- harmony::RunHarmony(object = GData, group.by.vars = "Batch", reduction.use = "pca",
                           dims.use = dims_use, assay.use = "SCT", reduction.save = "harmony",
                          verbose = TRUE)
GData <- FindNeighbors(GData, reduction = "harmony", dims = dims_use, k.param = 40)
GData <- FindClusters(GData, resolution = 0.2, cluster.name = "harmony_clusters")
GData <- RunUMAP(GData, reduction = "harmony", dims = dims_use, reduction.name = "umap.harmony", n.neighbors = 40, min.dist = 0.3, n.epochs = 500, seed.use = 12, n.threads = 1, verbose = TRUE)
invisible(gc())
ElbowPlot(GData, ndims = 100)


## FindAllMarkers annotation
Idents(GData) <- "harmony_clusters"
GData <- PrepSCTFindMarkers(GData)
markers <- FindAllMarkers(GData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
write.csv(markers, file = "find_all_markers_res0.2_dims22_1stQC.csv")

top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top_markers, file = "top20_find_all_markers_res0.2_dims22_1stQC.csv")


## Azimuth annotation
reference <- readRDS("pbmc_multimodal_2023.rds")

anchors <- FindTransferAnchors(
  reference = reference,
  query = GData,
  reference.assay = "SCT",
  query.assay = "SCT",
  normalization.method = "SCT",
  reference.reduction = "spca",
  k.filter = NA,
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

GData <- MapQuery(
  anchorset = anchors,
  query = GData,
  reference = reference,
  refdata = list(celltype.l1 = "celltype.l1"),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)


## scType annotation
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(GData[["SCT"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
GData_scaled <- if (seurat_package_v5) as.matrix(GData[["SCT"]]$scale.data) else as.matrix(GData[["SCT"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = GData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(GData@meta.data$harmony_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(GData@meta.data[GData@meta.data$harmony_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(GData@meta.data$harmony_clusters==cl)), 10)
}))
write.csv(cL_resutls, file = "sctype_scores_1stQC.csv")
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
write.csv(sctype_scores, file = "top_sctype_scores_1stQC.csv")


## Check UMAPS and annotate clusters 

DimPlot(GData, reduction = "umap.harmony", group.by = "Stimulation")
DimPlot(GData, reduction = "umap.harmony", group.by = "Batch")
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle")
DimPlot(GData, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE, raster = T)
DimPlot(GData, reduction = "umap.harmony", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T)
Idents(GData) <- "harmony_clusters"
Cluster_celltypes <- c("Memory T cells", "Cytotoxic T/NK cells", "Naive T cells", "B cells", "Unknown", "Myeloid cells",
                       "Plasma B cells")
names(Cluster_celltypes) <- levels(GData)
GData <- RenameIdents(GData, Cluster_celltypes)
DimPlot(GData, reduction = "umap.harmony", label = TRUE, pt.size = 0.5, label.size = 3, raster = FALSE)
GData$Cluster_celltypes <- Idents(GData)


## Save final RDS file post 1st QC filtering 
saveRDS(GData, file = "1stQC_harmony-dims22.rds")




