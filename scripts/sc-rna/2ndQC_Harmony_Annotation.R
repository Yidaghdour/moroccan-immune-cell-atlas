library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)
library(readxl)
library(presto)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(SeuratWrappers)
library(glmGamPoi)
library(harmony)
library(SingleCellExperiment)
library(scDblFinder)
options(future.globals.maxSize = 1e9)

setwd("/.../.../")
GData<-readRDS("1stQC_harmony-dims22.rds")

#### Clean and remove unwanted data in seurat object after the 1st leniant QC filtering round ####

GData@reductions <- list()
GData@graphs <- list()
GData@neighbors <- list()
GData@commands <- list()

## Keep only RNA assays
DefaultAssay(GData) <- "RNA"
keep_assays <- "RNA"
for (nm in names(GData@assays)) {
  if (!nm %in% keep_assays) GData[[nm]] <- NULL
}

## Ensure sparse matrices 
to_dgc <- function(m) if (!inherits(m, "dgCMatrix")) as(m, "dgCMatrix") else m
if ("RNA" %in% names(GData@assays)) {
  GData[["RNA"]]@counts <- to_dgc(GData[["RNA"]]@counts)
  if (length(GData[["RNA"]]@data)) GData[["RNA"]]@data <- to_dgc(GData[["RNA"]]@data)
}

## drop unwanted metadata columns 
GData@meta.data <- GData@meta.data[, -c(21:28)]



#### scDblFinder doublet removal per pool ####
set.seed(123)
GData$scDblFinder_class <- NA_character_
GData$scDblFinder_score <- NA_real_

pool_col    <- "Pool"               
cluster_col <- "harmony_clusters"   

pool_vec <- as.character(GData@meta.data[[pool_col]])
pools    <- sort(unique(pool_vec))

for (p in pools) {
  message("Running scDblFinder for pool: ", p)
  
  cells_p <- rownames(GData@meta.data)[pool_vec == p]
  obj_p   <- subset(GData, cells = cells_p)
  
  sce <- as.SingleCellExperiment(obj_p, assay = "RNA")
  
  # ensure cluster labels exist as a vector
  cl <- colData(sce)[[cluster_col]]
  sce <- scDblFinder(sce, clusters = cl)
  
  GData$scDblFinder_class[colnames(obj_p)] <- as.character(colData(sce)$scDblFinder.class)
  GData$scDblFinder_score[colnames(obj_p)] <- as.numeric(colData(sce)$scDblFinder.score)
  
  rm(obj_p, sce); gc()
}


table(GData$scDblFinder_class, useNA = "ifany")

with(GData@meta.data,
     table(harmony_clusters, scDblFinder_class))

GData$scDblFinder_singlet <- GData$scDblFinder_class == "singlet"
table(GData$scDblFinder_singlet, useNA = "ifany")


#### 2nd QC and stringent filtering per pool and per major cell type ####

md <- GData@meta.data
pool_col     <- "Pool"      
celltype_col <- "Cluster_celltypes" 

nmads_count <- 4   # stringent for counts
nmads_feat  <- 4   # stringent for features
nmads_mt    <- 3   # keep mito reasonably strict
min_n_group <- 10  # minimum cells required in a pool×celltype group to use group-specific thresholds

md$log10_nCount   <- log10(md$nCount_RNA + 1)
md$log10_nFeature <- log10(md$nFeature_RNA + 1)

pool <- md[[pool_col]]
ct   <- md[[celltype_col]]
grp  <- paste(pool, ct, sep = "||")

# compute per-group median/MAD thresholds
get_thr <- function(x, g, nmads=3) {
  split_x <- split(x, g)
  stats <- lapply(split_x, function(v) {
    m <- median(v, na.rm = TRUE)
    s <- mad(v, constant = 1, na.rm = TRUE)  # raw MAD (robust); nmads controls stringency
    c(med = m, mad = s, lower = m - nmads*s, upper = m + nmads*s, n = sum(!is.na(v)))
  })
  out <- do.call(rbind, stats)
  out <- data.frame(group = rownames(out), out, row.names = NULL, check.names = FALSE)
  out
}

# thresholds for pool×celltype groups
thr_gc_count <- get_thr(md$log10_nCount,   grp, nmads = nmads_count)
thr_gc_feat  <- get_thr(md$log10_nFeature, grp, nmads = nmads_feat)
thr_gc_mt    <- get_thr(md$percent.mt,     grp, nmads = nmads_mt)

# thresholds for pool-only fallback
thr_p_count  <- get_thr(md$log10_nCount,   pool, nmads = nmads_count)
thr_p_feat   <- get_thr(md$log10_nFeature, pool, nmads = nmads_feat)
thr_p_mt     <- get_thr(md$percent.mt,     pool, nmads = nmads_mt)

# function to map thresholds back to cells with fallback when group is too small
map_thr_with_fallback <- function(md, grp_vec, pool_vec, thr_group, thr_pool, field) {
  # field is one of: "lower", "upper"
  idx_g <- match(grp_vec,  thr_group$group)
  idx_p <- match(pool_vec, thr_pool$group)
  
  # group sizes for fallback decision
  n_g <- thr_group$n[idx_g]
  
  val_g <- thr_group[[field]][idx_g]
  val_p <- thr_pool[[field]][idx_p]
  
  # fallback to pool thresholds if group n < min_n_group or missing
  use_pool <- is.na(n_g) | n_g < min_n_group | is.na(val_g)
  ifelse(use_pool, val_p, val_g)
}

# attach per-cell thresholds
md$count_lower <- map_thr_with_fallback(md, grp, pool, thr_gc_count, thr_p_count, "lower")
md$count_upper <- map_thr_with_fallback(md, grp, pool, thr_gc_count, thr_p_count, "upper")

md$feat_lower  <- map_thr_with_fallback(md, grp, pool, thr_gc_feat,  thr_p_feat,  "lower")
md$feat_upper  <- map_thr_with_fallback(md, grp, pool, thr_gc_feat,  thr_p_feat,  "upper")

md$mt_upper    <- map_thr_with_fallback(md, grp, pool, thr_gc_mt,    thr_p_mt,    "upper")  # mt upper only

# flag cells to discard (more stringent pass)
md$discard_strict <- (md$log10_nCount   < md$count_lower) |
  (md$log10_nCount   > md$count_upper) |
  (md$log10_nFeature < md$feat_lower)  |
  (md$log10_nFeature > md$feat_upper)  |
  (md$percent.mt     > md$mt_upper)

# write flag back + subset using barcodes (robust)
GData$discard_strict <- md$discard_strict

table(GData$scDblFinder_class)
table(GData$discard_strict)
table(GData$scDblFinder_class, GData$discard_strict)

# Filter cells post after scDblFinder and 2nd stringent QC
GData <- subset(GData,subset = scDblFinder_class == "singlet")
GData <- subset(GData,subset = discard_strict == "FALSE")



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
dims_use <- 1:27
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
write.csv(markers, file = "find_all_markers_res0.2_dims27_2ndQC.csv")

top_markers <- markers %>%
  group_by(cluster) %>%
  filter(!grepl('RPL', gene)) %>%
  filter(!grepl('RPS', gene)) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top_markers, file = "top20_find_all_markers_res0.2_dims27_2ndQC.csv")


## Azimuth
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
  refdata = list(predicted.l1 = "celltype.l1", predicted.l2 = "celltype.l2"),
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
write.csv(cL_resutls, file = "sctype_scores_2ndQC.csv")
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
write.csv(sctype_scores, file = "top_sctype_scores_2ndQC.csv")


## Check UMAPS and annotate clusters 

DimPlot(GData, reduction = "umap.harmony", group.by = "Stimulation")
DimPlot(GData, reduction = "umap.harmony", group.by = "Batch")
DimPlot(GData, reduction = "umap.harmony", group.by = "Lifestyle")
DimPlot(GData, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE, raster = T)
DimPlot(GData, reduction = "umap.harmony", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE)
DimPlot(GData, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE, raster = T)
DimPlot(GData, reduction = "umap.harmony", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE)

DimPlot(GData, reduction = "umap.harmony", group.by = "harmony_clusters", label = T)
Idents(GData) <- "harmony_clusters"
Cluster_celltypes <- c("Cytotoxic T/NK cells", "Memory T cells", "Naive CD4 T cells", "B cells", "Naive CD8 T cells", "Unknown", "Monocytes",
                       "Dendritic cells")
names(Cluster_celltypes) <- levels(GData)
GData <- RenameIdents(GData, Cluster_celltypes)
DimPlot(subset(GData, subset = Cluster_celltypes != "Unknown"), reduction = "umap.harmony", label = TRUE, pt.size = 0.5, label.size = 4, raster = FALSE)
DimPlot(GData, reduction = "umap.harmony", group.by = "Cluster_celltypes", label = TRUE, pt.size = 0.5, label.size = 4, raster = FALSE)
GData$Cluster_celltypes <- as.character(Idents(GData)) 
table(Idents(GData))



## Save final RDS file post 1st QC filtering 
saveRDS(GData, file = "2ndQC_harmony-dims27.rds")


