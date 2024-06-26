#package install
options(repos=structure(c(CRAN="http://cran.mirrors.hoobly.com/")))
install.packages("modeltools")
install.packages("segmented")
install.packages("sn")
install.packages("flexmix")
install.packages("prabclus")
install.packages("diptest")
install.packages("mvtnorm")
install.packages("robustbase")
install.packages("trimcluster")
install.packages("kernlab")
install.packages("foreach")
install.packages("ModelMetrics")
install.packages("car")
install.packages("mclust")
install.packages("caTools")
install.packages("RColorBrewer")
install.packages("gtools")
install.packages("devtools")
install.packages("prabclus")
install.packages("Rmisc")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("http://s3-us-west-2.amazonaws.com/10x.files/code/cellrangerRkit-2.0.0.tar.gz", repos=NULL)
install.packages("cellrangerRkit")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("stringi")
library(devtools)
library(Matrix)
library(dplyr)
library(ggplot2)
library(stringi)
install_github("satijalab/seurat", force = TRUE)
library(Seurat)
# input and output pathway define
setwd("C:/Users/jimmy.fang/Documents")
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
# Read Data from cell ranger analysis
pbmc.data <- Read10X(data.dir = "C:/Users/jimmy.fang/Documents/test/")
# Create object of single cell analysis
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#Automatic estimate read count of mitochondria
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Violine plot to present the read count 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Filter lower quality cell
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Normalize data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
# Find Variant feature and ploting
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
# Seurat integration
features <- SelectIntegradtionFeatures(object.list = pbmc)
pbmc.anchors <- FindiIntegrationAnchors(object.list = pbmc, dim = 1:20, anchor.features= features)
pbmc <- IntegrateData(anchorset = pbmc.anchor, dim = 1:20)
# Scale data
pbmc <- ScaleData(pbmc, features = all.genes)
# PCA analysis (Cell feature)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# Three different plot shows the PCA result
VizDimLoadings(pbmc, dims = 1:5, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# PCA score differentiation(Top 15)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
# Cell clustering (dims define the top 13th cluster to find)
pbmc <- FindNeighbors(pbmc, dims = 1:13)
# Cell cluster resolution(0.4-1.2 can much easier to different cell)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
# UMAP/tSNE find the cluster 
pbmc <- RunUMAP(pbmc, dims = 1:13)
DimPlot(pbmc, reduction = "umap")
# Find marker gene of cell clustering(Top5)
#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#Violine plot different cell type of gene
VlnPlot(pbmc, features = c("BRCA1", "EGFR"))
VlnPlot(pbmc, features = c("CD4", "CD8"), slot = "counts", log = TRUE)
# Based on Umap data plot different gene of cluster
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
# Top10 DoHeatmap 
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# Define cell type indentity to the all cluster
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# Save RD5 as database
saveRDS(pbmc, file = "C:/Users/jimmy.fang/Documents/GSE_final.rds")
