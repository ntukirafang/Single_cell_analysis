# Install required dependencies
install.packages(c("ggplot2", "cowplot", "Matrix", "dplyr"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Seurat from Bioconductor
BiocManager::install("Seurat")
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
# Alternatively, install Seurat from CRAN
install.packages("Seurat")

# Load the Seurat package
library(Seurat)

# Install clustree
install.packages("clustree")
library("clustree")

# Install matrix to figure out %%
library(dplyr)
library(magrittr)

# input and output pathway define
directory<- "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/"

# Read Data from cell ranger analysis
Health1.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/Health/CD11b_ROS_neg_rep1")
Health2.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/Health/CD11b_ROS_neg_rep2")
EAE_ROS_neg1.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/EAE/EAE_CD11b_ROS_neg_rep1")
EAE_ROS_neg2.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/EAE/EAE_CD11b_ROS_neg_rep2")
EAE_ROS_pos1.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/EAE/EAE_CD11b_ROS_pos_rep1")
EAE_ROS_pos2.data <- Read10X(data.dir = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/EAE/EAE_CD11b_ROS_pos_rep2")

# Create object of single cell analysis
Health1 <- CreateSeuratObject(counts = Health1.data, project = "Health1", min.cells = 3, min.features = 200)
Health2 <- CreateSeuratObject(counts = Health2.data, project = "Health2", min.cells = 3, min.features = 200)
EAE_ROS_neg1 <- CreateSeuratObject(counts = EAE_ROS_neg1.data, project = "EAE_ROS_neg1", min.cells = 3, min.features = 200)
EAE_ROS_neg2 <- CreateSeuratObject(counts = EAE_ROS_neg2.data, project = "EAE_ROS_neg2", min.cells = 3, min.features = 200)
EAE_ROS_pos1 <- CreateSeuratObject(counts = EAE_ROS_pos1.data, project = "EAE_ROS_pos1", min.cells = 3, min.features = 200)
EAE_ROS_pos2 <- CreateSeuratObject(counts = EAE_ROS_pos2.data, project = "EAE_ROS_pos2", min.cells = 3, min.features = 200)

# Merge the datasets
EAE <- merge(Health1, y = list(Health2, EAE_ROS_neg1, EAE_ROS_neg2, EAE_ROS_pos1, EAE_ROS_pos2), 
             add.cell.ids = c("Health1", "Health2", "EAE_ROS_neg1", "EAE_ROS_neg2", "EAE_ROS_pos1", "EAE_ROS_pos2"), 
             project = "EAE")

# Calculate the percentage of mitochondrial genes
EAE[["percent.mt"]] <- PercentageFeatureSet(EAE, pattern = "^MT-")

# Plot QC metrics
VlnPlot(EAE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(EAE, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(EAE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter low-quality cells
EAE <- subset(EAE, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)

# Normalize data
EAE <- NormalizeData(EAE, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features and plot
EAE <- FindVariableFeatures(EAE, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(EAE), 10)
plot1 <- VariableFeaturePlot(EAE)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(EAE)

# Seurat integration
EAE.list <- SplitObject(EAE, split.by = "orig.ident")

# Filter genes for integration
filter_genes <- function(EAE) {
  EAE <- subset(EAE, features = all.genes)
  return(EAE)
}
EAE.list <- lapply(EAE.list, filter_genes)

# Select integration features
features <- SelectIntegrationFeatures(object.list = EAE.list)

# Scale data
EAE <- ScaleData(EAE, features = all.genes)

# PCA analysis
EAE <- RunPCA(EAE, features = VariableFeatures(object = EAE))
print(EAE[["pca"]], dims = 1:12, nfeatures = 20)
ElbowPlot(EAE)

# Find neighbors
EAE <- FindNeighbors(EAE, dims = 1:12)

# PCA visualization and clustering
VizDimLoadings(EAE, dims = 1:12, reduction = "pca")
DimPlot(EAE, reduction = "pca")
DimHeatmap(EAE, dims = 1:12, cells = 500, balanced = TRUE)

# JackStraw analysis
EAE <- JackStraw(EAE, num.replicate = 100)
EAE <- ScoreJackStraw(EAE, dims = 1:12)
JackStrawPlot(EAE, dims = 1:12)

# Find clusters
EAE <- FindNeighbors(EAE, dims = 1:12)
resolutions <- seq(0.1, 1.0, by = 0.1)
for (res in resolutions) {
  EAE <- FindClusters(EAE, resolution = res)
}
print(colnames(EAE@meta.data))
clustree(EAE)
EAE <- FindClusters(EAE, resolution = 0.8)

# Run tSNE and UMAP
EAE <- RunTSNE(EAE, dims = 1:12)
DimPlot(EAE, reduction = "tsne", group.by = "RNA_snn_res.0.8", label = TRUE, pt.size = 0.5)

EAE <- RunUMAP(EAE, dims = 1:12)
DimPlot(EAE, reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE, pt.size = 0.5)

# Find marker genes for clusters
EAE <- JoinLayers(EAE)
EAE.markers <- FindAllMarkers(EAE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
EAE.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

# Save markers to CSV file
write.csv(EAE.markers, file = file.path(directory, "EAE_all_marker.csv"), row.names = TRUE)

# Top10 DoHeatmap(log2 foldchange) 
EAE.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EAE, features = top10$gene) + NoLegend()
DoHeatmap(EAE, features = top10$gene)

# Top10 DoHeatmap(p value) 
EAE.markers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = p_val) -> top10
DoHeatmap(EAE, features = top10$gene) + NoLegend()
DoHeatmap(EAE, features = top10$gene)

# Find independent cluster marker(Not necessary)
#cluster0.markers <- FindMarkers(EAE, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Extract expression data for the specified genes along with cluster information
genes_of_interest_0 <- c("Spp1", "Elf1", "Fli1", "Foxp1", "Maf", "Mafb", "Mef2a", "Mef2c", "Rela", "Mef2c", "Runx1", "Tcf4")
expression_data_0 <- FetchData(EAE, vars = c(genes_of_interest_0, "seurat_clusters"))
genes_of_interest_1 <- c("Spp1", "Ets2", "Irf1", "Maff", "Max", "Nfkb1")
expression_data_1 <- FetchData(EAE, vars = c(genes_of_interest_1, "seurat_clusters"))
genes_of_interest_11 <- c("Spp1", "Maz", "Nfya")
expression_data_11 <- FetchData(EAE, vars = c(genes_of_interest_11, "seurat_clusters"))
genes_of_interest_12 <- c("Spp1", "Ctcf")
expression_data_12 <- FetchData(EAE, vars = c(genes_of_interest_12, "seurat_clusters"))
genes_of_interest_14 <- c("Spp1", "Lef1", "Nfatc1", "Rora")
expression_data_14 <- FetchData(EAE, vars = c(genes_of_interest_14, "seurat_clusters"))
genes_of_interest_2 <- c("Spp1", "Prdm1")
expression_data_2 <- FetchData(EAE, vars = c(genes_of_interest_2, "seurat_clusters"))
genes_of_interest_3 <- c("Spp1", "Nr4a2")
expression_data_3 <- FetchData(EAE, vars = c(genes_of_interest_3, "seurat_clusters"))
genes_of_interest_4 <- c("Spp1", "Bcl6", "Smad4", "Stat3", "Vdr")
expression_data_4 <- FetchData(EAE, vars = c(genes_of_interest_4, "seurat_clusters"))
genes_of_interest_5 <- c("Spp1", "Bach1", "Rara")
expression_data_5 <- FetchData(EAE, vars = c(genes_of_interest_5, "seurat_clusters"))
genes_of_interest_6 <- c("Spp1", "Rel", "Zfp467")
expression_data_6 <- FetchData(EAE, vars = c(genes_of_interest_6, "seurat_clusters"))
genes_of_interest_7 <- c("Spp1", "Mafg", "Stat2")
expression_data_7 <- FetchData(EAE, vars = c(genes_of_interest_7, "seurat_clusters"))


head(expression_data)


# Extract data for specific clusters, for example, cluster 1 and cluster 2
cluster0_data <- subset(expression_data_0, seurat_clusters == 0)
cluster1_data <- subset(expression_data_1, seurat_clusters == 1)
cluster2_data <- subset(expression_data_2, seurat_clusters == 2)
cluster3_data <- subset(expression_data_3, seurat_clusters == 3)
cluster4_data <- subset(expression_data_4, seurat_clusters == 4)
cluster5_data <- subset(expression_data_5, seurat_clusters == 5)
cluster6_data <- subset(expression_data_6, seurat_clusters == 6)
cluster7_data <- subset(expression_data_7, seurat_clusters == 7)
cluster11_data <- subset(expression_data_11, seurat_clusters == 11)
cluster12_data <- subset(expression_data_12, seurat_clusters == 12)
cluster14_data <- subset(expression_data_14, seurat_clusters == 14)

# Save the results to CSV files
write.csv(cluster0_data, file = file.path(directory, "cluster0_gene_expression.csv"), row.names = TRUE)
write.csv(cluster1_data, file = file.path(directory, "cluster1_gene_expression.csv"), row.names = TRUE)
write.csv(cluster2_data, file = file.path(directory, "cluster2_gene_expression.csv"), row.names = TRUE)
write.csv(cluster3_data, file = file.path(directory, "cluster3_gene_expression.csv"), row.names = TRUE)
write.csv(cluster4_data, file = file.path(directory, "cluster4_gene_expression.csv"), row.names = TRUE)
write.csv(cluster5_data, file = file.path(directory, "cluster5_gene_expression.csv"), row.names = TRUE)
write.csv(cluster6_data, file = file.path(directory, "cluster6_gene_expression.csv"), row.names = TRUE)
write.csv(cluster7_data, file = file.path(directory, "cluster7_gene_expression.csv"), row.names = TRUE)
write.csv(cluster11_data, file = file.path(directory, "cluster11_gene_expression.csv"), row.names = TRUE)
write.csv(cluster12_data, file = file.path(directory, "cluster12_gene_expression.csv"), row.names = TRUE)
write.csv(cluster14_data, file = file.path(directory, "cluster14_gene_expression.csv"), row.names = TRUE)

#Violine plot different cell type of gene(# marker_genes can adjust anytime)

#marker_genes = c("Spp1", "Elf1", "Fli1", "Foxp1", "Maf", "Mafb", "Mef2a", "Mef2c", "Rela", "Mef2c", "Runx1", "Tcf4")
#marker_genes = c("Spp1", "Ets2", "Irf1", "Maff", "Max", "Nfkb1")
#marker_genes = c("Spp1", "Maz", "Nfya")
#marker_genes = c("Spp1", "Ctcf")
#marker_genes = c("Spp1", "Lef1", "Nfatc1", "Rora")
marker_genes = c("Spp1", "Prdm1")
#marker_genes = c("Spp1", "Nr4a2")
#marker_genes = c("Spp1", "Bcl6", "Smad4", "Stat3", "Vdr")
#marker_genes = c("Spp1", "Bach1", "Rara")
#marker_genes = c("Spp1", "Rel", "Zfp467")
#marker_genes = c("Spp1", "Mafg", "Stat2")
VlnPlot(EAE, features =marker_genes)
VlnPlot(EAE, features =marker_genes, pt.size= 0)
FeaturePlot(EAE, features = marker_genes )
RidgePlot(EAE, features = marker_genes )

#Coexpression

FeaturePlot(object = EAE,
            features = c("Spp1", "Prdm1"),
            reduction = "umap", pt.size = 0.8, blend = TRUE, blend.threshold = 1, order = T)

# Define cell type indentity to the all cluster(Not yet to finish)
levels(EAE)
new.cluster.ids <- c("Mg_Health_1", "MP_EAE_1", "MP_EAE_2", "MP_EAE_3", "MP_EAE_4", "MP_EAE_5",
                     "Mg_Health_2", "Mg_EAE_1", "MP_EAE_6", "Mg_EAE_2", "DCs_EAE", "Mg_EAE_3", 
                     "MP_EAE_7", "Mg_EAE_4", "T_cells", "Mg_Health_3")
names(new.cluster.ids) <- levels(EAE)
EAE <- RenameIdents(EAE, new.cluster.ids)
DimPlot(EAE, reduction = "umap", label = TRUE, pt.size = 0.5)

# Save RD5 as database
saveRDS(EAE, file = "D:/admin/Desktop/SOP_PS/Duke_data/Reduced dims/SC_RNA/EAE_final.rds")
