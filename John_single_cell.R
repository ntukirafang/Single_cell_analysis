# Load libraries
library(Seurat)
library(dplyr)

# Step 1: Load Seurat objects
nxg_immune <- readRDS("NXG7_10_Named_reso.6_Seurat_integrated.rds")
gran_atlas <- readRDS("GranAtlas.rds")

# Subset immune cells from nxg_immune
nxg_immune_only <- subset(nxg_immune, idents = "immune")  # adjust this if immune cells are labeled differently

# Merge with gran_atlas
merged_obj <- merge(gran_atlas, nxg_immune_only, add.cell.ids = c("Gran", "NXG7_10"))

# Step 2: Rerun the Seurat pipeline
# Source the existing script or run your own analysis steps manually
# Assuming you are using a script like "complete_human_v10.R"
source("~/scripts/complete_human_v10.R")  # Update with actual path if needed

# Or, manually:
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- FindNeighbors(merged_obj, dims = 1:20)
merged_obj <- FindClusters(merged_obj, resolution = 0.6)
merged_obj <- RunUMAP(merged_obj, dims = 1:20)

# Step 3: Annotate the immune cells
# This assumes you have a method or signature list to label cell types
# Example (placeholder logic):
immune_obj <- subset(merged_obj, idents = "immune")  # again, adjust if needed
# Use marker genes or SingleR, scCATCH, etc. to annotate here

# Step 4: PCA on macrophages
# Subset macrophages (e.g., using cluster ID or marker gene)
macrophages <- subset(immune_obj, idents = "macrophage")  # update ident or metadata field as needed
macrophages <- RunPCA(macrophages, features = VariableFeatures(macrophages))
DimPlot(macrophages, reduction = "pca", group.by = "orig.ident")

# Step 5: Run DEG
# For example, comparing two conditions in macrophages
Idents(macrophages) <- "orig.ident"
deg_results <- FindMarkers(macrophages, ident.1 = "NXG7", ident.2 = "GranPat14", only.pos = FALSE)
write.csv(deg_results, "macrophage_DEGs_NXG7_vs_GranPat14.csv")
