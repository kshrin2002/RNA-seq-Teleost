library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(tidyr)
# Load data
pbmc.data <- Read10X(data.dir = "C:/Users/jayan/Downloads/41559_2021_1580_MOESM3_ESM/Supplemental_data/am_10x")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

featnames <- rownames(pbmc)

meta_path <- "C:/Users/jayan/Downloads/41559_2021_1580_MOESM3_ESM/Supplemental_data/2-object_metadata/Astyanax_mexicanus_meta_data.csv"

metadata <- read.csv(meta_path, row.names = 1, header = TRUE, stringsAsFactors = FALSE)

head(metadata)
dim(metadata)


# Add to Seurat object
pbmc <- AddMetaData(pbmc, metadata = metadata)


# Violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident", ncol = 2)


plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Normalize data
pbmc <- NormalizeData(pbmc)

# Find variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, set.seed = 0)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20, seed.use = 42)

DimPlot(pbmc, reduction = "umap", group.by = "Cluster", label = TRUE, repel = TRUE) 
pbmc.ast <- pbmc

saveRDS(pbmc.ast, file = "AstMex_Step1.RDS")
