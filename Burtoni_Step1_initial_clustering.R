library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(tidyr)
library(ggplot2)


load("C:/Users/jayan/Downloads/Burtoni_snseq_data_IMC/Burtoni_snseq_data_IMC/burtoni.snseq.combined.sct.RData/burtoni.snseq.combined.sct.RData")


meta_path_burtoni <- "C:/Users/jayan/Downloads/Burtoni_snseq_data_IMC/Burtoni_snseq_data_IMC/sctypemarkers.hypo.unknown.csv"
metadata_burtoni <- read.csv(meta_path_burtoni, row.names = 1, header = TRUE, stringsAsFactors = FALSE)

head(metadata_burtoni)
dim(metadata_burtoni)

# Add to Seurat object
burtoni.snseq.combined.sct <- AddMetaData(burtoni.snseq.combined.sct, metadata = metadata_burtoni)


# Skipping prep and regular clustering steps since it was already done, running it again to have umap in the script
hypo <- burtoni.snseq.combined.sct <- RunUMAP(burtoni.snseq.combined.sct, reduction = "pca", dims = 1:20, seed.use = 42)

DimPlot(burtoni.snseq.combined.sct, reduction = "umap", group.by = "sctypemarkers.hypo.broad", label = TRUE, repel = TRUE) 


# Save as Robj
hypo.ast <- hypo
saveRDS(hypo.ast, file = "C://Users//jayan//Desktop//RNA-seq-Teleost//Burtoni_cluster.rds")


