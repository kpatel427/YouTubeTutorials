# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/demo/singleCell_singleR/scripts")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Input Data 10X CellRanger .HDF5 format --------------
hdf5_obj <- Read10X_h5(filename = '../data/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5',
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering -----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)


# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# get reference data -----------
ref <- celldex::HumanPrimaryCellAtlasData()
View(as.data.frame(colData(ref)))

# expression values are log counts (log normalized counts)


# run SingleR (default mode) ---------
# default for SingleR is to perform annotation of each individual cell in the test dataset

pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')


# Annotation diagnostics ----------


# ...Based on the scores within cells -----------
pred
pred$scores

plotScoreHeatmap(pred)


# ...Based on deltas across cells ----------

plotDeltaDistribution(pred)




# ...Comparing to unsupervised clustering ------------

tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))




