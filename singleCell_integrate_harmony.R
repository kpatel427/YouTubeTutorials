# script to integrate across conditions using Harmony
# setwd("~/Desktop/demo/single_cell_integrate/")

# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)



# get data -------------------------
AvailableData()
# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")
str(ifnb)


# QC and filtering
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
View(ifnb@meta.data)
# explore QC

# filter
ifnb
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')


# run Harmony -----------
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

before|after
