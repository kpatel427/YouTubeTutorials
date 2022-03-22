# script to identify cluster identity -----------------
# Finding markers in every cluster
# Finding conserved markers 
# Finding markers DE between conditions

# setwd("~/Desktop/demo/single_cell_DEG/script")

set.seed(1234)

library(Seurat)
library(tidyverse)

# load data
ifnb_harmony <- readRDS('../ifnb_harmony.rds')
str(ifnb_harmony)
View(ifnb_harmony@meta.data)

# visualize data
clusters <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'stim')

condition|clusters

# findAll markers -----------------

FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')


# findConserved markers -------------

# Notes:
# slot depends on the type of the test used, 
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony) <- 'RNA'

DefaultAssay(ifnb_harmony)

markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                     ident.1 = 3,
                     grouping.var = 'stim')

head(markers_cluster3)

# let's visualize top features
FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')


# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))





# rename cluster 3 ident
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, `3` = 'CD16 Mono')

DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata
View(ifnb_harmony@meta.data)

# Settings cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual anotations (sometimes both)
# need to make sure each cell type forms a separate cluster

# setting Idents as Seurat annotations provided (also a sanity check!)
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)


# findMarkers between conditions ---------------------
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)
View(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(markers_cluster3)


FeaturePlot(ifnb_harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')



