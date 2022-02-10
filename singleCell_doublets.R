# filter out doublets: DoubletFinder
# setwd("~/Desktop/demo/single_cell_doublets/scripts")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

# create counts matrix
cts <- ReadMtx(mtx = '/Users/kr/Desktop/demo/single_cell_doublets/data/raw_feature_bc_matrix/matrix.mtx.gz',
        features = '/Users/kr/Desktop/demo/single_cell_doublets/data/raw_feature_bc_matrix/features.tsv.gz',
        cells = '/Users/kr/Desktop/demo/single_cell_doublets/data/raw_feature_bc_matrix/barcodes.tsv.gz')

cts[1:10,1:10]

# create Seurat object
pbmc.seurat <- CreateSeuratObject(counts = cts)
str(pbmc.seurat)

# QC and Filtering
# explore QC

pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')

pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)

pbmc.seurat
pbmc.seurat.filtered


# pre-process standard workflow
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)




## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, 
                                     PCs = 1:20, 
                                     pN = 0.25, 
                                     pK = pK, 
                                     nExp = nExp_poi.adj,
                                     reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")


# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)


