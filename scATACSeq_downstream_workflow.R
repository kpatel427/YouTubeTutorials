# script to perform differential peak accesibility analysis and visualize genomic regions using single-cell ATAC-Seq data
# Vignette: https://stuartlab.org/signac/articles/pbmc_vignette
# continued from: PART 1 link: https://youtu.be/yEKZJVjc5DY?si=cm0okOcJQMwkCvPo
# setwd("~/Desktop/demo/single_cell_ATACSeq")

# install packages
# remotes::install_github("stuart-lab/signac", ref="develop")
# install.packages("Matrix", type = "source")
# install.packages("irlba", type = "source")
# BiocManager::install("EnsDb.Hsapiens.v75")

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(SingleR)

# Pre-processed ATAC data ----------------------------------------------------
pbmc
DimPlot(object = pbmc, label = TRUE) + NoLegend()

# Create a gene activity matrix ------------------------------------------------
gene.activities <- GeneActivity(pbmc)
gene.activities[1:10,1:10]

# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc@assays
pbmc <- NormalizeData(object = pbmc,
              assay = 'RNA',
              normalization.method = 'LogNormalize',
              scale.factor = median(pbmc$nCount_RNA))


# to interpret ATAC-Seq clusters, visualizing activity of canonical marker genes
# assuming a general correspondence between gene body/promoter accessibility and gene expression which may not always be the case

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(pbmc, features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
            pt.size = 0.1,
            max.cutoff = 'q95',
            ncol = 3)


# Integrating with scRNA-Seq data. ---------------------------------------------
# Link to pre-processed RNA-Seq data: https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS('data/pbmc_10k_v3.rds')
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

View(pbmc_rna@meta.data)

# plot them before integrating
p1 <- DimPlot(pbmc, reduction = 'umap') + NoLegend() + ggtitle('scATAC-Seq')
p2 <- DimPlot(pbmc_rna, reduction = 'umap', group.by = 'celltype', repel = TRUE, label = TRUE) + ggtitle('scRNA-Seq') + NoLegend()

p1 | p2

# ** Should have the prior knowledge of cell types expected in your query dataset when using ref dataset


# ....Transfer Anchors by Seurat --------------
# Identify anchors

transfer.anchors <- FindTransferAnchors(reference = pbmc_rna,
                    query = pbmc,
                    reduction = 'pcaproject') # CCA is very slow


predicted.labels <- TransferData(anchorset = transfer.anchors,
             refdata = pbmc_rna$celltype,
             weight.reduction = pbmc[['lsi']],
             dims = 2:30)
head(predicted.labels)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
View(pbmc@meta.data)


plot1 <- DimPlot(pbmc, 
        reduction = 'umap',
        group.by = 'predicted.id',
        label = TRUE,
        repel = TRUE) + NoLegend() + ggtitle('scATAC-Seq')

plot2 <- DimPlot(pbmc_rna, 
                 reduction = 'umap',
                 group.by = 'celltype',
                 label = TRUE,
                 repel = TRUE) + NoLegend() + ggtitle('scRNA-Seq')

plot1 | plot2



# Finding differentially accessible peaks between cell types -------------------
Idents(pbmc) <- pbmc$predicted.id


# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'ATAC'

da_peaks <- FindMarkers(object = pbmc,
            ident.1 = 'CD4 Naive',
            ident.2 = 'CD14+ Monocytes',
            test.use = 'LR',
            latent.vars = 'nCount_ATAC')

head(da_peaks)

da_plot1 <- VlnPlot(object = pbmc,
        features = rownames(da_peaks)[1],
        pt.size = 0.1,
        idents = c('CD4 Naive','CD14+ Monocytes'))

da_plot2 <- FeaturePlot(object = pbmc,
            features = rownames(da_peaks)[1],
            pt.size = 0.1)

da_plot1 | da_plot2


# fold change between two groups of cells
fc <- FoldChange(object = pbmc, ident.1 = 'CD4 Naive', ident.2 = 'CD14+ Monocytes')
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE),]
head(fc)



# plotting genomic regions -----------------------------------


# set plotting order

levels(pbmc) <- unique(pbmc$predicted.id)

CoveragePlot(object = pbmc,
             region = rownames(da_peaks)[1],
             extend.upstream = 40000,
             extend.downstream = 20000)


# create interactive version of these plots?
CoverageBrowser(pbmc, region = 'CD8A')


