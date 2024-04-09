# script to process single-cell ATAC-Seq data
# Vignette: https://stuartlab.org/signac/articles/pbmc_vignette
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


# what is a fragment file? How is it generated?
# https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
frag.file <- read.delim('data/atac_v1_pbmc_10k_fragments.tsv.gz', header = F, nrows = 10)
head(frag.file)


# 1. Read in data -----------------

counts <- Read10X_h5('data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
counts[1:10,1:10]


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "data/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

str(chrom_assay)

metadata <- read.csv(file = 'data/atac_v1_pbmc_10k_singlecell.csv', header = T, row.names = 1)
View(metadata)


# create a seurat Object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)

str(pbmc)


# ....Adding Gene Annotation -------------------

pbmc@assays$ATAC@annotation
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)


# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))


# add the gene information to the object
Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation



# 2. Computing QC ---------------------

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

View(pbmc@meta.data)

# ....Visualizing QC --------------------

colnames(pbmc@meta.data)
a1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

a1 | a2

VlnPlot(object = pbmc, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
        pt.size = 0.1,
        ncol = 6)


# ....Filtering poor quality cells --------------------

pbmc <- subset(x = pbmc,
               subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 30000 &
                 pct_reads_in_peaks > 15 & 
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)



# 3. Normalization and linear dimensional reduction ------------------
pbmc <- RunTFIDF(pbmc) # normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # selecting top features
pbmc <- RunSVD(pbmc) # dimensionality reduction

DepthCor(pbmc)

# 4. Non-linear dimensional reduction and Clustering -------------------
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

DimPlot(object = pbmc, label = TRUE) + NoLegend()







