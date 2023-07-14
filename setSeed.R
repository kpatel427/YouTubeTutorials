# script to demonstrate set.seed
# setwd("~/Desktop/demo/set_seed/scripts")


# Example 1: -----------------------------------

# generate 10 random numbers
runif(10)

runif(10)

set.seed(123)
runif(10)

set.seed(123)
runif(10)

# Example 2: -----------------------------------

library(Seurat)
library(tidyverse)
library(gridExtra)

load('../data/pbmc.seurat.filtered.RData')

# before setting seed
before_setseed_run1 <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20, seed.use = NULL)
a1 <- DimPlot(before_setseed_run1, reduction = 'umap')

before_setseed_run2 <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20, seed.use = NULL)
a2 <- DimPlot(before_setseed_run2, reduction = 'umap')

grid.arrange(a1, a2, ncol = 2)


# After setting seed
after_setseed_run1 <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20, seed.use = 77)
b1 <- DimPlot(after_setseed_run1, reduction = 'umap')

after_setseed_run2 <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20, seed.use = 77)
b2 <- DimPlot(after_setseed_run2, reduction = 'umap')

grid.arrange(b1, b2, ncol = 2)



