# script to demonstrate reading single cell matrices in various format 
# and converting to seurat object 
# setwd("/Users/kr/Desktop/demo/InputFormat_singleCell")

# load libraries
library(Seurat)
library(SeuratDisk)

# .RDS format
rds_obj <- readRDS('ependymal_cells.rds')

# 10X CellRanger .HDF5 format 
hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
           use.names = TRUE,
           unique.features = TRUE)
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

# .mtx file
mtx_obj <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
        features = "raw_feature_bc_matrix/features.tsv.gz",
        cells = "raw_feature_bc_matrix/barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

# .loom files
loom_oj <- Connect(filename = "adult-hem-organs-10X-bone-marrow.loom", mode = 'r')
seurat_loom <- as.Seurat(loom_oj)

# .h5ad format 
# step 1: convert AnnData object to an h5Seurat file
Convert("adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

# step 2: Load h5Seurat file into a Seurat object 
seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.h5seurat")


