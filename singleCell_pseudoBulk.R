# script to perform pseudo-bulk DGA
# setwd("~/Desktop/demo/single_cell_DEG")

library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)


# get data
eh <- ExperimentHub()
query(eh, "Kang")

sce <- eh[["EH2259"]]
seu.obj <- as.Seurat(sce, data = NULL)
View(seu.obj@meta.data)


# QC and filtering
# explore QC


# get mito percent
seu.obj$mitoPercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
View(seu.obj@meta.data)

# filter
seu.filtered <- subset(seu.obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
         nCount_originalexp > 800 & 
         mitoPercent < 5 &
         multiplets == 'singlet')

seu.obj
seu.filtered

# run Seurat's standard workflow steps
seu.filtered <- NormalizeData(seu.filtered)
seu.filtered <- FindVariableFeatures(seu.filtered)
seu.filtered <- ScaleData(seu.filtered)
seu.filtered <- RunPCA(seu.filtered)
ElbowPlot(seu.filtered)
seu.filtered <- RunUMAP(seu.filtered, dims = 1:20)

# visualize 
cell_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'cell', label = TRUE)
cond_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'stim')

cell_plot|cond_plot

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(seu.filtered@meta.data)
seu.filtered$samples <- paste0(seu.filtered$stim, seu.filtered$ind)

DefaultAssay(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                    group.by = c("cell", "samples"),
                    assays = 'originalexp',
                    slot = "counts",
                    return.seurat = FALSE)

cts <- cts$originalexp

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                 f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_bcell <- cts.split.modified$`B cells`


# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('stim', samples), 'Stimulated', 'Control')) %>%
  column_to_rownames(var = 'samples')

# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_bcell,
                       colData = colData,
                       design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)



# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_Stimulated_vs_Control")
res



