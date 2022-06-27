# script to perform trajectory analysis
# https://www.nature.com/articles/s41467-019-10291-0
# setwd("~/Desktop/demo/monocle3")

set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)


# read in data
markers <- read.delim('ABC_Marker.txt', header = T) # gene metadata
metadata <- read.delim('ABC_Meta.txt', header = T) # cell metadata
expr <- read.delim('ABC_umi_matrix_7551_cells.csv', header = T, sep = ',') # expression matrix



# create seurat object ---------------
expr.t <- t(expr)
seu.obj <- CreateSeuratObject(counts = expr.t)
View(seu.obj@meta.data)
seu.obj@meta.data <- merge(seu.obj@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id')
View(seu.obj@meta.data)
seu.obj@meta.data <- seu.obj@meta.data %>% 
  column_to_rownames(var = 'Row.names')
seu.obj$mitopercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 &
                    nFeature_RNA > 500 &
                    mitopercent < 10)


# subset my seurat object - B cells

unique(seu.obj.filtered@meta.data$population)

Idents(seu.obj.filtered) <- seu.obj.filtered$population
b.seu <- subset(seu.obj.filtered, idents = "b")
b.seu
unique(b.seu@meta.data$redefined_cluster)

# pre-processing using seurat
b.seu <- NormalizeData(b.seu)
b.seu <- FindVariableFeatures(b.seu)
b.seu <- ScaleData(b.seu)
b.seu <- RunPCA(b.seu)
b.seu <- FindNeighbors(b.seu, dims = 1:30)
b.seu <- FindClusters(b.seu, resolution = 0.9)
b.seu <- RunUMAP(b.seu, dims = 1:30, n.neighbors = 50)

a1 <- DimPlot(b.seu, reduction = 'umap', group.by = 'redefined_cluster', label = T)
a2 <- DimPlot(b.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2


# MONOCLE3 WORKFLOW ---------------------
# monocle3 requires cell_data_set object
# convert seurat object to cell_data_set object for monocle3





# ...1 Convert to cell_data_set object ------------------------

cds <- as.cell_data_set(b.seu)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)



# ...2. Cluster cells (using clustering info from seurat's UMAP)---------------------------
# let's use the clustering information have

# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- b.seu@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- b.seu@reductions$umap@cell.embeddings



# plot

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
           color_cells_by = "redefined_cluster",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 5]))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(redefined_cluster, monocle3_pseudotime, median), fill = redefined_cluster)) +
  geom_boxplot()




# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

FeaturePlot(b.seu, features = c('E2F2', 'STMN1', 'CD52'))


# visualizing pseudotime in seurat

b.seu$pseudotime <- pseudotime(cds)
Idents(b.seu) <- b.seu$redefined_cluster
FeaturePlot(b.seu, features = "pseudotime", label = T)



