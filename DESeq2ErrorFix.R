# script to demonstrate how to fix: ncol(countData) == nrow(colData) is not TRUE
# setwd("~/Desktop/demo/DESeq2_tutorial/")

library(tidyverse)
library(DESeq2)

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('data/counts_data.csv')
head(counts_data)


# read in sample info
colData <- read.csv('data/sample_info.csv')


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))


# jumble rows in colData - recreate the error
colData <- colData[sample(1:nrow(colData)),]


# are they in the same order?
all(colnames(counts_data) == rownames(colData))
counts_data <- counts_data[, rownames(colData)]


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)
