# script to manipulate gene expression data (GSE183947)
# setwd("~/Desktop/demo/data_manipulation_R/scripts")

# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# read in the data ---------
dat <- read.csv(file = "../data/GSE183947_fpkm.csv")
dim(dat)

# get metadata --------
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
# Error: The size of the connection buffer (131072) was not large enough                                          0s
# to fit a complete line:
#   * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# select, mutate, rename ------------
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))


# looking at gene expression data ---------
head(dat)

# reshaping data - from wide to long--------
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)


# join dataframes = dat.long + metadata.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description")) 

# explore data ------
# filter, group_by, summarize and arrange 
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)

