# script to manipulate data
# setwd("~/Desktop/demo/data_manipulation_R/scripts")

#load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# 1. read in data --------------------
dat <- read.csv('../data/GSE183947_fpkm.csv')

# 2. Fetch metadata ----------
geo_id <- "GSE183947"

#gse <- getGEO(geo_id,GSEMatrix=TRUE)
# Error: The size of the connection buffer (131072) was not large enough                                          0s
# to fit a complete line:
#   * Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*1000)
gse <- getGEO(geo_id,GSEMatrix=TRUE)
gse[[1]]

metadata <- pData(phenoData(gse[[1]]))


# 3. subset ------------------
metadata <- metadata %>%
  select(1,2,10,11,17) %>% # select() columns
  rename(tissue=characteristics_ch1) %>% # rename() columns
  rename(metastasis=characteristics_ch1.1) %>%
  mutate(tissue = gsub('tissue:', '', tissue), # mutate() to change values in a column/create new variables
         metastasis = gsub('metastasis:', '', metastasis))


# 4. reshape ---------------
# wide format
dat[1:10,1:10]

# long
dat.long <- dat %>%
  rename(gene=X) %>%
  gather(key = 'samples', value = 'FPKM', -gene) # to convert wide format to long = facilitates to add more data to the dataframe

head(dat.long)


# 5. merge -------------
dat.long.metadata <- dat.long %>%
  left_join(., metadata, by = c("samples" = "description"))

# 6. Explore data -------------
# summary FPKM statistics for samples from breast tumor and normal tissues

dat.long.metadata %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean = mean(FPKM),
            median = median(FPKM),
            q1 = quantile(FPKM, 0.25),
            q2 = quantile(FPKM, 0.5),
            q3 = quantile(FPKM, 0.75))


# checking expression for housekeeping gene 
# A housekeeping gene ideally should be stable, expressed in the cells and tissues of interest 
# that do not show changes under the experimental conditions or disease state.
dat.long.metadata %>%
  filter(gene == 'GAPDH') %>%
  group_by(gene, tissue) %>%
  summarize(mean = mean(FPKM),
            median = median(FPKM),
            q1 = quantile(FPKM, 0.25),
            q2 = quantile(FPKM, 0.5),
            q3 = quantile(FPKM, 0.75))


# 7. find top 10 genes with highest -----------------
dat.long.metadata %>%
  filter(gene > 0) %>%
  group_by(gene) %>%
  mutate(mean_FPKM = round(mean(FPKM),2)) %>%
  select(gene, mean_FPKM) %>%
  distinct() %>%
  arrange(as.integer(-mean_FPKM)) %>%
  top_n(10)






