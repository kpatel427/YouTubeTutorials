# script to download data from TCGA using TCGAbiolinks
# setwd("~/Desktop/demo/TCGAbiolinks")


library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)


# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-BRCA')



# building a query
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
         data.category = 'Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-LL-A73Y-01A-11R-A33J-07', 'TCGA-E2-A1IU-01A-11R-A14D-07','TCGA-AO-A03U-01B-21R-A10J-07'))

getResults(query_TCGA)

# download data - GDCdownload
GDCdownload(query_TCGA)


# prepare data
tcga_brca_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, 'fpkm_unstrand')


# build a query to retrieve DNA methylation data --------------
query_methly <- GDCquery(project = 'TCGA-GBM',
         data.category = 'DNA Methylation',
         platform = 'Illumina Human Methylation 27',
         access = 'open',
         data.type = 'Methylation Beta Value',
         barcode = c('TCGA-19-0962-01B-01D-0521-05', 'TCGA-06-0137-01A-01D-0218-05'))

output_query_methyl <- getResults(query_methly)

GDCdownload(query_methly)


# plot probes showing differences in beta values between samples                                               
dna.meth <- GDCprepare(query_methly, summarizedExperiment = TRUE)
assay(dna.meth)  

idx <- dna.meth %>% 
  assay %>% 
  rowVars() %>% 
  order(decreasing = TRUE) %>% 
  head(10)

# plot
pheatmap(assay(dna.meth)[idx,])


# download and visualize mutation data from TCGA ----------------------
query_mutation <- GDCquery(project = 'TCGA-BRCA',
         data.category = 'Simple Nucleotide Variation',
         access = 'open',
         barcode = c('TCGA-LL-A73Y-01A-11D-A33E-09,TCGA-LL-A73Y-10B-01D-A33H-09',
                     'TCGA-E9-A1NH-01A-11D-A14G-09,TCGA-E9-A1NH-11A-33D-A14G-09'))

output_query_mutation <- getResults(query_mutation)

GDCdownload(query_mutation)

maf <- GDCprepare(query_mutation, summarizedExperiment = TRUE)

# maftools utils to read and create dashboard
maftools.input <- read.maf(maf)

plotmafSummary(maf = maftools.input,
               addStat = 'median',
               rmOutlier = TRUE,
               dashboard = TRUE)

# oncoprint
oncoplot(maf = maftools.input,
         top = 10,
         removeNonMutated = TRUE)






