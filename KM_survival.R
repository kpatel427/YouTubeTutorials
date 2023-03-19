# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")


library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# getting clinical data for TCGA-BRCA cohort -------------------
clinical_brca <- GDCquery_clinic("TCGA-BRCA")
any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(9,39,45)]


# looking at some variables associated with survival 
table(clinical_brca$vital_status)

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive",
                                         clinical_brca$days_to_last_follow_up,
                                         clinical_brca$days_to_death)





# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_brca_all = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

output_brca <- getResults(query_brca_all)
# get 20 primary tissue sample barcodes
tumor <- output_brca$cases[1:20]
# OR
tumor <- output_brca[output_brca$sample_type == "Primary Tumor", "cases"][1:20]
tumor

# # get gene expression data from 20 primary tumors 
query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_brca)

# get counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")
brca_matrix[1:10,1:10]


# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))


# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]


# Get data for TP53 gene and add gene metadata information to it -------------
brca_tp53 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "TP53")


# get median value
median_value <- median(brca_tp53$counts)

# denote which cases have higher or lower expression than median count
brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")

# Add clinical information to brca_tp53
brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x = 'case_id', by.y = 'submitter_id')


# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
fit
ggsurvplot(fit,
           data = brca_tp53,
           pval = T,
           risk.table = T)


fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)
