#!/bin/bash

# Script to filter and annotate variants
# This script is for demonstration purposes only


# directories
ref="/Users/kr/Desktop/demo/supporting_files/hg38/hg38.fa"
results="/Users/kr/Desktop/demo/variant_calling/results"

if false
then

# -------------------
# Filter Variants - GATK4
# -------------------

# Filter SNPs
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# Filter INDELS
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"





# Select Variants that PASS filters
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf


gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf


# to exclude variants that failed genotype filters
cat analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf




# -------------------
# Annotate Variants - GATK4 Funcotator
# -------------------

# Annotate using Funcotator
gatk Funcotator \
	--variant ${results}/analysis-ready-snps-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /Users/kr/Desktop/demo/tools/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF

gatk Funcotator \
	--variant ${results}/analysis-ready-indels-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /Users/kr/Desktop/demo/tools/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7.20200521g \
	--output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF


fi

# Extract fields from a VCF file to a tab-delimited table

gatk VariantsToTable \
	-V ${results}/analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/output_snps.table





