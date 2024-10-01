#!/bin/bash

# Script to call somatic variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
# This script is for demonstration purposes only


# set paths to directories and files
gatk_path=/Users/kr/Desktop/demo/tools/gatk/gatk
ref=/Users/kr/Desktop/demo/supporting_files/hg38/hg38.fa
project_dir=/Users/kr/Desktop/demo/somatic_mutect2
aligned_reads=$project_dir/aligned
reads=$project_dir/reads
results=$project_dir/results
mutect2_supporting_files=/Users/kr/Desktop/demo/supporting_files/mutect2_supporting_files



if false
then
echo "Download Mutect2 supporting files..."

################################################### Mutect2 files (TO BE DOWNLOADED ONLY ONCE) ##########################################################

# gnomAD
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz ~/Desktop/demo/supporting_files/mutect2_supporting_files
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi ~/Desktop/demo/supporting_files/mutect2_supporting_files

# PoN
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz ~/Desktop/demo/supporting_files/mutect2_supporting_files
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi ~/Desktop/demo/supporting_files/mutect2_supporting_files

# to create your own panel of normals: https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA

# intervals list
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list ~/Desktop/demo/supporting_files/mutect2_supporting_files


# ----------------------------------------------
# STEP 6: Call Somatic Variants - Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2)
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk mutect2 caller"


${gatk_path} Mutect2 -R ${ref} \
     -I ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam \
     -I ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam \
     -tumor HG008-T \
     -normal HG008-N-D \
     --germline-resource ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
     --panel-of-normals ${mutect2_supporting_files}/1000g_pon.hg38.vcf.gz \
     -O ${results}/HG008_somatic_variants_mutect2.vcf.gz \
     --f1r2-tar-gz ${results}/HG008_f1r2.tar.gz \


# ----------------------------------------------
# STEP 7: Estimate cross-sample contamination
# ----------------------------------------------

# GetPileupSummaries
# Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results are used with CalculateContamination.

# Intervals: wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list

echo "STEP 7: Estimate cross-sample contamination"

# tumor
${gatk_path} GetPileupSummaries \
    --java-options '-Xmx50G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/HG008-T_sorted_dedup_bqsr_reads.bam \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/HG008_T_getpileupsummaries.table

# normal
${gatk_path} GetPileupSummaries \
    --java-options '-Xmx50G' --tmp-dir ${project_dir}/tmp/ \
    -I ${aligned_reads}/HG008-N-D_sorted_dedup_bqsr_reads.bam  \
    -V ${mutect2_supporting_files}/af-only-gnomad.hg38.vcf.gz \
    -L ${mutect2_supporting_files}/exome_calling_regions.v1.1.interval_list \
    -O ${results}/HG008-N-D_getpileupsummaries.table



# Calculate contamination
${gatk_path} CalculateContamination \
    -I ${results}/HG008_T_getpileupsummaries.table \
    -matched ${results}/HG008-N-D_getpileupsummaries.table \
    -O ${results}/HG008_pair_calculatecontamination.table 


# ----------------------------------------------
# STEP 8: Estimate read orientation artifacts
# ----------------------------------------------

echo "STEP 8: Estimate read orientation artifacts"

# read orientation
${gatk_path} LearnReadOrientationModel \
    -I ${results}/HG008_f1r2.tar.gz \
    -O ${results}/read-orientation-model.tar.gz


# ----------------------------------------------
# STEP 9: Filter Variants Called By Mutect2
# ----------------------------------------------

echo "STEP 9: Filter Variants"
${gatk_path} FilterMutectCalls \
        -V ${results}/HG008_somatic_variants_mutect2.vcf.gz \
        -R ${ref} \
        --contamination-table ${results}/HG008_pair_calculatecontamination.table \
        --ob-priors ${results}/read-orientation-model.tar.gz \
        -O ${results}/HG008_somatic_variants_filtered_mutect2.vcf


fi

# ----------------------------------------------
# STEP 10: Annotate Variants - Funcotator
# ----------------------------------------------

echo "STEP 10: Annotate Variants"

# Annotate using Funcotator
${gatk_path} Funcotator \
    --variant ${results}/HG008_somatic_variants_filtered_mutect2.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path /Users/kr/Desktop/demo/tools/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.8.hg38.20230908s \
    --output ${results}/HG008_somatic_variants_functotated.vcf \
    --output-file-format VCF



