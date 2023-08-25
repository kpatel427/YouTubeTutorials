#!/bin/bash


pwd

# adjust paths to files based on your current working directory

# Run findMotifsGenome.pl on genomic regions
# To find line number for header line in the peak file
cat data/GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | grep -n ^chr
cat data/GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | tail -n +30 | cut -f 1-3,10 > data/GSE129314_RUNX1ChIP_peaks_homer.txt


# Run Homer
findMotifsGenome.pl data/GSE129314_RUNX1ChIP_peaks_homer.txt hg38 results/ -size 200

# size = size of the region for motif finding
# For recommendations to set this parameter, refer: http://homer.ucsd.edu/homer/ngs/peakMotifs.html



# find motif instances - method 1
findMotifsGenome.pl data/GSE129314_RUNX1ChIP_peaks_homer.txt hg38 results/ -find  results/knownResults/known41.motif > results/motif_instances.txt

# find motif instances - method 2
annotatePeaks.pl data/GSE129314_RUNX1ChIP_peaks_homer.txt hg38 -m results/knownResults/known41.motif > results/motif_instances_annotatePeak_method.txt