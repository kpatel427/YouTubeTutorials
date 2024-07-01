# script to convert hg19 coordinates to hg38
#setwd("~/Desktop/demo/liftOver")

library(rtracklayer)

coordinates <- read.delim('data/testhg19.bed', header = F)
names(coordinates) <- c('chromosome','start','end')

# specify coordinates to liftover
grObject <- GRanges(coordinates)

# import the chain file
chainObject <- import.chain("chain_files/hg19ToHg38.over.chain")

# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))
