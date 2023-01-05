# script to install R packages
# setwd("~/Desktop/demo/basicR")

library(readxl)

# 1. using CRAN
install.packages('readxl')
library(readxl)


# 2. install from local file
install.packages('/Users/kr/Desktop/demo/basicR/stringr_1.5.0.tar.gz', repos = NULL)
library(stringr)


# 3. install packages from github (public repositories)
# https://github.com/rstudio/shiny
library(shiny)
install.packages('remotes')
library(remotes)
remotes::install_github('rstudio/shiny')
library(shiny)


# 4. install from bioconductor
# we need a package called BiocManager which is available on CRAN
library(Biostrings)
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('Biostrings')

library(Biostrings)

dnastring=DNAString('AATCGAACTGG')
reverseComplement(dnastring)


# where are my packages saved?
.libPaths()


# how to remove packages?
remove.packages('readxl')

