### 4.1 GSVA Vignette Functional Enrichment Script of Leukaemia Data ###
library(GSEABase)
library(GSVA)
library(Biobase, genefilter, limma)
library(RColorBrewer)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVAdata")
library(GSVAdata)
setwd("L:/Richard B/Analysis/2019/May 2019/testing_the_GSVA_package/Input")

