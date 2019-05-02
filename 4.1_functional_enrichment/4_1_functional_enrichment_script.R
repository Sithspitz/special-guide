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

# Load C2BroadSets
## This is a GeneSetCollection object (gmt)
### Can use the getGmt("x") command to do this from gmt formatted text files, like 'test_enrich_Set_gmt_idea.txt' from the 'Input'
data("c2BroadSets")
View(c2BroadSets)


# Load Leukemia data
## This is very different data because is in a ExpressionSet object 
data(leukemia)
View(leukemia_eset)

# View the subtypes
head(pData(leukemia_eset))
table(leukemia_eset$subtype)




