### 4.2 GSVA Vignette Molecular Signature Idenitifcation Script of Brain Cell Data ###
mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter", "limma", "RColorBrewer", "GSVAdata")
lapply(mypackages, library, character.only = T)
setwd("L:/Richard B/Analysis/2019/May 2019/testing_the_GSVA_package/Input")

# Load Verhaak et al 2010 data
data(gbm_VerhaakEtAl)
head(featureNames(gbm_eset))
View(gbm_eset)
table(gbm_eset$subtype)

# Load Gene Expression Signatures
data(brainTxDbSets)
View(brainTxDbSets)
sapply(brainTxDbSets, length) # Convert to a vector or a matrix
lapply(brainTxDbSets, head) # Convert to a list

# GSVA
gbm_es <- gsva(gbm_eset, brainTxDbSets, mx.diff = F, verbose = F, parallel.sz = 1)
a <- heatmap(gbm_es@assayData[["exprs"]], Rowv = F, Colv = T)

