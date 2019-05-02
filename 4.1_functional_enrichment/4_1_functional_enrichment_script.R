### 4.1 GSVA Vignette Functional Enrichment Script of Leukaemia Data ###
mypackages <- c("GSEABase", "GSVA", "Biobase", "genefilter", "limma", "RColorBrewer", "GSVAdata")
lapply(mypackages, library, character.only = T)
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

# Clean up step
filtered_eset <- nsFilter(leukemia_eset, require.entrez = T, remove.dupEntrez = T, var.func = IQR, 
                          var.filter = T, var.cutoff = 0.5, filterByQuantile = T, feature.exclude = "^AFFX")
View(filtered_eset)

# Return all the information for these 4292 features
leukemia_filtered_eset <- filtered_eset$eset
View(leukemia_filtered_eset)

# GSVA
leukemia_es <- gsva(leukemia_filtered_eset, c2BroadSets, min.sz = 10, max.sz = 500, verbose = T)

