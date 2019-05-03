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

# Statistical Comparison
## Script below the break shows that are 34 differentially activated C2 pathways between MLL and ALL
### Second break script shows the 122 differentially expressed genes with min log fold change of 2 at 0.1 FDR
adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

design <- model.matrix(~ factor(leukemia_es$subtype))
colnames(design) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, design)
fit2 <- eBayes(fit)
allGeneSets <- topTable(fit2, coef = "MLLvsALL", number = Inf)
DEgeneSets <- topTable(fit2, coef = "MLLvsALL", number = Inf, p.value = adjPvalueCutoff, adjust = "BH")
res <- decideTests(fit2, p.value = adjPvalueCutoff)
summary(res)

fit3 <- lmFit(leukemia_filtered_eset, design)
fit4 <- eBayes(fit3)
allGenes <- topTable(fit4, coef = "MLLvsALL", number = Inf)
DEgenes <- topTable(fit4, coef = "MLLvsALL", number = Inf, p.value = adjPvalueCutoff, adjust = "BH",
                    lfc = logFCcutoff)
res2 <- decideTests(fit4, p.value = adjPvalueCutoff, lfc = logFCcutoff)
summary(res2)
