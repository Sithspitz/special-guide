#### RJB TEST GVSA WITH FAKE DATA ####
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") } else { 
    (BiocManager::install("GSVA", version = "3.8")) }
library("GSVA")
library("GSEABase")
library("Biobase", "genefilter", "limma", "RColorBrewer")
setwd("L:/Richard B/Analysis/2019/May 2019/testing_the_GSVA_package/Input")


## This is the non-S2 gene set enriched example ##
not_s2_enriched_1 <- read.csv("test_gene_set.csv", header = T, row.names = 1)
not_s2_enriched_2 <- data.matrix(not_s2_enriched_1) ## has to be a matrix
gmt_gene_set <- getGmt("test_enrich_Set_gmt_idea.txt") ## Using a gmt gene set

not_s2_enriched_output <- gsva(not_s2_enriched_2, gmt_gene_set, method = "gsva")

## This is the S2 gene set enriched example ##
s2_enriched_1 <- read.csv("test_gene_set_s2_enriched.csv", header = T, row.names = 1)
s2_enriched_2 <- data.matrix(s2_enriched_1)

s2_enriched_output <- gsva(s2_enriched_2, gmt_gene_set, method = "gsva")
s2_enriched_output_ssgsea <- gsva(s2_enriched_2, gmt_gene_set, method = "ssgsea") # Using another appropriate method, ssgsea which looks at gene set enrichment scores per sample

## This is the S2 high and S3 low gene set enriched example ##

s2_hi_s3_lo_1 <- read.csv("test_gene_set_s2_enriched_s3_low.csv", header = T, row.names = 1)
s2_hi_s3_lo_1 <- data.matrix(s2_hi_s3_lo_1)

s2_hi_s3_lo_output <- gsva(s2_hi_s3_lo_11, gmt_gene_set, method = "gsva")
s2_hi_s3_lo_output_ssgsea <- gsva(s2_hi_s3_lo_1, gmt_gene_set, method = "ssgsea")

# As for what to do after this at the moment, I am not sure..