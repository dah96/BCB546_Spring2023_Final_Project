# Some test, TL;DR

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(enrichplot)

# Read Files
df.ZIKV <- fread("GSE207347_A1B1_vs_A2B2_ZIKV_ribodiff_name.txt.gz")
df.ZIKV_DE <- fread("GSE207347_ZIKV_DESeq2_result_name.txt.gz")

df.ZIKV_DE.sub <- subset(df.ZIKV_DE, !is.na(padj))
idx <- df.ZIKV_DE.sub$log2FoldChange < 0 & df.ZIKV_DE.sub$padj < 0.05
df.sub <- df.ZIKV_DE.sub[idx,]


# ZIKV_DE Data

original_gene_list <- df.sub$log2FoldChange
names(original_gene_list) <- df.sub$ID
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
keytypes(org.Hs.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL",
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none") %>% pairwise_termsim()


emapplot(gse)

