#Author: K.T.Shreya
#Date: 10/08/2024
#Purpose: Gene set enrichment analysis of differentially expressed ion channels in patients with cancer

rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
#BiocManager::install("GseaVis")
#install.packages("MASS")
#BiocManager::install("ggpp")
library(GseaVis)
#devtools::install_github("junjunlab/GseaVis")
# load test data
#data(geneList, package="DOSE")
library(ggpp)
library(org.Hs.eg.db)
setwd('path\\to\\working\\directory')

df = read.csv("ic_tumor_mt.txt", header=TRUE, sep = '\t')
df = df[!is.na(df$logFC), ]
original_gene_list <- df$logFC
names(original_gene_list) <- df$genes
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# Convert gene symbols to ENTREZ IDs
gene_ids <- bitr(names(original_gene_list), fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

# Filter the gene list to match converted IDs
gene_list <- gene_list[names(gene_list) %in% gene_ids$SYMBOL]
names(gene_list) <- gene_ids$ENTREZID
# check
head(gene_list)


gse <- gseGO(geneList=gene_list, 
             ont ="CC", #Similar for BP and MF
             keyType = "ENTREZID", 
             nPerm = 10000, 
             minGSSize = 5, 
             maxGSSize = 300, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")

head(gse)
#cnetplot(gse)
dotplot(gse, showCategory=10, split=".sign",color = "p.adjust") + facet_grid(.~.sign)

