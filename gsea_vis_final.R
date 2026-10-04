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
library(ggplot2)
setwd('path\\to\\working\\directory')

df = read.csv("ic_tumor_mt.txt", header=TRUE, sep = '\t')
df = df[!is.na(df$logFC), ]
original_gene_list <- df$logFC
names(original_gene_list) <- df$genes
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

gene_ids <- bitr(
  names(gene_list),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_df <- data.frame(
  SYMBOL = names(gene_list),
  logFC = as.numeric(gene_list),
  stringsAsFactors = FALSE
)

gene_df <- merge(
  gene_df,
  gene_ids,
  by = "SYMBOL"
)

gene_df <- gene_df[!duplicated(gene_df$ENTREZID), ]

gene_list <- gene_df$logFC
names(gene_list) <- gene_df$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)
head(gene_list)


gse <- gseGO(
    geneList = gene_list,
    ont = "MF",
    keyType = "ENTREZID",
    nPerm = 10000,
    minGSSize = 5,
    maxGSSize = 300,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = org.Hs.eg.db,
    pAdjustMethod = "none"
)

head(gse)


p <- dotplot(
  gse,
  showCategory = 10,
  split = ".sign",
  color = "pvalue"
) +
  facet_grid(. ~ .sign) +
  theme_bw(base_size = 18) +
  theme(
    axis.text.y = element_text(
      size = 28,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 18,
      face = "bold"
    ),
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      size = 15,
      face = "bold"
    )
  )


