#Author: K.T.Shreya
#Date: 29/04/2024
#Purpose: Estimate differentially expressed ion channels in RNA-Seq fpkm datasets of patients with cancer (limma) post removin batch effects
rm(list = ls())
setwd('path\\to\\working\\directory')
#Importing libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#install.packages('limma')
library(limma)
library("ggplot2")
library("sva")
library ("afffy")
#install.packages("plot3D")
library(ggExtra)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages("volcano3D")
#BiocManager::install("pcplot")

library(mixOmics)

# Read the data
input_file = "tumor_tcga_met500_gtex_combined.txt"
file <- read.table(file = input_file, sep = '\t', header = TRUE,row.names = 1)

# Check dimensions
dim(file)
head(file)

allmisscols <- apply(file,2, function(x)all(is.na(x)));  
colswithallmiss <-names(allmisscols[allmisscols>0]);    
print("the columns with all values missing");    
print(colswithallmiss);

# Remove rows with all NAs
file <- file[rowSums(is.na(file)) != ncol(file), ]

# Ensure no NAs in the data
file <- na.omit(file)
dim(file)

constant_columns <- which(apply(file, 2, var) == 0)
if (length(constant_columns) > 0) {
  print("Constant columns detected:")
  print(constant_columns)
  file <- file[, -constant_columns] # Remove constant columns
} else {
  print("No constant columns detected")}

# Create batch vector
batch <- c(rep("GTEX",279), rep("TCGA", 571), rep("MET500", 1))

pca_res <- prcomp(t(file), scale. = TRUE)
pca_data <- data.frame(pca_res$x)
pca_data$Batch <- batch

ggplot(pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot to Check for Batch Effects",
       x = "Principal Component 1",
       y = "Principal Component 2")



# Check the length of batch and dimensions of file
length(batch)
dim(file)

#Design matrix
sample = factor(c(rep("Normal",338), rep("Tumor",504), rep("Metastatic", 9)))
design.mat = model.matrix(~0+sample)
colnames(design.mat) = levels(sample)
design.mat

corrected_fpkm <- removeBatchEffect(file, batch=batch, design=design.mat)
write.table(corrected_fpkm, file = 'tumor_corrected_fpkm.txt', sep = '\t', row.names=TRUE, col.names=TRUE)

pca_res <- prcomp(t(corrected_fpkm), scale. = TRUE)
pca_data <- data.frame(pca_res$x)
pca_data$Batch <- batch

ggplot(pca_data, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot to Check for Batch Effects",
       x = "Principal Component 1",
       y = "Principal Component 2")

#Contrast matrix
contrast.mat = makeContrasts(TVsN = Tumor-Normal, 
                             MVsN = Metastatic-Normal,
                             MVsT = Metastatic-Tumor, levels = design.mat)
contrast.mat

#Fit Bayes method
fit = lmFit(corrected_fpkm, design.mat)
fit2 = contrasts.fit(fit, contrast.mat)
fit3 = eBayes(fit2)
fit3

deg = topTable(fit3,coef = 'TVsN', number = nrow(file))

#Upregulated and downregulated gene identification
deg$diffexpressed = "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
deg$diffexpressed[deg$adj.P.Val<0.05 &  deg$logFC > 0.6] = "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
deg$diffexpressed[deg$adj.P.Val<0.05 & deg$logFC < -0.6] = "DOWN"


deg <- cbind(rownames(deg), data.frame(deg, row.names=NULL))
deg
colnames(deg)[1] = "genes"
deg

#Export results to csv
write.table(deg, file = 'tumor_tn_differential.txt', sep = '\t', row.names=FALSE, col.names=TRUE)



