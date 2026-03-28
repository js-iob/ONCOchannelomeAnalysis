#Author: K.T.Shreya
#Date: 10/12/2024
#Purpose: Estimate co-expressed ion channels from RNA-seq datasets of patients with cancer (WGCNA)
rm(list = ls())
#Import libraries
library(WGCNA)
#install.packages("flashClust")
library(flashClust)
setwd('path\\to\\working\\directory')
#options(stringsAsFactors = T)

#Data pre-processing
dat = read.csv("ic_tumor_corrected_fpkm.txt", sep = '\t', header = TRUE)

dat = t(dat)
write.table(dat, file = 'transposed_ic_tumor_mt_corrected_fpkm.txt', sep = '\t', row.names=TRUE, col.names=FALSE)
dat = read.table("transposed_ic_tumor_mt_corrected_fpkm.txt", sep = '\t', header = TRUE,row.names = 1)
dat

gsg = goodSamplesGenes(dat, verbose = 3)
gsg
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dat)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dat)[!gsg$goodSamples], collapse = ", ")))
  dat = dat[gsg$goodSamples, gsg$goodGenes]
}
dim(dat)

#Soft-thresholding parameter
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, networkType = "signed")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence_HNSC_MT"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity_tumor_MT"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")


#Adjacency matrix
softPower = 12
adjacency = adjacency(dat, power = softPower, type = "signed")
dissTOM = 1-TOMsimilarity(adjacency, TOMType="signed")
geneTree  = flashClust(as.dist(dissTOM), method="average")

#Heirarchical clustering
par(mfrow=c(1,1))
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 4, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicMods
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
colors = table(dynamicColors)
colors

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "hnsc_MT")


par(mfrow=c(1,1))
traitdata = read.csv("tumor_mt_trait.txt", sep = '\t',header = TRUE)
patient = rownames(dat)
traitRows = match(patient, traitdata$Sample_ID)
datTraits = traitdata[traitRows, -1]
#datTraits
names(datTraits)
nGenes = ncol(dat)
nSamples = nrow(dat)
MEs0 = moduleEigengenes(dat, dynamicColors)$eigengenes
MEs1 = orderMEs(MEs0)
MEs1
write.table(MEs1, file = "hnsc_mt_ME.txt", sep="\t")
moduleTraitCor = cor(MEs1, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs1),
               ySymbols = NULL,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1), 
               main = paste("Module- Trait relationships - Tumor"), ylab = "Gene Expression-Based Modules")

par(mfrow = c(1,2))
metastatic = as.data.frame(datTraits$Metastatic)
names(metastatic) = "metastatic"
modNames = substring(names(MEs1), 3)
modNames
geneModuleMembership = as.data.frame(cor(dat, MEs1, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(dat, metastatic, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(metastatic), sep="")
names(GSPvalue) = paste("p.GS.", names(metastatic), sep="")
module = "blue"
column = match(module, modNames)
moduleGenes = dynamicColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Metastatic",
                   main = paste("Tumor\nModule membership vs. gene significance\n"),
                   cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, col = module)
tumor_mn_blue = names(dat)[dynamicColors=="blue"]
write.table(tumor_mn_blue, file = "tumor_mt_blue_mod.txt", sep="\t")



#Correlation with transcription factor

tf_expression_file <- read.table('tumor_mt_blue_tf_exp.txt', sep='\t',header = TRUE, row.names = 1)
tf_expression = t(tf_expression_file)
dim(tf_expression)
module_eigengene <- moduleEigengenes(dat, dynamicColors)$eigengenes
dim(module_eigengene)

blue_module_eigengene <- module_eigengene[, "MEblue"]
blue_module_eigengene_df <- data.frame(row.names = rownames(module_eigengene),
                                       MEblue = blue_module_eigengene)
#correlation_results <- cor(blue_module_eigengene, tf_expression, method = "pearson")
#correlation_results

# Create a data frame to store both correlation values and p-values
correlation_results_df <- data.frame(Correlation = numeric(ncol(tf_expression)), 
                                     P_value = numeric(ncol(tf_expression)), 
                                     row.names = colnames(tf_expression))  # Use colnames as row names for transcription factors

# Loop over each transcription factor to calculate correlation and p-value
for (i in 1:ncol(tf_expression)) {
  correlation_test <- cor.test(blue_module_eigengene, tf_expression[, i], method = "pearson")
  
  # Store the correlation and p-value
  correlation_results_df[i, "Correlation"] <- correlation_test$estimate
  correlation_results_df[i, "P_value"] <- correlation_test$p.value
}

# View the results
correlation_results_df

write.table(correlation_results_df, file = "tumor_mt_blue_tf_corr.txt", sep="\t")
