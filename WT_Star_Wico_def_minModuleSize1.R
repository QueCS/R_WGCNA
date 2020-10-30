###################################################
# MAKE SURE YOU WORKING DIRECTORY IS SET PROPERLY #
#       RNAseq DATA MUST BE IN THIS DIRECTORY     #
###################################################

# Load mandatory libraries
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
library(xlsx)
library(readxl)

# The following setting is important, do not omit
options(stringsAsFactors = FALSE)

# Import RNASeq data
rawData = read.table("all_star-wico_vs_Petunia_axillaris_v1.6.2_genomeHiCassembly.DESeq2_normalized_counts.txt", comment="", header=TRUE)
dim(rawData)

# Transpose the data (invert rows and columns)
transposedRawData = as.data.frame(t(rawData))
dim(transposedRawData)

# Get rid of genes and samples with too many missing values
# minNSamples = minimal number of samples you want in a gene (here we have 30 samples and we want all genes without all samples present to be discarded)
gsg = goodSamplesGenes(transposedRawData, minNSamples = 30, verbose = 3)
gsg$allOK
cleanData = transposedRawData[gsg$goodSamples, gsg$goodGenes]

# Plot and export (.png) a Cluster Dendrogram showing obvious sample outliers
sampleTree = hclust(dist(cleanData), method = "average")
plot(sampleTree)
png("WT_Star_Wico_def_minModuleSize1_Sample Cluster Dendrogram.png", width = 800, height = 600, units = "px")
plot(sampleTree)
dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(transposedRawData, powerVector = powers, verbose = 5)
cex1 = 0.9

# Plot the Mean connectivity as a function of the SFTP
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Thresholding Power",
ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Export the same plot in .png
png("WT_Star_Wico_def_minModuleSize1_Mean Connectivity.png", width = 800, height = 600, units = "px")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Thresholding Power",
ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Construct the gene network and identify modules
# "maxBlockSize" parameter should be set at the highest your machine can handle considering it's available RAM
# A 4GB workstation should handle up to 8000-10000 probes
# A 8GB workstation should handle up to 12000-14000 probes
# A 16GB workstation should handle up to 20000-22000 probes
# A 32GB workstation should handle up to 30000-32000 probes
# Calculus (256GB WS in 2020) should handle 62000 probes
# "power" parameter should be set using the Mean connectivity plot
net = blockwiseModules(transposedRawData, maxBlockSize = 62000, power = 10,
TOMType = "unsigned", minModuleSize = 1,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE)

# See how many modules where identified
table(net$colors)

# Store valuable data from the TOM
moduleLabels = net$colors
geneTree = net$dendrograms
MEs = net$MEs
# Convert labels to colors for plotting
moduleColors = labels2colors(moduleLabels)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(geneTree[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Export the same plot in .png
png("WT_Star_Wico_def_minModuleSize1_Gene Cluster Dendrogram.png", width = 800, height = 600, units = "px")
plotDendroAndColors(geneTree, moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Plot the eigengene dendrogram
datME = moduleEigengenes(cleanData, moduleLabels)$eigengenes
signif(cor(datME, use = "p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering Tree based on the module eigengene")
# Export the same plot in .png
png("WT_Star_Wico_def_minModuleSize1_Clustering Tree of the Module Eigengene.png", width = 800, height = 600, units = "px")
plot(hclustdatME, main="Clustering Tree based on the module eigengenes")
dev.off()

# Plot for each module the gene expression levels and the eigengen expression levels in the different samples (red > green)
# Store the plots in a .pdf file
pdf("WT_Star_Wico_def_minModuleSize1_Gene and Eigengene expression level.pdf")
which.module="0"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="1"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="2"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="3"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="4"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="5"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="6"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="7"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="8"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="9"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="10"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="11"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="12"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="13"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="14"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="15"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="16"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="17"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="18"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="19"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="20"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="21"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="22"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="23"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="24"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="25"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="26"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="27"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="28"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="29"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="30"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="31"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="32"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="33"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="34"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="35"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="36"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="37"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="38"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="39"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="40"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="41"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="42"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="43"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="44"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="45"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="46"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="47"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="48"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="49"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="50"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="51"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="52"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="53"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="54"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="55"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="56"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="57"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="58"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="59"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="60"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="61"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="62"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="63"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="64"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="65"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="66"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="67"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="68"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="69"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="70"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="71"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="72"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="73"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="74"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="75"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="76"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="77"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="78"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,rlabels=F,rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4.7, 0.2, 1.1))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")
dev.off()

# Store and export the Eigengenes expression levels per Modules
tMEs = t(MEs)
MEsExport = as.data.frame(tMEs)
write.xlsx2(MEsExport, "WT_Star_Wico_def_minModuleSize1_Eigengene expression level.xlsx")
# Store and export the Modules Labels
moduleLabelsExport = as.data.frame(moduleLabels)
moduleLabelsExport = cbind(rownames(moduleLabelsExport), data.frame(moduleLabels, row.names = NULL))
names(moduleLabelsExport)[1] = "PeaxiGene"
names(moduleLabelsExport)[2] = "Module"
write.xlsx2(moduleLabelsExport, row.names = FALSE, "WT_Star_Wico_def_minModuleSize1_Module membership.xlsx")

# Merge Annotation & Module membership
mergedAnnotations <- read_excel("mergedAnnotations.xlsx")
completeData = merge(moduleLabelsExport, mergedAnnotations, by = "PeaxiGene")
write.xlsx2(completeData, row.names = FALSE, "WT_Star_Wico_def_minModuleSize1_Final output.xlsx")

# Store an .RData image of the working space
save.image(file = "WT_Star_Wico_def_minModuleSize1.RData")
