###################################################
# MAKE SURE YOU WORKING DIRECTORY IS SET PROPERLY #
#    ALL DATA FILES MUST BE IN THIS DIRECTORY     #
#  EXPORTED FILES WILL END UP IN THIS DIRECTORY   #
###################################################
#     SET THOSE UP BEFORE RUNNING THE SCRIPT      #
###################################################
# Exported files prefix (filePrefix_File.extension)
filePrefix = "publicationDataTer01_gsg_power15_mergeCutHeight030"
# Width and height of .png exports (px)
imageWidth = 1600
imageHeight = 1200
# Set you number of samples, this setting is used to get rid of genes not present in all samples in data preprocessing
# Manually overwrite it when the gsg function is called if you want to keep genes with missing values
sampleNb = 30

# Load mandatory libraries
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
library(xlsx)
library(readxl)

###################################################
# The following setting is important, do not omit #
###################################################
options(stringsAsFactors = FALSE)

# Import RNAseq data
rawData = read.table("publicationData.txt", comment="", header=TRUE)
dim(rawData)

# Get rid of unwanted columns
rawData = rawData[grep(pattern = "IC", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg1.2", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg3_Sep", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg3_Pet", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg1_Bud", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg2_Bud", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg1_Bud", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Stg3_Bud", colnames(rawData), invert = TRUE)]
rawData = rawData[grep(pattern = "Veg.Mrstm", colnames(rawData), invert = TRUE)]

# Get rid of rows containing columns with 0 reads
rawData = rawData[apply(rawData, 1, function(row) all(row != 0)),]

# Transpose the data (invert rows and columns)
transposedRawData = as.data.frame(t(rawData))
dim(transposedRawData)

# Get rid of genes and samples with too many missing values
gsg = goodSamplesGenes(transposedRawData, minNSamples = sampleNb, verbose = 3)
cleanData = transposedRawData[gsg$goodSamples, gsg$goodGenes]

# Plot and export (.png) a Cluster Dendrogram showing obvious sample outliers
sampleTree = hclust(dist(cleanData), method = "average")
plot(sampleTree)
png("Samples Cluster Dendrogram.png", width = imageWidth, height = imageHeight, units = "px")
plot(sampleTree, main = paste("Samples Cluster Dendrogram"))
dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(cleanData, powerVector = powers, verbose = 5)
cex1 = 0.9

# Plot the Mean Connectivity = f(sft)
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Thresholding Power",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Export the same plot in .png
png("Mean Connectivity.png", width = imageWidth, height = imageHeight, units = "px")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Thresholding Power",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Display sof thresholding powers (sft) values to chose the power fitting your dataset the best
# Usually chose the first power that has an SFT.R.sq > 0.90 (abline function in the plot)
sft

# Plot the Scale Independence : signed R² = f(sft)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.90, col="red")
# Export the same plot in .png
png("Scale Independence.png", width = imageWidth, height = imageHeight, units = "px")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale Independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")
abline(h=0.90, col="red")
dev.off()

# Construct the gene network and identify modules
# "maxBlockSize" parameter should be set at the highest your workstation can handle considering it's available RAM
# A 4GB workstation should handle up to 8000-10000 probes
# A 8GB workstation should handle up to 12000-14000 probes
# A 16GB workstation should handle up to 20000-22000 probes
# A 32GB workstation should handle up to 30000-32000 probes
# Calculus (256GB workstation in 2020) should handle 62000 probes
# "power" parameter should be set using sft table, Mean Connectivity and Scale Independence plots
# "mergeCutHeight" parameter should be set specifically to your dataset
# According to Peter Langfelder, 50-100 samples works well with 0.20 to 0.25, for 30 samples I find that 0.30 is fine
# The bigger mergeCutHeight is the less modules you will obtain, you need to balance it in order to have biologically coherent modules
# "networkType" is to chose wether you want positively and negatively corelated genes to be grouped in the same modules (unsigned) or not (signed)
# "minModuleSize" sets how small a module can be
net = blockwiseModules(cleanData, maxBlockSize = 62000, power = 15,
networkType = "signed", TOMType = "signed", minModuleSize = 1,
reassignThreshold = 0, mergeCutHeight = 0.30,
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
png(file = sprintf("%s_Genes Cluster Dendrogram.png", filePrefix), width = imageWidth, height = imageHeight, units = "px")
plotDendroAndColors(geneTree[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Plot the eigengene dendrogram
datME = moduleEigengenes(cleanData, moduleLabels)$eigengenes
signif(cor(datME, use = "p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
par(mfrow=c(1,1))
plot(hclustdatME, main="Module Eigengenes Clustering Tree")
# Export the same plot in .png
png(file = sprintf("%s_Module Eigengenes Clustering Tree.png", filePrefix), width = imageWidth, height = imageHeight, units = "px")
plot(hclustdatME, main="Module Eigengenes Clustering Tree")
dev.off()

# Plot for each module the gene expression levels and the eigengen expression levels in the different samples (red > green)
# Store the plots in a .pdf file
pdf(file = sprintf("%s_Genes and Eigengenes expression levels.pdf", filePrefix))
which.module="0"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="1"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="2"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="3"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="4"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="5"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="6"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="7"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="8"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="9"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="10"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="11"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="12"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="13"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="14"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="15"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="16"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="17"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="18"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="19"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="20"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="21"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="22"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="23"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="24"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="25"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="26"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="27"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="28"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="29"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="30"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="31"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="32"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="33"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="34"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="35"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="36"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="37"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="38"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="39"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="40"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="41"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="42"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="43"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="44"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")

which.module="45"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 7, 2))
plotMat(t(scale(cleanData[,moduleLabels==which.module])), nrgcols=100,clabels=colnames(rawData),rcols=which.module, main=which.module, cex.main=2)
par(mar=c(5, 4, 0.2, 0.5))
barplot(ME, col=0, main="", cex.main=2, ylab="eigengene expr. lvl.",xlab="Samples")
dev.off()

# Store and export the Eigengenes expression levels per Modules
tMEs = t(MEs)
MEsExport = as.data.frame(tMEs)
write.xlsx2(MEsExport, file = sprintf("%s_Eigengene expression level.xlsx", filePrefix))
# Store and export the Modules Labels
moduleLabelsExport = as.data.frame(moduleLabels)
moduleLabelsExport = cbind(rownames(moduleLabelsExport), data.frame(moduleLabels, row.names = NULL))
names(moduleLabelsExport)[1] = "PeaxiGene"
names(moduleLabelsExport)[2] = "Module"
write.xlsx2(moduleLabelsExport, row.names = FALSE, file = sprintf("%s_Module membership.xlsx", filePrefix))

# Merge mergedAnnotation & Module membership
mergedAnnotation = read_xlsx("mergedAnnotation.xlsx")
mergedData1 = merge(moduleLabelsExport, mergedAnnotation, by = "PeaxiGene")
write.xlsx2(mergedData1, row.names = FALSE, file = sprintf("%s_mergedAnnotation Output.xlsx", filePrefix))

# Merge manualAnnotation & Module membership
manualAnnotation = read_xlsx("manualAnnotation.xlsx")
mergedData2 = merge(moduleLabelsExport, manualAnnotation, by = "PeaxiGene")
write.xlsx2(mergedData2, row.names = FALSE, file = sprintf("%s_manualAnnotation Output.xlsx", filePrefix))

# Store an .RData image of the working space
save.image(file = sprintf("%s_Image.RData", filePrefix))
