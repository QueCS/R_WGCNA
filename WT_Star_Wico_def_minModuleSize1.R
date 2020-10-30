###################################################
# MAKE SURE YOU WORKING DIRECTORY IS SET PROPERLY #
#       RNAseq DATA MUST BE IN THIS DIRECTORY     #
###################################################

#Load mandatory libraries
library(fastcluster)
library(dynamicTreeCut)
library(WGCNA)
library(xlsx)

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
png("WT_Star_Wico_def_minModuleSize1_Samples Cluster Dendrogram.png", width = 800, height = 600, units = "px")
plot(sampleTree)
dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(transposedRawData, powerVector = powers, verbose = 5)
cex1 = 0.9

# Plot the Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
ylab="Mean Connectivity", type="n",
main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# Export the same plot in .png
png("WT_Star_Wico_def_minModuleSize1_Mean Connectivity.png", width = 800, height = 600, units = "px")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
ylab="Mean Connectivity", type="n",
main = paste("Mean Connectivity"))
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



# Store a .RData image of the working space
save.image(file = "WT_Star_Wico_def_minModuleSize1.RData")
