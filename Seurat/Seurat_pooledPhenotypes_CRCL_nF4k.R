# Choose a project name
projectName = "pooledPhenotypes_CRCL_nF4k"
# Width and height of .png exports (px)
imageWidth = 800
imageHeight = 600

# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(limma)
library(xlsx)
library(readxl)

# Import annotation data
mergedAnnotation = read_xlsx("mergedAnnotation.xlsx")
manualAnnotation = read_xlsx("manualAnnotation.xlsx")

# Import the .h5 dataset
data = Read10X_h5("filtered_feature_bc_matrix.h5")
# Create the Seurat object
data = CreateSeuratObject(counts = data, project = projectName, min.cells = 3, min.features = 200)
data

# Visualize QC metrics as violin plots
QCPlot = VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
QCPlot
# Export
png(file = sprintf("%s_QC.png", projectName), width = imageWidth, height = imageHeight, units = "px")
QCPlot
dev.off()

# Visualize Count-Feature relationships
countFeaturePlot = FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
countFeaturePlot
#Export
png(file = sprintf("%s_Count-Feature.png", projectName), width = imageWidth, height = imageHeight, units = "px")
countFeaturePlot
dev.off()

# Create a subset of the object
# nFeature_RNA > X & nFeature_RNA < Y with X and Y chosen based on QC plots
X = 200
Y = 4000
data = subset(data, subset = nFeature_RNA > X & nFeature_RNA < Y)

# Normalize data
data = NormalizeData(data)

# Keep and show the 2000 highest variable features for downstream analysis
data = FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(data), 10)
plot1 = VariableFeaturePlot(data)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2
# Export
png(file = sprintf("%s_Variable Features_NOID.png", projectName), width = imageWidth, height = imageHeight, units = "px")
plot1
dev.off()
png(file = sprintf("%s_Variable Features_ID.png", projectName), width = imageWidth, height = imageHeight, units = "px")
plot2
dev.off()

# Linear transformation scaling & PCA computing
allGenes = rownames(data)
data = ScaleData(data, features = allGenes)
data = RunPCA(data, features = VariableFeatures(object = data))

# Show how many PC should be used later on
# JackStraw
dataJS = JackStraw(data, num.replicate = 100, dims = 20)
dataJS = ScoreJackStraw(dataJS, dims = 1:20)
JSPlot = JackStrawPlot(dataJS, dims = 1:20)
JSPlot
# Export
png(file = sprintf("%s_JackStraw Plot.png", projectName), width = imageWidth, height = imageHeight, units = "px")
JSPlot
dev.off()
# ElbowPlot
ElbowPlot = ElbowPlot(data, ndims = 20, reduction = "pca")
# Export
png(file = sprintf("%s_ElbowPlot.png", projectName), width = imageWidth, height = imageHeight, units = "px")
ElbowPlot
dev.off()

# Cluster the cells and plot the clusters
# dims (1:20) must be chosen based on JackStraw and Elbow plots
data = FindNeighbors(data, dims = 1:20)
data = FindClusters(data, resolution = 0.5)
data = RunUMAP(data, dims = 1:20)
clusterPlot = DimPlot(data, reduction = "umap")
clusterPlot
png(file = sprintf("%s_Clusters.png", projectName), width = imageWidth, height = imageHeight, units = "px")
clusterPlot
dev.off()

# Find markers for each cluster
cluster0Markers = FindMarkers(data, ident.1 = 0, min.pct = 0.25)
cluster1Markers = FindMarkers(data, ident.1 = 1, min.pct = 0.25)
cluster2Markers = FindMarkers(data, ident.1 = 2, min.pct = 0.25)
cluster3Markers = FindMarkers(data, ident.1 = 3, min.pct = 0.25)
cluster4Markers = FindMarkers(data, ident.1 = 4, min.pct = 0.25)
cluster5Markers = FindMarkers(data, ident.1 = 5, min.pct = 0.25)
cluster6Markers = FindMarkers(data, ident.1 = 6, min.pct = 0.25)
cluster7Markers = FindMarkers(data, ident.1 = 7, min.pct = 0.25)
cluster8Markers = FindMarkers(data, ident.1 = 8, min.pct = 0.25)
cluster9Markers = FindMarkers(data, ident.1 = 9, min.pct = 0.25)
cluster10Markers = FindMarkers(data, ident.1 = 10, min.pct = 0.25)
cluster11Markers = FindMarkers(data, ident.1 = 11, min.pct = 0.25)
# cluster12Markers = FindMarkers(data, ident.1 = 12, min.pct = 0.25)
# cluster13Markers = FindMarkers(data, ident.1 = 13, min.pct = 0.25)
# cluster14Markers = FindMarkers(data, ident.1 = 14, min.pct = 0.25)
# cluster15Markers = FindMarkers(data, ident.1 = 15, min.pct = 0.25)
# cluster16Markers = FindMarkers(data, ident.1 = 16, min.pct = 0.25)
# cluster17Markers = FindMarkers(data, ident.1 = 17, min.pct = 0.25)
# cluster18Markers = FindMarkers(data, ident.1 = 18, min.pct = 0.25)
# cluster19Markers = FindMarkers(data, ident.1 = 19, min.pct = 0.25)
# Move rownames to the first column of the tables
cluster0Markers = cbind(rownames(cluster0Markers), data.frame(cluster0Markers, row.names = NULL))
cluster1Markers = cbind(rownames(cluster1Markers), data.frame(cluster1Markers, row.names = NULL))
cluster2Markers = cbind(rownames(cluster2Markers), data.frame(cluster2Markers, row.names = NULL))
cluster3Markers = cbind(rownames(cluster3Markers), data.frame(cluster3Markers, row.names = NULL))
cluster4Markers = cbind(rownames(cluster4Markers), data.frame(cluster4Markers, row.names = NULL))
cluster5Markers = cbind(rownames(cluster5Markers), data.frame(cluster5Markers, row.names = NULL))
cluster6Markers = cbind(rownames(cluster6Markers), data.frame(cluster6Markers, row.names = NULL))
cluster7Markers = cbind(rownames(cluster7Markers), data.frame(cluster7Markers, row.names = NULL))
cluster8Markers = cbind(rownames(cluster8Markers), data.frame(cluster8Markers, row.names = NULL))
cluster9Markers = cbind(rownames(cluster9Markers), data.frame(cluster9Markers, row.names = NULL))
cluster10Markers = cbind(rownames(cluster10Markers), data.frame(cluster10Markers, row.names = NULL))
cluster11Markers = cbind(rownames(cluster11Markers), data.frame(cluster11Markers, row.names = NULL))
# cluster12Markers = cbind(rownames(cluster12Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster13Markers = cbind(rownames(cluster13Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster14Markers = cbind(rownames(cluster14Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster15Markers = cbind(rownames(cluster15Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster16Markers = cbind(rownames(cluster16Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster17Markers = cbind(rownames(cluster17Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster18Markers = cbind(rownames(cluster18Markers), data.frame(cluster12Markers, row.names = NULL))
# cluster19Markers = cbind(rownames(cluster19Markers), data.frame(cluster12Markers, row.names = NULL))
# Give a name "PeaxiGene" to the first column of the tables
names(cluster0Markers)[1] = "PeaxiGene"
names(cluster1Markers)[1] = "PeaxiGene"
names(cluster2Markers)[1] = "PeaxiGene"
names(cluster3Markers)[1] = "PeaxiGene"
names(cluster4Markers)[1] = "PeaxiGene"
names(cluster5Markers)[1] = "PeaxiGene"
names(cluster6Markers)[1] = "PeaxiGene"
names(cluster7Markers)[1] = "PeaxiGene"
names(cluster8Markers)[1] = "PeaxiGene"
names(cluster9Markers)[1] = "PeaxiGene"
names(cluster10Markers)[1] = "PeaxiGene"
names(cluster11Markers)[1] = "PeaxiGene"
# names(cluster12Markers)[1] = "PeaxiGene"
# names(cluster13Markers)[1] = "PeaxiGene"
# names(cluster14Markers)[1] = "PeaxiGene"
# names(cluster15Markers)[1] = "PeaxiGene"
# names(cluster16Markers)[1] = "PeaxiGene"
# names(cluster17Markers)[1] = "PeaxiGene"
# names(cluster18Markers)[1] = "PeaxiGene"
# names(cluster19Markers)[1] = "PeaxiGene"

# Add annotation data to the cluster tables
mergedAnnotationCluster0Markers = merge(cluster0Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster0Markers = merge(cluster0Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster1Markers = merge(cluster1Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster1Markers = merge(cluster1Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster2Markers = merge(cluster2Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster2Markers = merge(cluster2Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster3Markers = merge(cluster3Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster3Markers = merge(cluster3Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster4Markers = merge(cluster4Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster4Markers = merge(cluster4Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster5Markers = merge(cluster5Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster5Markers = merge(cluster5Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster6Markers = merge(cluster6Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster6Markers = merge(cluster6Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster7Markers = merge(cluster7Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster7Markers = merge(cluster7Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster8Markers = merge(cluster8Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster8Markers = merge(cluster8Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster9Markers = merge(cluster9Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster9Markers = merge(cluster9Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster10Markers = merge(cluster10Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster10Markers = merge(cluster10Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster11Markers = merge(cluster11Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster11Markers = merge(cluster11Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster12Markers = merge(cluster12Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster12Markers = merge(cluster12Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster13Markers = merge(cluster13Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster13Markers = merge(cluster13Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster14Markers = merge(cluster14Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster14Markers = merge(cluster14Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster15Markers = merge(cluster15Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster15Markers = merge(cluster15Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster16Markers = merge(cluster16Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster16Markers = merge(cluster16Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster17Markers = merge(cluster17Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster17Markers = merge(cluster17Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster18Markers = merge(cluster18Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster18Markers = merge(cluster18Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster19Markers = merge(cluster19Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster19Markers = merge(cluster19Markers, manualAnnotation, by = "PeaxiGene")
# Export
write.xlsx2(mergedAnnotationCluster0Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.0")
write.xlsx2(manualAnnotationCluster0Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.0", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster1Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.1")
write.xlsx2(manualAnnotationCluster1Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.1", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster2Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.2", append = TRUE)
write.xlsx2(manualAnnotationCluster2Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.2", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster3Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.3", append = TRUE)
write.xlsx2(manualAnnotationCluster3Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.3", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster4Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.4", append = TRUE)
write.xlsx2(manualAnnotationCluster4Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.4", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster5Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.5", append = TRUE)
write.xlsx2(manualAnnotationCluster5Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.5", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster6Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.6", append = TRUE)
write.xlsx2(manualAnnotationCluster6Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.6", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster7Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.7", append = TRUE)
write.xlsx2(manualAnnotationCluster7Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.7", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster8Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.8", append = TRUE)
write.xlsx2(manualAnnotationCluster8Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.8", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster9Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.9", append = TRUE)
write.xlsx2(manualAnnotationCluster9Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.9", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster10Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.10", append = TRUE)
write.xlsx2(manualAnnotationCluster10Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.10", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster11Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.11", append = TRUE)
write.xlsx2(manualAnnotationCluster11Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.11", append = TRUE)
gc()
# write.xlsx2(mergedAnnotationCluster12Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.12", append = TRUE)
# write.xlsx2(manualAnnotationCluster12Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.12", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster13Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.13", append = TRUE)
# write.xlsx2(manualAnnotationCluster13Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.13", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster14Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.14", append = TRUE)
# write.xlsx2(manualAnnotationCluster14Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.14", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster15Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.15", append = TRUE)
# write.xlsx2(manualAnnotationCluster15Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.15", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster16Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.16", append = TRUE)
# write.xlsx2(manualAnnotationCluster16Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.16", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster17Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.17", append = TRUE)
# write.xlsx2(manualAnnotationCluster17Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.17", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster18Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.18", append = TRUE)
# write.xlsx2(manualAnnotationCluster18Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.18", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster19Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "mrgC.19", append = TRUE)
# write.xlsx2(manualAnnotationCluster19Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), sheetName = "manC.19", append = TRUE)
# gc()

# Find all markers of each cluster
dataMarkers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Find the top 10 markers of each cluster
top10 = dataMarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# Generate an expression heatmap of the top 10 markers of each cluster
top10HeatMap = DoHeatmap(data, features = top10$gene) + NoLegend()
top10HeatMap
# Export
png(file = sprintf("%s_Top 10 markers Heatmap.png", projectName), width = 4000, height = 3000, units = "px")
top10HeatMap
dev.off()
# Generate a dataframe of the top 10 markers of each cluster
top10df = as.data.frame(top10)
names(top10df)[7] = "PeaxiGene"
# Add annotation data to the dataframe
mergedAnnotationTop10Markers = merge(top10df, mergedAnnotation, by = "PeaxiGene")
manualAnnotationTop10Markers = merge(top10df, manualAnnotation, by = "PeaxiGene")
# Export
# mergedAnnotation
write.xlsx2(mergedAnnotationTop10Markers, file = sprintf("%s_Top 10 markers.xlsx", projectName), sheetName = "mergedAnnotation")
gc()
# manualAnnotation
write.xlsx2(manualAnnotationTop10Markers, file = sprintf("%s_Top 10 markers.xlsx", projectName), sheetName = "manualAnnotation", append = TRUE)
gc()

# Store an .RData image of the working space
save.image(file = sprintf("%s_Image.RData", projectName))
