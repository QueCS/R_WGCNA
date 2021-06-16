# Choose a project name
projectName = "pooledPhenotypes_CRCL_SCT_001"
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
library(glmGamPoi)

# Import annotation data
mergedAnnotation = read_xlsx("mergedAnnotation.xlsx")
manualAnnotation = read_xlsx("manualAnnotation.xlsx")

# Import the .h5 dataset
rawData = Read10X_h5("filtered_feature_bc_matrix.h5")

# Create the Seurat object
seuratObjectData = CreateSeuratObject(counts = rawData)
# Perform SCTransform (normalization and variance stabilization by regularized negative binomial regression)
sctData = SCTransform(seuratObjectData, method = "glmGamPoi")
# Perform PCA
pcaData <- RunPCA(sctData, verbose = FALSE)
# Perform UMPA (Uniform Manifold Approximation and Projection dimensional reduction)
umapData = RunUMAP(pcaData, dims = 1:30)
# Find neighbors (computes the k.param nearest neighbors)
nghData = FindNeighbors(umapData, dims = 1:30)
# Find clusters
clstData = FindClusters(nghData)
# Plot the clustering
Clustering = DimPlot(clstData, label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
# Export
pdf(file = sprintf("%s_Clustering.pdf", projectName), width = 14.8, height = 11.1)
Clustering
dev.off()

# Find cluster markers (all) & Format the df
cluster0Markers = FindMarkers(clstData, ident.1 = 0, min.pct = 0.25)
cluster0Markers = cbind(rownames(cluster0Markers), data.frame(cluster0Markers, row.names = NULL))
names(cluster0Markers)[1] = "PeaxiGene"
cluster1Markers = FindMarkers(clstData, ident.1 = 1, min.pct = 0.25)
cluster1Markers = cbind(rownames(cluster1Markers), data.frame(cluster1Markers, row.names = NULL))
names(cluster1Markers)[1] = "PeaxiGene"
cluster2Markers = FindMarkers(clstData, ident.1 = 2, min.pct = 0.25)
cluster2Markers = cbind(rownames(cluster2Markers), data.frame(cluster2Markers, row.names = NULL))
names(cluster2Markers)[1] = "PeaxiGene"
cluster3Markers = FindMarkers(clstData, ident.1 = 3, min.pct = 0.25)
cluster3Markers = cbind(rownames(cluster3Markers), data.frame(cluster3Markers, row.names = NULL))
names(cluster3Markers)[1] = "PeaxiGene"
cluster4Markers = FindMarkers(clstData, ident.1 = 4, min.pct = 0.25)
cluster4Markers = cbind(rownames(cluster4Markers), data.frame(cluster4Markers, row.names = NULL))
names(cluster4Markers)[1] = "PeaxiGene"
cluster5Markers = FindMarkers(clstData, ident.1 = 5, min.pct = 0.25)
cluster5Markers = cbind(rownames(cluster5Markers), data.frame(cluster5Markers, row.names = NULL))
names(cluster5Markers)[1] = "PeaxiGene"
cluster6Markers = FindMarkers(clstData, ident.1 = 6, min.pct = 0.25)
cluster6Markers = cbind(rownames(cluster6Markers), data.frame(cluster6Markers, row.names = NULL))
names(cluster6Markers)[1] = "PeaxiGene"
cluster7Markers = FindMarkers(clstData, ident.1 = 7, min.pct = 0.25)
cluster7Markers = cbind(rownames(cluster7Markers), data.frame(cluster7Markers, row.names = NULL))
names(cluster7Markers)[1] = "PeaxiGene"
cluster8Markers = FindMarkers(clstData, ident.1 = 8, min.pct = 0.25)
cluster8Markers = cbind(rownames(cluster8Markers), data.frame(cluster8Markers, row.names = NULL))
names(cluster8Markers)[1] = "PeaxiGene"
cluster9Markers = FindMarkers(clstData, ident.1 = 9, min.pct = 0.25)
cluster9Markers = cbind(rownames(cluster9Markers), data.frame(cluster9Markers, row.names = NULL))
names(cluster9Markers)[1] = "PeaxiGene"
cluster10Markers = FindMarkers(clstData, ident.1 = 10, min.pct = 0.25)
cluster10Markers = cbind(rownames(cluster10Markers), data.frame(cluster10Markers, row.names = NULL))
names(cluster10Markers)[1] = "PeaxiGene"
cluster11Markers = FindMarkers(clstData, ident.1 = 11, min.pct = 0.25)
cluster11Markers = cbind(rownames(cluster11Markers), data.frame(cluster11Markers, row.names = NULL))
names(cluster11Markers)[1] = "PeaxiGene"
cluster12Markers = FindMarkers(clstData, ident.1 = 12, min.pct = 0.25)
cluster12Markers = cbind(rownames(cluster12Markers), data.frame(cluster12Markers, row.names = NULL))
names(cluster12Markers)[1] = "PeaxiGene"
cluster13Markers = FindMarkers(clstData, ident.1 = 13, min.pct = 0.25)
cluster13Markers = cbind(rownames(cluster13Markers), data.frame(cluster13Markers, row.names = NULL))
names(cluster13Markers)[1] = "PeaxiGene"
cluster14Markers = FindMarkers(clstData, ident.1 = 14, min.pct = 0.25)
cluster14Markers = cbind(rownames(cluster14Markers), data.frame(cluster14Markers, row.names = NULL))
names(cluster14Markers)[1] = "PeaxiGene"
cluster15Markers = FindMarkers(clstData, ident.1 = 15, min.pct = 0.25)
cluster15Markers = cbind(rownames(cluster15Markers), data.frame(cluster15Markers, row.names = NULL))
names(cluster15Markers)[1] = "PeaxiGene"
cluster16Markers = FindMarkers(clstData, ident.1 = 16, min.pct = 0.25)
cluster16Markers = cbind(rownames(cluster16Markers), data.frame(cluster16Markers, row.names = NULL))
names(cluster16Markers)[1] = "PeaxiGene"
# cluster17Markers = FindMarkers(clstData, ident.1 = 17, min.pct = 0.25)
# cluster17Markers = cbind(rownames(cluster17Markers), data.frame(cluster17Markers, row.names = NULL))
# names(cluster17Markers)[1] = "PeaxiGene"
# cluster18Markers = FindMarkers(clstData, ident.1 = 18, min.pct = 0.25)
# cluster18Markers = cbind(rownames(cluster18Markers), data.frame(cluster18Markers, row.names = NULL))
# names(cluster18Markers)[1] = "PeaxiGene"
# cluster19Markers = FindMarkers(clstData, ident.1 = 19, min.pct = 0.25)
# cluster19Markers = cbind(rownames(cluster19Markers), data.frame(cluster19Markers, row.names = NULL))
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
mergedAnnotationCluster12Markers = merge(cluster12Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster12Markers = merge(cluster12Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster13Markers = merge(cluster13Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster13Markers = merge(cluster13Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster14Markers = merge(cluster14Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster14Markers = merge(cluster14Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster15Markers = merge(cluster15Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster15Markers = merge(cluster15Markers, manualAnnotation, by = "PeaxiGene")
mergedAnnotationCluster16Markers = merge(cluster16Markers, mergedAnnotation, by = "PeaxiGene")
manualAnnotationCluster16Markers = merge(cluster16Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster17Markers = merge(cluster17Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster17Markers = merge(cluster17Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster18Markers = merge(cluster18Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster18Markers = merge(cluster18Markers, manualAnnotation, by = "PeaxiGene")
# mergedAnnotationCluster19Markers = merge(cluster19Markers, mergedAnnotation, by = "PeaxiGene")
# manualAnnotationCluster19Markers = merge(cluster19Markers, manualAnnotation, by = "PeaxiGene")
# Export
write.xlsx2(mergedAnnotationCluster0Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.0")
write.xlsx2(manualAnnotationCluster0Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.0", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster1Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.1", append = TRUE)
write.xlsx2(manualAnnotationCluster1Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.1", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster2Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.2", append = TRUE)
write.xlsx2(manualAnnotationCluster2Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.2", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster3Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.3", append = TRUE)
write.xlsx2(manualAnnotationCluster3Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.3", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster4Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.4", append = TRUE)
write.xlsx2(manualAnnotationCluster4Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.4", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster5Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.5", append = TRUE)
write.xlsx2(manualAnnotationCluster5Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.5", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster6Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.6", append = TRUE)
write.xlsx2(manualAnnotationCluster6Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.6", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster7Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.7", append = TRUE)
write.xlsx2(manualAnnotationCluster7Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.7", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster8Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.8", append = TRUE)
write.xlsx2(manualAnnotationCluster8Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.8", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster9Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.9", append = TRUE)
write.xlsx2(manualAnnotationCluster9Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.9", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster10Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.10", append = TRUE)
write.xlsx2(manualAnnotationCluster10Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.10", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster11Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.11", append = TRUE)
write.xlsx2(manualAnnotationCluster11Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.11", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster12Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.12", append = TRUE)
write.xlsx2(manualAnnotationCluster12Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.12", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster13Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.13", append = TRUE)
write.xlsx2(manualAnnotationCluster13Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.13", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster14Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.14", append = TRUE)
write.xlsx2(manualAnnotationCluster14Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.14", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster15Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.15", append = TRUE)
write.xlsx2(manualAnnotationCluster15Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.15", append = TRUE)
gc()
write.xlsx2(mergedAnnotationCluster16Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.16", append = TRUE)
write.xlsx2(manualAnnotationCluster16Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.16", append = TRUE)
gc()
# write.xlsx2(mergedAnnotationCluster17Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.17", append = TRUE)
# write.xlsx2(manualAnnotationCluster17Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.17", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster18Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.18", append = TRUE)
# write.xlsx2(manualAnnotationCluster18Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.18", append = TRUE)
# gc()
# write.xlsx2(mergedAnnotationCluster19Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "mrgC.19", append = TRUE)
# write.xlsx2(manualAnnotationCluster19Markers, file = sprintf("%s_Cluster markers.xlsx", projectName), row.names = FALSE, sheetName = "manC.19", append = TRUE)
# gc()

# Find all markers of each cluster
dataMarkers = FindAllMarkers(clstData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Find the top 12 markers of each cluster
top12 = dataMarkers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_log2FC)
names(top12)[7] = "PeaxiGene"

# Find the top 12 markers for each cluster
cluster0Top12Markers = top12[1:12,]
cluster1Top12Markers = top12[13:24,]
cluster2Top12Markers = top12[25:36,]
cluster3Top12Markers = top12[37:48,]
cluster4Top12Markers = top12[49:60,]
cluster5Top12Markers = top12[61:72,]
cluster6Top12Markers = top12[73:84,]
cluster7Top12Markers = top12[85:96,]
cluster8Top12Markers = top12[97:108,]
cluster9Top12Markers = top12[109:120,]
cluster10Top12Markers = top12[121:132,]
cluster11Top12Markers = top12[133:144,]
cluster12Top12Markers = top12[145:156,]
cluster13Top12Markers = top12[157:168,]
cluster14Top12Markers = top12[169:180,]
cluster15Top12Markers = top12[181:192,]
cluster16Top12Markers = top12[193:210,]
# cluster17Top12Markers = top12[211:222,]
# cluster18Top12Markers = top12[223:234,]
# cluster19Top12Markers = top12[235:246,]
# cluster20Top12Markers = top12[247:258,]

# .pdf export of Feature & Violin Plots
pdf(file = sprintf("%s_Top 12 Markers.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  cluster0Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster0Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster1Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster1Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster2Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster2Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster3Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster3Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster4Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster4Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster5Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster5Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster6Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster6Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster7Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster7Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster8Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster8Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster9Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster9Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster10Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster10Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster11Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster11Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster12Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster12Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster13Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster13Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster14Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster14Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster15Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster15Top12Markers$PeaxiGene)
FeaturePlot(clstData, features =  cluster16Top12Markers$PeaxiGene, label = TRUE)
VlnPlot(clstData, features =  cluster16Top12Markers$PeaxiGene)
# FeaturePlot(clstData, features =  cluster17Top12Markers$PeaxiGene, label = TRUE)
# VlnPlot(clstData, features =  cluster17Top12Markers$PeaxiGene)
# FeaturePlot(clstData, features =  cluster18Top12Markers$PeaxiGene, label = TRUE)
# VlnPlot(clstData, features =  cluster18Top12Markers$PeaxiGene)
# FeaturePlot(clstData, features =  cluster19Top12Markers$PeaxiGene, label = TRUE)
# VlnPlot(clstData, features =  cluster19Top12Markers$PeaxiGene)
# FeaturePlot(clstData, features =  cluster20Top12Markers$PeaxiGene, label = TRUE)
# VlnPlot(clstData, features =  cluster20Top12Markers$PeaxiGene)
dev.off()

# .pdf exports of genes of interest
# AN1
pdf(file = sprintf("%s_AN1.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf00338g00912", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf00338g00912")  + NoLegend()
dev.off()
# AN2
pdf(file = sprintf("%s_AN2.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf00118g00310", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf00118g00310")  + NoLegend()
dev.off()
# DEF
pdf(file = sprintf("%s_DEF.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf01084g00119", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf01084g00119")  + NoLegend()
dev.off()
# TM6
pdf(file = sprintf("%s_TM6.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf00083g01812", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf00083g01812")  + NoLegend()
dev.off()
# GLO1
pdf(file = sprintf("%s_GLO1.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf00922g00026", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf00922g00026")  + NoLegend()
dev.off()
# GLO2
pdf(file = sprintf("%s_GLO2.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  "Peaxi162Scf00591g00074", label = TRUE, pt.size = 1.3, label.size = 5) + NoLegend()
VlnPlot(clstData, features =  "Peaxi162Scf00591g00074")  + NoLegend()
dev.off()

#.pdf exports of genes groups of interest
# Petal identity
pdf(file = sprintf("%s_Petal identity.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf01084g00119", "Peaxi162Scf00922g00026", "Peaxi162Scf00591g00074", "Peaxi162Scf00481g00721", "Peaxi162Scf00020g02338",
                                    "Peaxi162Scf00021g00015", "Peaxi162Scf00017g03246", "Peaxi162Scf00128g01417", "Peaxi162Scf00175g00918", "Peaxi162Scf00365g00212"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf01084g00119", "Peaxi162Scf00922g00026", "Peaxi162Scf00591g00074", "Peaxi162Scf00481g00721", "Peaxi162Scf00020g02338",
                                "Peaxi162Scf00021g00015", "Peaxi162Scf00017g03246", "Peaxi162Scf00128g01417", "Peaxi162Scf00175g00918", "Peaxi162Scf00365g00212")) + NoLegend()
dev.off()
# Conical cells & striations
pdf(file = sprintf("%s_Conical cells & striations.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf00166g01118", "Peaxi162Scf00332g00746", "Peaxi162Scf01269g00017", "Peaxi162Scf00073g01424"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00166g01118", "Peaxi162Scf00332g00746", "Peaxi162Scf01269g00017", "Peaxi162Scf00073g01424")) + NoLegend()
dev.off()
# Epidermis
pdf(file = sprintf("%s_Epidermis.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf00688g00432", "Peaxi162Scf00079g00092", "Peaxi162Scf00950g00058", "Peaxi162Scf00007g01020", "Peaxi162Scf00932g00044",
                                    "Peaxi162Scf00262g00121", "Peaxi162Scf00387g00614", "Peaxi162Scf00021g00922", "Peaxi162Scf00171g00618", "Peaxi162Scf00456g00213"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00688g00432", "Peaxi162Scf00079g00092", "Peaxi162Scf00950g00058", "Peaxi162Scf00007g01020", "Peaxi162Scf00932g00044",
                                "Peaxi162Scf00262g00121", "Peaxi162Scf00387g00614", "Peaxi162Scf00021g00922", "Peaxi162Scf00171g00618", "Peaxi162Scf00456g00213")) + NoLegend()
FeaturePlot(clstData, features =  c("Peaxi162Scf00233g00069", "Peaxi162Scf00915g00122", "Peaxi162Scf00040g00326", "Peaxi162Scf00040g00320", "Peaxi162Scf00643g00645",
                                    "Peaxi162Scf00213g00121", "Peaxi162Scf00015g00043", "Peaxi162Scf00063g02324", "Peaxi162Scf00922g00110"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00233g00069", "Peaxi162Scf00915g00122", "Peaxi162Scf00040g00326", "Peaxi162Scf00040g00320", "Peaxi162Scf00643g00645",
                                "Peaxi162Scf00213g00121", "Peaxi162Scf00015g00043", "Peaxi162Scf00063g02324", "Peaxi162Scf00922g00110")) + NoLegend()
dev.off()
# Adaxial polarity
pdf(file = sprintf("%s_Adaxial polarity.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf00654g00024", "Peaxi162Scf00038g00093", "Peaxi162Scf00267g00059", "Peaxi162Scf00270g00419"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00654g00024", "Peaxi162Scf00038g00093", "Peaxi162Scf00267g00059", "Peaxi162Scf00270g00419")) + NoLegend()
dev.off()
# Abaxial polarity
pdf(file = sprintf("%s_Abaxial polarity.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf00078g00177", "Peaxi162Scf00045g00221", "Peaxi162Scf00185g00095", "Peaxi162Scf00424g00110", "Peaxi162Scf00170g00063", "Peaxi162Scf00008g00498"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00078g00177", "Peaxi162Scf00045g00221", "Peaxi162Scf00185g00095", "Peaxi162Scf00424g00110", "Peaxi162Scf00170g00063", "Peaxi162Scf00008g00498")) + NoLegend()
dev.off()
# Pigmentation
pdf(file = sprintf("%s_Pigmentation.pdf", projectName), width = 14.8, height = 11.1)
FeaturePlot(clstData, features =  c("Peaxi162Scf00338g00912", "Peaxi162Scf00118g00310", "Peaxi162Scf00472g00077", "Peaxi162Scf00349g00057", "Peaxi162Scf00177g00620",
                                    "Peaxi162Scf00366g00630", "Peaxi162Scf00620g00533", "Peaxi162Scf00163g00081", "Peaxi162Scf00487g00064"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00338g00912", "Peaxi162Scf00118g00310", "Peaxi162Scf00472g00077", "Peaxi162Scf00349g00057", "Peaxi162Scf00177g00620",
                                "Peaxi162Scf00366g00630", "Peaxi162Scf00620g00533", "Peaxi162Scf00163g00081", "Peaxi162Scf00487g00064")) + NoLegend()
FeaturePlot(clstData, features =  c("Peaxi162Scf00378g00113", "Peaxi162Scf00518g00430", "Peaxi162Scf00089g00427", "Peaxi162Scf00316g00055", "Peaxi162Scf00713g00038"), label = TRUE) + NoLegend()
VlnPlot(clstData, features =  c("Peaxi162Scf00378g00113", "Peaxi162Scf00518g00430", "Peaxi162Scf00089g00427", "Peaxi162Scf00316g00055", "Peaxi162Scf00713g00038")) + NoLegend()
dev.off()

# Store an .RData image of the working space
save.image(file = sprintf("%s_Image.RData", projectName))
