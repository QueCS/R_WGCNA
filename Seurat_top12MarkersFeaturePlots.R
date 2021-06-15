# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(limma)
library(xlsx)
library(readxl)

# Load .RData image

# Keep top 12 markers for each cluster
cluster0Top12Markers = head(cluster0Markers, 12)
cluster1Top12Markers = head(cluster1Markers, 12)
cluster2Top12Markers = head(cluster2Markers, 12)
cluster3Top12Markers = head(cluster3Markers, 12)
cluster4Top12Markers = head(cluster4Markers, 12)
cluster5Top12Markers = head(cluster5Markers, 12)
cluster6Top12Markers = head(cluster6Markers, 12)
cluster7Top12Markers = head(cluster7Markers, 12)
cluster8Top12Markers = head(cluster8Markers, 12)
cluster9Top12Markers = head(cluster9Markers, 12)
cluster10Top12Markers = head(cluster10Markers, 12)
# cluster11Top12Markers = head(cluster11Markers, 12)
# cluster12Top12Markers = head(cluster12Markers, 12)
# cluster13Top12Markers = head(cluster13Markers, 12)
# cluster14Top12Markers = head(cluster14Markers, 12)
# cluster15Top12Markers = head(cluster15Markers, 12)
# cluster16Top12Markers = head(cluster16Markers, 12)
# cluster17Top12Markers = head(cluster17Markers, 12)
# cluster18Top12Markers = head(cluster18Markers, 12)
# cluster19Top12Markers = head(cluster19Markers, 12)
# cluster20Top12Markers = head(cluster20Markers, 12)

# FeaturePlot & Export
cluster0FeaturePlot = FeaturePlot(data, features =  cluster0Top12Markers$PeaxiGene, label = TRUE)
cluster0FeaturePlot
png(file = sprintf("%s_cluster0Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster0FeaturePlot
dev.off()
cluster1FeaturePlot = FeaturePlot(data, features =  cluster1Top12Markers$PeaxiGene, label = TRUE)
cluster1FeaturePlot
png(file = sprintf("%s_cluster1Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster1FeaturePlot
dev.off()
cluster2FeaturePlot = FeaturePlot(data, features =  cluster2Top12Markers$PeaxiGene, label = TRUE)
cluster2FeaturePlot
png(file = sprintf("%s_cluster2Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster2FeaturePlot
dev.off()
cluster3FeaturePlot = FeaturePlot(data, features =  cluster3Top12Markers$PeaxiGene, label = TRUE)
cluster3FeaturePlot
png(file = sprintf("%s_cluster3Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster3FeaturePlot
dev.off()
cluster4FeaturePlot = FeaturePlot(data, features =  cluster4Top12Markers$PeaxiGene, label = TRUE)
cluster4FeaturePlot
png(file = sprintf("%s_cluster4Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster4FeaturePlot
dev.off()
cluster5FeaturePlot = FeaturePlot(data, features =  cluster5Top12Markers$PeaxiGene, label = TRUE)
cluster5FeaturePlot
png(file = sprintf("%s_cluster5Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster5FeaturePlot
dev.off()
cluster6FeaturePlot = FeaturePlot(data, features =  cluster6Top12Markers$PeaxiGene, label = TRUE)
cluster6FeaturePlot
png(file = sprintf("%s_cluster6Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster6FeaturePlot
dev.off()
cluster7FeaturePlot = FeaturePlot(data, features =  cluster7Top12Markers$PeaxiGene, label = TRUE)
cluster7FeaturePlot
png(file = sprintf("%s_cluster7Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster7FeaturePlot
dev.off()
cluster8FeaturePlot = FeaturePlot(data, features =  cluster8Top12Markers$PeaxiGene, label = TRUE)
cluster8FeaturePlot
png(file = sprintf("%s_cluster8Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster8FeaturePlot
dev.off()
cluster9FeaturePlot = FeaturePlot(data, features =  cluster9Top12Markers$PeaxiGene, label = TRUE)
cluster9FeaturePlot
png(file = sprintf("%s_cluster9Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster9FeaturePlot
dev.off()
cluster10FeaturePlot = FeaturePlot(data, features =  cluster10Top12Markers$PeaxiGene, label = TRUE)
cluster10FeaturePlot
png(file = sprintf("%s_cluster10Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
cluster10FeaturePlot
dev.off()
# cluster11FeaturePlot = FeaturePlot(data, features =  cluster11Top12Markers$PeaxiGene, label = TRUE)
# cluster11FeaturePlot
# png(file = sprintf("%s_cluster11Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster11FeaturePlot
# dev.off()
# cluster12FeaturePlot = FeaturePlot(data, features =  cluster12Top12Markers$PeaxiGene, label = TRUE)
# cluster12FeaturePlot
# png(file = sprintf("%s_cluster12Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster12FeaturePlot
# dev.off()
# cluster13FeaturePlot = FeaturePlot(data, features =  cluster13Top12Markers$PeaxiGene, label = TRUE)
# cluster13FeaturePlot
# png(file = sprintf("%s_cluster13Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster13FeaturePlot
# dev.off()
# cluster14FeaturePlot = FeaturePlot(data, features =  cluster14Top12Markers$PeaxiGene, label = TRUE)
# cluster14FeaturePlot
# png(file = sprintf("%s_cluster14Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster14FeaturePlot
# dev.off()
# cluster15FeaturePlot = FeaturePlot(data, features =  cluster15Top12Markers$PeaxiGene, label = TRUE)
# cluster15FeaturePlot
# png(file = sprintf("%s_cluster15Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster15FeaturePlot
# dev.off()
# cluster16FeaturePlot = FeaturePlot(data, features =  cluster16Top12Markers$PeaxiGene, label = TRUE)
# cluster16FeaturePlot
# png(file = sprintf("%s_cluster16Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster16FeaturePlot
# dev.off()
# cluster17FeaturePlot = FeaturePlot(data, features =  cluster17Top12Markers$PeaxiGene, label = TRUE)
# cluster17FeaturePlot
# png(file = sprintf("%s_cluster17Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster17FeaturePlot
# dev.off()
# cluster18FeaturePlot = FeaturePlot(data, features =  cluster18Top12Markers$PeaxiGene, label = TRUE)
# cluster18FeaturePlot
# png(file = sprintf("%s_cluster18Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster18FeaturePlot
# dev.off()
# cluster19FeaturePlot = FeaturePlot(data, features =  cluster19Top12Markers$PeaxiGene, label = TRUE)
# cluster19FeaturePlot
# png(file = sprintf("%s_cluster19Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster19FeaturePlot
# dev.off()
# cluster20FeaturePlot = FeaturePlot(data, features =  cluster20Top12Markers$PeaxiGene, label = TRUE)
# cluster20FeaturePlot
# png(file = sprintf("%s_cluster20Top12Markers.png", projectName), width = 1600, height = 1200, units = "px")
# cluster20FeaturePlot
# dev.off()
