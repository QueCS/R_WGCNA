# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(limma)
library(xlsx)
library(readxl)

# Load .RData image

# Violin plot for AN1
violinPlot = VlnPlot(data, features = c("Peaxi162Scf00338g00912"))
violinPlot
# Export
png(file = sprintf("%s_AN1.png", projectName), width = 1600, height = 1200, units = "px")
violinPlot
dev.off()

# Violin plot for AN2
violinPlot = VlnPlot(data, features = c("Peaxi162Scf00118g00310"))
violinPlot
# Export
png(file = sprintf("%s_AN2.png", projectName), width = 1600, height = 1200, units = "px")
violinPlot
dev.off()

# Violin plot for PhML1
violinPlot = VlnPlot(data, features = c("Peaxi162Scf00079g00092", "Peaxi162Scf00262g00121", "Peaxi162Scf00950g00058"))
violinPlot
# Export
png(file = sprintf("%s_PhML1.png", projectName), width = 1600, height = 1200, units = "px")
violinPlot
dev.off()

# Violin plot for KANADI
violinPlot = VlnPlot(data, features = c("Peaxi162Scf00078g00177", "Peaxi162Scf00045g00221", "Peaxi162Scf00185g00095", "Peaxi162Scf00424g00110", "Peaxi162Scf00170g00063", "Peaxi162Scf00008g0049"))
violinPlot
# Export
png(file = sprintf("%s_KANADI.png", projectName), width = 1600, height = 1200, units = "px")
violinPlot
dev.off()
