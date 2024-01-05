parabiosis.metadata <- as.data.frame(multimodal_gd_gex@meta.data)
parabiosis.metadata <- subset(parabiosis.metadata, parabiosis.metadata$parabiosis %in% c("circulating","resident"))

# FIGURE 6B

DimPlot(multimodal_gd_gex, group.by = "parabiosis",reduction = "umap",cols = c("#00AFBB", "#E7B800"), cells = rownames(parabiosis.metadata))

# FIGURE 6C

DimPlot(multimodal_gd_gex, group.by = "organ",reduction = "umap",cols = c('#12CBC4','#FFC312','#833471','#1289A7','#D980FA','#ED4C67'), cells = rownames(parabiosis.metadata))

# FIGURE 6D

# download gd_parabiosis_milo.rds from here: https://drive.google.com/file/d/1rX6hu2aHXnn5WKn25AbzarOxx6W25oKy/view?usp=sharing

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

gd_parabiosis_milo.rds <- readRDS("gd_parabiosis_milo.rds")
plotNhoodGraphDA(gd_parabiosis_milo.rds, da_results, layout="UMAP",alpha=0.1) 

# FIGURE 6E

# download da_results from here: https://drive.google.com/file/d/1uWUlp_rZgD1-d_8C3ji_a-FxSUuF29VA/view?usp=sharing

da_results <- readRDS("da_results_parabiosis_milo.rds")
da_results <- annotateNhoods(gd_parabiosis_milo.rds, da_results, coldata_col = "organ")
da_results$organ <- ifelse(da_results$organ_fraction < 0.6, "all_others", da_results$organ)
plotDAbeeswarm(da_results, group.by = "organ") 
