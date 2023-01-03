multimodal_gd_gex <- readRDS("/data/gruen/sagar2/R/multimodal_gd/combined_analysis/gd_holygrail_filtered_final.rds")

library(Seurat)
library(ggplot2)
library(tidyverse)

custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266'
)

set.seed(123)
clust.col = sample(rainbow(10))

custom_colors$discrete <- c(colors_dutch, clust.col)
DimPlot(multimodal_gd_gex, reduction = "umap",cols = custom_colors$discrete)

