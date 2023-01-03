# download gd_holygrail_filtered_final.rds from here: 
multimodal_gd_gex <- readRDS("gd_holygrail_filtered_final.rds")

library(Seurat)
library(ggplot2)
library(tidyverse)

# FIGURE 1C

custom_colors <- list()
colors_dutch <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67','#A3CB38','#1289A7','#D980FA','#B53471','#EE5A24','#009432','#0652DD','#9980FA','#833471','#EA2027','#006266')
set.seed(123)
clust.col = sample(rainbow(10))
custom_colors$discrete <- c(colors_dutch, clust.col)

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(multimodal_gd_gex@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(multimodal_gd_gex@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = multimodal_gd_gex@meta.data$seurat_clusters
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

bind_cols(multimodal_gd_gex@meta.data, as.data.frame(multimodal_gd_gex@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 7,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.8,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank())+
  scale_color_manual(
    name = 'Cluster', values = custom_colors$discrete,
    guide = guide_legend(override.aes = list(size = 2))
  )

# FIGURE 1D

custom_colors <- list()
custom_colors$discrete <- c('#12CBC4','#FFC312','#833471','#1289A7',"grey",'#A3CB38','#D980FA','#ED4C67')
DimPlot(multimodal_gd_gex, group.by = "organ",reduction = "umap",cols = custom_colors$discrete)
