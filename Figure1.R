# download gd_holygrail_filtered_final.rds from here: 
multimodal_gd_gex <- readRDS("gd_holygrail_filtered_final.rds")

library(Seurat)
library(ggplot2)
library(tidyverse)

# FIGURE 1B

library(tidyr)
multimodal_gd_gex@meta.data$organ <- multimodal_gd_gex@meta.data$organ %>% replace_na('pooled')

types <- as.data.frame(multimodal_gd_gex@meta.data)
types2 <- types$organ
counts <- as.data.frame.matrix(table(types2,types$organ))
counts <- as.data.frame(t(colSums(counts)))
counts <- as.data.frame(t(counts))
colnames(counts) <- c("no_of_cells")
counts$organs <- rownames(counts)
counts <- counts[order(counts$no_of_cells, decreasing = T),]
barplot(counts$no_of_cells, main="No of cells", horiz=TRUE, names.arg=rownames(counts),las=2, xlim = c(0,25000),xaxt='n', col = c("grey", '#FFC312','#12CBC4','#833471','#ED4C67','#A3CB38','#1289A7','#D980FA'), border = NA)
axis(side=1, at=seq(0, 25000, by=2500), las=2)

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


# FIGURE 1E

metadata.wo.pooled <- subset(multimodal_gd_gex@meta.data, multimodal_gd_gex@meta.data$organ %in% c("LI","liver","LN","lung","SI","skin","spleen"))

table_clusters_by_organ <- metadata.wo.pooled %>%
  dplyr::rename('cluster' = 'seurat_clusters') %>%
  group_by(cluster, organ) %>%
  summarize(count = n()) %>%
  spread(organ, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  select(c('cluster', 'total_cell_count', everything())) %>%
  arrange(factor(cluster, levels = levels(metadata.wo.pooled$seurat_clusters)))

knitr::kable(table_clusters_by_organ)

temp_labels <- metadata.wo.pooled %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename('cluster' = seurat_clusters)

custom_colors$discrete <- c('#12CBC4','#FFC312','#833471','#1289A7','#A3CB38','#D980FA','#ED4C67')
table_clusters_by_organ %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'cluster') %>%
  mutate(cluster = factor(cluster, levels = levels(multimodal_gd_gex@meta.data$seurat_clusters))) %>%
  ggplot(aes(cluster, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity') +
  scale_fill_manual(name = 'Sample', values = custom_colors$discrete) +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

# FIGURE 1F

avg.exp <- AverageExpression(multimodal_gd_gex,features = VariableFeatures(multimodal_gd_gex), group.by = "organ")
avg.exp <- as.data.frame(avg.exp$RNA)
avg.exp$pooled <- NULL

library(factoextra)
res.pca <- prcomp(t(avg.exp))
fviz_eig(res.pca)

dat <- as.data.frame(res.pca$x[,c(1:2)])
dat$organ <- rownames(dat)
library(ggrepel)
ggplot(dat, aes(x = PC1, y = PC2)) + geom_point(aes(color = organ), size=3) +
  geom_text_repel(aes(label = organ,  color = organ), size = 5)+
  scale_color_manual(values = c('#12CBC4','#FFC312','#833471','#1289A7','#A3CB38','#D980FA','#ED4C67')) +
  theme_bw() + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")

