# download gd_with_vdj_final.rds from here: https://drive.google.com/file/d/19mnl5Ytigr98xKFvj1jPqFQcksvOQQcf/view?usp=sharing

gd.vdj <- readRDS("gd_with_vdj_final.rds")

# FIGURE 5A 

custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266'
)

set.seed(123)
clust.col = sample(rainbow(10))

gd.vdj@meta.data <- subset(gd.vdj@meta.data, gd.vdj@meta.data$seurat_clusters %in% c(0:18))
custom_colors$discrete <- c(colors_dutch, clust.col)

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gd.vdj@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gd.vdj@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gd.vdj@meta.data$seurat_clusters
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

p <- bind_cols(gd.vdj@meta.data, as.data.frame(gd.vdj@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.2) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 6,
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
    guide = guide_legend(override.aes = list(size = 1))
  )

p

# FIGURE 5B

DimPlot(gd.vdj, group.by = 'celltypes', pt.size = 0.1, label = F, label.box = F, cols = c('#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471'))

# FIGURE 5C

DimPlot(gd.vdj, group.by = 'organ', pt.size = 0.1, label = F, label.box = F, cols = c('#FFC312','#833471','#1289A7','#A3CB38','#ED4C67'))

