# download combined_commonpeaks_hm_integrated_only_atac.rds from here:

gd_atac <- readRDS("combined_commonpeaks_hm_integrated_only_atac.rds")

library(tidyverse)
library(scran)
library(patchwork)
library(viridis)
library(ggforce)
library(gghalves)
library(ggridges)
library(scDblFinder)
library(SingleR)
library(wesanderson)
library(pheatmap)

# figure 5a

custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266'
)

set.seed(123)
clust.col = sample(rainbow(5))
custom_colors$discrete <- c(colors_dutch, clust.col)

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gd_atac@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gd_atac@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gd_atac@meta.data$seurat_clusters
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

bind_cols(gd_atac@meta.data, as.data.frame(gd_atac@reductions$umap@cell.embeddings)) %>%
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


# fugure 5b

DimPlot(gd_atac, group.by = 'holygrail.celltype', pt.size = 0.1, cols = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471'))

# figure 5c

da_peaks_all <- readRDS("/data/gruen/sagar2/R/gd_multiome_scRNAseq_scATACseq/da_peaks_all_clusters.rds")
avg.exp <- AverageExpression(gd_atac, group.by = "ident", slot="data", return.seurat = F, features = unique(da_peaks_all$gene))
avg.exp <- as.data.frame(avg.exp$ATAC)
avg.exp[avg.exp > 5] <- 5

genes <- ClosestFeature(gd_atac, regions = rownames(avg.exp))
genes$heatmap_id <- paste(genes$query_region, genes$gene_name, sep = "_")
rownames(genes) <- genes$query_region
heatmap_dataframe <- merge(genes, avg.exp, by = "row.names")
rownames(heatmap_dataframe) <- heatmap_dataframe$heatmap_id
heatmap_dataframe <- heatmap_dataframe[,c(11:27)]
heatmap_dataframe2 <- heatmap_dataframe[, c(1,6,11,3,4,8,5,2,10,13,14,9,12,16,17,15,7)]
pheatmap(heatmap_dataframe2, color = wes_palette("Zissou1", 10, type = "continuous"),show_rownames = F, fontsize_row = 1, cluster_rows = T, cluster_cols = T)

# figure 5d

factor(Idents(gd_atac), levels= c(0,5,10,2,3,7,4,1,9,12,13,8,11,15,16,14,6))
Idents(gd_atac) <- factor(Idents(gd_atac), levels= c(0,5,10,2,3,7,4,1,9,12,13, 8,11,15,16,14,6))
color.atac <- custom_colors$discrete[c(0,5,10,2,3,7,4,1,9,12,13,8,11,15,16,14,6)+1]

p <- CoveragePlot(
  object = gd_atac,
  region = c("Cd24a"),
  extend.upstream = 2000,
  extend.downstream = 2000
)

p & scale_fill_manual(values = color.atac)

# figure 5e

p <- CoveragePlot(
  object = gd_atac,
  region = c("Ifng"),
  extend.upstream = 2000,
  extend.downstream = 2000
)

p & scale_fill_manual(values = color.atac)

# figure 5f

p <- CoveragePlot(
  object = gd_atac,
  region = c("Rorc"),
  extend.upstream = 2000,
  extend.downstream = 2000
)

p & scale_fill_manual(values = color.atac)

