# download gd_holygrail_filtered_final.rds from here: https://drive.google.com/file/d/16SeXuyOkRm3y33TqxwIcWewmBY_Nbpfw/view?usp=share_link

multimodal_gd_gex <- readRDS("gd_holygrail_filtered_final.rds")

# figure 2a

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

custom_colors2 <- list()
custom_colors2$discrete <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')


clusters <- levels(multimodal_gd_gex@meta.data$seurat_clusters)
cell_types <- sort(unique(multimodal_gd_gex@meta.data$celltypes))

color_assignments <- setNames(
  c(custom_colors$discrete[1:length(clusters)], custom_colors2$discrete[1:length(cell_types)]),
  c(clusters,cell_types)
)

library(ggforce)
data <- multimodal_gd_gex@meta.data %>%
  dplyr::rename(cell_type = celltypes) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
  group_by(seurat_clusters, cell_type) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = c(clusters,cell_types))
  )

data_labels <- tibble(
  group = c(
    rep('seurat_clusters', length(clusters)),
    rep('cell_type', length(cell_types))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'seurat_clusters', 1, 0),
    nudge_x = ifelse(group == 'seurat_clusters', -0.1, 0.1)
  )

ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = seurat_clusters), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Cluster','Cell type')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

# figure 2b

custom_colors2 <- list()
custom_colors2$discrete <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')
DimPlot(multimodal_gd_gex, group.by = "celltypes",reduction = "umap",cols = custom_colors2$discrete)

# figure 2c

avg.exp <- AverageExpression(multimodal_gd_gex,features = VariableFeatures(multimodal_gd_gex), group.by = "celltypes")
avg.exp <- as.data.frame(avg.exp$RNA)
avg.exp$pooled <- NULL

library(factoextra)
res.pca <- prcomp(t(avg.exp))
fviz_eig(res.pca)


dat <- as.data.frame(res.pca$x[,c(1:2)])
dat$celltype <- rownames(dat)
library(ggrepel)
ggplot(dat, aes(x = PC1, y = PC2)) + geom_point(aes(color = celltype), size=3) +
  geom_text_repel(aes(label = celltype,  color = celltype), size = 5)+
  scale_color_manual(values = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')) +
  theme_bw() + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")


# figure 2d

Idents(multimodal_gd_gex) <- "seurat_clusters"
cluster.order <- c(0,7,11,2,8,18,14,12,4,1,9,16,20,6,10,15,17,3,5,19,13,21)
factor(Idents(multimodal_gd_gex), levels= c(0,7,11,2,8,18,14,12,4,1,9,16,20,6,10,15,17,3,5,19,13,21))
Idents(multimodal_gd_gex) <- factor(Idents(multimodal_gd_gex), levels= c(0,7,11,2,8,18,14,12,4,1,9,16,20,6,10,15,17,3,5,19,13,21))
DotPlot(multimodal_gd_gex, features = c("Sell","S1pr1","Tcf7","Lef1","Cd8b1","Myb","Sox4","Eomes","Ly6c2","Ifit1","Ifit3","Irf7","Itga4","Zbtb20","Prdm16","Klri1","Cd40lg","Ifng","Klrb1c","Cd160","Fgl2","Cd8a","Gzma","Gzmb","Itgae","Cd200r2","Themis","Stat3","Stat4","Kit","Ccr9","Nr4a1","Cd69","Mki67","Zeb2","Klrg1","Tbx21","Rorc","Il17a","Cd163l1","Il23r","Zbtb16","Ctla2a","Gem","Areg","Mest","Odc1","Rora"),cols = c("RdYlBu")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() +theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  text = element_text(size = 10),
  panel.grid.major.x = element_line(color = "grey80"),
  panel.grid.major.y = element_line(color = "grey80") 
  
)


# figure 2e

metadata.organ <- subset(multimodal_gd_gex@meta.data, multimodal_gd_gex@meta.data$organ %in% c("LI","liver","SI","lung","spleen","skin","LN"))
clusters <- unique(metadata.organ$organ)
cell_types <- sort(unique(metadata.organ$celltypes))

color_assignments <- setNames(
  c('#12CBC4','#FFC312','#A3CB38','#1289A7','#ED4C67','#D980FA','#833471', custom_colors2$discrete[1:length(cell_types)]),
  c(clusters,cell_types)
)

library(ggforce)
data <- metadata.organ %>%
  dplyr::rename(cell_type = celltypes) %>%
  dplyr::mutate(cell_type = factor(cell_type, levels = cell_types)) %>%
  group_by(organ, cell_type) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = c(clusters,cell_types))
  )

data_labels <- tibble(
  group = c(
    rep('organ', length(clusters)),
    rep('cell_type', length(cell_types))
  )
) %>%
  mutate(
    hjust = ifelse(group == 'organ', 1, 0),
    nudge_x = ifelse(group == 'organ', -0.1, 0.1)
  )

ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = organ), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +
  geom_text(
    aes(y = n, split = y), stat = 'parallel_sets_axes', fontface = 'bold',
    hjust = data_labels$hjust, nudge_x = data_labels$nudge_x
  ) +
  scale_x_discrete(labels = c('Organ','Cell type')) +
  scale_fill_manual(values = color_assignments) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(face = 'bold', colour = 'black', size = 15),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
