# please download gd_holygrail_final_oct20_citeseq_batch.rds here: https://drive.google.com/file/d/1WrHVbtPh_hk4YffVOoVXXJrzYkWP0dc-/view?usp=share_link

gd_citeseq <- readRDS("gd_holygrail_final_oct20_citeseq_batch.rds")

# FIGURE 3B

DimPlot(gd_citeseq, group.by = "celltypes",reduction = "umap",cols = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471'))

# FIGURE 3C
DotPlot(gd_citeseq, features = c("CD62L-CITESEQ",
                                               "CD5-CITESEQ",
                                               "CD24-CITESEQ",
                                               "CD27-CITESEQ",
                                               "CD122-CITESEQ",
                                               "NK-1.1-CITESEQ",
                                               "CD69-CITESEQ",
                                               "CD8a-CITESEQ",
                                               "KIT-CITESEQ",
                                               "KLRG1-CITESEQ",
                                               "CD44-CITESEQ",
                                               "CD25-CITESEQ",
                                               "ICOS-CITESEQ",
                                               "IL7R-CITESEQ",
                                               "SCA1-CITESEQ",
                                               "NKp46-CITESEQ",
                                               "Vg2-CITESEQ",
                                               "CCR6-CITESEQ",
                                               "CCR2-CITESEQ",
                                               "CCR7-CITESEQ",
                                               "CD4-CITESEQ",
                                               "TCRbeta-CITESEQ",
                                               "CD19-CITESEQ",
                                               "CSF1R-CITESEQ",
                                               "CD11c-CITESEQ",
                                               "GR-1-CITESEQ"),cols = c("RdYlBu")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() +theme(
                                                 panel.background = element_blank(),
                                                 panel.border = element_rect(fill = NA),
                                                 text = element_text(size = 10),
                                                 panel.grid.major.x = element_line(color = "grey80"),
                                                 panel.grid.major.y = element_line(color = "grey80") 
                                                 
                                               )

# FIGURE 3D

gd_citeseq_reduction <- gd_citeseq
DefaultAssay(gd_citeseq_reduction) <- "ADT"
gd_citeseq_reduction <- RunPCA(gd_citeseq_reduction, features = rownames(gd_citeseq_reduction)[c(1:4,6,9,10,12,14:26)], reduction.name = "pca_adt", reduction.key = "pca_adt_", verbose = FALSE)
gd_citeseq_reduction <- RunUMAP(gd_citeseq_reduction, dims = 1:20, reduction = "pca_adt", reduction.key = "adtUMAP_", reduction.name = "umap_adt")
DimPlot(gd_citeseq_reduction, reduction = "umap_adt", pt.size = 0.1, group.by = "celltypes", cols = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')) 

# FIGURE 3E
Idents(multimodal_gd_gex) <- "celltypes"
factor(Idents(multimodal_gd_gex), levels= c("others" ,"Areg+","Rorc+","Proliferating","Klrg1+ effector","Gzmb+ IELs","Cd160+","Sell+Ly6c+","Sell+Ly6c-"))
Idents(multimodal_gd_gex) <- factor(Idents(multimodal_gd_gex), levels= c("others" ,"Areg+","Rorc+","Proliferating","Klrg1+ effector","Gzmb+ IELs","Cd160+","Sell+Ly6c+","Sell+Ly6c-"))

DotPlot(multimodal_gd_gex, features = c("Sell","Ly6c2","Cd160","Gzmb","Klrg1","Mki67","Rorc","Areg"),cols = c("RdYlBu")) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() +theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  text = element_text(size = 10),
  panel.grid.major.x = element_line(color = "grey80"),
  panel.grid.major.y = element_line(color = "grey80") 
  
)

# FIGURE 3F

library(dplyr)
table_organ_by_celltypes <- multimodal_gd_gex@meta.data %>%
  group_by(organ, celltypes) %>%
  summarize(count = n()) %>%
  spread(celltypes, count, fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('organ', 'total_cell_count', everything())) %>%
  arrange(factor(organ, levels = levels(multimodal_gd_gex@meta.data$organ)))

x <- c("LN","spleen","liver","lung","skin","LI","SI","pooled")

table_organ_by_celltypes <- table_organ_by_celltypes %>% mutate(organ =  factor(organ, levels = x)) %>% arrange(organ) 
temp_labels <- multimodal_gd_gex@meta.data %>% group_by(organ) %>% tally()
temp_labels <-  temp_labels %>% mutate(organ =  factor(organ, levels = x)) %>% arrange(organ)

table_organ_by_celltypes %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'organ') %>%
  ggplot(aes(organ, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = organ, y = Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'celltypes', values = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

