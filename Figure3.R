# please download gd_holygrail_final_oct20_citeseq_batch.rds here: https://drive.google.com/file/d/1WrHVbtPh_hk4YffVOoVXXJrzYkWP0dc-/view?usp=share_link

gd_citeseq <- readRDS("gd_holygrail_final_oct20_citeseq_batch.rds")

# figure 3b

DimPlot(gd_citeseq, group.by = "celltypes",reduction = "umap",cols = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471'))

# figure 3c
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

# figure 3d
FeaturePlot(gd_citeseq, features = c("Klrb1c"), min.cutoff = "q9",cols = c("grey90","darkgreen"))
FeaturePlot(gd_citeseq, features = c("NK-1.1-CITESEQ"), min.cutoff = "q10",cols = c("grey90","purple4"))

# figure 3e
FeaturePlot(gd_citeseq, features = c("Vg2-CITESEQ"), min.cutoff = "q10",cols = c("grey90","purple4"))

# figure 3f

gd_citeseq_reduction <- gd_citeseq
DefaultAssay(gd_citeseq_reduction) <- "ADT"
gd_citeseq_reduction <- RunPCA(gd_citeseq_reduction, features = rownames(gd_citeseq_reduction)[c(1:4,6,9,10,12,14:26)], reduction.name = "pca_adt", reduction.key = "pca_adt_", verbose = FALSE)
gd_citeseq_reduction <- RunUMAP(gd_citeseq_reduction, dims = 1:20, reduction = "pca_adt", reduction.key = "adtUMAP_", reduction.name = "umap_adt")
DimPlot(gd_citeseq_reduction, reduction = "umap_adt", pt.size = 0.1, group.by = "celltypes", cols = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','grey','#ED4C67','#1289A7','#D980FA','#B53471')) 


