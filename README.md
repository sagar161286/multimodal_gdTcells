# multimodal_gdTcells

du Halgouet, A., Bruder, K., Peltokangas, N. et al. Multimodal profiling reveals site-specific adaptation and tissue residency hallmarks of γδ T cells across organs in mice. Nat Immunol (2024). https://doi.org/10.1038/s41590-023-01710-y

This repository contains the Seurat R object files and the scripts used to generate the figures in the above mentioned manuscript published in Nature Immunology. 

Please contact me at sagar@uniklinik-freiburg.de for more details.

Please download Seurat R object for gene expression analysis used in Figures 1 and 2 here:
https://drive.google.com/file/d/16SeXuyOkRm3y33TqxwIcWewmBY_Nbpfw/view?usp=share_link
Afterwards, please run R scripts Figure1.R and Figure2.R to reproduce figures in panels 1 and 2 in the paper.

Please download the Seurat R object for gene expression as well as TotalSeq analysis used in Figure 3 here:
https://drive.google.com/file/d/1WrHVbtPh_hk4YffVOoVXXJrzYkWP0dc-/view?usp=share_link
Afterwards, please run R script Figure3.R to reproduce figures in panel 3 in the paper.

Please download the Seurat R object for chromatin accessibility analysis used in Figure 4 here:
https://drive.google.com/file/d/1LewftPzazsR1Ozy5rp-he2syG18Gr5L1/view?usp=share_link
Afterwards, please run R script Figure4.R to reproduce figures in panel 4 in the paper.

Please download the Seurat R object for TCR repertoire analysis used in Figure 5 here:
https://drive.google.com/file/d/19mnl5Ytigr98xKFvj1jPqFQcksvOQQcf/view?usp=sharing
Afterwards, please run R script Figure5.R to reproduce figures in panel 4 in the paper.

The filtered contig annotations file detailing the TCR repertoire of each cell can be download from here:
https://drive.google.com/file/d/1oMa9cqgvseVPE6Uo5w9oOKZmg4vdINFC/view?usp=sharing 
This file should be used in conjunction with the Figure 5 Seurat R object previously downloaded from https://drive.google.com/file/d/19mnl5Ytigr98xKFvj1jPqFQcksvOQQcf/view?usp=sharing

Please download the milo R object for Milo analysis used in Figure 6 here:
https://drive.google.com/file/d/1rX6hu2aHXnn5WKn25AbzarOxx6W25oKy/view?usp=sharing
Also download the differential abundance testing results here:
https://drive.google.com/file/d/1uWUlp_rZgD1-d_8C3ji_a-FxSUuF29VA/view?usp=sharing
Afterwards, please run R script Figure6.R to reproduce figures in panel 6 in the paper.

ABSTRACT

γδ T cells perform heterogeneous functions in homeostasis and disease across tissues. However, it is unclear whether these roles correspond to distinct γδ subsets or to a homogeneous population of cells exerting context-dependent functions. Here, by cross-organ multimodal single-cell profiling, we reveal that various mouse tissues harbor unique site-adapted γδ subsets. Epidermal and intestinal intraepithelial γδ T cells are transcriptionally homogeneous and exhibit epigenetic hallmarks of functional diversity. Through parabiosis experiments, we uncovered cellular states associated with cytotoxicity, innate-like rapid interferon-γ production and tissue repair functions displaying tissue residency hallmarks. Notably, our observations add nuance to the link between interleukin-17-producing γδ T cells and tissue residency. Moreover, transcriptional programs associated with tissue-resident γδ T cells are analogous to those of CD8+ tissue-resident memory T cells. Altogether, this study provides a multimodal landscape of tissue-adapted γδ T cells, revealing heterogeneity, lineage relationships and their tissue residency program.

Special thanks goes to: https://romanhaa.github.io/projects/scrnaseq_workflow/ for the beautiful codes.

