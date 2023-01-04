gd_atac <- readRDS("/data/gruen/sagar2/R/gd_multiome_scRNAseq_scATACseq/combined_commonpeaks_hm_integrated_only_atac.rds")

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)

DefaultAssay(gd_atac) <- "ATAC"

main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- which(as.character(seqnames(granges(gd_atac))) %in% main.chroms)
gd_atac[["ATAC"]] <- subset(gd_atac[["ATAC"]], features = rownames(gd_atac[["ATAC"]])[keep.peaks])

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
gd_atac <- AddMotifs(
  object = gd_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)


# Run chromVAR

gd_atac <- RunChromVAR(
  object = gd_atac,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(gd_atac) <- 'chromvar'

# figure 6a

factor(Idents(gd_atac), levels= c(0,5,10,2,3,7,4,1,9,12,13,8,11,15,16,14,6))
Idents(gd_atac) <- factor(Idents(gd_atac), levels= c(0,5,10,2,3,7,4,1,9,12,13, 8,11,15,16,14,6))

dotplot.motif <-  c("MA0523.1","MA0768.1","MA1421.1","MA1522.1","MA0753.2","MA0769.2","MA1653.1","MA0442.2","MA1578.1","MA0688.1",
                    "MA0689.1","MA0690.1","MA0800.1","MA0801.1","MA0802.1","MA0803.1","MA0805.1","MA0806.1","MA1567.1","MA0807.1",
                    "MA0640.2","MA0473.3","MA0598.3","MA0761.2","MA0062.3","MA0750.2","MA0076.2","MA0764.2","MA0645.1","MA0765.2",
                    "MA0051.1","MA0160.1","MA0517.1","MA0652.1","MA0653.1","MA0017.2","MA0772.1","MA1111.1","MA1418.1","MA1419.1",
                    "MA0625.1","MA0624.1","MA1525.1","MA1420.1","MA1649.1","MA1547.1","MA1573.1","MA0493.1","MA1511.1","MA1515.1",
                    "MA0071.1","MA1150.1","MA1151.1","MA0742.1","MA1516.1","MA1517.1","MA1532.1","MA0857.1","MA0038.2","MA0483.1",
                    "MA0482.2","MA0842.2","MA0037.3","MA0117.2","MA0107.1","MA0101.1","MA0002.2","MA0511.2","MA0684.2","MA1112.2",
                    "MA0476.1","MA0488.1","MA0496.3","MA1633.1","MA1622.1")

DotPlot(gd_atac, features = dotplot.motif,cols = c("RdYlBu"), cluster.idents = F) + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + RotatedAxis() + scale_x_discrete( labels=ConvertMotifID(gd_atac,id = dotplot.motif, assay = "ATAC")) + theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA),
  text = element_text(size = 10),
  panel.grid.major.x = element_line(color = "grey80"),
  panel.grid.major.y = element_line(color = "grey80") 
  
)

# figures 6b-6f

MotifPlot(gd_atac,motifs = c("MA0769.2"),assay = "ATAC")
MotifPlot(gd_atac,motifs = c("MA0690.1"),assay = "ATAC")
MotifPlot(gd_atac,motifs = c("MA0653.1"),assay = "ATAC")
MotifPlot(gd_atac,motifs = c("MA0476.1"),assay = "ATAC")
MotifPlot(gd_atac,motifs = c("MA1151.1"),assay = "ATAC")

FeaturePlot(gd_atac,features = c("MA0769.2","MA0690.1","MA0653.1","MA0476.1","MA1151.1"),min.cutoff = 'q10',max.cutoff = 'q90',pt.size = 0.1, cols = c("grey","red4"))

