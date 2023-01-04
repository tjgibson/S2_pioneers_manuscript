# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)
library(RColorBrewer)


# define input files ===========================================================
Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]
Zld_WT_RNAseq_bw <- snakemake@input[["Zld_WT_RNAseq_bw"]]
Zld_Zld_RNAseq_bw <- snakemake@input[["Zld_Zld_RNAseq_bw"]]

Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
Grh_Grh_ATAC_bw <- snakemake@input[["Grh_Grh_ATAC_bw"]]
Grh_WT_RNAseq_bw <- snakemake@input[["Grh_WT_RNAseq_bw"]]
Grh_Grh_RNAseq_bw <- snakemake@input[["Grh_Grh_RNAseq_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

# # create blank layout for plot ===============================================
pdf(snakemake@output[[1]], useDingbats = FALSE)
pageCreate(width = 18, height = 12.5, default.units = "cm", showGuides = TRUE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors

# model ========================================================================
  PostScriptTrace(
    "manuscript/figures/models/system_schematic.eps",
    outfilename = "manuscript/figures/models/system_schematic.eps.xml",
    charpath = FALSE,
    charpos = FALSE
  )
schematic <- readPicture("manuscript/figures/models/system_schematic.eps.xml")
# grid.picture(schematic)
# schematic_grob <- pictureGrob(schematic, exp = 0)
schematic_gtree <- grid.grabExpr(grid.picture(schematic))

plotGG(schematic_gtree, x = 0.75, y = 0.5, width = 6, height = 2.5, default.units = "cm")



plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = 0.5, y = 0.5, just = "center", default.units = "cm"
)

plotGG(schematic_grob, x = 0.75, y = 0.5, width = 6, height = 2.5, default.units = "cm")


# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)




# close graphics device ========================================================
dev.off()




