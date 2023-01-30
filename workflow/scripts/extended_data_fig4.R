# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(rtracklayer)
library(EBImage)

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
# Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
# Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]
# Zld_WT_RNAseq_bw <- snakemake@input[["Zld_WT_RNAseq_bw"]]
# Zld_Zld_RNAseq_bw <- snakemake@input[["Zld_Zld_RNAseq_bw"]]
# 
# Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
# Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
# Grh_Grh_ATAC_bw <- snakemake@input[["Grh_Grh_ATAC_bw"]]
# Grh_WT_RNAseq_bw <- snakemake@input[["Grh_WT_RNAseq_bw"]]
# Grh_Grh_RNAseq_bw <- snakemake@input[["Grh_Grh_RNAseq_bw"]]
# 
# zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
# grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

# zld_ATAC_results_fn <- "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results_annotated.tsv"
# grh_ATAC_results_fn  <- "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results_annotated.tsv"


twi_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Twi_RNAseq_S2-Twi-vs-S2-WT_results_annotated.tsv"
twi_ATAC_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Twi_ATACseq_S2-Twi-vs-S2-WT-40uM_results.tsv"

twi_ChIP_peaks_fn <- "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP.narrowPeak"


# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig4.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 12, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
twi_color <- "#00A08A"
  

# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot
blot_image <- readImage("data/immunoblot_raw_images/2021-9-13_Twi-vs-embryo/Twi-vs-embryos.tif")

# rotate image
blot_image <- blot_image |> 
  rotate(89, bg.col = "white") |> 
  flop()

# crop image
blot_image <- blot_image[1335:3074,830:1211]

blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 10

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 1.5,
  width = plot_width,
  height = plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
  )

# add labels to western blot
plotSegments(
  x0 = ref_x + 1.1, y0 = ref_y + 0.25, x1 = ref_x + 4.7, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 4.9, y0 = ref_y + 0.25, x1 = ref_x + 8.3, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 8.6, y0 = ref_y + 0.25, x1 = ref_x + 10.6, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

# panel label
plotText(
  label = "S2 Twi", params = large_text_params, fontface = "bold",
  x = ref_x + 2.75, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "S2 HA-Twi", params = large_text_params, fontface = "bold",
  x = ref_x + 6.5, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "3-4H embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 9.5, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 0.75, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "n embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 1.5, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "5", params = large_text_params, 
  x = ref_x + 2.25, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "10", params = large_text_params, 
  x = ref_x + 3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "20", params = large_text_params, 
  x = ref_x + 3.75, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "40", params = large_text_params,
  x = ref_x + 4.5, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 5.25, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "5", params = large_text_params, 
  x = ref_x + 6, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "10", params = large_text_params, 
  x = ref_x + 6.65, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "20", params = large_text_params, 
  x = ref_x + 7.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "40", params = large_text_params, 
  x = ref_x + 8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "5", params = large_text_params, 
  x = ref_x + 8.75, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)


plotText(
  label = "10", params = large_text_params, 
  x = ref_x + 9.6, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "20", params = large_text_params, 
  x = ref_x + 10.45, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)



# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 12
ref_y <- 0.5

# read in ATAC-seq results
twi_ATAC_results <- read_tsv(twi_ATAC_results_fn)

# annotate ATAC peaks as overlapping a ChIP-seq peak
twi_ChIP_peaks <- import(twi_ChIP_peaks_fn)

twi_ATAC_gr <- twi_ATAC_results |> 
  makeGRangesFromDataFrame()

twi_ATAC_results$has_ChIP_peak <- FALSE
overlaps <- findOverlaps(twi_ATAC_gr, twi_ChIP_peaks)@from
twi_ATAC_results$has_ChIP_peak[overlaps] <- TRUE

# define categories for coloring points
twi_ATAC_results <- twi_ATAC_results |> 
  mutate(diff_class = case_when(
    is_diff & has_ChIP_peak ~ "twi bound",
    is_diff & !has_ChIP_peak ~ "not twi bound",
    !is_diff ~ "ns"
  ))

# volcano plot
b_plot <- twi_ATAC_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not twi bound", "twi bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("grey", "black", twi_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# place plots on plotGardener page
plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "ATAC-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel c ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5


# read in RNAseq_results
twi_RNAseq_results <- read_tsv(twi_RNAseq_results_fn)


# define categories for coloring points
twi_RNAseq_results <- twi_RNAseq_results |>
  mutate(
    diff_class = case_when(
      is_diff & has_ChIP_peak ~ "Twi bound",
      is_diff & !has_ChIP_peak ~ "not Twi bound",
      !is_diff ~ "ns"
    )
  )

# volcano plot
c_plot <- twi_RNAseq_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Twi bound", "Twi bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = diff_class
  )) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey", "black", twi_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# place plots on plotGardener page
plotGG(
  plot = c_plot,
  x = (ref_x),
  y = (ref_y),
  width = 4.5,
  height = 4,
  just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "c",
  params = panel_label_params,
  fontface = "bold",
  x = ref_x,
  y = ref_y,
  just = "bottom",
  default.units = "cm"
)

plotText(
  label = "RNA-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)


# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 5

# panel label
plotText(
  label = "d",
  params = panel_label_params,
  fontface = "bold",
  x = ref_x,
  y = ref_y,
  just = "bottom",
  default.units = "cm"
)

test_genes <- twi_RNAseq_results |>
  dplyr::filter(is_diff, log2FoldChange > 0) |>
  pull(gene_id)


ego <- enrichGO(
  gene          = test_genes,
  # universe      = background_genes,
  OrgDb         = org.Dm.eg.db,
  ont           = "BP",
  keyType = "FLYBASE",
  qvalueCutoff = 0.05,
  readable      = TRUE
)


# generate GO enrichment plot
c_plot <- barplot(ego, showCategory = 10) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = c_plot,
  x = (ref_x),
  y = (ref_y),
  width = 10,
  height = 4,
  just = c("left", "top"),
  default.units = "cm"
)

# close graphics device ========================================================
dev.off()




