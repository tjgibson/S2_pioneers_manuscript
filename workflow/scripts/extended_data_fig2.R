# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]

S2_Zld_CR_0H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Zld_aZld_0H_accessble_zscore.bw"
S2_Zld_CR_4H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Zld_aZld_4H_accessble_zscore.bw"
S2_Zld_CR_12H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Zld_aZld_12H_accessble_zscore.bw"
S2_Zld_CR_24H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Zld_aZld_24H_accessble_zscore.bw"
S2_Zld_CR_48H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Zld_aZld_48H_accessble_zscore.bw"

S2_Grh_CR_0H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Grh_aGrh_0H_accessble_zscore.bw"
S2_Grh_CR_4H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Grh_aGrh_4H_accessble_zscore.bw"
S2_Grh_CR_12H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Grh_aGrh_12H_accessble_zscore.bw"
S2_Grh_CR_24H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Grh_aGrh_24H_accessble_zscore.bw"
S2_Grh_CR_48H_bw <- "data/2021-03_CR/results/bigwigs/z_score_normalized/S2-Grh_aGrh_48H_accessble_zscore.bw"

zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"

# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/extended_data_fig2.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 18.5, default.units = "cm", showGuides = TRUE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")


# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500


# impor ChIP classes ===========================================================
zld_ChIP_classes <- zld_ChIP_classes_fn |> 
  read_tsv()

grh_ChIP_classes <- grh_ChIP_classes_fn |> 
  read_tsv()

# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# placeholder for western blot
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 8, height = 3, default.units = "cm", just = c("top","left"))

# panel b ======================================================================
# reference points for positioning figure components
ref_x <- 9
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# placeholder for western blot
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 8, height = 3, default.units = "cm", just = c("top","left"))

# Panel C ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 4.5 

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Zld_CR_0H =  S2_Zld_CR_0H_bw,
  S2_Zld_CR_4H = S2_Zld_CR_4H_bw,
  S2_Zld_CR_12H = S2_Zld_CR_12H_bw,
  S2_Zld_CR_24H = S2_Zld_CR_24H_bw,
  S2_Zld_CR_48H = S2_Zld_CR_48H_bw
)


regions <- zld_ChIP_classes %>%
  filter(class != "i") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$class, 
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

h_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = h_hm,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 5.25, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_4_heatmap_body_2_1")
grid.xaxis(at = c(0, 0.5, 1), label = c(paste0("-",hm_upstream / 1000, "KB"), "peak center", paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.375), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 Grh ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class II", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 1.5), just = c("right", "center"), default.units = "cm"
)

plotText(
  label = "class III", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 3), just = c("right","center"), default.units = "cm"
)


# close graphics device ========================================================
dev.off()




