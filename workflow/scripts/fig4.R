# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
# library(grImport)
suppressPackageStartupMessages(library(grid))
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# zld_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")
# grh_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")
# twi_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")
# 
# S2_ZLD_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# embryo_Zld_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-nc14_aZld.bw"
# brain_Zld_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/brain_aZld.bw"
# S2_H3K27me3_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE151983_S2_aH3K27me3_IP.bw"
# S2_H3K9me3_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw"
# 
# S2_Grh_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# embryo_Grh_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw"
# wing_disc_Grh_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/wing-disc_aGrh.bw"
# 
# S2_Twi_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# embryo_Twi_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw"

zld_tissue_occupancy_fn <- snakemake@input[["zld_tissue_occupancy"]]
grh_tissue_occupancy_fn <- snakemake@input[["grh_tissue_occupancy"]]
twi_tissue_occupancy_fn <- snakemake@input[["twi_tissue_occupancy"]]

S2_ZLD_ChIP_bw <-   snakemake@input[["S2_ZLD_ChIP_bw"]]
embryo_Zld_ChIP_bw <-  snakemake@input[["embryo_Zld_ChIP_bw"]]
brain_Zld_ChIP_bw <-  snakemake@input[["brain_Zld_ChIP_bw"]]
S2_H3K27me3_bw <-  snakemake@input[["S2_H3K27me3_bw"]]
S2_H3K9me3_bw <-  snakemake@input[["S2_H3K9me3_bw"]]

S2_Grh_ChIP_bw <-   snakemake@input[["S2_Grh_ChIP_bw"]]
embryo_Grh_ChIP_bw <-  snakemake@input[["embryo_Grh_ChIP_bw"]]
wing_disc_Grh_ChIP_bw <-  snakemake@input[["wing_disc_Grh_ChIP_bw"]]

S2_Twi_ChIP_bw <-  snakemake@input[["S2_Twi_ChIP_bw"]]
embryo_Twi_ChIP_bw <-  snakemake@input[["embryo_Twi_ChIP_bw"]]

S2_Zld_ChIP_DMSO_bw <- snakemake@input[["S2_Zld_ChIP_DMSO_bw"]]
S2_Zld_ChIP_taz_bw  <- snakemake@input[["S2_Zld_ChIP_taz_bw"]]
S2_Zld_H3K27me3_DMSO_bw <-  snakemake@input[["S2_Zld_H3K27me3_DMSO_bw"]]
S2_Zld_H3K27me3_taz_bw <-  snakemake@input[["S2_Zld_H3K27me3_taz_bw"]]

S2_Grh_ChIP_DMSO_bw <-   snakemake@input[["S2_Grh_ChIP_DMSO_bw"]]
S2_Grh_ChIP_taz_bw <-  snakemake@input[["S2_Grh_ChIP_taz_bw"]]
S2_Grh_H3K27me3_DMSO_bw <-  snakemake@input[["S2_Grh_H3K27me3_DMSO_bw"]]
S2_Grh_H3K27me3_taz_bw <-  snakemake@input[["S2_Grh_H3K27me3_taz_bw"]]

S2_Twi_ChIP_DMSO_bw <-   snakemake@input[["S2_Twi_ChIP_DMSO_bw"]]
S2_Twi_ChIP_taz_bw <-  snakemake@input[["S2_Twi_ChIP_taz_bw"]]
S2_Twi_H3K27me3_DMSO_bw <-  snakemake@input[["S2_Twi_H3K27me3_DMSO_bw"]]
S2_Twi_H3K27me3_taz_bw <-  snakemake@input[["S2_Twi_H3K27me3_taz_bw"]]

## create blank layout for plot =================================================
# define figure dimensions in cm
fig_width <-  18
fig_height <- 12.5

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
grh_color <- "#F98400"
twi_color <- "#00A08A"

zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")
twi_heatmap_colors <- brewer.pal(9, "GnBu")

# # set genome browser height
# gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  1000
hm_downstream <-  1000

# read in tissue occupancy classes =============================================
zld_tissue_occupancy <- zld_tissue_occupancy_fn %>% 
  read_tsv()

grh_tissue_occupancy <- grh_tissue_occupancy_fn %>% 
  read_tsv()

twi_tissue_occupancy <- twi_tissue_occupancy_fn %>% 
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


bw <- c(
  S2_ZLD_ChIP =  S2_ZLD_ChIP_bw,
  embryo_Zld_ChIP = embryo_Zld_ChIP_bw,
  brain_Zld_ChIP = brain_Zld_ChIP_bw,
  S2_H3K27me3 = S2_H3K27me3_bw,
  S2_H3K9me3 = S2_H3K9me3_bw 
)


regions <- zld_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  individual_scales = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE, 
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.2), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 0.6), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.3), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.58), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.9), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.5), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Zld", "\n", "S2 cells"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.68), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Zld", "\n", "embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 1.7), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Zld", "\n", "NSC"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.71), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K27me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.74), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K9me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.8), y = (ref_y), just = c("center"), default.units = "cm"
)

# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Grh_ChIP =   S2_Grh_ChIP_bw,
  embryo_Grh_ChIP =  embryo_Grh_ChIP_bw,
  wing_disc_Grh_ChIP =  wing_disc_Grh_ChIP_bw,
  S2_H3K27me3 = S2_H3K27me3_bw,
  S2_H3K9me3 = S2_H3K9me3_bw 
)


regions <- grh_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  individual_scales = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_6_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.2), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 0.5), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.35), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.8), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.3), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.7), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Grh", "\n", "S2 cells"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.68), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Grh", "\n", "embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 1.7), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Grh", "\n", "wing disc"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.71), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K27me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.74), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K9me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.8), y = (ref_y), just = c("center"), default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 0.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Twi_ChIP =   S2_Twi_ChIP_bw,
  embryo_Twi_ChIP =  embryo_Twi_ChIP_bw,
  S2_H3K27me3 = S2_H3K27me3_bw,
  S2_H3K9me3 = S2_H3K9me3_bw 
)


regions <- twi_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  individual_scales = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_11_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.4), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.25), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 3.25), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.35), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.55), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.75), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 5), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Twi", "\n", "S2 cells"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.75), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Twi", "\n", "embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.15), y = (ref_y), just = c("center"), default.units = "cm"
)


plotText(
  label = "H3K27me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.5), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K9me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 6.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)



bw <- c(
  S2_Zld_ChIP_DMSO =  S2_Zld_ChIP_DMSO_bw,
  S2_Zld_ChIP_taz = S2_Zld_ChIP_taz_bw,
  S2_Zld_H3K27me3_DMSO = S2_Zld_H3K27me3_DMSO_bw,
  S2_Zld_H3K27me3_taz = S2_Zld_H3K27me3_taz_bw 
)


regions <- zld_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

group <- c(1,1,2,2)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$class, 
  order_by_samples = 1,
  scale_group = group,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE, 
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_15_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.2), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 0.6), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.3), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.58), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.9), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.5), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Zld ChIP", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.75), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Zld ChIP", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)


# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 6.5

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Grh_ChIP_DMSO =  S2_Grh_ChIP_DMSO_bw,
  S2_Grh_ChIP_taz = S2_Grh_ChIP_taz_bw,
  S2_Grh_H3K27me3_DMSO = S2_Grh_H3K27me3_DMSO_bw,
  S2_Grh_H3K27me3_taz = S2_Grh_H3K27me3_taz_bw
)


regions <- grh_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

group <- c(1,1,2,2)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_split = regions$class, 
  order_by_samples = 1,
  scale_group = group,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE, 
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_19_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.2), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 0.5), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.35), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.8), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.3), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.7), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4), just = c("center"), default.units = "cm"
)


plotText(
  label = paste0("Grh ChIP", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.75), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Grh ChIP", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)


# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <-6.5

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Twi_ChIP_DMSO =  S2_Twi_ChIP_DMSO_bw,
  S2_Twi_ChIP_taz = S2_Twi_ChIP_taz_bw,
  S2_Twi_H3K27me3_DMSO = S2_Twi_H3K27me3_DMSO_bw,
  S2_Twi_H3K27me3_taz = S2_Twi_H3K27me3_taz_bw
)


regions <- twi_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

group <- c(1,1,2,2)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_split = regions$class, 
  order_by_samples = 1,
  scale_group = group,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE, 
  return_heatmap_list = TRUE
)


b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_23_heatmap_body_6_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 0.4), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.25), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 3.25), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.35), just = c("center"), default.units = "cm"
)

plotText(
  label = "IV", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.55), just = c("center"), default.units = "cm"
)

plotText(
  label = "V", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 4.75), just = c("center"), default.units = "cm"
)

plotText(
  label = "VI", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 5), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Twi ChIP", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.75), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("Twi ChIP", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "DMSO"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("H3K27me3", "\n", "taz"), params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)




dev.off()
