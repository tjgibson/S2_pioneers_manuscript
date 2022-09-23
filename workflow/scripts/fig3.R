# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")
# function to get range for browser tracks -------------------------------------
signal_range <- function(x, extra = 0.05) {
  min <- floor(min(x))
  max <- ceiling(max(x))
  range <- max - min
  margin <- range * extra
  output <- c(min - margin, max + margin)
  return(output)
}

# define input files ===========================================================
class_I_bed_fn <- c(
  zld_class_I = "results/ChIP_peak_classes/zld_class_I.bed",
  grh_class_I = "results/ChIP_peak_classes/grh_class_I.bed",
  twi_class_I = "results/ChIP_peak_classes/twi_class_I.bed"
)

ns_sites_bed_fn <- "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed"

S2_Zld_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
S2_Grh_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
S2_Twi_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
H3K27ac_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw"
Nej_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw"
H3K4me1_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw"
H3K4me3_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw"
H2AV_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw"

## create blank layout for plot =================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pageCreate(width = 12.0, height = 12.5, default.units = "cm", showGuides = TRUE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
twi_color <- "#00A08A"
twi_heatmap_colors <- c("white","#00A08A", "#154734")

zld_color <- "#5BBCD6"
grh_color <- "#F98400"

# reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_brower_label <- 1.75
x_offset_browser <- 2.5

# set genome browser height
gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500

# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate euler plot of overlap between Zld, Grh and Twi class I sites
nonspecific_sites <- class_I_bed_fn %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()

euler_fit <- nonspecific_sites %>% 
  dplyr::select(6:8) %>% 
  as.matrix() %>% 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(zld_color, grh_color, twi_color))


# place euler plot on page
plotGG(
  plot = euler_plot,
  x = (ref_x), y = (ref_y + 0.25),
  width = 2.5, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Twi class I", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = "top", default.units = "cm", fontcolor = twi_color
)

plotText(
  label = paste0("Zld", "\n", "class I"), params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.5), just = "center", default.units = "cm", fontcolor = zld_color
)

plotText(
  label = paste0("Grh", "\n", "class I"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y + 2.5), just = "center", default.units = "cm", fontcolor = grh_color
)

# Panel B ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 3.75 

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in overlapping class I sites
ns_sites <- rtracklayer::import(ns_sites_bed_fn)

# generate heatmap
bw <- c(
  S2_Zld_ChIP =  S2_Zld_ChIP_bw,
  S2_Grh_ChIP =  S2_Grh_ChIP_bw,
  S2_Twi_ChIP =  S2_Twi_ChIP_bw,
  H3K27ac = H3K27ac_bw,
  Nej = Nej_bw,
  H3K4me1 =  H3K4me1_bw,
  H3K4me3 = H3K4me3_bw,
  H2AV = H2AV_bw
)

hm_upstream <- 500
hm_downstream <- 500

regions <- ns_sites 

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = brewer.pal(9, "Greys"), 
  order_by_samples = 1:3, 
  individual_scales = TRUE,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 8, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_1_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = paste(length(regions), "overlapping", "\n", "class I sites"), params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.5), just = c("center"), default.units = "cm", rot = 90
)

plotText(
  label = "Zld", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.70625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Grh", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.71875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Twi", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.73125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K27ac", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.74375), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "CBP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.75625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me1", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.76875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 6.78125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H2AV", params = small_text_params, fontface = "bold",
  x = (ref_x + 7.79375), y = (ref_y), just = c("center"), default.units = "cm"
)

# Panel C ======================================================================
# panel label
ref_x <- 3.5
ref_y <- 0.5 

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# dev.off()
