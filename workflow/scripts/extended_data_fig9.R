# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(EBImage)
library(rtracklayer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]

zld_titration_blot_aZld <- "data/immunoblot_raw_images/2022-7-12_titration/Zld-titration_aZld.tif"

grh_titration_blot_aGrh <- "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aGrh_short.tif"
grh_titration_blot_aTub <- "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aTub.tif"

twi_titration_blot_aTwi <- "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTwi.tif"
twi_titration_blot_aTub <- "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTub.tif"

zld_tissue_classes <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv") |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |> 
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |> 
  mutate(class = replace(class, class == "repressed_other", "class VI"))

grh_tissue_classes <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv") |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |> 
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |> 
  mutate(class = replace(class, class == "repressed_other", "class VI"))

twi_tissue_classes <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv") |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |> 
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |> 
  mutate(class = replace(class, class == "repressed_other", "class VI"))


# open graphics device =========================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig9.pdf")
# create blank layout for plot =================================================
pageCreate(width = 18.3, height = 12.5, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
grh_color <- "#F98400"
twi_color <- "#00A08A"
H3K27me3_color <- "gray40"      
        
zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")
twi_heatmap_colors <- brewer.pal(9, "GnBu")


# reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_browser_label <- 1
x_offset_browser <- 1.5

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
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot
blot_image <- readImage("data/immunoblot_raw_images/2022-7-12_titration/Zld-titration_aZld.tif")

# rotate image
blot_image <- blot_image |> 
  rotate(85, bg.col = "white") |> 
  flop()

# crop image
blot_image <- blot_image[977:1565,665:1178]

# adjust brightness and contrast
blot_image <- (blot_image - 0.1)

# get aspect ratio
blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]

# place blot on page
plot_width <- 4

plotRaster(
  blot_image,
  x = ref_x + 1.5,
  y = ref_y + 0.5,
  width = plot_width,
  height =  plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)


# label lanes
plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 0.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 1.9, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "500", params = large_text_params, 
  x = ref_x + 2.8, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "1000", params = large_text_params, 
  x = ref_x + 3.7, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "1500", params = large_text_params, 
  x = ref_x + 4.6, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)



plotText(
  label = "Zld", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
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

# anti-Grh
# read in tiff of western blot
grh_blot_image <- readImage("data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aGrh_short.tif")

# rotate image
grh_blot_image <- grh_blot_image |> 
  rotate(90, bg.col = "white") |> 
  flop()

# crop image
grh_blot_image <- grh_blot_image[1024:1638,833:1110]

# adjust brightness and contrast
grh_blot_image <- (grh_blot_image - 0.3)

# anti-Tub
# read in tiff of western blot
tub_blot_image <- readImage("data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aTub.tif")

# rotate image
tub_blot_image <- tub_blot_image |> 
  rotate(90, bg.col = "white") |> 
  flop()

# crop image
tub_blot_image <- tub_blot_image[1191:1805,940:1043]

# get aspect ratio
grh_blot_dim <- dim(grh_blot_image)
grh_blot_aspect_ratio <- grh_blot_dim[2] / grh_blot_dim[1]

tub_blot_dim <- dim(tub_blot_image)
tub_blot_aspect_ratio <- tub_blot_dim[2] / tub_blot_dim[1]



# check that images were cropped to the same width
if (grh_blot_dim[1] != tub_blot_dim[1]) {stop("Grh and anti-tub blot are different widths")}


# place blot on page
plot_width <- 4

plotRaster(
  grh_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 0.5,
  width = plot_width,
  height =  plot_width * grh_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  tub_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 2.6,
  width = plot_width,
  height =  plot_width * tub_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)


# label lanes
plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 0.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 2, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "25", params = large_text_params, 
  x = ref_x + 2.9, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "100", params = large_text_params, 
  x = ref_x + 3.8, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "400", params = large_text_params, 
  x = ref_x + 4.7, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)



plotText(
  label = "Grh", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 1.1, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "tubulin", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 3, just = c("right","center"), default.units = "cm"
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

# anti-Twi
# read in tiff of western blot
twi_blot_image <- readImage("data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTwi.tif")

# rotate image
twi_blot_image <- twi_blot_image |> 
  rotate(90, bg.col = "white") |> 
  flop()

# crop image
twi_blot_image <- twi_blot_image[466:1103,1043:1215]

# adjust brightness and contrast
twi_blot_image <- (twi_blot_image - 0.2)

# anti-Tub
# read in tiff of western blot
tub_blot_image <- readImage("data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTub.tif")

# rotate image
tub_blot_image <- tub_blot_image |> 
  rotate(87, bg.col = "white") |> 
  flop()

# crop image
tub_blot_image <- tub_blot_image[979:1616,1120:1233]

# get aspect ratio
twi_blot_dim <- dim(twi_blot_image)
twi_blot_aspect_ratio <- twi_blot_dim[2] / twi_blot_dim[1]

tub_blot_dim <- dim(tub_blot_image)
tub_blot_aspect_ratio <- tub_blot_dim[2] / tub_blot_dim[1]



# check that images were cropped to the same width
if (twi_blot_dim[1] != tub_blot_dim[1]) {stop("Twi and anti-tub blot are different widths")}


# place blot on page
plot_width <- 4

plotRaster(
  twi_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 0.5,
  width = plot_width,
  height =  plot_width * twi_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  tub_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 2,
  width = plot_width,
  height =  plot_width * tub_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)


# label lanes
plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 0.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 2.25, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "10", params = large_text_params, 
  x = ref_x + 3.15, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "40", params = large_text_params, 
  x = ref_x + 3.95, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "160", params = large_text_params, 
  x = ref_x + 4.85, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)



plotText(
  label = "HA-Twi", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 1, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "tubulin", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 2.25, just = c("right","center"), default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5



class_binding <- GRangesList(
  zld_tissue_classes = makeGRangesFromDataFrame(zld_tissue_classes),
  `0` = import("ChIPseq/results/peaks/merged/narrow/S2-Zld-0uM_aZld_IP_peaks.narrowPeak"),
  `500` = import("ChIPseq/results/peaks/merged/narrow/S2-Zld-500uM_aZld_IP_peaks.narrowPeak"),
  `1000` = import("ChIPseq/results/peaks/merged/narrow/S2-Zld-1000uM_aZld_IP_peaks.narrowPeak"),
  `1500` = import("ChIPseq/results/peaks/merged/narrow/S2-Zld-1500uM_aZld_IP_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))


d_plot <- zld_tissue_classes |> 
  left_join(class_binding) |> 
  pivot_longer(`0`:`1500`, names_to = "CuSO4_conc", values_to = "has_ChIP_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","500", "1000", "1500"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_bound =  100 *mean(has_ChIP_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_bound, group = class)) +
  geom_bar(stat = "identity", fill = zld_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent bound") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = d_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 5



class_binding <- GRangesList(
  grh_tissue_classes = makeGRangesFromDataFrame(grh_tissue_classes),
  `0` = import("ChIPseq/results/peaks/merged/narrow/S2-Grh-0uM_aGrh_IP_peaks.narrowPeak"),
  `25` = import("ChIPseq/results/peaks/merged/narrow/S2-Grh-25uM_aGrh_IP_peaks.narrowPeak"),
  `100` = import("ChIPseq/results/peaks/merged/narrow/S2-Grh-100uM_aGrh_IP_peaks.narrowPeak"),
  `400` = import("ChIPseq/results/peaks/merged/narrow/S2-Grh-400uM_aGrh_IP_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))

 
e_plot <- grh_tissue_classes |> 
  left_join(class_binding) |> 
  pivot_longer(`0`:`400`, names_to = "CuSO4_conc", values_to = "has_ChIP_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","25", "100", "400"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_bound =  100 *mean(has_ChIP_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_bound, group = class)) +
  geom_bar(stat = "identity", fill = grh_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent bound") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = e_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 5



class_binding <- GRangesList(
  twi_tissue_classes = makeGRangesFromDataFrame(twi_tissue_classes),
  `0` = import("ChIPseq/results/peaks/merged/narrow/S2-HA-Twi-0uM_aHA_IP_peaks.narrowPeak"),
  `10` = import("ChIPseq/results/peaks/merged/narrow/S2-HA-Twi-10uM_aHA_IP_peaks.narrowPeak"),
  `40` = import("ChIPseq/results/peaks/merged/narrow/S2-HA-Twi-40uM_aHA_IP_peaks.narrowPeak"),
  `160` = import("ChIPseq/results/peaks/merged/narrow/S2-HA-Twi-160uM_aHA_IP_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))


f_plot <- twi_tissue_classes |> 
  left_join(class_binding) |> 
  
  pivot_longer(`0`:`160`, names_to = "CuSO4_conc", values_to = "has_ChIP_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","10", "40", "160"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_bound =  100 *mean(has_ChIP_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_bound, group = class)) +
  geom_bar(stat = "identity", fill = twi_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent bound") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = f_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 9


class_accessibility <- GRangesList(
  zld_tissue_classes = makeGRangesFromDataFrame(zld_tissue_classes),
  `0` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Zld_0uM_peaks.narrowPeak"),
  `500` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Zld_500uM_peaks.narrowPeak"),
  `1000` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Zld_1000uM_peaks.narrowPeak"),
  `1500` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Zld_1500uM_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))



g_plot <- zld_tissue_classes |> 
  left_join(class_accessibility) |> 
  pivot_longer(`0`:`1500`, names_to = "CuSO4_conc", values_to = "has_ATAC_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","500", "1000", "1500"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_accessible =  100 *mean(has_ATAC_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_accessible, group = class)) +
  geom_bar(stat = "identity", fill = zld_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent accessible") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = g_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel H ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 9


class_accessibility <- GRangesList(
  grh_tissue_classes = makeGRangesFromDataFrame(grh_tissue_classes),
  `0` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Grh_0uM_peaks.narrowPeak"),
  `25` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Grh_25uM_peaks.narrowPeak"),
  `100` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Grh_100uM_peaks.narrowPeak"),
  `400` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-Grh_400uM_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))


h_plot <- grh_tissue_classes |> 
  left_join(class_accessibility) |> 
  pivot_longer(`0`:`400`, names_to = "CuSO4_conc", values_to = "has_ATAC_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","25", "100", "400"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_accessible =  100 *mean(has_ATAC_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_accessible, group = class)) +
  geom_bar(stat = "identity", fill = grh_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent accessible") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = h_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel I ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 9

class_accessibility <- GRangesList(
  twi_tissue_classes = makeGRangesFromDataFrame(twi_tissue_classes),
  `0` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-HA-Twi_0uM_peaks.narrowPeak"),
  `10` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-HA-Twi_10uM_peaks.narrowPeak"),
  `40` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-HA-Twi_40uM_peaks.narrowPeak"),
  `160` = import("ATACseq/results/peaks/merged_by_sample/titration_S2-HA-Twi_160uM_peaks.narrowPeak")
)  |> 
  peak_overlap_table() |> 
  dplyr::select(-c(width,strand))


i_plot <- twi_tissue_classes |> 
  left_join(class_accessibility) |> 
  pivot_longer(`0`:`160`, names_to = "CuSO4_conc", values_to = "has_ATAC_peak") |> 
  mutate(CuSO4_conc = factor(CuSO4_conc, levels = c("0","10", "40", "160"))) |> 
  group_by(class, CuSO4_conc) |> 
  summarise(percent_accessible =  100 *mean(has_ATAC_peak)) |> 
  ggplot(aes(x = CuSO4_conc, y = percent_accessible, group = class)) +
  geom_bar(stat = "identity", fill = twi_color) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab("percent accessible") +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # geom_line() +
  facet_grid(cols = vars(class))



plotGG(
  plot = i_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5,  height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)
# close graphics device ========================================================
dev.off()
