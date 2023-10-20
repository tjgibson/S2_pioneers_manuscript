# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(plyranges)
library(rtracklayer)
library(EBImage)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# define input files explicitly for interactively testing script
# zld_DBD_blot <- "data/immunoblot_raw_images/2021-10-14_Zld-DBD/Zld-Grh-DBDs.tif"
# 
# grh_DBD_blot <- "data/immunoblot_raw_images/2021-08-19_Grh-DBD/Grh-DBD.tif"
# 
# zld_FL_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# zld_DBD_peaks_fn <- "ChIPseq/results/peaks/filtered/S2-Zld-DBD_aZld_IP.narrowPeak"
# 
# grh_FL_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# grh_DBD_peaks_fn <- "ChIPseq/results/peaks/filtered/S2-Grh-DBD_aGrh_IP.narrowPeak"
# 
# zld_FL_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# zld_DBD_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-DBD_aZld_IP.bw"
# 
# grh_FL_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# grh_DBD_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-DBD_aGrh_IP.bw"
# 
# `S2-Zld-FL_aZld_GFP` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_GFP_RGB_470X 520M (GFP).tif"
# `S2-Zld-FL_aZld_DAPI` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-Zld-FL_aZld_BF` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_BF_RGB_BF.tif"
# 
# `S2-Zld-DBD_aZld_GFP` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_GFP_RGB_470X 520M (GFP).tif"
# `S2-Zld-DBD_aZld_DAPI` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-Zld-DBD_aZld_BF` <- "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_BF_RGB_BF.tif"
# 
# `S2-WT_aZld_GFP` <- "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_GFP_RGB_470X 520M (GFP).tif"
# `S2-WT_aZld_DAPI` <- "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-WT_aZld_BF` <- "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_BF_RGB_BF.tif"
# 
# `S2-Grh-FL_aGrh_GFP` <- "data/DBD_immunostaining_raw_images/20211027 Grh-FL_GFP_RGB_470X 520M (GFP).tif"
# `S2-Grh-FL_aGrh_DAPI` <- "data/DBD_immunostaining_raw_images/20211027 Grh-FL_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-Grh-FL_aGrh_BF` <- "data/DBD_immunostaining_raw_images/20211027 Grh-FL_BF_RGB_BF.tif"
# 
# `S2-Grh-DBD_aGrh_GFP` <- "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_GFP_RGB_470X 520M (GFP).tif"
# `S2-Grh-DBD_aGrh_DAPI` <- "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-Grh-DBD_aGrh_BF` <- "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_BF_RGB_BF.tif"
# 
# `S2-WT_aGrh_GFP` <- "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_GFP_RGB_470X 520M (GFP).tif"
# `S2-WT_aGrh_DAPI` <-  "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_DAPI_RGB_395X 445M (DAPI).tif"
# `S2-WT_aGrh_BF` <- "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_BF_RGB_BF.tif"


# get input files from snakemake
zld_DBD_blot <- snakemake@input[["zld_DBD_blot"]]

grh_DBD_blot <- snakemake@input[["grh_DBD_blot"]]

zld_FL_ChIP_classes_fn <- snakemake@input[["zld_FL_ChIP_classes_fn"]]
zld_DBD_peaks_fn <- snakemake@input[["zld_DBD_peaks_fn"]]

grh_FL_ChIP_classes_fn <- snakemake@input[["grh_FL_ChIP_classes_fn"]]
grh_DBD_peaks_fn <- snakemake@input[["grh_DBD_peaks_fn"]]

zld_FL_bw <-  snakemake@input[["zld_FL_bw"]]
zld_DBD_bw <-  snakemake@input[["zld_DBD_bw"]]

grh_FL_bw <- snakemake@input[["grh_FL_bw"]]
grh_DBD_bw <- snakemake@input[["grh_DBD_bw"]]

`S2-Zld-FL_aZld_GFP` <- snakemake@input[["S2_Zld_FL_aZld_GFP"]]
`S2-Zld-FL_aZld_DAPI` <- snakemake@input[["S2_Zld_FL_aZld_DAPI"]]
`S2-Zld-FL_aZld_BF` <- snakemake@input[["S2_Zld_FL_aZld_BF"]]

`S2-Zld-DBD_aZld_GFP` <- snakemake@input[["S2_Zld_DBD_aZld_GFP"]]
`S2-Zld-DBD_aZld_DAPI` <- snakemake@input[["S2_Zld_DBD_aZld_DAPI"]]
`S2-Zld-DBD_aZld_BF` <- snakemake@input[["S2_Zld_DBD_aZld_BF"]]

`S2-WT_aZld_GFP` <- snakemake@input[["S2_WT_aZld_GFP"]]
`S2-WT_aZld_DAPI` <- snakemake@input[["S2_WT_aZld_DAPI"]]
`S2-WT_aZld_BF` <- snakemake@input[["S2_WT_aZld_BF"]]

`S2-Grh-FL_aGrh_GFP` <- snakemake@input[["S2_Grh_FL_aGrh_GFP"]]
`S2-Grh-FL_aGrh_DAPI` <- snakemake@input[["S2_Grh_FL_aGrh_DAPI"]]
`S2-Grh-FL_aGrh_BF` <- snakemake@input[["S2_Grh_FL_aGrh_BF"]]

`S2-Grh-DBD_aGrh_GFP` <- snakemake@input[["S2_Grh_DBD_aGrh_GFP"]]
`S2-Grh-DBD_aGrh_DAPI` <- snakemake@input[["S2_Grh_DBD_aGrh_DAPI"]]
`S2-Grh-DBD_aGrh_BF` <- snakemake@input[["S2_Grh_DBD_aGrh_BF"]]

`S2-WT_aGrh_GFP` <- snakemake@input[["S2_WT_aGrh_GFP"]]
`S2-WT_aGrh_DAPI` <-  snakemake@input[["S2_WT_aGrh_DAPI"]]
`S2-WT_aGrh_BF` <- snakemake@input[["S2_WT_aGrh_BF"]]

# impor ChIP classes ===========================================================
zld_ChIP_classes <- zld_FL_ChIP_classes_fn |> 
  read_tsv()

grh_ChIP_classes <- grh_FL_ChIP_classes_fn |> 
  read_tsv()

# create blank layout for plot =================================================
# define figure dimensions in cm
fig_width <-  18
fig_height <- 17

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig10.pdf", useDingbats = FALSE, height = 8.5)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

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
blot_image <- readImage(zld_DBD_blot)

# rotate image
blot_image <- blot_image |> 
  rotate(88, bg.col = "white") |> 
  flop()

# crop image
blot_image <- blot_image[150:1356,237:1271]

# adjust brightness
blot_image <- blot_image - 0.3

# get blot aspect ratio
blot_dim <- dim(blot_image)
IF_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 1,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# add labels to western blot
plotSegments(
  x0 = ref_x + 1.5, y0 = ref_y + 0.25, x1 = ref_x + 5.8, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 6, y0 = ref_y + 0.25, x1 = ref_x + 7.5, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotText(
  label = "S2 Zld DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 3.5, y = ref_y, just = "center", default.units = "cm"
)

plotText(
  label = "S2 Zld FL", params = large_text_params, fontface = "bold",
  x = ref_x + 6.75, y = ref_y, just = "center", default.units = "cm"
)


plotText(
  label = "FL", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 1.35, just = "right", default.units = "cm"
)

plotText(
  label = "DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 6.2, just = "right", default.units = "cm"
)

plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 0.75, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 1.85, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "50", params = large_text_params, 
  x = ref_x + 2.6, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "100", params = large_text_params, fontcolor = "firebrick3",
  x = ref_x + 3.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "200", params = large_text_params, 
  x = ref_x + 4.1, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "400", params = large_text_params, 
  x = ref_x + 4.85, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "800", params = large_text_params, 
  x = ref_x + 5.5, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 6.25, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1000", params = large_text_params, 
  x = ref_x + 7.1, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)




# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 9.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot
blot_image <- readImage(grh_DBD_blot)

# rotate image
blot_image <- blot_image |> 
  rotate(88, bg.col = "white")

# crop image
blot_image <- blot_image[3537:4583,540:1350]

# get blot aspect ratio
blot_dim <- dim(blot_image)
IF_aspect_ratio <- blot_dim[2] / blot_dim[1]


# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 1,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# add labels to western blot
plotSegments(
  x0 = ref_x + 1.5, y0 = ref_y + 0.25, x1 = ref_x + 2.25, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 2.5, y0 = ref_y + 0.25, x1 = ref_x + 7.25, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotText(
  label = "S2 Grh FL", params = large_text_params, fontface = "bold",
  x = ref_x + 1.75, y = ref_y, just = "center", default.units = "cm"
)

plotText(
  label = "S2 Grh DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 4.6, y = ref_y, just = "center", default.units = "cm"
)


plotText(
  label = "FL", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 1.6, just = "right", default.units = "cm"
)

plotText(
  label = "DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 6, just = "right", default.units = "cm"
)

plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 0.75, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "100", params = large_text_params, 
  x = ref_x + 1.85, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, 
  x = ref_x + 2.75, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "5", params = large_text_params, 
  x = ref_x + 3.6, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "10", params = large_text_params, 
  x = ref_x + 4.5, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "20", params = large_text_params, fontcolor = "firebrick3",
  x = ref_x + 5.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "40", params = large_text_params, 
  x = ref_x + 6.3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "80", params = large_text_params, 
  x = ref_x + 7.2, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 8.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# S2-Zld-FL --------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 681:839,
  y = 1737:1895
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.2),
  DAPI = c(0,0.6),
  BF = c(0,0.7)
)
# read input images
images <- list(
  GFP = readImage(`S2-Zld-FL_aZld_GFP`),
  DAPI = readImage(`S2-Zld-FL_aZld_DAPI`),
  BF =  readImage(`S2-Zld-FL_aZld_BF`)
) 

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# S2-Zld-DBD -------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 606:764,
  y = 205:363
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.15),
  DAPI = c(0,0.4),
  BF = c(0,0.6)
)
# read input images
images <- list(
  GFP = readImage(`S2-Zld-DBD_aZld_GFP`),
  DAPI = readImage(`S2-Zld-DBD_aZld_DAPI`),
  BF =  readImage(`S2-Zld-DBD_aZld_BF`)
)

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# S2-WT aZld -------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 792:956,
  y = 863:1027
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.15),
  DAPI = c(0,0.4),
  BF = c(0,0.6)
)
# read input images
images <- list(
  GFP = readImage(`S2-WT_aZld_GFP`),
  DAPI = readImage(`S2-WT_aZld_DAPI`),
  BF =  readImage(`S2-WT_aZld_BF`)
)

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# add labels -------------------------------------------------------------------
plotText(
  label = "anti-Zld", params = large_text_params, fontface = "bold",
  x = ref_x + 2, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "DAPI", params = large_text_params, fontface = "bold",
  x = ref_x + 4.5, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "brightfield", params = large_text_params, fontface = "bold",
  x = ref_x + 7, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "Zld FL", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 1.5, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "Zld DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 4, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "WT", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 6.5, just = c("right","center"), default.units = "cm"
)


# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 9.5
ref_y <- 8.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# S2-Grh-FL --------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 970:1366,
  y = 642:1038
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.12),
  DAPI = c(0,1),
  BF = c(0,0.6)
)
# read input images
images <- list(
  GFP = readImage(`S2-Grh-FL_aGrh_GFP`),
  DAPI = readImage(`S2-Grh-FL_aGrh_DAPI`),
  BF =  readImage(`S2-Grh-FL_aGrh_BF`)
) 

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# S2-Grh-DBD -------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 1593:1717,
  y = 1723:1847
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.15),
  DAPI = c(0,0.25),
  BF = c(0,0.65)
)
# read input images
images <- list(
  GFP = readImage(`S2-Grh-DBD_aGrh_GFP`),
  DAPI = readImage(`S2-Grh-DBD_aGrh_DAPI`),
  BF =  readImage(`S2-Grh-DBD_aGrh_BF`)
)

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 3,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# S2-WT aGrh -------------------------------------------------------------------
# parameters for image processing

# coordinates for cell of interest
coordinates <- list(
  x = 257:411,
  y = 1419:1573
)

# signal intensity range for histogram normalization
intensity_range <- list(
  GFP = c(0,0.15),
  DAPI = c(0,0.25),
  BF = c(0,0.65)
)
# read input images
images <- list(
  GFP = readImage(`S2-WT_aGrh_GFP`),
  DAPI = readImage(`S2-WT_aGrh_DAPI`),
  BF =  readImage(`S2-WT_aGrh_BF`)
)

# process images
images <- process_IF_images(x = images, coordinates = coordinates, normalize = TRUE, intensity_range = intensity_range)

# get aspect ratio
IF_dim <- dim(images$GFP)
IF_aspect_ratio <- IF_dim[2] / IF_dim[1]



# place images on page
plot_width <- 2

plotRaster(
  images$GFP,
  x = ref_x + 1,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$DAPI,
  x = ref_x + 3.5,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  images$BF,
  x = ref_x + 6,
  y = ref_y + 5.5,
  width = plot_width,
  height = plot_width * IF_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# add labels -------------------------------------------------------------------
plotText(
  label = "anti-Grh", params = large_text_params, fontface = "bold",
  x = ref_x + 2, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "DAPI", params = large_text_params, fontface = "bold",
  x = ref_x + 4.5, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "brightfield", params = large_text_params, fontface = "bold",
  x = ref_x + 7, y = ref_y + 0.25, just = "center", default.units = "cm"
)

plotText(
  label = "Grh FL", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 1.5, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "Grh DBD", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 4, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "WT", params = large_text_params, fontface = "bold",
  x = ref_x + 0.75, y = ref_y + 6.5, just = c("right","center"), default.units = "cm"
)


# close graphics device ========================================================
dev.off()
