# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(EBImage)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]

zld_titration_blot_aZld <- "data/immunoblot_raw_images/2022-7-12_titration/Zld-titration_aZld.tif"

grh_titration_blot_aGrh <- "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aGrh_short.tif"
grh_titration_blot_aTub <- "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aTub.tif"

twi_titration_blot_aTwi <- "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTwi.tif"
twi_titration_blot_aTub <- "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTub.tif"

# open graphics device =========================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig9.pdf")
# create blank layout for plot =================================================
pageCreate(width = 18.3, height = 12, default.units = "cm", showGuides = FALSE)

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


# close graphics device ========================================================
dev.off()
