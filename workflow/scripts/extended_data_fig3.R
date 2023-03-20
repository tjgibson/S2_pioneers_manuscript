# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(EBImage)

source("workflow/scripts/plot_heatmap.R")

# define input files ===========================================================
# specify input files explicitly for testing script interactively
# S2_Zld_CR_0H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_0H_small.bw"
# S2_Zld_CR_4H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_4H_small.bw"
# S2_Zld_CR_12H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_12H_small.bw"
# S2_Zld_CR_24H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_24H_small.bw"
# S2_Zld_CR_48H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_48H_small.bw"
# 
# S2_Grh_CR_0H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_0H_small.bw"
# S2_Grh_CR_4H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_4H_small.bw"
# S2_Grh_CR_12H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_12H_small.bw"
# S2_Grh_CR_24H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_24H_small.bw"
# S2_Grh_CR_48H_bw <- "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_48H_small.bw"
# 
# 
# zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# 
# zld_blot_image <- "data/immunoblot_raw_images/2018-11-27_timecourse/zld_timecourse.tif"
# grh_blot_image <- "data/immunoblot_raw_images/2018-11-27_timecourse/grh_timecourse.tif"

# get input files from snakemake
S2_Zld_CR_0H_bw <- snakemake@input[["S2_Zld_CR_0H_bw"]]
S2_Zld_CR_4H_bw <- snakemake@input[["S2_Zld_CR_4H_bw"]]
S2_Zld_CR_12H_bw <- snakemake@input[["S2_Zld_CR_12H_bw"]]
S2_Zld_CR_24H_bw <- snakemake@input[["S2_Zld_CR_24H_bw"]]
S2_Zld_CR_48H_bw <- snakemake@input[["S2_Zld_CR_48H_bw"]]

S2_Grh_CR_0H_bw <- snakemake@input[["S2_Grh_CR_0H_bw"]]
S2_Grh_CR_4H_bw <- snakemake@input[["S2_Grh_CR_4H_bw"]]
S2_Grh_CR_12H_bw <- snakemake@input[["S2_Grh_CR_12H_bw"]]
S2_Grh_CR_24H_bw <- snakemake@input[["S2_Grh_CR_24H_bw"]]
S2_Grh_CR_48H_bw <- snakemake@input[["S2_Grh_CR_48H_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes_fn"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes_fn"]]

zld_blot_image <- snakemake@input[["zld_blot_image"]]
grh_blot_image <- snakemake@input[["grh_blot_image"]]

# create blank layout for plot ===============================================
fig_width <-  18
fig_height <- 14

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig3.pdf", useDingbats = FALSE)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

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

# read in tiff of western blot in grayscale
blot_image <- readImage(zld_blot_image) |> 
  channel("gray")

# rotate image
blot_image <- blot_image |>
  rotate(1, bg.col = "white")

# adjust brightness and contrast
blot_image <- (blot_image * 1.1 - 0.3)

# crop image
blot_image <- blot_image[501:776,393:518]

# get blot aspect ratio
blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotText(
  label = "0H", params = large_text_params,
  x = ref_x + 1.8, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "1H", params = large_text_params,
  x = ref_x + 2.9, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "4H", params = large_text_params, 
  x = ref_x + 4, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "12H", params = large_text_params,
  x = ref_x + 5.1, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "24H", params = large_text_params,
  x = ref_x + 6.2, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "48H", params = large_text_params,
  x = ref_x + 7.3, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

# add arrowheads to blot
plotPolygon(x = c(1, 1.25, 1), y = c(1.4,1.5,1.6), default.units = "cm", fill = "black")

plotPolygon(x = c(1, 1.25, 1), y = c(3.6,3.7,3.8), default.units = "cm", fill = "grey")


# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 9
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot in grayscale
blot_image <- readImage(grh_blot_image) |> 
  channel("gray")

# rotate image
blot_image <- blot_image |>
  rotate(1, bg.col = "white")

# adjust brightness and contrast
blot_image <- (blot_image - 0.1)

# crop image
blot_image <- blot_image[849:1131,577:744]

# get blot aspect ratio
blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 0.5,
  width = plot_width,
  height = plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotText(
  label = "0H", params = large_text_params,
  x = ref_x + 1.8, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "1H", params = large_text_params,
  x = ref_x + 2.9, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "4H", params = large_text_params, 
  x = ref_x + 4, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "12H", params = large_text_params,
  x = ref_x + 5, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "24H", params = large_text_params,
  x = ref_x + 6.1, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "48H", params = large_text_params,
  x = ref_x + 7.1, y = ref_y + 0.25, just = c("center"), default.units = "cm"
)

# add arrowheads to blot
plotPolygon(x = c(9.5, 9.75, 9.5), y = c(1.7,1.8,1.9), default.units = "cm", fill = "black")

plotPolygon(x = c(9.5, 9.75, 9.5), y = c(4.6,4.7,4.8), default.units = "cm", fill = "grey")


# Panel C ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 6

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


# define regions to plot
regions <- zld_ChIP_classes |>
  mutate(class = str_to_upper(class)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.1,3.5)
plot_colors <- c(
  `I` = "#BDD7E7",
  `II` = "#6BAED6",
  `III` = "#2171B5"
)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        # legend.key.size = unit(2, 'mm'),
        # legend.title = element_blank(),
        # legend.position = "bottom",
        # legend.margin=margin(-5,-5,-5,-5),
        # legend.box.margin=margin(0,0,0,0)
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        # legend.key.size = unit(2, 'mm'),
        # legend.title = element_blank(),
        # legend.position = "bottom",
        # legend.margin=margin(-5,-5,-5,-5),
        # legend.box.margin=margin(0,0,0,0)
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 3.75), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_3 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_3,
  x = (ref_x + 7.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_4,
  x = (ref_x + 10.75), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_5 <- plot_average(bw[5], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_5,
  x = (ref_x + 14.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "CUT&RUN signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.75), just = "center", default.units = "cm", rot = 90
)

plotText(
  label = "0H", params = large_text_params,
  x = (ref_x + 1.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "4H", params = large_text_params,
  x = (ref_x + 5.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "12H", params = large_text_params,
  x = (ref_x + 8.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "24H", params = large_text_params,
  x = (ref_x + 12.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "48H", params = large_text_params,
  x = (ref_x + 15.75), y = (ref_y), just = "center", default.units = "cm"
)

# add legend for metaplots
plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 16.8, y = 6.5, width = 1, height = 1,
  just = c("center", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)


# Panel D ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 9.75

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_Grh_CR_0H = S2_Grh_CR_0H_bw,
  S2_Grh_CR_4H = S2_Grh_CR_4H_bw,
  S2_Grh_CR_4H = S2_Grh_CR_12H_bw, 
  S2_Grh_CR_24H = S2_Grh_CR_24H_bw,
  S2_Grh_CR_48H = S2_Grh_CR_48H_bw
)


# define regions to plot
regions <- grh_ChIP_classes |>
  mutate(class = str_to_upper(class)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.1,3)
plot_colors <- c(
  `I` = "#FDBE85",
  `II` = "#FD8D3C",
  `III` = "#D94701"
)
  
metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        # legend.key.size = unit(2, 'mm'),
        # legend.title = element_blank(),
        # legend.position = "bottom",
        # legend.margin=margin(-5,-5,-5,-5),
        # legend.box.margin=margin(0,0,0,0)
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        # legend.key.size = unit(2, 'mm'),
        # legend.title = element_blank(),
        # legend.position = "bottom",
        # legend.margin=margin(-5,-5,-5,-5),
        # legend.box.margin=margin(0,0,0,0)
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 3.75), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_3 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_3,
  x = (ref_x + 7.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_4,
  x = (ref_x + 10.75), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_5 <- plot_average(bw[5], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_5,
  x = (ref_x + 14.25), y = (ref_y + 0.25),
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "CUT&RUN signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.75), just = "center", default.units = "cm", rot = 90
)

plotText(
  label = "0H", params = large_text_params,
  x = (ref_x + 1.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "4H", params = large_text_params,
  x = (ref_x + 5.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "12H", params = large_text_params,
  x = (ref_x + 8.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "24H", params = large_text_params,
  x = (ref_x + 12.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "48H", params = large_text_params,
  x = (ref_x + 15.75), y = (ref_y), just = "center", default.units = "cm"
)

# add legend for metaplots
plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 16.8, y = 10.25, width = 1, height = 1,
  just = c("center", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

# close graphics device ========================================================
dev.off()




