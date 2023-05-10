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
# define input files explicitly for interactively testing script
# CR_spikeIn_counts_fn <-"histone_CUTandRUN/results/scaling_factors/epiCypher_barcode_counts.tsv"
# 
# taz_blot_aH3K27me3_image <- "data/immunoblot_raw_images/2021-06-29_taz/anti-H3K27me3_3.tif"
# taz_blot_aTub_image <- "data/immunoblot_raw_images/2021-06-29_taz/anti-tubulin_2.tif"

# get input files from snakemake
CR_spikeIn_counts_fn <-snakemake@input[["CR_spikeIn_counts_fn"]]
taz_blot_aH3K27me3_image <- snakemake@input[["taz_blot_aH3K27me3_image"]]
taz_blot_aTub_image <- snakemake@input[["taz_blot_aTub_image"]]


# create blank layout for plot ===============================================
fig_width <-  12.5
fig_height <- 6

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig8.pdf", useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

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


# panel A ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 0.5

plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import images of western blot of H3K27me3 and tubulin in DMSO or control cells
# read in tiff of western blot
H3K27me3_blot_image <- readImage(taz_blot_aH3K27me3_image)

# rotate image
H3K27me3_blot_image <- H3K27me3_blot_image |> 
  rotate(91, bg.col = "white") |> 
  flop()

# crop image
H3K27me3_blot_image <- H3K27me3_blot_image[1536:2188,880:1102]

# anti-Tubulin
tub_blot_image <- readImage(taz_blot_aTub_image)

# rotate image
tub_blot_image <- tub_blot_image |> 
  rotate(91, bg.col = "white") |> 
  flop()

# crop image
tub_blot_image <- tub_blot_image[3280:3932,1184:1406]

# adjust brightness and contrast
tub_blot_image <- (tub_blot_image - 0.3) 

# get aspect ratio
H3K27me3_blot_dim <- dim(H3K27me3_blot_image)
H3K27me3_blot_aspect_ratio <- H3K27me3_blot_dim[2] / H3K27me3_blot_dim[1]

tub_blot_dim <- dim(tub_blot_image)
tub_blot_aspect_ratio <- tub_blot_dim[2] / tub_blot_dim[1]


# check that images were cropped to the same width
if (H3K27me3_blot_dim[1] != tub_blot_dim[1]) {stop("K27me3 and anti-tub blot are different widths")}


# place blot on page
plot_width <- 4

plotRaster(
  H3K27me3_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 0.5,
  width = plot_width,
  height =  plot_width * H3K27me3_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotRaster(
  tub_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 2.25,
  width = plot_width,
  height =  plot_width * tub_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

plotSegments(
  x0 = ref_x + 1.8, y0 = ref_y + 0.25, x1 = ref_x + 3.4, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 3.6, y0 = ref_y + 0.25, x1 = ref_x + 5.1, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

# panel label
plotText(
  label = "DMSO", params = large_text_params, fontface = "bold",
  x = ref_x + 2.7, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "taz", params = large_text_params, fontface = "bold",
  x = ref_x + 4.25, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "H3K27me3", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "tubulin", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 2.8, just = c("right","center"), default.units = "cm"
)

# panel B ==================================================================
# panel label
ref_x <- 6.5
ref_y <- 0.5

plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in barcode counts
CR_spikeIn_counts <- read_tsv(CR_spikeIn_counts_fn)

AB_specificity <- CR_spikeIn_counts |>
  group_by(sample_name, target) |>
  summarise(barcode_count = sum(count))

barcode_totals <- CR_spikeIn_counts |>
  group_by(sample_name) |>
  summarise(total_count = sum(count))

AB_specificity <- AB_specificity |>
  left_join(barcode_totals) |>
  mutate(barcode_percent = barcode_count / total_count * 100)


sample_order <- c(
  "S2-Zld_DMSO_aH3K27me3_rep1",
  "S2-Zld_DMSO_aH3K27me3_rep2",
  "S2-Grh_DMSO_aH3K27me3_rep1",
  "S2-Grh_DMSO_aH3K27me3_rep2",
  "S2-HA-Twi_DMSO_aH3K27me3_rep1",
  "S2-HA-Twi_DMSO_aH3K27me3_rep2",
  
  "S2-Zld_Taz_aH3K27me3_rep1",
  "S2-Zld_Taz_aH3K27me3_rep2",
  "S2-Grh_Taz_aH3K27me3_rep1",
  "S2-Grh_Taz_aH3K27me3_rep2",
  "S2-HA-Twi_Taz_aH3K27me3_rep1",
  "S2-HA-Twi_Taz_aH3K27me3_rep2",
  
  "S2-Zld_DMSO_IgG_rep1",
  "S2-Zld_DMSO_IgG_rep2",
  "S2-Grh_DMSO_IgG_rep1",
  "S2-Grh_DMSO_IgG_rep2",
  "S2-HA-Twi_DMSO_IgG_rep1",
  "S2-HA-Twi_DMSO_IgG_rep2"
)

sample_labels <- sample_order |>
  tibble() |>
  separate(sample_order, into = c("cell_line", "treatment", "antibody", "replicate"), sep = "_") |>
  unite("label", cell_line, replicate) |>
  pull(label)

target_order <- c(
  "Unmodified",
  "H3K4me1",
  "H3K4me2",
  "H3K4me3",
  "H3K9me1",
  "H3K9me2",
  "H3K9me3",
  "H3K27me1",
  "H3K27me2",
  "H3K27me3",
  "H3K36me1",
  "H3K36me2",
  "H3K36me3",
  "H4K20me1",
  "H4K20me2",
  "H4K20me3"
)


# generate plot
b_plot <- AB_specificity |>
  mutate(sample_name = factor(sample_name, levels = sample_order)) |>
  mutate(sample_name = fct_rev(sample_name)) |>
  mutate(target = factor(target, levels = target_order)) |>
  ggplot(aes(x = target, y = sample_name, fill = barcode_percent)) + geom_tile() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    text = element_text(size = small_text_params$fontsize),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10)
  ) +
  scale_y_discrete(labels = sample_labels)

# place heatmap on page
plotGG(
  plot = b_plot,
  x = (ref_x + 1),
  y = (ref_y),
  width = 4.5,
  height = 4,
  just = c("left", "top"),
  default.units = "cm"
)

# add sample labels to heatmap
plotSegments(
  x0 = (ref_x + 1), y0 = (ref_y + 0.25), x1 = (ref_x + 1), y1 = (ref_y + 0.95),
  default.units = "cm"
)

plotSegments(
  x0 = (ref_x + 1), y0 = (ref_y + 1.05), x1 = (ref_x + 1), y1 = (ref_y + 1.85),
  default.units = "cm"
)

plotSegments(
  x0 = (ref_x + 0.25), y0 = (ref_y + 0.25), x1 = (ref_x + 0.25), y1 = (ref_y + 1.85),
  default.units = "cm"
)


plotSegments(
  x0 = (ref_x + 0.25), y0 = (ref_y + 1.9), x1 = (ref_x + 0.25), y1 = (ref_y + 2.7),
  default.units = "cm"
)

plotText(
  label = "DMSO",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x + 0.8,
  y = ref_y + 0.6,
  default.units = "cm"
)

plotText(
  label = "taz",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x + 0.8,
  y = ref_y + 1.4,
  default.units = "cm"
)

plotText(
  label = "anti-H3K27me3",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x,
  y = ref_y + 1,
  default.units = "cm"
)

plotText(
  label = "IgG",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x,
  y = ref_y + 2.2,
  default.units = "cm"
)

# close graphics device ========================================================
dev.off()

