# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# library(grImport)
library(grid)
library(RColorBrewer)


source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# define input files explictly for interactive testing
# Twi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# Twi_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw"
# Twi_Twi_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Twi_40uM_small.bw"
# Zld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# Grh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# twi_ChIP_classes_fn <- "results/ChIP_peak_classes/twi_ChIP_classes.tsv"
# zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# twi_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Twi_RNAseq_S2-Twi-vs-S2-WT_results.tsv"
# zld_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results.tsv"
# grh_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results.tsv"

# get input files from snakemake
Twi_ChIP_bw <- snakemake@input[["Twi_ChIP_bw"]]
Twi_WT_ATAC_bw <- snakemake@input[["Twi_WT_ATAC_bw"]]
Twi_Twi_ATAC_bw <- snakemake@input[["Twi_Twi_ATAC_bw"]]

Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]

twi_ChIP_classes_fn <- snakemake@input[["twi_ChIP_classes"]]
zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

twi_RNAseq_results_fn <- snakemake@input[["twi_RNAseq_results"]]
zld_RNAseq_results_fn <- snakemake@input[["zld_RNAseq_results"]]
grh_RNAseq_results_fn <- snakemake@input[["grh_RNAseq_results"]]

# create blank layout for plot =================================================
# define figure dimensions in cm
fig_width <-  8.8
fig_height <- 14.5

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 7)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
twi_color <- "#00A08A"
# twi_heatmap_colors <- c("white","#00A08A", "#154734")
twi_heatmap_colors <- brewer.pal(9, "GnBu")

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


# class I genome browser tracks ------------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# class label
plotText(
  label = "class I", params = large_text_params, fontface = "bold",
  x = (ref_x + x_offset_class_label), y = (ref_y + 0.75), just = c("center"), default.units = "cm",
  rot = 90
)

# class I browser track labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Twi ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)


# plot Twi genome browser for class I
region <- pgParams(
  chrom = "chr2R",
  chromstart = 11466356, chromend = 11472917,
  assembly = "dm6"
)

`S2-Twi_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Twi_WT_ATAC_bw, 
                           params = region
)

`S2-Twi_ATAC` <- readBigwig(file = Twi_Twi_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Twi_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Twi_ATAC`$score))

s1 <- plotSignal(
  data = `S2-Twi_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Twi_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

annoGenomeLabel(
  plot = s3, scale = "bp", fontsize = small_text_params$fontsize,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.5), just = c("left", "top"),
  default.units = "cm"
)

# plotGenes(
#   params = region, assembly = "dm6",
#   x = 2.75, y = 2.25, height = 0.4, width = 3,
#   just = c("left", "center"),
#   default.units = "cm"
# )



# class II genome browser tracks -----------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 2.5

# class label
plotText(
  label = "class II", params = large_text_params, fontface = "bold",
  x = (ref_x + x_offset_class_label), y = (ref_y + 0.75), just = c("center"), default.units = "cm",
  rot = 90
)

# class II browser track labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Twi ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

# plot Twi genome browser for class II
region <- pgParams(
  chrom = "chr2L",
  chromstart = 16408217, chromend = 16421641,
  assembly = "dm6"
)

`S2-Twi_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Twi_WT_ATAC_bw, 
                           params = region
)

`S2-Twi_ATAC` <- readBigwig(file = Twi_Twi_ATAC_bw, 
                            params = region
)

ChIP_range <- c(-0.5,3)
ATAC_range <- c(-0.5, 10)


s1 <- plotSignal(
  data = `S2-Twi_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)


annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Twi_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

annoGenomeLabel(
  plot = s3, scale = "bp", fontsize = small_text_params$fontsize,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.5), just = c("left", "top"),
  default.units = "cm"
)

annoHighlight(
  plot = s1,
  chrom = "chr2L",
  chromstart = 16414631,
  chromend = 16415188,
  y = (ref_y), height = 1.5, just = c("left", "top"),
  default.units = "cm"
)

# plotGenes(
#   params = region, assembly = "dm6",
#   x = 2.75, y = 2.25, height = 0.4, width = 3,
#   just = c("left", "center"),
#   default.units = "cm"
# )


# class III genome browser tracks ----------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 4.5

# class label
plotText(
  label = "class III", params = large_text_params, fontface = "bold",
  x = (ref_x + x_offset_class_label), y = (ref_y + 0.75), just = c("center"), default.units = "cm",
  rot = 90
)

# class III browser track labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Twi ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)



# plot Zld genome browser for class III
region <- pgParams(
  chrom = "chrX",
  chromstart = 10390659, chromend = 10417540,
  assembly = "dm6"
)

`S2-Twi_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Twi_WT_ATAC_bw, 
                           params = region
)

`S2-Twi_ATAC` <- readBigwig(file = Twi_Twi_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Twi_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Twi_ATAC`$score))

s1 <- plotSignal(
  data = `S2-Twi_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s3 <- plotSignal(
  data = `S2-Twi_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_x + 5.25), width = 5.5, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

annoGenomeLabel(
  plot = s3, scale = "bp", fontsize = small_text_params$fontsize,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.5), just = c("left", "top"),
  default.units = "cm"
)

# plotGenes(
#   params = region, assembly = "dm6",
#   x = (ref_x + 2.25), y = (ref_x + 5.75), height = 0.4, width = 3,
#   just = c("left", "center"),
#   default.units = "cm"
# )


# Panel B ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 6.5 

# panel label

plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in ChIP classes
twi_chip_classes <- read_tsv(twi_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

# generate pie charts
twi_class_plot <- twi_chip_classes  |> 
  # mutate(zygotic_class = factor(zygotic_class, levels = c("I", "II", "III"))) |>
  group_by(class) |> summarise(n = n()) |>
  ggplot(aes(x='', y = n, fill = class)) + 
  geom_bar(stat = "identity", position = "fill") + 
  coord_polar("y", direction = 1, start = pi / 2) + 
  # geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 0.35 * small_text_params$fontsize) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#DFF3F0", "#5FC3B5", "#154734"))


# place pie charts on page
plotGG(
  plot = twi_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add labels
plotText(
  label = paste("class I", "\n", "n =", sum(twi_chip_classes$class == "I")), fontsize = small_text_params$fontsize,
  x = (ref_x + 1.5), y = (ref_y + 1), just = "center", default.units = "cm"
)

plotText(
  label = paste("class II", "\n", "n =", sum(twi_chip_classes$class == "II")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2), y = (ref_y + 2), just = "center", default.units = "cm"
)

plotText(
  label = paste("class III", "\n", "n =", sum(twi_chip_classes$class == "III")), fontsize = small_text_params$fontsize,
  x = (ref_x + 3.2), y = (ref_y + 1.5), just = "center", default.units = "cm"
)


# Panel C ======================================================================
# panel label
ref_x <- 4.5
ref_y <- 6.5 

# generate chart
twi_feature_plot <- twi_chip_classes  |> 
  ggplot(aes(x=class, fill = feature)) + 
  geom_bar(position = "fill") + 
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_manual(values=c("#DFF3F0", "#5FC3B5")) +
  ylab("proportion")


# place chart on page
plotGG(
  plot = twi_feature_plot,
  x = (ref_x), y = ref_y,
  width = 4, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# panel label

plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel D ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 9.5


# panel label

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# generate heatmap
bw <- c(
  S2_Twi_ChIP =  Twi_ChIP_bw,
  S2_WT_ATAC = Twi_WT_ATAC_bw,
  S2_Twi_ATAC = Twi_Twi_ATAC_bw
)

groups <- c(1,2,2)

regions <- twi_chip_classes |>
  filter(class != "I") |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_split = regions$class, order_by_samples = 1, 
  scale_group = groups,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)


d_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = d_hm,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 7.25, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_2_1")
grid.xaxis(at = c(0, 0.5, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), "peak center", paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.325), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 Twi ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 6.875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class II", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 2.5), just = c("right", "center"), default.units = "cm"
)

plotText(
  label = "class III", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 4.1), just = c("right","center"), default.units = "cm"
)


# close graphics device ========================================================
dev.off()

