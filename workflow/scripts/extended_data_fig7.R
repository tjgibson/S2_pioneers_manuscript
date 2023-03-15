# setup ------------------------------------------------------------------------
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(eulerr)
library(RColorBrewer)
library(plotgardener)

source("workflow/scripts/utils.R")
source("workflow/scripts/plot_heatmap.R")
source("../useful_functions/peak_processing.R")

# define input files ===========================================================
# define input files explicitly for interactively testing script
# S2_Zld_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Zld_aZld_IP.bed"
# nc14_Zld_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/embryo-nc14_aZld_summits.bed"
# brain_Zld_ChIP_fn <-  "../S2_pioneers/results/bed_files/peaks/brain_aZld.bed"
# 
# Zld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# Zld_nc14_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/individual/embryo-nc14_aZld.bw"
# 
# S2_Grh_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP.bed"
# embryo_Grh_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/embryo-15-16H_aGrh_summits.bed"
# wing_Grh_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/wing-disc_aGrh_summits.bed"
# 
# Grh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# Grh_embryo_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw"
# 
# S2_twi_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP.bed"
# embryo_twi_ChIP_fn <- "published_ChIPseq/results/peaks/merged/narrow/embryo-1-3H_aTwi_summits.bed"
# 
# Twi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# Twi_embryo_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw"
# 
# H3K27me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/individual/GSE151983_S2_aH3K27me3_IP.bw"
# H3K9me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw"

# get input files from snakemake
S2_Zld_ChIP_fn <- snakemake@input[["S2_Zld_ChIP_fn"]]
nc14_Zld_ChIP_fn <- snakemake@input[["nc14_Zld_ChIP_fn"]]
brain_Zld_ChIP_fn <-  snakemake@input[["brain_Zld_ChIP_fn"]]

Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
Zld_nc14_ChIP_bw <- snakemake@input[["Zld_nc14_ChIP_bw"]]

S2_Grh_ChIP_fn <- snakemake@input[["S2_Grh_ChIP_fn"]]
embryo_Grh_ChIP_fn <- snakemake@input[["embryo_Grh_ChIP_fn"]]
wing_Grh_ChIP_fn <- snakemake@input[["wing_Grh_ChIP_fn"]]

Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
Grh_embryo_ChIP_bw <- snakemake@input[["Grh_embryo_ChIP_bw"]]

S2_twi_ChIP_fn <- snakemake@input[["S2_twi_ChIP_fn"]]
embryo_twi_ChIP_fn <- snakemake@input[["embryo_twi_ChIP_fn"]]

Twi_ChIP_bw <- snakemake@input[["Twi_ChIP_bw"]]
Twi_embryo_ChIP_bw <- snakemake@input[["Twi_embryo_ChIP_bw"]]

H3K27me3_ChIP_bw <- snakemake@input[["H3K27me3_ChIP_bw"]]
H3K9me3_ChIP_bw <- snakemake@input[["H3K9me3_ChIP_bw"]]


# create blank layout for plot ===============================================
fig_width <-  18
fig_height <- 18

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig7.pdf", useDingbats = FALSE)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
zld_heatmap_colors <- brewer.pal(9, "Blues")

grh_color <- "#F98400"
grh_heatmap_colors <- brewer.pal(9, "Oranges")

twi_color <- "#00A08A"
twi_heatmap_colors <- brewer.pal(9, "GnBu")

H3K27me3_color <- "gray40"      
  
  # reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_browser_label <- 1
x_offset_browser <- 1.5

# set genome browser height
gb_height <- 0.3
  
# panel A ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
peak_overlaps <- 
  c(
    S2_Zld_ChIP = S2_Zld_ChIP_fn,
    nc14_Zld_ChIP = nc14_Zld_ChIP_fn,
    brain_Zld_ChIP = brain_Zld_ChIP_fn
  ) |> 
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  
  
  GRangesList() |> 
  
  peak_overlap_table() |> 
  mutate(overlap_group = case_when(
    S2_Zld_ChIP & !nc14_Zld_ChIP & !brain_Zld_ChIP ~ 1,
    !S2_Zld_ChIP & nc14_Zld_ChIP & !brain_Zld_ChIP ~ 2,
    !S2_Zld_ChIP & !nc14_Zld_ChIP & brain_Zld_ChIP ~ 3,
    S2_Zld_ChIP & nc14_Zld_ChIP & !brain_Zld_ChIP ~ 4,
    S2_Zld_ChIP & !nc14_Zld_ChIP & brain_Zld_ChIP ~ 5,
    !S2_Zld_ChIP & nc14_Zld_ChIP & brain_Zld_ChIP ~ 6,
    S2_Zld_ChIP & nc14_Zld_ChIP & brain_Zld_ChIP ~ 7,
  ))


# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:8) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(zld_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.5), y = (ref_y),
  width = 4, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Zld NSC", params = small_text_params, fontface = "bold",
  x = (ref_x + 2), y = (ref_y), just = "top", default.units = "cm"
)

plotText(
  label = "Zld S2", params = small_text_params, fontface = "bold",
  x = (ref_x + 1), y = (ref_y + 2.5), just = "center", default.units = "cm"
)

plotText(
  label = paste0("Zld nc14 embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.5), y = (ref_y + 2.3), just = "center", default.units = "cm"
)



# panel B ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 6.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
# generate overlap table and annotate ChIP classes


peak_overlaps <- 
  c(
    S2_Grh_ChIP = S2_Grh_ChIP_fn,
    embryo_Grh_ChIP = embryo_Grh_ChIP_fn,
    wing_Grh_ChIP = wing_Grh_ChIP_fn
  ) |> 
  
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  GRangesList() |> 
  peak_overlap_table() |> 
  
  mutate(overlap_group = case_when(
    S2_Grh_ChIP & !embryo_Grh_ChIP & !wing_Grh_ChIP ~ 1,
    !S2_Grh_ChIP & embryo_Grh_ChIP & !wing_Grh_ChIP ~ 2,
    !S2_Grh_ChIP & !embryo_Grh_ChIP & wing_Grh_ChIP ~ 3,
    S2_Grh_ChIP & embryo_Grh_ChIP & !wing_Grh_ChIP ~ 4,
    S2_Grh_ChIP & !embryo_Grh_ChIP & wing_Grh_ChIP ~ 5,
    !S2_Grh_ChIP & embryo_Grh_ChIP & wing_Grh_ChIP ~ 6,
    S2_Grh_ChIP & embryo_Grh_ChIP & wing_Grh_ChIP ~ 7,
  ))



# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:8) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(grh_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.25), y = (ref_y),
  width = 4, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Grh wing disc", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.5), y = (ref_y + 0.25), just = "top", default.units = "cm"
)

plotText(
  label = paste0("Grh S2"), params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 2), just = "center", default.units = "cm"
)

plotText(
  label = paste0("Grh","\n", "15-16H embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y + 1.5), just = "center", default.units = "cm"
)


# panel C ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 12.5
ref_y <- 0.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
# generate overlap table and annotate ChIP classes


peak_overlaps <- 
  c(
    S2_twi_ChIP = S2_twi_ChIP_fn,
    embryo_twi_ChIP = embryo_twi_ChIP_fn
  ) |> 
  
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  GRangesList() |> 
  peak_overlap_table() |> 
  
  mutate(overlap_group = case_when(
    S2_twi_ChIP & !embryo_twi_ChIP ~ 1,
    !S2_twi_ChIP & embryo_twi_ChIP ~ 2,
    S2_twi_ChIP & embryo_twi_ChIP ~ 3
  ))



# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:7) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(twi_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.25), y = (ref_y),
  width = 4, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = paste0("Twi S2"), params = small_text_params, fontface = "bold",
  x = (ref_x + 1.6), y = (ref_y + 1.8), just = "center", default.units = "cm"
)

plotText(
  label = paste0("Twi", "\n", "2-3H embryo"), params = small_text_params, fontface = "bold",
  x = (ref_x + 4.5), y = (ref_y + 1.2), just = "center", default.units = "cm"
)

# panel D ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 4.5

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class IV", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chr3R",
  chromstart = 30824716, chromend = 30888891,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr3R",
  chromstart = 30846356, chromend = 30858282,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1.5),
  default.units = "cm",
  lty = 2
)



plotText(
  label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                              params = zoomRegion
)

`Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Zld_S2_ChIP$score, Zld_nc14_ChIP$score))

s2 <- plotSignal(
  data = `Zld_nc14_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Zld_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)


# panel E ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 8.5

plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class V", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chr3L",
  chromstart = 11796001, chromend = 12064216,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr3L",
  chromstart = 11919547, chromend = 11940115,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                              params = zoomRegion
)

`Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Zld_S2_ChIP$score, Zld_nc14_ChIP$score))

s2 <- plotSignal(
  data = `Zld_nc14_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Zld_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)

# panel F ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 12.5

plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class VI", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chrX",
  chromstart = 393315, chromend = 402092,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- c(-1.2,8)

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chrX",
  chromstart = 393449, chromend = 400413,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                              params = zoomRegion
)

`Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Zld_S2_ChIP$score, Zld_nc14_ChIP$score))

s2 <- plotSignal(
  data = `Zld_nc14_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Zld_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)


# panel G ==================================================================
# panel label
ref_x <- 6.5
ref_y <- 4.5

plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class IV", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chr3L",
  chromstart = 3232675, chromend = 3324202,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr3L",
  chromstart = 3254519, chromend = 3272949,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Grh_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Grh_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)

# panel H ==================================================================
# panel label
ref_x <- 6.5
ref_y <- 8.5

plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class V", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chr2L",
  chromstart = 4443392, chromend = 4669912,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr2L",
  chromstart = 4562557, chromend = 4626155,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Grh_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Grh_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)


# panel I ==================================================================
# panel label
ref_x <- 6.5
ref_y <- 12.5

plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class VI", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

# region candidates:
# chr3R:22,748,753-22,807,827
# chrX:13,216,271-13,269,214

region <- pgParams(
  chrom = "chrX",
  chromstart = 13181251, chromend = 13295616,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- c(-1.1,8)

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chrX",
  chromstart = 13218840, chromend = 13238073,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Grh_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Grh_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)

# panel J ==================================================================
# panel label
ref_x <- 12.5
ref_y <- 4.5

plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class IV", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

region <- pgParams(
  chrom = "chr2R",
  chromstart = 12753782, chromend = 12796150,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr2R",
  chromstart = 12773275, chromend = 12791478,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Twi_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Twi_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)

# panel K ==================================================================
# panel label
ref_x <- 12.5
ref_y <- 8.5

plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class V", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)

# candidate regions
# chrX:18,315,069-18,367,911
# chrX:17891901-18095637
# chr3L:13,423,435-14,067,450

region <- pgParams(
  chrom = "chr3L",
  chromstart = 13423435, chromend = 14067450,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr3L",
  chromstart = 13658008, chromend = 13683890,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Twi_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Twi_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)


# panel L ==================================================================
# panel label
ref_x <- 12.5
ref_y <- 12.5

plotText(
  label = "l", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class VI", params = large_text_params, fontface = "bold",
  x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
)


# add broad H3K27me3 region
plotText(
  label = paste("H3K27me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.25),
  just = c("right", "center"),
  default.units = "cm"
)

plotText(
  label = paste("H3K9me3"),
  params = small_text_params,
  fontface = "bold",
  x = (ref_x + x_offset_browser_label),
  y = (ref_y + 0.75),
  just = c("right", "center"),
  default.units = "cm"
)


region <- pgParams(
  chrom = "chr3L",
  chromstart = 10650450, chromend = 10684646,
  assembly = "dm6"
)

`H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                              params = region
)

`H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                             params = region
)

ChIP_range <- c(-1.1,8)

s1 <- plotSignal(
  data = `H3K27me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `H3K9me3_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


# add zoomed in view of ChIP
zoomRegion <- pgParams(
  chrom = "chr3L",
  chromstart = 10667170, chromend = 10675448,
  assembly = "dm6"
)

annoZoomLines(
  plot = s2, params = zoomRegion,
  y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 1.4), extend = c(1, 1),
  default.units = "cm",
  lty = 2
)

plotText(
  label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)



`Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                params = zoomRegion
)

`Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                            params = zoomRegion
)



ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))

s2 <- plotSignal(
  data = `Twi_embryo_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s3 <- plotSignal(
  data = `Twi_S2_ChIP`, params = zoomRegion,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
  range = ChIP_range
)

annoYaxis(
  plot = s3, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



plotGenes(
  params = zoomRegion, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
  just = c("left", "top"),
  default.units = "cm"
)


dev.off()
