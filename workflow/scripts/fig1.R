# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")

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
Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]
Zld_WT_RNAseq_bw <- snakemake@input[["Zld_WT_RNAseq_bw"]]
Zld_Zld_RNAseq_bw <- snakemake@input[["Zld_Zld_RNAseq_bw"]]

Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
Grh_Grh_ATAC_bw <- snakemake@input[["Grh_Grh_ATAC_bw"]]
Grh_WT_RNAseq_bw <- snakemake@input[["Grh_WT_RNAseq_bw"]]
Grh_Grh_RNAseq_bw <- snakemake@input[["Grh_Grh_RNAseq_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

# # create blank layout for plot ===============================================
pdf(snakemake@output[[1]], useDingbats = FALSE)
# tiff("results/pub_figs/2022_drafts/fig1.tiff", width = 4.392, height = 3, units = "in", res = 300)
pageCreate(width = 18.3, height = 12.5, default.units = "cm", showGuides = FALSE)

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

# reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_brower_label <- 1.75
x_offset_browser <- 2.5

# set genome browser height
gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500

# model ========================================================================
# PostScriptTrace(
#   "results/pub_figs/2022_drafts/system_schematic.eps",
#   outfilename = "results/pub_figs/2022_drafts/system_schematic.eps.xml",
#   charpath = FALSE,
#   charpos = FALSE
#   )
# schematic <- readPicture("results/pub_figs/2022_drafts/system_schematic.eps.xml")
# # grid.picture(schematic)
# # schematic_grob <- pictureGrob(schematic, exp = 0)
# schematic_gtree <- grid.grabExpr(grid.picture(schematic))
# 
# plotGG(schematic_gtree, x = 0.75, y = 0.5, width = 6, height = 2.5, default.units = "cm")
# 
# 
# 
# plotText(
#   label = "a", params = panel_label_params, fontface = "bold",
#   x = 0.5, y = 0.5, just = "center", default.units = "cm"
# )
# 
# plotGG(schematic_grob, x = 0.75, y = 0.5, width = 6, height = 2.5, default.units = "cm")


# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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


# browser track labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Zld ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

# browser tracks
region <- pgParams(
  chrom = "chr2R",
  chromstart = 11472922, chromend = 11478341,
  assembly = "dm6"
)

`S2-Zld_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                           params = region
)

`S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Zld_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Zld_ATAC`$score))

s1 <- plotSignal(
  data = `S2-Zld_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s3 <- plotSignal(
  data = `S2-Zld_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
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
#   x = (ref_x + x_offset_browser), y = (ref_y + 1.75), height = gb_height, width = 3, fontsize = small_text_params$fontsize,
#   just = c("left", "center"),
#   fill = c("gray50", "gray50"),
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

# browser track labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Zld ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

# plot Zld genome browser for class II
region <- pgParams(
  chrom = "chr3L",
  chromstart = 14964675, chromend = 14968087,
  assembly = "dm6"
)

`S2-Zld_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                           params = region
)

`S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Zld_ChIP`$score)
ATAC_range <- c(-1,3)


s1 <- plotSignal(
  data = `S2-Zld_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Zld_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
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
#   x = (ref_x + 2.25), y = (ref_y + 3.75), height = gb_height, width = 3,
#   fill = c("gray50", "gray50"),
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


# browser track labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)


plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Zld ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)



# plot Zld genome browser for class III
region <- pgParams(
  chrom = "chr2R",
  chromstart = 24345003, chromend = 24347154,
  assembly = "dm6"
)

`S2-Zld_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                           params = region
)

`S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Zld_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Zld_ATAC`$score))

s1 <- plotSignal(
  data = `S2-Zld_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_x + 4.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_x + 4.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Zld_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_x + 5.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoGenomeLabel(
  plot = s3, scale = "bp", fontsize = small_text_params$fontsize,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.5), just = c("left", "top"),
  default.units = "cm"
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

# plotGenes(
#   params = region, assembly = "dm6",
#   x = (ref_x + x_offset_browser), y = (ref_x + 5.75), height = gb_height, width = 3,
#   fill = c("gray50", "gray50"),
#   just = c("left", "center"),
#   default.units = "cm"
# )


# Panel B ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 6.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# class I genome browser tracks ------------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 6.5

# class label
plotText(
  label = "class I", params = large_text_params, fontface = "bold",
  x = (ref_x + x_offset_class_label), y = (ref_y + 1), just = c("center"), default.units = "cm",
  rot = 90
)

# class I browser track labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Grh ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)


# plot Grh genome browser for class I
region <- pgParams(
  chrom = "chr3R",
  chromstart = 21442922, chromend = 21450542,
  assembly = "dm6"
)

`S2-Grh_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                           params = region
)

`S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Grh_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))

s1 <- plotSignal(
  data = `S2-Grh_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range =ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Grh_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
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
#   x = (ref_x + x_offset_browser), y = (ref_y + 1.75), height = gb_height, width = 3,
#     fill = c("gray50", "gray50"),
#   just = c("left", "center"),
#   default.units = "cm"
# )


# class II genome browser tracks -----------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 8.5

# class label
plotText(
  label = "class II", params = large_text_params, fontface = "bold",
  x = (ref_x +x_offset_class_label), y = (ref_y + 0.75), just = c("center"), default.units = "cm",
  rot = 90
)

# class II browser track labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Grh ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

# plot Grh genome browser for class II
region <- pgParams(
  chrom = "chr2R",
  chromstart = 13923261, chromend = 13927608,
  assembly = "dm6"
)

`S2-Grh_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                           params = region
)

`S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Grh_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))


s1 <- plotSignal(
  data = `S2-Grh_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Grh_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
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
#   x = (ref_x + 2.25), y = (ref_y + 1.75), height = gb_height, width = 3,
#   fill = c("gray50", "gray50"),
#   just = c("left", "center"),
#   default.units = "cm"
# )

annoHighlight(
  plot = s1,
  chrom = "chr2R",
  chromstart = 13924529,
  chromend = 13925169,
  y = (ref_y), height = 1.5, just = c("left", "top"),
  default.units = "cm"
)


# class III genome browser tracks ----------------------------------------------
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 10.5

# class label
plotText(
  label = "class III", params = large_text_params, fontface = "bold",
  x = (ref_x + x_offset_class_label), y = (ref_y + 0.75), just = c("center"), default.units = "cm",
  rot = 90
)

# class III browser track labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Grh ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

# plot Grh genome browser for class III
region <- pgParams(
  chrom = "chr2R",
  chromstart = 15828720, chromend = 15843156,
  assembly = "dm6"
)

`S2-Grh_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                           params = region
)

`S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(`S2-Grh_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))


s1 <- plotSignal(
  data = `S2-Grh_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-Grh_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
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
#   x = (ref_x + 2.25), y = (ref_y + 1.75), height = gb_height, width = 3,
#   fill = c("gray50", "gray50"),
#   just = c("left", "center"),
#   default.units = "cm"
# )


# Panel C ======================================================================
# panel label
ref_x <- 6.5
ref_y <- 0.5 

# panel label

plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in ChIP classes
zld_chip_classes <- read_tsv(zld_ChIP_classes_fn)

# generate pie charts
zld_class_plot <- zld_chip_classes  %>% 
  # mutate(zygotic_class = factor(zygotic_class, levels = c("i", "ii", "iii"))) %>%
  group_by(class) %>% summarise(n = n()) %>%
  ggplot(aes(x='', y = n, fill = class)) + 
  geom_bar(stat = "identity", position = "fill") + 
  coord_polar("y", direction = 1, start = pi / 2) + 
  # geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 0.35 * small_text_params$fontsize) +
  theme_void() + 
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Blues")


# place pie charts on page
plotGG(
  plot = zld_class_plot,
  x = (ref_x), y = ref_y,
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels
plotText(
  label = paste("class I", "\n", "n =", sum(zld_chip_classes$class == "i")), fontsize = small_text_params$fontsize,
  x = (ref_x + 1), y = (ref_y + 0.5), just = "center", default.units = "cm"
)

plotText(
  label = paste("class II", "\n", "n =", sum(zld_chip_classes$class == "ii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1.5), just = "center", default.units = "cm"
)

plotText(
  label = paste("class III", "\n", "n =", sum(zld_chip_classes$class == "iii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1), just = "center", default.units = "cm"
)


# Panel D ======================================================================
# panel label
ref_x <- 9.5
ref_y <- 0.5 


# generate chart
zld_feature_plot <- zld_chip_classes  %>% 
  ggplot(aes(x=class, fill = feature)) + 
  geom_bar(position = "fill") + 
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_brewer(name = NULL, palette = "Blues") +
  ylab("proportion")


# place chart on page
plotGG(
  plot = zld_feature_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel label

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel E ======================================================================
# panel label
ref_x <- 6.5
ref_y <- 2.75 

# panel label

plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define parameters for heatmaps
bw <- c(
  S2_ZLD_ChIP =  Zld_ChIP_bw,
  S2_WT_ATAC = Zld_WT_ATAC_bw,
  S2_ZLD_ATAC = Zld_Zld_ATAC_bw
)

groups <- c(1,2,2)

regions <- zld_chip_classes %>%
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$class, order_by_samples = 1, 
  scale_group = groups,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

e_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = e_hm,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 5.25, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_2_1")
grid.xaxis(at = c(0, 0.5, 1), label = c(paste0("-",hm_upstream / 1000, "KB"), "peak center", paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.375), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 Zld ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class II", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 1.5), just = c("right", "center"), default.units = "cm"
)

plotText(
  label = "class III", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 2.9), just = c("right","center"), default.units = "cm"
)


# Panel F ======================================================================
# panel label
ref_x <- 6.5
ref_y <- 6.5 

# panel label

plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in ChIP classes
grh_chip_classes <- read_tsv(grh_ChIP_classes_fn)

# generate pie charts
grh_class_plot <- grh_chip_classes  %>% 
  # mutate(zygotic_class = factor(zygotic_class, levels = c("i", "ii", "iii"))) %>%
  group_by(class) %>% summarise(n = n()) %>%
  ggplot(aes(x='', y = n, fill = class, label = n)) + 
  geom_bar(stat = "identity", position = "fill") + 
  coord_polar("y", direction = 1, start = pi / 2) + 
  theme_void() + 
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Oranges") 

# place pie charts on page
plotGG(
  plot = grh_class_plot,
  x = ref_x, y = ref_y,
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm",
)

# add labels
plotText(
  label = paste("class I", "\n", "n =", sum(grh_chip_classes$class == "i")), fontsize = small_text_params$fontsize,
  x = (ref_x + 1), y = (ref_y + 0.5), just = "center", default.units = "cm"
)

plotText(
  label = paste("class II", "\n", "n =", sum(grh_chip_classes$class == "ii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1.75), just = "center", default.units = "cm"
)

plotText(
  label = paste("class III", "\n", "n =", sum(grh_chip_classes$class == "iii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1), just = "center", default.units = "cm"
)


# Panel G ======================================================================
# panel label
ref_x <- 9.5
ref_y <- 6.5 


# generate chart
grh_feature_plot <- grh_chip_classes  %>% 
  ggplot(aes(x=class, fill = feature)) + 
  geom_bar(position = "fill") + 
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_brewer(name = NULL, palette = "Oranges") +
  ylab("proportion")


# place chart on page
plotGG(
  plot = grh_feature_plot,
  x = ref_x, y = ref_y,
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel label

plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel H ======================================================================
# panel label
ref_x <- 6.5
ref_y <- 8.75 

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  S2_ZLD_ChIP =  Grh_ChIP_bw,
  S2_WT_ATAC = Grh_WT_ATAC_bw,
  S2_ZLD_ATAC = Grh_Grh_ATAC_bw
)

groups <- c(1,2,2)

regions <- grh_chip_classes %>%
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_split = regions$class, order_by_samples = 1, 
  scale_group = groups,
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


# Panel I ======================================================================
# panel label
ref_x <- 13
ref_y <- 0.5 
# reference points for positioning figure components
x_offset_brower_label <- 1.25
x_offset_browser <- 2.0

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)
# 
# i_plot <- zld_chip_classes %>% 
#   ggplot(aes(x = class, y = atac_FC)) +
#   geom_violin(fill = zld_color, lwd = 0.1) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, lwd = 0.1) +
#   ylab(expression(ATAC~log[2]~fold~change) ) +
#   theme_classic(base_size = small_text_params$fontsize)
# 
# plotGG(
#   plot = i_plot,
#   x = (ref_x + 0.25), y = (ref_y + 0.5),
#   width = 2.5, height = 2, just = c("left", "top"),
#   default.units = "cm"
# )

# browser track labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Zld ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT RNA-seq", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Zld RNA-seq", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)


# browser tracks
region <- pgParams(
  chrom = "chr3R",
  chromstart = 12644310, chromend = 12653775,
  assembly = "dm6"
)

`S2-Zld_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                           params = region
)

`S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                            params = region
)

`S2-WT_RNAseq` <- readBigwig(file = Zld_WT_RNAseq_bw, 
                             params = region
)

`S2-Zld_RNAseq` <- readBigwig(file = Zld_Zld_RNAseq_bw, 
                              params = region
)

ChIP_range <- signal_range(`S2-Zld_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Zld_ATAC`$score))
RNAseq_range <- signal_range(c(`S2-WT_RNAseq`$score, `S2-Zld_RNAseq`$score))

s1 <- plotSignal(
  data = `S2-Zld_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s3 <- plotSignal(
  data = `S2-Zld_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s4 <- plotSignal(
  data = `S2-WT_RNAseq`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = RNAseq_range
)

annoYaxis(
  plot = s4, at = round(RNAseq_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s5 <- plotSignal(
  data = `S2-Zld_RNAseq`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = RNAseq_range
)

annoYaxis(
  plot = s5, at = round(RNAseq_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

g1 <- plotGenes(
  params = region, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 0.6, width = 2.75, fontsize = small_text_params$fontsize,
  fill = c("gray50", "gray50"),
  just = c("left", "top"),
  default.units = "cm"
)

annoHighlight(
  plot = s1,
  chrom = "chr3R",
  chromstart = 12645060,
  chromend = 12650125,
  y = (ref_y), height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# annoGenomeLabel(
#   plot = g1, scale = "bp", fontsize = small_text_params$fontsize,
#   x = (ref_x + x_offset_browser), y = (ref_y + 2), just = c("left", "top"),
#   default.units = "cm"
# )



# Panel J ======================================================================
# panel label
ref_x <- 13
ref_y <- 3.6 



n_sites <- zld_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% 
  group_by(feature, class, RNAseq_is_diff) %>% 
  summarise(n = n())

n_diff <-  zld_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% pull(RNAseq_is_diff) %>% sum()

percent_diff <- n_diff / nrow(zld_chip_classes)


j_plot <- zld_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% 
  ggplot(aes(x = class, fill = RNAseq_is_diff)) + 
  geom_bar(position = "fill") + 
  geom_hline(yintercept = percent_diff, lty = 2, lwd = 0.5, color = "gray50") +
  geom_text(data = n_sites, aes(label = n, y = n), position = position_fill(vjust=0.5), size = small_text_params$fontsize * 0.35) +
  facet_grid(~feature) + 
  scale_fill_brewer(palette = "Blues") + 
  labs(fill = "RNA-seq DE") +
  ylab("proportion") +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))


# j_plot <- zld_chip_classes %>% 
#   ggplot(aes(x = class, y = RNAseq_FC)) +
#   geom_violin(fill = zld_color, lwd = 0.1) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, lwd = 0.1) +
#   ylab(expression(RNA-seq~log[2]~fold~change) ) +
#   theme_classic(base_size = small_text_params$fontsize)
# 
plotGG(
  plot = j_plot,
  x = (ref_x), y = (ref_y),
  width = 5, height = 2.65, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel K ======================================================================
# panel label
ref_x <- 13
ref_y <- 6.5 
# reference points for positioning figure components
x_offset_brower_label <- 1.25
x_offset_browser <- 2.0

# panel label
plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# browser track labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Grh ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT RNA-seq", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 Grh RNA-seq", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_brower_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
)


# browser tracks
region <- pgParams(
  chrom = "chr3R",
  chromstart = 11766394, chromend = 11783534,
  assembly = "dm6"
)

`S2-Grh_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                           params = region
)

`S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                            params = region
)

`S2-WT_RNAseq` <- readBigwig(file = Grh_WT_RNAseq_bw, 
                             params = region
)

`S2-Grh_RNAseq` <- readBigwig(file = Grh_Grh_RNAseq_bw, 
                              params = region
)

ChIP_range <- signal_range(`S2-Grh_ChIP`$score)
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))
RNAseq_range <- signal_range(c(`S2-WT_RNAseq`$score, `S2-Grh_RNAseq`$score))

s1 <- plotSignal(
  data = `S2-Grh_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s2, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s3 <- plotSignal(
  data = `S2-Grh_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s4 <- plotSignal(
  data = `S2-WT_RNAseq`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = RNAseq_range
)

annoYaxis(
  plot = s4, at = round(RNAseq_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s5 <- plotSignal(
  data = `S2-Grh_RNAseq`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 2.75, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = RNAseq_range
)

annoYaxis(
  plot = s5, at = round(RNAseq_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

g1 <- plotGenes(
  params = region, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 0.6, width = 2.75, fontsize = small_text_params$fontsize, 
  fill = c("gray50", "gray50"),
  just = c("left", "top"),
  default.units = "cm"
)

annoHighlight(
  plot = s1,
  chrom = "chr3R",
  chromstart = 11778754,
  chromend = 11780187,
  y = (ref_y), height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# annoGenomeLabel(
#   plot = g1, scale = "bp", fontsize = small_text_params$fontsize,
#   x = (ref_x + x_offset_browser), y = (ref_y + 2), just = c("left", "top"),
#   default.units = "cm"
# )



# Panel L ======================================================================
# panel label
ref_x <- 13
ref_y <- 9.6

# panel label
plotText(
  label = "l", params = panel_label_params, fontface = "bold",
  x = ref_x, y = (ref_y), just = "bottom", default.units = "cm"
)

n_sites <- grh_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% 
  group_by(feature, class, RNAseq_is_diff) %>% 
  summarise(n = n())

n_diff <-  grh_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% pull(RNAseq_is_diff) %>% sum()

percent_diff <- n_diff / nrow(grh_chip_classes)


l_plot <- grh_chip_classes %>% 
  replace_na(list(RNAseq_is_diff = FALSE)) %>% 
  ggplot(aes(x = class, fill = RNAseq_is_diff)) + 
  geom_bar(position = "fill") + 
  geom_hline(yintercept = percent_diff, lty = 2, lwd = 0.5, color = "gray50") +
  geom_text(data = n_sites, aes(label = n, y = n), position = position_fill(vjust=0.5), size = small_text_params$fontsize * 0.35) +
  facet_grid(~feature) + 
  scale_fill_brewer(palette = "Oranges") + 
  labs(fill = "RNA-seq DE") +
  ylab("proportion") +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = l_plot,
  x = (ref_x), y = (ref_y),
  width = 5, height = 2.65, just = c("left", "top"),
  default.units = "cm"
)
# l_plot <- grh_chip_classes %>% 
#   ggplot(aes(x = class, y = RNAseq_FC)) +
#   geom_violin(fill = grh_color, lwd = 0.1) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, lwd = 0.1) +
#   ylab(expression(RNA-seq~log[2]~fold~change) ) +
#   theme_classic(base_size = small_text_params$fontsize)
# 
# plotGG(
#   plot = l_plot,
#   x = (ref_x + 0.25), y = (ref_y + 0.5),
#   width = 2.5, height = 2, just = c("left", "top"),
#   default.units = "cm"
# )


# close graphics device ========================================================
dev.off()




