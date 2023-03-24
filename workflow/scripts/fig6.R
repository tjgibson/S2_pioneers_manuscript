# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
library(grImport)
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# define input files explicitly for interactive testing
# Zld_FL_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# Zld_DBD_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-DBD_aZld_IP.bw"
# Zld_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw"
# Zld_Zld_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw"
# 
# Grh_FL_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# Grh_DBD_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-DBD_aGrh_IP.bw"
# Grh_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw"
# Grh_Grh_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw"
# 
# zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# 
# zld_FL_atac_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Zld_ATACseq_S2-Zld-FL-vs-S2-WT_results.tsv"
# zld_DBD_atac_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Zld-DBD_ATACseq_S2-Zld-DBD-vs-S2-WT_results.tsv"
# grh_FL_atac_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Grh_ATACseq_S2-Grh-FL-vs-S2-WT_results.tsv"
# grh_DBD_atac_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Grh-DBD_ATACseq_S2-Grh-DBD-vs-S2-WT_results.tsv"

# get input files from snakemake
Zld_FL_ChIP_bw <- snakemake@input[["Zld_FL_ChIP_bw"]]
Zld_DBD_ChIP_bw <- snakemake@input[["Zld_DBD_ChIP_bw"]]
Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]

Grh_FL_ChIP_bw <- snakemake@input[["Grh_FL_ChIP_bw"]]
Grh_DBD_ChIP_bw <- snakemake@input[["Grh_DBD_ChIP_bw"]]
Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
Grh_Grh_ATAC_bw <- snakemake@input[["Grh_Grh_ATAC_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

zld_FL_atac_results_fn <- snakemake@input[["zld_FL_atac_results"]]
zld_DBD_atac_results_fn <- snakemake@input[["zld_DBD_atac_results"]]
grh_FL_atac_results_fn <- snakemake@input[["grh_FL_atac_results"]]
grh_DBD_atac_results_fn <- snakemake@input[["grh_DBD_atac_results"]]
  

## create blank layout for plot =================================================
# define figure dimensions in cm
fig_width <-  18
fig_height <- 14

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
zld_heatmap_colors <- brewer.pal(9, "Blues")

grh_color <- "#F98400"
grh_heatmap_colors <- brewer.pal(9, "Oranges")

# reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_browser_label <- 1.5
x_offset_browser <- 2.25

# set genome browser height
gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500

# panel A ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 0.5

plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Zelda broswer tracks 
plotText(
  label = "Zld FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "Zld DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "ZLD FL ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)

region <- pgParams(
  chrom = "chr3R",
  chromstart = 6750296, chromend = 6761204,
  assembly = "dm6"
)

`Zld_FL_ChIP` <- readBigwig(file = Zld_FL_ChIP_bw, 
                            params = region
)

`Zld_DBD_ChIP` <- readBigwig(file = Zld_DBD_ChIP_bw, 
                            params = region
)

`S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                           params = region
)

`S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(c(Zld_FL_ChIP$score, Zld_DBD_ChIP$score))
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Zld_ATAC`$score))

s1 <- plotSignal(
  data = `Zld_FL_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s2 <- plotSignal(
  data = `Zld_DBD_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s3 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)



s4 <- plotSignal(
  data = `S2-Zld_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
  range = ATAC_range
)

annoYaxis(
  plot = s4, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

plotGenes(
  params = region, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2), height = 0.75, width = 6, fontsize = small_text_params$fontsize,
  fill = c("gray50", "gray50"),
  fontcolor = c("gray50", "gray50"),
  just = c("left", "top"),
  default.units = "cm"
)

annoHighlight(
  plot = s1,
  chrom = "chr3R",
  chromstart = 6755358,
  chromend = 6755908,
  y = (ref_y), height = 2, just = c("left", "top"),
  default.units = "cm"
)


# panel B ======================================================================
# panel label
ref_x <- 9.25
ref_y <- 0.5

plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Zelda broswer tracks 
plotText(
  label = "Grh FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 0.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "Grh DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
)

plotText(
  label = "Grh FL ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
)

region <- pgParams(
  chrom = "chr2R",
  chromstart = 15830447, chromend = 15853445,
  assembly = "dm6"
)

`Grh_FL_ChIP` <- readBigwig(file = Grh_FL_ChIP_bw, 
                            params = region
)

`Grh_DBD_ChIP` <- readBigwig(file = Grh_DBD_ChIP_bw, 
                             params = region
)

`S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                           params = region
)

`S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                            params = region
)

ChIP_range <- signal_range(c(Grh_FL_ChIP$score, Grh_DBD_ChIP$score))
ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))

s1 <- plotSignal(
  data = `Grh_FL_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s1, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s2 <- plotSignal(
  data = `Grh_DBD_ChIP`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ChIP_range
)

annoYaxis(
  plot = s2, at = round(ChIP_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

s3 <- plotSignal(
  data = `S2-WT_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s3, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)


s4 <- plotSignal(
  data = `S2-Grh_ATAC`, params = region,
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 6, height = gb_height,
  just = c("left", "center"), default.units = "cm",
  linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
  range = ATAC_range
)

annoYaxis(
  plot = s4, at = round(ATAC_range, 1),
  axisLine = TRUE, fontsize = 5, lwd = 0.5, 
)

plotGenes(
  params = region, assembly = "dm6",
  x = (ref_x + x_offset_browser), y = (ref_y + 2), height = 0.75, width = 6, fontsize = small_text_params$fontsize,
  fill = c("gray50", "gray50"),
  fontcolor = c("gray50", "gray50"),
  just = c("left", "top"),
  default.units = "cm"
)

annoHighlight(
  plot = s1,
  chrom = "chr2R",
  chromstart = 15834134,
  chromend = 15836230,
  y = (ref_y), height = 2, just = c("left", "top"),
  default.units = "cm"
)


# Panel C ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 3.5 

# panel label

plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

bw <- c(
  ZLD_FL_ChIP =  Zld_FL_ChIP_bw,
  Zld_DBD_ChIP = Zld_DBD_ChIP_bw
)


zld_chip_classes <- read_tsv(zld_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

regions <- zld_chip_classes |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

plot_range <- c(0,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.15),
        legend.position = "none", 
        axis.title.y = element_blank(),
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2.5, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2.5, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add labels to metaplots
plotText(
  label = "Zld FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.5), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "Zld DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.25), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "z-score-normalized signal", params = small_text_params, rot = 90,
  x = (ref_x + 0.25), y = (ref_y + 1.25), just = c("center","center"), default.units = "cm"
)


# add legend
plot_colors <- c(
  I = "#DEEBF7",
  II = "#9ECAE1",
  III = "#3182BD"
)


plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 3.5, y = 6.25, width =3, height = 0.5,
  just = c("center", "center"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "h"
)

# Panel D ======================================================================
# panel label
ref_x <- 6
ref_y <- 3.5 


# import ATAC-seq results
zld_FL_atac_results <- read_tsv(zld_FL_atac_results_fn)
zld_DBD_atac_results <- read_tsv(zld_DBD_atac_results_fn)


plot_range <- c(-log10(zld_FL_atac_results$padj),-log10(zld_DBD_atac_results$padj)) |> na.omit() |> range()

d_plot <-
  zld_DBD_atac_results |> 
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = is_diff)) +
  geom_point(alpha = 0.7, size = 0.1) +
  theme_classic(base_size = 5) +
  theme( legend.position = "none") +
  scale_color_manual(values = c("grey", "black")) +
  ylim(plot_range) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))

plotGG(
  plot = d_plot,
  x = (ref_x), y = (ref_y + 0.4),
  width = 2.75, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel label

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)



# Panel E ======================================================================
# panel label
ref_x <- 9.25
ref_y <- 3.5 

# panel label

plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

bw <- c(
  Grh_FL_ChIP =  Grh_FL_ChIP_bw,
  Grh_DBD_ChIP = Grh_DBD_ChIP_bw
)


grh_chip_classes <- read_tsv(grh_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

regions <- grh_chip_classes |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

plot_range <- c(0,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2.5, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2.5, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add labels to metaplots
plotText(
  label = "Grh FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.5), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "Grh DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.25), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "z-score-normalized signal", params = small_text_params, rot = 90,
  x = (ref_x + 0.25), y = (ref_y + 1.25), just = c("center","center"), default.units = "cm"
)


# add legend
plot_colors <- c(
  I = "#FEE6CE",
  II = "#FDAE6B",
  III = "#E6550D"
)


plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 12.25, y = 6.25, width =3, height = 0.5,
  just = c("center", "center"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "h"
)

# Panel F ======================================================================
# panel label
ref_x <- 14.75
ref_y <- 3.5 

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


grh_FL_atac_results <- read_tsv(grh_FL_atac_results_fn)
grh_DBD_atac_results <- read_tsv(grh_DBD_atac_results_fn)

plot_range <- c(-log10(grh_FL_atac_results$padj),-log10(grh_DBD_atac_results$padj)) |> na.omit() |> range()

f_plot <-
  grh_DBD_atac_results |> 
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = is_diff)) +
  geom_point(alpha = 0.7, size = 0.1) +
  theme_classic(base_size = 5) +
  theme( legend.position = "none") +
  scale_color_manual(values = c("grey", "black")) +
  ylim(plot_range) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))


plotGG(
  plot = f_plot,
  x = (ref_x), y = (ref_y + 0.4),
  width = 2.75, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 6.75

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# # get model for EPS file
# PostScriptTrace(
#   "manuscript/figures/models/model.eps",
#   outfilename = "manuscript/figures/models/model.eps.xml",
#   charpath = FALSE,
#   charpos = FALSE
# )
# schematic <- readPicture("manuscript/figures/models/model.eps.xml")
# schematic_gtree <- grid.grabExpr(grid.picture(schematic))
# 
# 
# # plot model on page
# plotGG(schematic_gtree, x = ref_x, y = ref_y, width = 17, height = 6.5, default.units = "cm")

# placeholder for model
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 17, height = 6.5, default.units = "cm", just = c("top","left"))

# close graphics device ========================================================
dev.off()
