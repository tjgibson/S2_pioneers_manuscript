# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
# library(grImport)
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
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
fig_width <-  16
fig_height <- 6

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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 2), height = 0.75, width = 5, fontsize = small_text_params$fontsize,
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
ref_x <- 8.5
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 5, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 2), height = 0.75, width = 5, fontsize = small_text_params$fontsize,
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


zld_chip_classes <- read_tsv(zld_ChIP_classes_fn)

regions <- zld_chip_classes %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

plot_range <- c(0,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 2.5), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to metaplots
plotText(
  label = "Zld FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "Zld DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.5), y = (ref_y), just = c("center","center"), default.units = "cm"
)

# Panel D ======================================================================
# panel label
ref_x <- 5.5
ref_y <- 3.5 

# panel label

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import ATAC-seq results
zld_FL_atac_results <- read_tsv(zld_FL_atac_results_fn)
zld_DBD_atac_results <- read_tsv(zld_DBD_atac_results_fn)


zld_FL_range <- -log10(zld_FL_atac_results$padj) %>% na.omit() %>% range()

d_plot <-
  zld_DBD_atac_results %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = is_diff)) +
  geom_point(alpha = 0.7, size = 0.1) +
  theme_classic(base_size = 5) +
  theme( legend.position = "none") +
  scale_color_manual(values = c("grey", "black")) +
  ylim(zld_FL_range)

plotGG(
  plot = d_plot,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)



# Panel E ======================================================================
# panel label
ref_x <- 8.5
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


grh_chip_classes <- read_tsv(grh_ChIP_classes_fn)

regions <- grh_chip_classes %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

plot_range <- c(0,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 2.5), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to metaplots
plotText(
  label = "Grh FL ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = c("center","center"), default.units = "cm"
)

plotText(
  label = "Grh DBD ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.5), y = (ref_y), just = c("center","center"), default.units = "cm"
)


# Panel F ======================================================================
# panel label
ref_x <- 13.5
ref_y <- 3.5 

# panel label

plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


grh_FL_atac_results <- read_tsv(grh_FL_atac_results_fn)
grh_DBD_atac_results <- read_tsv(grh_DBD_atac_results_fn)

grh_FL_range <- -log10(grh_FL_atac_results$padj) %>% na.omit() %>% range()

f_plot <-
  grh_DBD_atac_results %>% 
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = is_diff)) +
  geom_point(alpha = 0.7, size = 0.1) +
  theme_classic(base_size = 5) +
  theme( legend.position = "none") +
  scale_color_manual(values = c("grey", "black")) +
  ylim(grh_FL_range)

plotGG(
  plot = f_plot,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)




# close graphics device ========================================================
dev.off()
