# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)

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
pdf(snakemake@output[[1]], useDingbats = FALSE)
pageCreate(width = 18.3, height = 12.5, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
twi_color <- "#00A08A"
twi_heatmap_colors <- c("white","#00A08A", "#154734")

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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3, height = gb_height,
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
  x = (ref_x + x_offset_browser), y = (ref_x + 5.25), width = 3, height = gb_height,
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
ref_x <- 6.5
ref_y <- 0.5 

# panel label

plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in ChIP classes
twi_chip_classes <- read_tsv(twi_ChIP_classes_fn)

# generate pie charts
twi_class_plot <- twi_chip_classes  %>% 
  # mutate(zygotic_class = factor(zygotic_class, levels = c("i", "ii", "iii"))) %>%
  group_by(class) %>% summarise(n = n()) %>%
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
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels
plotText(
  label = paste("class I", "\n", "n =", sum(twi_chip_classes$class == "i")), fontsize = small_text_params$fontsize,
  x = (ref_x + 1), y = (ref_y + 0.5), just = "center", default.units = "cm"
)

plotText(
  label = paste("class II", "\n", "n =", sum(twi_chip_classes$class == "ii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1.5), just = "center", default.units = "cm"
)

plotText(
  label = paste("class III", "\n", "n =", sum(twi_chip_classes$class == "iii")), fontsize = small_text_params$fontsize,
  x = (ref_x + 2.2), y = (ref_y + 1), just = "center", default.units = "cm"
)


# Panel C ======================================================================
# panel label
ref_x <- 9.5
ref_y <- 0.5 

# generate chart
twi_feature_plot <- twi_chip_classes  %>% 
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
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel label

plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel D ======================================================================
# panel label
ref_x <- 6.5
ref_y <- 2.75 


# panel label

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

bw <- c(
  S2_Twi_ChIP =  Twi_ChIP_bw,
  S2_WT_ATAC = Twi_WT_ATAC_bw,
  S2_Twi_ATAC = Twi_Twi_ATAC_bw
)

groups <- c(1,2,2)

regions <- twi_chip_classes %>%
  filter(class != "i") %>% 
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


plotGG(
  plot = d_hm,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 5.25, height = 3, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_2_1")
grid.xaxis(at = c(0, 0.5, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), "peak center", paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.375), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "S2 Twi ATAC", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class II", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 2), just = c("right", "center"), default.units = "cm"
)

plotText(
  label = "class III", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.5), y = (ref_y + 3.15), just = c("right","center"), default.units = "cm"
)



# # Panel E ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 6.75

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# e_plot <- twi_chip_classes %>%
#   ggplot(aes(x = class, y = atac_FC)) +
#   geom_violin(fill = twi_color, lwd = 0.1) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, lwd = 0.1) +
#   ylab(expression(ATAC~log[2]~fold~change) ) +
#   theme_classic(base_size = small_text_params$fontsize)
# 
# plotGG(
#   plot = e_plot,
#   x = ref_x, y = ref_y,
#   width = 2.5, height = 2, just = c("left", "top"),
#   default.units = "cm"
# )

zld_chip_classes <- read_tsv(zld_ChIP_classes_fn)
grh_chip_classes <- read_tsv(grh_ChIP_classes_fn)

zld_closed <- zld_chip_classes %>% 
  filter(class != "i") %>% 
  group_by(class) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  add_column(factor = "Zld")

grh_closed <- grh_chip_classes %>% 
  filter(class != "i") %>% 
  group_by(class) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  add_column(factor = "Grh")

twi_closed <- twi_chip_classes %>% 
  filter(class != "i") %>% 
  group_by(class) %>% 
  summarise(n = n()) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  add_column(factor = "Twi")

p_data <- bind_rows(zld_closed, grh_closed, twi_closed) %>% 
  filter(class == "iii")

e_plot <- p_data %>% 
  mutate(factor = factor(factor, levels = c("Twi", "Grh", "Zld"))) %>% 
  ggplot(aes(x = factor, y = percent, fill = factor)) + geom_bar(stat = "identity") +
  theme_classic(base_size = 5) + 
  theme(
        legend.position = "none",
        axis.title.x = element_blank()) +
  scale_fill_manual(values=c(twi_color, grh_color, zld_color)) +
  ylab("% closed sites opened")

plotGG(
  plot = e_plot,
  x = ref_x, y = (ref_y + 0.25),
  width = 2.5, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# # Panel F ======================================================================
# panel label
ref_x <- 3.25
ref_y <- 6.75

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# f_plot <- twi_chip_classes %>%
#   ggplot(aes(x = class, y = RNAseq_FC)) +
#   geom_violin(fill = twi_color, lwd = 0.1) +
#   geom_boxplot(width = 0.1, outlier.shape = NA, lwd = 0.1) +
#   ylab(expression(RNA-seq~log[2]~fold~change) ) +
#   theme_classic(base_size = small_text_params$fontsize)
# 
# plotGG(
#   plot = f_plot,
#   x = ref_x, y = ref_y,
#   width = 2.5, height = 2, just = c("left", "top"),
#   default.units = "cm"
# )
# 

# define bigwig files to use
bw <- c(
  Twi_ChIP_bw,
  Grh_ChIP_bw,
  Zld_ChIP_bw
)

# define regions to use for metaplots
twi_regions <- twi_chip_classes %>% 
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

grh_regions <- grh_chip_classes %>% 
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

zld_regions <- zld_chip_classes %>% 
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



metaplot_1 <- plot_average(bw[1], regions = twi_regions, row_split = twi_regions$class, line_width = 0.2) +
  scale_color_manual(values=c("#5FC3B5", "#154734")) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
        )

plotGG(
  plot = metaplot_1,
  x = ref_x, y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

metaplot_2 <- plot_average(bw[2], regions = grh_regions, row_split = grh_regions$class, line_width = 0.2) +
  scale_color_manual(values=c("#FDAE6B", "#E6550D")) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  )

plotGG(
  plot = metaplot_2,
  x = ref_x + 2.5, y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

metaplot_3 <- plot_average(bw[3], regions = zld_regions, row_split = zld_regions$class, line_width = 0.2) +
  scale_color_manual(values=c("#9ECAE1", "#3182BD")) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0)
  )

plotGG(
  plot = metaplot_3,
  x = ref_x + 5, y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add titles to metaplots
plotText(
  label = "Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.6), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 6.1), y = (ref_y), just = c("center"), default.units = "cm"
)




# # Panel G ======================================================================
# panel label
ref_x <- 10.75
ref_y <- 6.75

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# import RNA-seq results and reformat for plotting
twi_RNAseq_results <- read_tsv(twi_RNAseq_results_fn) %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  add_column(factor = "Twi") %>% 
  left_join(dplyr::select(twi_chip_classes, gene_id, peak_id), by = "gene_id") %>% 
  mutate(bound_by_factor = !is.na(peak_id)) %>% 
  dplyr::select(gene_id, factor,bound_by_factor) %>% 
  distinct(.keep_all = TRUE)

zld_RNAseq_results <- read_tsv(zld_RNAseq_results_fn) %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  add_column(factor = "Zld") %>% 
  left_join(dplyr::select(zld_chip_classes, gene_id, peak_id), by = "gene_id") %>% 
  mutate(bound_by_factor = !is.na(peak_id)) %>% 
  dplyr::select(gene_id, factor,bound_by_factor) %>% 
  distinct(.keep_all = TRUE)

grh_RNAseq_results <- read_tsv(grh_RNAseq_results_fn) %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  add_column(factor = "Grh") %>% 
  left_join(dplyr::select(grh_chip_classes, gene_id, peak_id), by = "gene_id") %>% 
  mutate(bound_by_factor = !is.na(peak_id)) %>% 
  dplyr::select(gene_id, factor,bound_by_factor) %>% 
  distinct(.keep_all = TRUE)

p_data <- bind_rows(twi_RNAseq_results, zld_RNAseq_results, grh_RNAseq_results)

h_plot <- p_data %>% 
  ggplot(aes(x=factor, fill = bound_by_factor)) + geom_bar(position = "stack") +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_fill_manual(values=c("#DFF3F0", "#5FC3B5")) +
  ylab("n upregulated genes") +
  theme(legend.key.size = unit(2, 'mm'), 
        axis.title.x = element_blank(),
        # legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0),
        plot.margin = margin(0,0,0,0))

plotGG(
  plot = h_plot,
  x = ref_x, y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)


# close graphics device ========================================================
dev.off()