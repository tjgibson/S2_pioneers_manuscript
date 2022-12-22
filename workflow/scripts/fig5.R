# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
# library(grImport)
suppressPackageStartupMessages(library(grid))
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# zld_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")
# grh_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")
# twi_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")
# 
# zld_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-0uM_aZld_IP.bw"
# zld_500uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-500uM_aZld_IP.bw"
# zld_1000uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-1000uM_aZld_IP.bw"
# zld_1500uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-1500uM_aZld_IP.bw"
# 
# grh_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-0uM_aGrh_IP.bw"
# grh_25uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-25uM_aGrh_IP.bw"
# grh_100uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-100uM_aGrh_IP.bw"
# grh_400uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-400uM_aGrh_IP.bw"
# 
# twi_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-0uM_aHA_IP.bw"
# twi_10uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-10uM_aHA_IP.bw"
# twi_40uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-40uM_aHA_IP.bw"
# twi_160uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-160uM_aHA_IP.bw"
# 
# zld_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_0uM_small.bw"
# zld_500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_500uM_small.bw"
# zld_1000uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_1000uM_small.bw"
# zld_1500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_1500uM_small.bw"
# 
# grh_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_0uM_small.bw"
# grh_25uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_25uM_small.bw"
# grh_100uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_100uM_small.bw"
# grh_400uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_400uM_small.bw"
# 
# twi_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_0uM_small.bw"
# twi_10uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_10uM_small.bw"
# twi_40uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_40uM_small.bw"
# twi_160uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_160uM_small.bw"
# 
# zld_titration_classes_fn <- "results/ChIP_titration_classes/zld_titration_classes.tsv"
# grh_titration_classes_fn <- "results/ChIP_titration_classes/grh_titration_classes.tsv"
# twi_titration_classes_fn <- "results/ChIP_titration_classes/twi_titration_classes.tsv"

zld_tissue_occupancy_fn <- snakemake@input[["zld_tissue_occupancy"]]
grh_tissue_occupancy_fn <- snakemake@input[["grh_tissue_occupancy"]]
twi_tissue_occupancy_fn <- snakemake@input[["twi_tissue_occupancy"]]

zld_0uM_ChIP_bw <- snakemake@input[["zld_0uM_ChIP_bw"]]
zld_500uM_ChIP_bw <- snakemake@input[["zld_500uM_ChIP_bw"]]
zld_1000uM_ChIP_bw <- snakemake@input[["zld_1000uM_ChIP_bw"]]
zld_1500uM_ChIP_bw <- snakemake@input[["zld_1500uM_ChIP_bw"]]

grh_0uM_ChIP_bw <- snakemake@input[["grh_0uM_ChIP_bw"]]
grh_25uM_ChIP_bw <- snakemake@input[["grh_25uM_ChIP_bw"]]
grh_100uM_ChIP_bw <- snakemake@input[["grh_100uM_ChIP_bw"]]
grh_400uM_ChIP_bw <- snakemake@input[["grh_400uM_ChIP_bw"]]

twi_0uM_ChIP_bw <- snakemake@input[["twi_0uM_ChIP_bw"]]
twi_10uM_ChIP_bw <- snakemake@input[["twi_10uM_ChIP_bw"]]
twi_40uM_ChIP_bw <- snakemake@input[["twi_40uM_ChIP_bw"]]
twi_160uM_ChIP_bw <- snakemake@input[["twi_160uM_ChIP_bw"]]

zld_0uM_ATAC_bw <- snakemake@input[["zld_0uM_ATAC_bw"]]
zld_500uM_ATAC_bw <- snakemake@input[["zld_500uM_ATAC_bw"]]
zld_1000uM_ATAC_bw <- snakemake@input[["zld_1000uM_ATAC_bw"]]
zld_1500uM_ATAC_bw <- snakemake@input[["zld_1500uM_ATAC_bw"]]

grh_0uM_ATAC_bw <- snakemake@input[["grh_0uM_ATAC_bw"]]
grh_25uM_ATAC_bw <- snakemake@input[["grh_25uM_ATAC_bw"]]
grh_100uM_ATAC_bw <- snakemake@input[["grh_100uM_ATAC_bw"]]
grh_400uM_ATAC_bw <- snakemake@input[["grh_400uM_ATAC_bw"]]

twi_0uM_ATAC_bw <- snakemake@input[["twi_0uM_ATAC_bw"]]
twi_10uM_ATAC_bw <- snakemake@input[["twi_10uM_ATAC_bw"]]
twi_40uM_ATAC_bw <- snakemake@input[["twi_40uM_ATAC_bw"]]
twi_160uM_ATAC_bw <- snakemake@input[["twi_160uM_ATAC_bw"]]

zld_titration_classes_fn <- snakemake@input[["zld_titration_classes_fn"]]
grh_titration_classes_fn <- snakemake@input[["grh_titration_classes_fn"]]
twi_titration_classes_fn <- snakemake@input[["twi_titration_classes_fn"]]

## create blank layout for plot =================================================
pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/fig5.pdf")
pageCreate(width = 18.0, height = 18, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
grh_color <- "#F98400"
twi_color <- "#00A08A"

zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")
twi_heatmap_colors <- brewer.pal(9, "GnBu")

# # set genome browser height
# gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500

# read in tissue occupancy classes =============================================
zld_tissue_occupancy <- zld_tissue_occupancy_fn %>% 
  read_tsv()

grh_tissue_occupancy <- grh_tissue_occupancy_fn %>% 
  read_tsv()

twi_tissue_occupancy <- twi_tissue_occupancy_fn %>% 
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

# define bigwig files
bw <- c(
  zld_0uM_ChIP = zld_0uM_ChIP_bw,
  zld_500uM_ChIP = zld_500uM_ChIP_bw,
  zld_1000uM_ChIP = zld_1000uM_ChIP_bw, 
  zld_1500uM_ChIP = zld_1500uM_ChIP_bw
)


# define regions to plot
regions <- zld_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.1,3)
plot_colors <- c(
  `class I` = "#DEEBF7",
  `class II` = "#9ECAE1",
  `class III` = "#3182BD",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
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
  width = 2, height = 2, just = c("left", "top"),
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
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
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
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
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
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ChIP-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

plotText(
  label = "0 µM", params = large_text_params,
  x = (ref_x + 1.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "500 µM", params = large_text_params,
  x = (ref_x + 3.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "1000 µM", params = large_text_params,
  x = (ref_x + 6.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "1500 µM", params = large_text_params,
  x = (ref_x + 8.75), y = (ref_y), just = "center", default.units = "cm"
)


# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 3

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  zld_0uM_ATAC = zld_0uM_ATAC_bw,
  zld_500uM_ATAC = zld_500uM_ATAC_bw,
  zld_1000uM_ATAC = zld_1000uM_ATAC_bw, 
  zld_1500uM_ATAC = zld_1500uM_ATAC_bw
)


regions <- zld_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.3,3)

plot_colors <- c(
  `class I` = "#DEEBF7",
  `class II` = "#9ECAE1",
  `class III` = "#3182BD",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
)

metaplot_5 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_5,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_6 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_6,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_7 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_7,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_8 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_8,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ATAC-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

# legend for A and B ===========================================================
plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 10.25, y = 1.5, width =1, height = 2.5,
  just = c("left", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 0.5


# generate plot
zld_class_plot <- zld_titration_classes_fn |> 
  read_tsv() |> 
    ggplot(aes(x = CuSO4, fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_brewer(palette = "Blues") +
  ylab("n peaks")

# place chart on page
plotGG(
  plot = zld_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel D ======================================================================
library(eulerr)

# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 3

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define input files
input_files <- c(
  zld_0uM = "ChIPseq/results/peaks/filtered/S2-Zld-0uM_aZld_IP.narrowPeak",
  zld_500uM = "ChIPseq/results/peaks/filtered/S2-Zld-500uM_aZld_IP.narrowPeak",
  zld_1000uM = "ChIPseq/results/peaks/filtered/S2-Zld-1000uM_aZld_IP.narrowPeak",
  zld_1500uM = "ChIPseq/results/peaks/filtered/S2-Zld-1500uM_aZld_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


zld_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = zld_euler_plot,
  x = (ref_x), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5.5

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define bigwig files
bw <- c(
  grh_0uM_ChIP = grh_0uM_ChIP_bw,
  grh_25uM_ChIP = grh_25uM_ChIP_bw,
  grh_100uM_ChIP = grh_100uM_ChIP_bw, 
  grh_400uM_ChIP = grh_400uM_ChIP_bw
)


# define regions to plot
regions <- grh_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.3,8)
plot_colors <- c(
  `class I` = "#FEE6CE",
  `class II` = "#FDAE6B",
  `class III` = "#E6550D",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
)

metaplot_9 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_9,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_10 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_10,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

 # plot metaplot 3
metaplot_11 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_11,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_12 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_12,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ChIP-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

plotText(
  label = "0 µM", params = large_text_params,
  x = (ref_x + 1.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "25 µM", params = large_text_params,
  x = (ref_x + 3.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "100 µM", params = large_text_params,
  x = (ref_x + 6.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "400 µM", params = large_text_params,
  x = (ref_x + 8.75), y = (ref_y), just = "center", default.units = "cm"
)


# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 8

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  grh_0uM_ATAC = grh_0uM_ATAC_bw,
  grh_25uM_ATAC = grh_25uM_ATAC_bw,
  grh_100uM_ATAC = grh_100uM_ATAC_bw, 
  grh_400uM_ATAC = grh_400uM_ATAC_bw
)


regions <- grh_tissue_occupancy %>%
  # filter(class == "class II" | class == "class III") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# plot metaplot 1
plot_range <- c(-0.3,6)
plot_colors <- c(
  `class I` = "#FEE6CE",
  `class II` = "#FDAE6B",
  `class III` = "#E6550D",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
)

metaplot_13 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_13,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_14 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_14,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_15 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_15,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_16 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_16,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ATAC-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

# legend for E and F ===========================================================
plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 10.25, y = 6.5, width =1, height = 2.5,
  just = c("left", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 5.5


# generate plot
grh_class_plot <- grh_titration_classes_fn |> 
  read_tsv() |> 
  ggplot(aes(x = as.factor(CuSO4), fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_brewer(palette = "Oranges") +
  ylab("n peaks")

# place chart on page
plotGG(
  plot = grh_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# panel H ======================================================================

# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 8

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define input files
input_files <- c(
  grh_0uM = "ChIPseq/results/peaks/filtered/S2-Grh-0uM_aGrh_IP.narrowPeak",
  grh_25uM = "ChIPseq/results/peaks/filtered/S2-Grh-25uM_aGrh_IP.narrowPeak",
  grh_100uM = "ChIPseq/results/peaks/filtered/S2-Grh-100uM_aGrh_IP.narrowPeak",
  grh_400uM = "ChIPseq/results/peaks/filtered/S2-Grh-400uM_aGrh_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


grh_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = grh_euler_plot,
  x = (ref_x), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel I ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 10.5

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define bigwig files
bw <- c(
  twi_0uM_ChIP = twi_0uM_ChIP_bw,
  twi_10uM_ChIP = twi_10uM_ChIP_bw,
  twi_40uM_ChIP = twi_40uM_ChIP_bw, 
  twi_160uM_ChIP = twi_160uM_ChIP_bw
)


# define regions to plot
regions <- twi_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# plot metaplot 1
plot_range <- c(-0.3,8)
plot_colors <- c(
  `class I` =  "#E0F3DB",
  `class II` = "#A8DDB5",
  `class III` = "#43A2CA",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
)

metaplot_17 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_17,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_18 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_18,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_19 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_19,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_20 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
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
  plot = metaplot_20,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ChIP-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

plotText(
  label = "0 µM", params = large_text_params,
  x = (ref_x + 1.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "10 µM", params = large_text_params,
  x = (ref_x + 3.75), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "40 µM", params = large_text_params,
  x = (ref_x + 6.25), y = (ref_y), just = "center", default.units = "cm"
)

plotText(
  label = "160 µM", params = large_text_params,
  x = (ref_x + 8.75), y = (ref_y), just = "center", default.units = "cm"
)

# panel J ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 13

# panel label
plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  twi_0uM_ATAC = twi_0uM_ATAC_bw,
  twi_10uM_ATAC = twi_10uM_ATAC_bw,
  twi_40uM_ATAC = twi_40uM_ATAC_bw, 
  twi_160uM_ATAC = twi_160uM_ATAC_bw
)


regions <- twi_tissue_occupancy %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


plot_range <- c(-0.3,3)
plot_colors <- c(
  `class I` =  "#E0F3DB",
  `class II` = "#A8DDB5",
  `class III` = "#43A2CA",
  repressed_H3K27me3 = "#D9D9D9",
  repressed_H3K9me3 = "#969696",
  repressed_other = "#737373"
)

metaplot_21 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_21,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 2
metaplot_22 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
        
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_22,
  x = (ref_x + 2.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 3
metaplot_23 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_23,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_24 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "none"
  ) +
  ylim(plot_range)

plotGG(
  plot = metaplot_24,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# add labels to plots
plotText(
  label = "ATAC-seq signal", params = small_text_params,
  x = ref_x, y = (ref_y + 1.25), just = "center", default.units = "cm", rot = 90
)

# legend for I and J ===========================================================
plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = 10.25, y = 11.5, width =1, height = 2.5,
  just = c("left", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

# panel K ======================================================================
# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 10.5


# generate plot
twi_class_plot <- twi_titration_classes_fn |> 
  read_tsv() |> 
  ggplot(aes(x = as.factor(CuSO4), fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_brewer(palette = "GnBu") +
  ylab("n peaks")

# place chart on page
plotGG(
  plot = twi_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel L ======================================================================

# reference points for positioning figure components
ref_x <- 13.5
ref_y <- 13

# panel label
plotText(
  label = "l", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define input files
input_files <- c(
  twi_0uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-0uM_aHA_IP.narrowPeak",
  twi_10uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-10uM_aHA_IP.narrowPeak",
  twi_40uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-40uM_aHA_IP.narrowPeak",
  twi_160uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-160uM_aHA_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


twi_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = twi_euler_plot,
  x = (ref_x), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

dev.off()
