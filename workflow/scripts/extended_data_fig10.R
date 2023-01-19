# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(plyranges)
library(rtracklayer)

source("workflow/scripts/plot_heatmap.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]

zld_FL_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
zld_DBD_peaks_fn <- "ChIPseq/results/peaks/filtered/S2-Zld-DBD_aZld_IP.narrowPeak"

grh_FL_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
grh_DBD_peaks_fn <- "ChIPseq/results/peaks/filtered/S2-Grh-DBD_aGrh_IP.narrowPeak"

zld_FL_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
zld_DBD_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-DBD_aZld_IP.bw"

grh_FL_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
grh_DBD_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-DBD_aGrh_IP.bw"


# impor ChIP classes ===========================================================
zld_ChIP_classes <- zld_FL_ChIP_classes_fn |> 
  read_tsv()

grh_ChIP_classes <- grh_FL_ChIP_classes_fn |> 
  read_tsv()

# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig10.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 12, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")


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

zld_FL_peaks <- read_tsv(zld_FL_ChIP_classes_fn) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


regions <- zld_FL_peaks |> 
  # plyranges::filter(class != "i") |> 
  plyranges::filter(n_motifs > 0)

bw <- c(
  Zld_FL = zld_FL_bw,
  Zld_DBD = zld_DBD_bw
)


# plot metaplot 1
plot_colors <- c(
  i = "#DEEBF7",
  ii = "#9ECAE1",
  iii = "#3182BD"
)


plot_range <- c(0,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_manual(values = plot_colors) +
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
  scale_color_manual(values = plot_colors) +
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





dev.off()
