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
zld_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")
grh_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")
twi_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")

zld_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-0uM_aZld_IP.bw"
zld_500uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-500uM_aZld_IP.bw"
zld_1000uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-1000uM_aZld_IP.bw"
zld_1500uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-1500uM_aZld_IP.bw"

grh_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-0uM_aGrh_IP.bw"
grh_25uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-25uM_aGrh_IP.bw"
grh_100uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-100uM_aGrh_IP.bw"
grh_400uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-400uM_aGrh_IP.bw"

twi_0uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-0uM_aHA_IP.bw"
twi_10uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-10uM_aHA_IP.bw"
twi_40uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-40uM_aHA_IP.bw"
twi_160uM_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi-160uM_aHA_IP.bw"

zld_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_0uM_small.bw"
zld_500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_500uM_small.bw"
zld_1000uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_1000uM_small.bw"
zld_1500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Zld_1500uM_small.bw"

grh_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_0uM_small.bw"
grh_500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_25uM_small.bw"
grh_1000uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_100uM_small.bw"
grh_1500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-Grh_400uM_small.bw"

twi_0uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_0uM_small.bw"
twi_500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_10uM_small.bw"
twi_1000uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_40uM_small.bw"
twi_1500uM_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/titration_S2-HA-Twi_160uM_small.bw"


## create blank layout for plot =================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/fig5.pdf")
pageCreate(width = 18.0, height = 12, default.units = "cm", showGuides = TRUE)

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

# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


bw <- c(
  zld_0uM_ChIP = zld_0uM_ChIP_bw,
  zld_500uM_ChIP = zld_500uM_ChIP_bw,
  zld_1000uM_ChIP = zld_1000uM_ChIP_bw, 
  zld_1500uM_ChIP = zld_1500uM_ChIP_bw
)


regions <- zld_tissue_occupancy %>%
  filter(class != "class I") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions,
  upstream = hm_upstream, downstream = hm_downstream,
  colors  = zld_heatmap_colors,
  row_split = regions$class,
  # order_by_samples = 1,
  individual_scales = FALSE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# # add axes to heatmaps
# seekViewport(name = "matrix_1_heatmap_body_6_1")
# grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
# seekViewport(name = "page")
# 
# # add heatmap labels
# # plotText(
# #   label = paste(length(regions), "overlapping", "\n", "class I sites"), params = small_text_params, fontface = "bold",
# #   x = (ref_x), y = (ref_y + 1.5), just = c("center"), default.units = "cm", rot = 90
# # )
# 
# plotText(
#   label = paste0("Zld", "\n", "S2 cells"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 0.68), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = paste0("Zld", "\n", "embryo"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 1.7), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = paste0("Zld", "\n", "brain"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 2.71), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = "H3K27me3", params = small_text_params, fontface = "bold",
#   x = (ref_x + 3.74), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = "H3K9me3", params = small_text_params, fontface = "bold",
#   x = (ref_x + 4.8), y = (ref_y), just = c("center"), default.units = "cm"
# )

# plot metaplot 1
plot_range <- c(-0.1,3)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
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
  scale_color_brewer(palette = "Blues") +
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
  scale_color_brewer(palette = "Blues") +
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
  plot = metaplot_3,
  x = (ref_x + 5.25), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# plot metaplot 4
metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
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
  plot = metaplot_4,
  x = (ref_x + 7.75), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 4

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
  # filter(class == "class II" | class == "class III") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# hm <- plot_heatmap_minimal(
#   bw, regions, 
#   upstream = hm_upstream, downstream = hm_downstream, 
#   colors  = zld_heatmap_colors, 
#   row_split = regions$class, 
#   # order_by_samples = 1, 
#   individual_scales = FALSE,
#   use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
#   border = "black",
#   border_gp = gpar(lwd = 0.5),
#   show_heatmap_legend = FALSE, 
#   return_heatmap_list = TRUE
# )
# 
# b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))
# 
# # place heatmap on page
# plotGG(
#   plot = b_hm,
#   x = (ref_x + 0.25), y = (ref_y + 0.25),
#   width = 5, height = 5, just = c("left", "top"),
#   default.units = "cm"
# )
# 
# # add axes to heatmaps
# seekViewport(name = "matrix_1_heatmap_body_6_1")
# grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
# seekViewport(name = "page")
# 
# # add heatmap labels
# # plotText(
# #   label = paste(length(regions), "overlapping", "\n", "class I sites"), params = small_text_params, fontface = "bold",
# #   x = (ref_x), y = (ref_y + 1.5), just = c("center"), default.units = "cm", rot = 90
# # )
# 
# plotText(
#   label = paste0("Zld", "\n", "S2 cells"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 0.68), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = paste0("Zld", "\n", "embryo"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 1.7), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = paste0("Zld", "\n", "brain"), params = small_text_params, fontface = "bold",
#   x = (ref_x + 2.71), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = "H3K27me3", params = small_text_params, fontface = "bold",
#   x = (ref_x + 3.74), y = (ref_y), just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = "H3K9me3", params = small_text_params, fontface = "bold",
#   x = (ref_x + 4.8), y = (ref_y), just = c("center"), default.units = "cm"
# )

# plot metaplot 1
plot_range <- c(-0.1,3)

metaplot_5 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
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
  scale_color_brewer(palette = "Blues") +
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
  scale_color_brewer(palette = "Blues") +
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
  scale_color_brewer(palette = "Blues") +
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



dev.off()
