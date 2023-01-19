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

S2_Zld_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
S2_Grh_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
S2_Twi_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
MNase_seq_bw <-  "data/published_datasets/Chereji_2019/GSE128689_RAW/GSM3934479_S2_exp2_seq1_40min.bw"

zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
twi_ChIP_classes_fn <- "results/ChIP_peak_classes/twi_ChIP_classes.tsv"

zld_motif_instances_fn <- "results/motif_instances/zld_motifs.bed"
grh_motif_instances_fn <- "results/motif_instances/grh_motifs.bed"
twi_motif_instances_fn <- "results/motif_instances/twi_motifs.bed"

# impor ChIP classes ===========================================================
zld_ChIP_classes <- zld_ChIP_classes_fn |> 
  read_tsv()

grh_ChIP_classes <- grh_ChIP_classes_fn |> 
  read_tsv()

twi_ChIP_classes <- twi_ChIP_classes_fn |> 
  read_tsv()


# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig6.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 12, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_heatmap_colors <- brewer.pal(9, "Blues")
grh_heatmap_colors <- brewer.pal(9, "Oranges")
twi_heatmap_colors <- brewer.pal(9, "GnBu")


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

zld_motif_instances <- import(zld_motif_instances_fn) 

zld_bound_motifs <- zld_motif_instances |> 
  plyranges::join_overlap_left(makeGRangesFromDataFrame(zld_ChIP_classes, keep.extra.columns = TRUE)) |> 
  plyranges::filter(!is.na(class))

# generate heatmap
bw <- c(
  # S2_Zld_ChIP =  S2_Zld_ChIP_bw,
  MNase = MNase_seq_bw
)

regions <- zld_bound_motifs


# plot metaplot 1
# plot_range <- c(-0.1,3)
plot_colors <- c(
  i = "#DEEBF7",
  ii = "#9ECAE1",
  iii = "#3182BD"
)

a_plot <- plot_average(bw, regions = regions, row_split = regions$class, line_width = 0.5, upstream = hm_upstream, downstream = hm_downstream) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        # axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        # legend.position = "none"
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0)
        
  )

a_plot <- a_plot +
  geom_vline(xintercept = (max(a_plot$data$position) / 2), lty = 2, color = "gray")

plotGG(
  plot = a_plot,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 4, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)


# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

grh_motif_instances <- import(grh_motif_instances_fn) 

grh_bound_motifs <- grh_motif_instances |> 
  plyranges::join_overlap_left(makeGRangesFromDataFrame(grh_ChIP_classes, keep.extra.columns = TRUE)) |> 
  plyranges::filter(!is.na(class))

# generate heatmap
bw <- c(
  MNase = MNase_seq_bw
)

regions <- grh_bound_motifs


# plot metaplot 1
# plot_range <- c(-0.1,3)
plot_colors <- c(
  i = "#FEE6CE",
  ii = "#FDAE6B",
  iii = "#E6550D"
)

b_plot <- plot_average(bw, regions = regions, row_split = regions$class, line_width = 0.5, upstream = hm_upstream, downstream = hm_downstream) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        # axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        # legend.position = "none"
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0)
        
  )

b_plot <- b_plot +
  geom_vline(xintercept = (max(b_plot$data$position) / 2), lty = 2, color = "gray")

plotGG(
  plot = b_plot,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 4, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 10.5
ref_y <- 0.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

twi_motif_instances <- import(twi_motif_instances_fn) 

twi_bound_motifs <- twi_motif_instances |> 
  plyranges::join_overlap_left(makeGRangesFromDataFrame(twi_ChIP_classes, keep.extra.columns = TRUE)) |> 
  plyranges::filter(!is.na(class))

# generate heatmap
bw <- c(
  MNase = MNase_seq_bw
)

regions <- twi_bound_motifs


# plot metaplot 1
# plot_range <- c(-0.1,3)
plot_colors <- c(
  i =  "#E0F3DB",
  ii = "#A8DDB5",
  iii = "#43A2CA"
)

c_plot <- plot_average(bw, regions = regions, row_split = regions$class, line_width = 0.5, upstream = hm_upstream, downstream = hm_downstream) +
  # scale_color_brewer(palette = "Blues") +
  scale_color_manual(values = plot_colors) +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        # axis.title.y = element_blank(),
        plot.margin = margin(0,0,0,0),
        # legend.position = "none"
        legend.key.size = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(0,0,0,0)
        
  )

c_plot <- c_plot +
  geom_vline(xintercept = (max(c_plot$data$position) / 2), lty = 2, color = "gray")

plotGG(
  plot = c_plot,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 4, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 4.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

hm_upstream <- 250
hm_downstream <- 250

bw <- c(
  # Zld_ChIP = S2_Zld_ChIP_bw,
  MNase = MNase_seq_bw
)



regions <- zld_bound_motifs %>%
  plyranges::filter(class == "ii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

d_hm_1 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = d_hm_1,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

regions <- zld_bound_motifs %>%
  plyranges::filter(class == "iii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

d_hm_2 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = d_hm_2,
  x = (ref_x + 0.75), y = (ref_y + 3.5),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class II", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 1.5, just = "bottom", default.units = "cm", rot = 90,
)


plotText(
  label = "class III", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 5, just = "bottom", default.units = "cm", rot = 90,
)


# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 4.5

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

hm_upstream <- 250
hm_downstream <- 250

bw <- c(
  MNase = MNase_seq_bw
)



regions <- grh_bound_motifs %>%
  plyranges::filter(class == "ii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

e_hm_1 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = e_hm_1,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

regions <- grh_bound_motifs %>%
  plyranges::filter(class == "iii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

e_hm_2 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = e_hm_2,
  x = (ref_x + 0.75), y = (ref_y + 3.5),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class II", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 1.5, just = "bottom", default.units = "cm", rot = 90,
)


plotText(
  label = "class III", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 5, just = "bottom", default.units = "cm", rot = 90,
)


# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 10.5
ref_y <- 4.5

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

hm_upstream <- 250
hm_downstream <- 250

bw <- c(
  MNase = MNase_seq_bw
)



regions <- twi_bound_motifs %>%
  plyranges::filter(class == "ii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

f_hm_1 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = f_hm_1,
  x = (ref_x + 0.75), y = (ref_y + 0.25),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

regions <- twi_bound_motifs %>%
  plyranges::filter(class == "iii")

cl_data <- coverage_matrix(
  file = MNase_seq_bw, 
  regions, 
  upstream = hm_upstream, 
  downstream = hm_downstream)

cl_data[is.na(cl_data)] <- 0
hc <- hclust(dist(cl_data))

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_order = hc$order,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
  
)

f_hm_2 <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))


# place heatmap on page
plotGG(
  plot = f_hm_2,
  x = (ref_x + 0.75), y = (ref_y + 3.5),
  width = 3.5, height = 3, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "MNase-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "class II", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 1.5, just = "bottom", default.units = "cm", rot = 90,
)


plotText(
  label = "class III", params = large_text_params, fontface = "bold",
  x = ref_x + 0.5, y = ref_y + 5, just = "bottom", default.units = "cm", rot = 90,
)



dev.off()
