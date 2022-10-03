# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# library(grImport)
library(grid)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# function to get range for browser tracks -------------------------------------
signal_range <- function(x, extra = 0.05) {
  min <- floor(min(x))
  max <- ceiling(max(x))
  range <- max - min
  margin <- range * extra
  output <- c(min - margin, max + margin)
  return(output)
}

# function to extract legend from ggplot ---------------------------------------
get_legend <- function(plot){ 
  tmp <- ggplot_gtable(ggplot_build(plot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

# define input files ===========================================================
class_I_bed_fn <- c(
  zld_class_I = "results/ChIP_peak_classes/zld_class_I.bed",
  grh_class_I = "results/ChIP_peak_classes/grh_class_I.bed",
  twi_class_I = "results/ChIP_peak_classes/twi_class_I.bed"
)

ns_sites_bed_fn <- "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed"

S2_Zld_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
S2_Grh_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
S2_Twi_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
H3K27ac_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw"
Nej_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw"
H3K4me1_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw"
H3K4me3_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw"
H2AV_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw"

zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
twi_ChIP_classes_fn <- "results/ChIP_peak_classes/twi_ChIP_classes.tsv"

zld_motifs_fn <- "results/motif_instances/zld_motifs.tsv"
grh_motifs_fn <- "results/motif_instances/grh_motifs.tsv"
twi_motifs_fn <- "results/motif_instances/twi_motifs.tsv"

zld_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed"
grh_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed"
twi_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed"

## create blank layout for plot =================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/fig3.pdf")
pageCreate(width = 12.0, height = 12.5, default.units = "cm", showGuides = TRUE)

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

# generate euler plot of overlap between Zld, Grh and Twi class I sites
nonspecific_sites <- class_I_bed_fn |>
  map(rtracklayer::import) |>
  GRangesList() |>
  peak_overlap_table()

euler_fit <- nonspecific_sites |> 
  dplyr::select(6:8) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(zld_color, grh_color, twi_color))


# place euler plot on page
plotGG(
  plot = euler_plot,
  x = (ref_x), y = (ref_y + 0.25),
  width = 2.5, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Twi class I", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = "top", default.units = "cm", fontcolor = twi_color
)

plotText(
  label = paste0("Zld", "\n", "class I"), params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 2.5), just = "center", default.units = "cm", fontcolor = zld_color
)

plotText(
  label = paste0("Grh", "\n", "class I"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y + 2.5), just = "center", default.units = "cm", fontcolor = grh_color
)

# Panel B ======================================================================
ref_x <- 3.5
ref_y <- 0.5 

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in overlapping class I sites
ns_sites <- rtracklayer::import(ns_sites_bed_fn)

# generate heatmap
bw <- c(
  S2_Zld_ChIP =  S2_Zld_ChIP_bw,
  S2_Grh_ChIP =  S2_Grh_ChIP_bw,
  S2_Twi_ChIP =  S2_Twi_ChIP_bw,
  H3K27ac = H3K27ac_bw,
  Nej = Nej_bw,
  H3K4me1 =  H3K4me1_bw,
  H3K4me3 = H3K4me3_bw,
  H2AV = H2AV_bw
)

hm_upstream <- 500
hm_downstream <- 500

regions <- ns_sites 

hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = brewer.pal(9, "Greys"), 
  order_by_samples = 1:3, 
  individual_scales = TRUE,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
)

b_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = b_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 8, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_1_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = paste(length(regions), "overlapping", "\n", "class I sites"), params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 1.5), just = c("center"), default.units = "cm", rot = 90
)

plotText(
  label = "Zld", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.70625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Grh", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.71875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Twi", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.73125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K27ac", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.74375), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "CBP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.75625), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me1", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.76875), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 6.78125), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H2AV", params = small_text_params, fontface = "bold",
  x = (ref_x + 7.79375), y = (ref_y), just = c("center"), default.units = "cm"
)

# Panel C ======================================================================
# panel label
ref_x <- 0.5
ref_y <- 3.75 


# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import ChIP classes
zld_ChIP_classes <- read_tsv(zld_ChIP_classes_fn)
grh_ChIP_classes <- read_tsv(grh_ChIP_classes_fn)
twi_ChIP_classes <- read_tsv(twi_ChIP_classes_fn)

# get background motif frequency
zld_motifs_gr <- 
  zld_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

zld_ATAC_peaks <- 
  zld_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
zld_ATAC_peaks$class <- "all_ATAC_peaks"

zld_ATAC_peaks$n_motifs <- countOverlaps(zld_ATAC_peaks, zld_motifs_gr)

zld_merge_peaks <- zld_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)

grh_motifs_gr <- 
  grh_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

grh_ATAC_peaks <- 
  grh_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
grh_ATAC_peaks$class <- "all_ATAC_peaks"

grh_ATAC_peaks$n_motifs <- countOverlaps(grh_ATAC_peaks, grh_motifs_gr)

grh_merge_peaks <- grh_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)

twi_motifs_gr <- 
  twi_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

twi_ATAC_peaks <- 
  twi_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
twi_ATAC_peaks$class <- "all_ATAC_peaks"

twi_ATAC_peaks$n_motifs <- countOverlaps(twi_ATAC_peaks, twi_motifs_gr)

twi_merge_peaks <- twi_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)


# generate heatmap ggplots and extract legends
hm_limits <- c(0,100)

zld_plot_base <- zld_ChIP_classes |> 
  bind_rows(zld_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Zld motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=hm_limits, breaks = hm_limits, name = NULL, guide = guide_colourbar(direction = "horizontal")) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text = element_blank(),
        axis.title = element_blank())

zld_legend <- get_legend(zld_plot_base)

zld_plot <- 
  zld_plot_base +
  guides(fill = "none")

grh_plot_base <- grh_ChIP_classes |> 
  bind_rows(grh_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Grh motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Oranges", direction = 1, limits=hm_limits, breaks = hm_limits, name = NULL, guide = guide_colourbar(direction = "horizontal")) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text = element_blank(),
        axis.title = element_blank())

grh_legend <- get_legend(grh_plot_base)

grh_plot <- 
  grh_plot_base +
  guides(fill = "none")

twi_plot_base <- twi_ChIP_classes |> 
  bind_rows(twi_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Twi motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(
    palette = "GnBu", 
    direction = 1, 
    limits=hm_limits, 
    breaks = hm_limits, 
    name = NULL, guide = 
      guide_colourbar(direction = "horizontal")
    ) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text = element_blank(),
        axis.title = element_blank())

twi_legend <- get_legend(twi_plot_base)

twi_plot <- 
  twi_plot_base +
  guides(fill = "none")


# place heatmaps on page
plotGG(
  plot = zld_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 3.5, height = 0.5, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = grh_plot,
  x = ref_x + 0.25, y = ref_y + 0.3,
  width = 3.5, height = 0.5, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = twi_plot,
  x = ref_x + 0.25, y = ref_y + 0.6,
  width = 3.5, height = 0.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels to heatmaps
plotText(
  label = "Zld motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.22), just = c("right"), default.units = "cm"
)

plotText(
  label = "Grh motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.52), just = c("right"), default.units = "cm"
)

plotText(
  label = "Twi motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.82), just = c("right"), default.units = "cm"
)

plotText(
  label = "all accessible", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.8), y = (ref_y + 1), just = c("right", "top"), default.units = "cm", rot = 45
)

plotText(
  label = "class I", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.6), y = (ref_y + 1), just = c("right", "top"), default.units = "cm", rot = 45
)

plotText(
  label = "class II", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.4), y = (ref_y + 1), just = c("right", "top"), default.units = "cm", rot = 45
)

plotText(
  label = "class III", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.2), y = (ref_y + 1), just = c("right", "top"), default.units = "cm", rot = 45
)

# place legends on page
plotGG(
  plot = zld_legend,
  width = 0.5, height = 0.2,
  x = ref_x + 0.75, y = ref_y + 1.75,
  default.units = "cm"
)

plotGG(
  plot = grh_legend,
  width = 0.5, height = 0.2,
  x = ref_x + 1.75, y = ref_y + 1.75,
  default.units = "cm"
)

plotGG(
  plot = twi_legend,
  width = 0.5, height = 0.2,
  x = ref_x + 2.75, y = ref_y + 1.75,
  default.units = "cm"
)


# # add logos plots next to heatmap
# library(universalmotif)
# 
# zld_motif <- readRDS("published_ChIPseq/results/motifs/embryo-nc14_aZld/embryo-nc14_aZld_PWM.rds")
# zld_logo <- 
#   view_motifs(zld_motif) + 
#   theme(
#     axis.text.x = element_blank(), 
#     axis.text.y = element_blank(), 
#     axis.line.y = element_blank(), 
#     axis.ticks.y = element_blank(), 
#     axis.title.y = element_blank(),
#     plot.margin = margin(0, 0, 0, 0, "cm")
#     )
# 
# plotGG(
#   plot = zld_logo,
#   x = ref_x + 0.5, y = ref_y + 0.4,
#   width = 0.5, height = 0.5, just = c("right", "center"),
#   default.units = "cm"
# )

# Panel D ======================================================================
# set y limits for all n_motifs boxplots
n_motifs_ylim <- c(0,8)

# panel label
ref_x <- 4.25
ref_y <- 3.75 

# boxplot
zld_n_motifs_plot <- zld_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim)
  

plotGG(
  plot = zld_n_motifs_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel E ======================================================================
# panel label
ref_x <- 6
ref_y <- 3.75 

# boxplot
zld_n_motifs_plot <- zld_ChIP_classes |>
  ggplot(aes(x = class, y = average_motif_score)) + 
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize)


plotGG(
  plot = zld_n_motifs_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# boxplot
grh_n_motifs_plot <- grh_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = grh_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim)

plotGG(
  plot = grh_n_motifs_plot,
  width = 2, height = 1.5,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel F ======================================================================
# panel label
ref_x <- 9.5
ref_y <- 3.75 

# boxplot
twi_n_motifs_plot <- twi_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = twi_color, outlier.shape = NA, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim)

plotGG(
  plot = twi_n_motifs_plot,
  width = 2, height = 1.5,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)



# 
# zld_ChIP_classes |> 
#   ggplot(aes(x = class, y = average_motif_score)) + geom_boxplot()
# 
# grh_ChIP_classes |> 
#   ggplot(aes(x = class, y = n_motifs)) + geom_boxplot()
# 
# grh_ChIP_classes |> 
#   ggplot(aes(x = class, y = average_motif_score)) + geom_boxplot()
# 
# twi_ChIP_classes |> 
#   ggplot(aes(x = class, y = n_motifs)) + geom_boxplot()
# 
# twi_ChIP_classes |> 
#   ggplot(aes(x = class, y = average_motif_score)) + geom_boxplot()

dev.off()
