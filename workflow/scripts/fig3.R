# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
# library(grImport)
suppressPackageStartupMessages(library(grid))

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")


# define input files ===========================================================
class_I_bed_fn <- c(
  zld_class_I = snakemake@input[["zld_class_I_bed_fn"]],
  grh_class_I = snakemake@input[["grh_class_I_bed_fn"]],
  twi_class_I = snakemake@input[["twi_class_I_bed_fn"]]
)

ns_sites_bed_fn <- snakemake@input[["ns_sites_bed_fn"]]

S2_Zld_ChIP_bw <-  snakemake@input[["S2_Zld_ChIP_bw"]]
S2_Grh_ChIP_bw <-  snakemake@input[["S2_Grh_ChIP_bw"]]
S2_Twi_ChIP_bw <-  snakemake@input[["S2_Twi_ChIP_bw"]]
H3K27ac_bw <- snakemake@input[["H3K27ac_bw"]]
Nej_bw <- snakemake@input[["Nej_bw"]]
H3K4me1_bw <-  snakemake@input[["H3K4me1_bw"]]
H3K4me3_bw <- snakemake@input[["H3K4me3_bw"]]
H2AV_bw <- snakemake@input[["H2AV_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes_fn"]]
grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes_fn"]]
twi_ChIP_classes_fn <- snakemake@input[["twi_ChIP_classes_fn"]]

zld_motifs_fn <- snakemake@input[["zld_motifs_fn"]]
grh_motifs_fn <- snakemake@input[["grh_motifs_fn"]]
twi_motifs_fn <- snakemake@input[["twi_motifs_fn"]]

zld_WT_atac_fn <- snakemake@input[["zld_WT_atac_fn"]]
grh_WT_atac_fn <- snakemake@input[["grh_WT_atac_fn"]]
twi_WT_atac_fn <- snakemake@input[["twi_WT_atac_fn"]]

# class_I_bed_fn <- c(
#   zld_class_I = "results/ChIP_peak_classes/zld_class_I.bed",
#   grh_class_I = "results/ChIP_peak_classes/grh_class_I.bed",
#   twi_class_I = "results/ChIP_peak_classes/twi_class_I.bed"
# )
# 
# ns_sites_bed_fn <- "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed"
# 
# S2_Zld_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# S2_Grh_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# S2_Twi_ChIP_bw <-  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# H3K27ac_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw"
# Nej_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw"
# H3K4me1_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw"
# H3K4me3_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw"
# H2AV_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw"
# 
# zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# twi_ChIP_classes_fn <- "results/ChIP_peak_classes/twi_ChIP_classes.tsv"
# 
# zld_motifs_fn <- "results/motif_instances/zld_motifs.tsv"
# grh_motifs_fn <- "results/motif_instances/grh_motifs.tsv"
# twi_motifs_fn <- "results/motif_instances/twi_motifs.tsv"
# 
# zld_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed"
# grh_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed"
# twi_WT_atac_fn <- "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed"

## create blank layout for plot =================================================
# define figure dimensions in cm
fig_width <-  18
fig_height <- 12

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = TRUE)

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
  x = (ref_x + 0.65), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Grh", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.68), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Twi", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.71), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K27ac", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.74), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "CBP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.77), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me1", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.80), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H3K4me3", params = small_text_params, fontface = "bold",
  x = (ref_x + 6.83), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "H2AV", params = small_text_params, fontface = "bold",
  x = (ref_x + 7.86), y = (ref_y), just = c("center"), default.units = "cm"
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
zld_ChIP_classes <- read_tsv(zld_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

# get background motif frequency
zld_motifs_gr <- 
  zld_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

zld_ATAC_peaks <- 
  zld_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
zld_ATAC_peaks$class <- "all ATAC peaks"

zld_ATAC_peaks$n_motifs <- countOverlaps(zld_ATAC_peaks, zld_motifs_gr)

zld_merge_peaks <- zld_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)

# generate heatmap ggplots and extract legends
hm_limits <- c(0,100)

zld_plot <- zld_ChIP_classes |> 
  bind_rows(zld_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Zld motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=hm_limits, breaks = hm_limits, name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.7,-1.2), 
        legend.direction = "horizontal"
        )



# place heatmap on page
plotGG(
  plot = zld_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 3.25, height = 1.5, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "Zld motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.2), just = c("right","center"), default.units = "cm"
)



# add logos plots next to heatmap
library(universalmotif)

zld_motif <- readRDS("published_ChIPseq/results/motifs/embryo-nc14_aZld/embryo-nc14_aZld_PWM.rds")
zld_logo <-
  view_motifs(zld_motif) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm")
    )

plotGG(
  plot = zld_logo,
  x = ref_x, y = ref_y + 0.6,
  width = 0.9, height = 0.5, just = c("center", "center"),
  default.units = "cm"
)

# Panel D ======================================================================
# panel label
ref_x <- 4.5
ref_y <- 3.75 


# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import ChIP classes
grh_ChIP_classes <- read_tsv(grh_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

# get background motif frequency
grh_motifs_gr <- 
  grh_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

grh_ATAC_peaks <- 
  grh_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
grh_ATAC_peaks$class <- "all ATAC peaks"

grh_ATAC_peaks$n_motifs <- countOverlaps(grh_ATAC_peaks, grh_motifs_gr)

grh_merge_peaks <- grh_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)


grh_plot <- grh_ChIP_classes |> 
  bind_rows(grh_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "grh motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Oranges", direction = 1, limits=hm_limits, breaks = hm_limits, name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.7,-1.2), 
        legend.direction = "horizontal"
  )



# place heatmap on page
plotGG(
  plot = grh_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 3.25, height = 1.5, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "Grh motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.2), just = c("right","center"), default.units = "cm"
)



# add logos plots next to heatmap
grh_motif <- readRDS("published_ChIPseq/results/motifs/embryo-15-16H_aGrh/embryo-15-16H_aGrh_PWM.rds")
grh_logo <-
  view_motifs(grh_motif) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )

plotGG(
  plot = grh_logo,
  x = ref_x, y = ref_y + 0.6,
  width = 0.9, height = 0.5, just = c("center", "center"),
  default.units = "cm"
)

# Panel E ======================================================================
# panel label
ref_x <- 8.5
ref_y <- 3.75 


# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import ChIP classes
twi_ChIP_classes <- read_tsv(twi_ChIP_classes_fn) |> 
  mutate(class = str_to_upper(class))

# get background motif frequency
twi_motifs_gr <- 
  twi_motifs_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

twi_ATAC_peaks <- 
  twi_WT_atac_fn |> 
  rtracklayer::import() |> 
  resize(width = 201, fix = "center")
twi_ATAC_peaks$class <- "all ATAC peaks"

twi_ATAC_peaks$n_motifs <- countOverlaps(twi_ATAC_peaks, twi_motifs_gr)

twi_merge_peaks <- twi_ATAC_peaks |> 
  as.data.frame() |> 
  dplyr::select(-c(name, score, width, strand)) |> 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end)


twi_plot <- twi_ChIP_classes |> 
  bind_rows(twi_merge_peaks) |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "twi motif") |> 
  
  ggplot(aes(x=class, y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "GnBu", direction = 1, limits=hm_limits, breaks = hm_limits, name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.7,-1.2), 
        legend.direction = "horizontal"
  )



# place heatmap on page
plotGG(
  plot = twi_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 3.25, height = 1.5, just = c("left", "top"),
  default.units = "cm"
)

plotText(
  label = "Twi motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.4), y = (ref_y + 0.2), just = c("right","center"), default.units = "cm"
)



# add logos plots next to heatmap
twi_motif <- readRDS("published_ChIPseq/results/motifs/embryo-1-3H_aTwi/embryo-1-3H_aTwi_PWM.rds")
twi_logo <-
  view_motifs(twi_motif) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )

plotGG(
  plot = twi_logo,
  x = ref_x, y = ref_y + 0.6,
  width = 0.9, height = 0.5, just = c("center", "center"),
  default.units = "cm"
)


# Panel F ======================================================================
# set y limits for all n_motifs boxplots
n_motifs_ylim <- c(0,8)

# panel label
ref_x <- 0.5
ref_y <- 5.25 

# boxplot
zld_n_motifs_plot <- zld_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim) +
  ylab("n motifs")
  

plotGG(
  plot = zld_n_motifs_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel G ======================================================================
# panel label
ref_x <- 2.25
ref_y <- 5.25 

# boxplot
zld_motif_score_plot <- zld_ChIP_classes |>
  ggplot(aes(x = class, y = average_motif_score)) + 
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  ylab("average motif score")



plotGG(
  plot = zld_motif_score_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel H ======================================================================
# panel label
ref_x <- 4.5
ref_y <- 5.25 

# boxplot
grh_n_motifs_plot <- grh_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = grh_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim) +
  ylab("n motifs")



plotGG(
  plot = grh_n_motifs_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel I ======================================================================
# panel label
ref_x <- 6.25
ref_y <- 5.25 

# boxplot
grh_motif_score_plot <- grh_ChIP_classes |>
  ggplot(aes(x = class, y = average_motif_score)) + 
  geom_boxplot(fill = grh_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  ylab("average motif score")



plotGG(
  plot = grh_motif_score_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel J ======================================================================
# panel label
ref_x <- 8.5
ref_y <- 5.25 

# boxplot
twi_n_motifs_plot <- twi_ChIP_classes |>
  ggplot(aes(x = class, y = n_motifs)) + 
  geom_boxplot(fill = twi_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  scale_y_continuous(breaks=seq(n_motifs_ylim[1],n_motifs_ylim[2],1), limits = n_motifs_ylim) +
  ylab("n motifs")


plotGG(
  plot = twi_n_motifs_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# Panel K ======================================================================
# panel label
ref_x <- 10.25
ref_y <- 5.25 

# boxplot
twi_motif_score_plot <- twi_ChIP_classes |>
  ggplot(aes(x = class, y = average_motif_score)) + 
  geom_boxplot(fill = twi_color, outlier.size = 0.01, lwd = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  ylab("average motif score")


plotGG(
  plot = twi_motif_score_plot,
  width = 1.5, height = 2,
  x = ref_x, y = ref_y,
  default.units = "cm"
)

# panel label
plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)





dev.off()
