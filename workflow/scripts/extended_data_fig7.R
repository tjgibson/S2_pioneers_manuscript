# setup ------------------------------------------------------------------------
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(eulerr)
library(RColorBrewer)
library(plotgardener)

source("workflow/scripts/utils.R")
source("workflow/scripts/plot_heatmap.R")
source("../useful_functions/peak_processing.R")

# define input files ===========================================================
# Zld input files
S2_Zld_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Zld_aZld_IP.bed"
nc14_Zld_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/embryo-nc14_aZld_summits.bed"
brain_Zld_ChIP_fn <-  "../S2_pioneers/results/bed_files/peaks/brain_aZld.bed"

S2_ZLD_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
embryo_Zld_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-nc14_aZld.bw"
brain_Zld_ChIP_bw <-  "../S2_pioneers/results/bigwigs/normalized/Zld_NB_merged_zscore.bw"
zld_motif_bw <-  "results/motif_instances/zld_motifs_pvalue.bw"

zld_motif_instances_fn <- "results/motif_instances/zld_motifs.tsv"

# Grh input files
S2_Grh_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP.bed"
embryo_Grh_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/embryo-15-16H_aGrh_summits.bed"
wing_Grh_ChIP_fn <- "published_ChIPseq/results/peaks/individual/narrow/wing-disc_aGrh_summits.bed"

S2_Grh_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
embryo_Grh_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw"
wing_disc_Grh_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/wing-disc_aGrh.bw"
grh_motif_bw <-  "results/motif_instances/grh_motifs_pvalue.bw"

grh_motif_instances_fn <- "results/motif_instances/grh_motifs.tsv"

# Twi input files
S2_twi_ChIP_fn <- "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP.bed"
embryo_twi_ChIP_fn <- "published_ChIPseq/results/peaks/merged/narrow/embryo-1-3H_aTwi_summits.bed"

S2_twi_ChIP_bw <-   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
embryo_twi_ChIP_bw <-  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw"
twi_motif_bw <-  "results/motif_instances/twi_motifs_pvalue.bw"

twi_motif_instances_fn <- "results/motif_instances/twi_motifs.tsv"

# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig7.pdf", useDingbats = FALSE)
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

twi_color <- "#00A08A"
twi_heatmap_colors <- brewer.pal(9, "GnBu")

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500

# panel A ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
peak_overlaps <- 
  c(
    S2_Zld_ChIP = S2_Zld_ChIP_fn,
    nc14_Zld_ChIP = nc14_Zld_ChIP_fn,
    brain_Zld_ChIP = brain_Zld_ChIP_fn
  ) |> 
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  
  
  GRangesList() |> 
  
  peak_overlap_table() |> 
  mutate(overlap_group = case_when(
    S2_Zld_ChIP & !nc14_Zld_ChIP & !brain_Zld_ChIP ~ 1,
    !S2_Zld_ChIP & nc14_Zld_ChIP & !brain_Zld_ChIP ~ 2,
    !S2_Zld_ChIP & !nc14_Zld_ChIP & brain_Zld_ChIP ~ 3,
    S2_Zld_ChIP & nc14_Zld_ChIP & !brain_Zld_ChIP ~ 4,
    S2_Zld_ChIP & !nc14_Zld_ChIP & brain_Zld_ChIP ~ 5,
    !S2_Zld_ChIP & nc14_Zld_ChIP & brain_Zld_ChIP ~ 6,
    S2_Zld_ChIP & nc14_Zld_ChIP & brain_Zld_ChIP ~ 7,
  ))


# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:8) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(zld_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.25), y = (ref_y),
  width = 2.5, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Brain Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = "top", default.units = "cm"
)

plotText(
  label = paste0("S2 Zld ChIP"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.25), y = (ref_y + 2.25), just = "center", default.units = "cm"
)

plotText(
  label = paste0("nc14 embryo", "\n", "Zld Chip"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y + 2.5), just = "center", default.units = "cm"
)


# panel B ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 3.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


# plot heatmap of overlaps for Zld
bw <- c(
  S2_ZLD_ChIP =  S2_ZLD_ChIP_bw,
  embryo_Zld_ChIP = embryo_Zld_ChIP_bw,
  brain_Zld_ChIP = brain_Zld_ChIP_bw,
  zld_motif = zld_motif_bw
)

group <- c(1,1,1,2)

regions <- peak_overlaps |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$overlap_group,
  scale_group = group,
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
  width = 6, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_7_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "S2 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.9), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "nc14 Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Brain Zld ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Zld motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.6), y = (ref_y), just = c("center"), default.units = "cm"
)


# panel C ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 0.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# quantify n motifs per peak ---------------------------------------------------
peak_overlaps$name <- paste0("peak_", seq(nrow(peak_overlaps)))
peaks_gr <- makeGRangesFromDataFrame(peak_overlaps, keep.extra.columns = TRUE)

motifs_gr <- read_tsv(zld_motif_instances_fn) %>% 
  rename(motif_score = score) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

peak_overlaps$n_motifs <- countOverlaps(peaks_gr, motifs_gr)

# get average motif score per peak ---------------------------------------------
peak_motif_scores <- peaks_gr %>% 
  plyranges::join_overlap_left(motifs_gr) %>% 
  as.data.frame() %>% 
  group_by(name) %>% 
  summarise(average_motif_score = mean(motif_score)) 

peak_overlaps <- peak_overlaps %>% 
  left_join(peak_motif_scores, by = "name")



zld_motif_percent_hm <- 
  peak_overlaps |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(overlap_group) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Zld motif") |> 
  ggplot(aes(x= as.factor(overlap_group), y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits=c(0,100), breaks = c(0,100), name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        # axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank()
  )



# place heatmap on page
plotGG(
  plot = zld_motif_percent_hm,
  x = ref_x + 0.25, y = ref_y,
  width = 6, height = 0.75, just = c("left", "top"),
  default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 1.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

d_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = n_motifs)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = zld_color) +
  theme_classic(base_size = small_text_params$fontsize)

plotGG(
  plot = d_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel E ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 13
ref_y <- 1.5

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


e_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = average_motif_score)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = zld_color) +
  theme_classic(base_size = small_text_params$fontsize)


plotGG(
  plot = e_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)


# panel F ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 0.5
ref_y <- 3.5

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
# generate overlap table and annotate ChIP classes


peak_overlaps <- 
  c(
    S2_Grh_ChIP = S2_Grh_ChIP_fn,
    embryo_Grh_ChIP = embryo_Grh_ChIP_fn,
    wing_Grh_ChIP = wing_Grh_ChIP_fn
  ) |> 
  
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  GRangesList() |> 
  peak_overlap_table() |> 
  
  mutate(overlap_group = case_when(
    S2_Grh_ChIP & !embryo_Grh_ChIP & !wing_Grh_ChIP ~ 1,
    !S2_Grh_ChIP & embryo_Grh_ChIP & !wing_Grh_ChIP ~ 2,
    !S2_Grh_ChIP & !embryo_Grh_ChIP & wing_Grh_ChIP ~ 3,
    S2_Grh_ChIP & embryo_Grh_ChIP & !wing_Grh_ChIP ~ 4,
    S2_Grh_ChIP & !embryo_Grh_ChIP & wing_Grh_ChIP ~ 5,
    !S2_Grh_ChIP & embryo_Grh_ChIP & wing_Grh_ChIP ~ 6,
    S2_Grh_ChIP & embryo_Grh_ChIP & wing_Grh_ChIP ~ 7,
  ))



# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:8) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(grh_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.25), y = (ref_y),
  width = 2.5, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = "Wing disc Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.25), y = (ref_y), just = "top", default.units = "cm"
)

plotText(
  label = paste0("S2 Grh ChIP"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.25), y = (ref_y + 2.25), just = "center", default.units = "cm"
)

plotText(
  label = paste0("15-16H embryo", "\n", "Grh Chip"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y + 2.5), just = "center", default.units = "cm"
)


# panel G ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 3.5
ref_y <- 3.5

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


# plot heatmap of overlaps for Zld
bw <- c(
  S2_Grh_ChIP =  S2_Grh_ChIP_bw,
  embryo_Grh_ChIP = embryo_Grh_ChIP_bw,
  wing_disc_Grh_ChIP = wing_disc_Grh_ChIP_bw,
  grh_motif = grh_motif_bw
)

group <- c(1,1,1,2)

regions <- peak_overlaps |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_split = regions$overlap_group,
  scale_group = group,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
)

g_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = g_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 6, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_7_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "S2 Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 0.9), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "embryo Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Wing disc Grh ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "Grh motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.6), y = (ref_y), just = c("center"), default.units = "cm"
)


# panel H ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 3.5

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# quantify n motifs per peak ---------------------------------------------------
peak_overlaps$name <- paste0("peak_", seq(nrow(peak_overlaps)))
peaks_gr <- makeGRangesFromDataFrame(peak_overlaps, keep.extra.columns = TRUE)

motifs_gr <- read_tsv(grh_motif_instances_fn) %>% 
  rename(motif_score = score) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

peak_overlaps$n_motifs <- countOverlaps(peaks_gr, motifs_gr)

# get average motif score per peak ---------------------------------------------
peak_motif_scores <- peaks_gr %>% 
  plyranges::join_overlap_left(motifs_gr) %>% 
  as.data.frame() %>% 
  group_by(name) %>% 
  summarise(average_motif_score = mean(motif_score)) 

peak_overlaps <- peak_overlaps %>% 
  left_join(peak_motif_scores, by = "name")



grh_motif_percent_hm <- 
  peak_overlaps |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(overlap_group) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Grh motif") |> 
  ggplot(aes(x= as.factor(overlap_group), y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "Oranges", direction = 1, limits=c(0,100), breaks = c(0,100), name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        # axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank()
  )



# place heatmap on page
plotGG(
  plot = grh_motif_percent_hm,
  x = ref_x + 0.25, y = ref_y,
  width = 6, height = 0.75, just = c("left", "top"),
  default.units = "cm"
)

# panel I ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 4.5

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

i_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = n_motifs)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = grh_color) +
  theme_classic(base_size = small_text_params$fontsize)

plotGG(
  plot = i_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel j ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 13
ref_y <- 4.5

# panel label
plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


i_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = average_motif_score)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = grh_color) +
  theme_classic(base_size = small_text_params$fontsize)


plotGG(
  plot = i_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel K ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 0.5
ref_y <- 6.5

# panel label
plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# generate overlap table and annotate ChIP classes -----------------------------
# generate overlap table and annotate ChIP classes


peak_overlaps <- 
  c(
    S2_twi_ChIP = S2_twi_ChIP_fn,
    embryo_twi_ChIP = embryo_twi_ChIP_fn
  ) |> 
  
  map(rtracklayer::import) |> 
  map(resize, width = 201, fix = "center") |> 
  GRangesList() |> 
  peak_overlap_table() |> 
  
  mutate(overlap_group = case_when(
    S2_twi_ChIP & !embryo_twi_ChIP ~ 1,
    !S2_twi_ChIP & embryo_twi_ChIP ~ 2,
    S2_twi_ChIP & embryo_twi_ChIP ~ 3
  ))



# plot overlaps
euler_fit <- peak_overlaps |> 
  dplyr::select(6:7) |> 
  as.matrix() |> 
  eulerr::euler()

euler_plot <- plot(euler_fit,
                   quantities = list(fontsize = small_text_params$fontsize), labels = NULL,
                   fills = c(twi_color, "white", "grey"))

plotGG(
  plot = euler_plot,
  x = (ref_x + 0.25), y = (ref_y),
  width = 2.5, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# add labels for different groups
plotText(
  label = paste0("S2 Twi ChIP"), params = small_text_params, fontface = "bold",
  x = (ref_x + 0.25), y = (ref_y + 0.25), just = "center", default.units = "cm"
)

plotText(
  label = paste0("2-3H embryo", "\n", "Twi Chip"), params = small_text_params, fontface = "bold",
  x = (ref_x + 2.5), y = (ref_y + 0.25), just = "center", default.units = "cm"
)


# panel l ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 3.5
ref_y <- 6.5

# panel label
plotText(
  label = "l", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


# plot heatmap of overlaps for Twi
bw <- c(
  S2_twi_ChIP =  S2_twi_ChIP_bw,
  embryo_twi_ChIP = embryo_twi_ChIP_bw,
  twi_motif = twi_motif_bw
)

group <- c(1,1,2)

regions <- peak_overlaps |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_split = regions$overlap_group,
  scale_group = group,
  return_heatmap_list = TRUE,
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5)
)

l_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap on page
plotGG(
  plot = l_hm,
  x = (ref_x + 0.25), y = (ref_y + 0.25),
  width = 6, height = 2.25, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_7_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# add heatmap labels
plotText(
  label = "S2 Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 1.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "embryo Twi ChIP", params = small_text_params, fontface = "bold",
  x = (ref_x + 3.2), y = (ref_y), just = c("center"), default.units = "cm"
)


plotText(
  label = "Twi motif", params = small_text_params, fontface = "bold",
  x = (ref_x + 5.3), y = (ref_y), just = c("center"), default.units = "cm"
)


# panel M ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 6.5

# panel label
plotText(
  label = "m", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

# quantify n motifs per peak ---------------------------------------------------
peak_overlaps$name <- paste0("peak_", seq(nrow(peak_overlaps)))
peaks_gr <- makeGRangesFromDataFrame(peak_overlaps, keep.extra.columns = TRUE)

motifs_gr <- read_tsv(twi_motif_instances_fn) %>% 
  rename(motif_score = score) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

peak_overlaps$n_motifs <- countOverlaps(peaks_gr, motifs_gr)

# get average motif score per peak ---------------------------------------------
peak_motif_scores <- peaks_gr %>% 
  plyranges::join_overlap_left(motifs_gr) %>% 
  as.data.frame() %>% 
  group_by(name) %>% 
  summarise(average_motif_score = mean(motif_score)) 

peak_overlaps <- peak_overlaps %>% 
  left_join(peak_motif_scores, by = "name")



twi_motif_percent_hm <- 
  peak_overlaps |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(overlap_group) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  add_column(motif_name = "Grh motif") |> 
  ggplot(aes(x= as.factor(overlap_group), y = motif_name, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "GnBu", direction = 1, limits=c(0,100), breaks = c(0,100), name = NULL) +
  geom_text(aes(label = percent_with_motif), size = small_text_params$fontsize * 0.35) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(1, 'mm'), 
        # axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text.y = element_blank(),
        axis.title = element_blank()
  )



# place heatmap on page
plotGG(
  plot = twi_motif_percent_hm,
  x = ref_x + 0.25, y = ref_y,
  width = 6, height = 0.75, just = c("left", "top"),
  default.units = "cm"
)

# panel N ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 10
ref_y <- 7.5

# panel label
plotText(
  label = "n", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)

n_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = n_motifs)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = twi_color) +
  theme_classic(base_size = small_text_params$fontsize)

plotGG(
  plot = n_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel O ======================================================================
# reference points for positioning figure components ---------------------------
ref_x <- 13
ref_y <- 7.5

# panel label
plotText(
  label = "o", params = panel_label_params, fontface = "bold",
  x = (ref_x), y = (ref_y), just = "bottom", default.units = "cm"
)


o_plot <- 
  peak_overlaps |> 
  ggplot(aes(x = as.factor(overlap_group), y = average_motif_score)) +
  geom_boxplot(outlier.size = 0.01, lwd = 0.1, fill = twi_color) +
  theme_classic(base_size = small_text_params$fontsize)


plotGG(
  plot = o_plot,
  x = ref_x + 0.25, y = ref_y,
  width = 2.75, height = 2, just = c("left", "top"),
  default.units = "cm"
)


dev.off()
