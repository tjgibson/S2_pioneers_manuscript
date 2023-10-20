# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(RColorBrewer)
library(EBImage)
library(rtracklayer)
library(memes)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# define input files explicitly for interactively testing script
# CR_spikeIn_counts_fn <-"histone_CUTandRUN/results/scaling_factors/epiCypher_barcode_counts.tsv"
# 
# taz_blot_aH3K27me3_image <- "data/immunoblot_raw_images/2021-06-29_taz/anti-H3K27me3_3.tif"
# taz_blot_aTub_image <- "data/immunoblot_raw_images/2021-06-29_taz/anti-tubulin_2.tif"
# 
# bw <- c(
#   H3K27ac = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw",
#   Nej = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw",
#   H3K4me1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw",
#   H3K4me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw",
#   H4K16ac = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE37865_aH4K16ac.bw",
#   H2AV = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw",
# 
#   Rpb3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Rpb3_IP.bw",
#   `PolII-pS2` = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_PolII_phosphoSer2_IP.bw",
#   H3K36me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3K36me3_IP.bw",
#   SSRP1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_SSRP1_IP.bw",
#   Spt16 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Spt16_IP.bw",
# 
#   H3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3_IP.bw",
#   H1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE127227_aH1.bw",
# 
#   Pho = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE84502_aPho_IP.bw",
#   Ez = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_Ez_IP.bw",
#   Pc = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Pc_IP.bw",
#   Psc = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Psc_IP.bw",
#   Ph = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Ph_IP.bw",
#   dRing = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_dRing_IP.bw",
#   H3K27me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE151983_S2_aH3K27me3_IP.bw",
# 
#   H3K9me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw",
#   HP1a = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE56101_aHP1a.bw",
# 
#   M1BP = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_M1BP_IP.bw",
#   GAF = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_GAGA_IP.bw",
#   BEAF32 =
#     "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE52887_aBEAF32_IP.bw",
# 
#   CTCF = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCTCF.bw",
#   CP190 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCP190_IP.bw",
#   Su_hw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aSu-Hw_IP.bw",
#   Mod_mdg4 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aMod-mdg4.bw"
# 
# )
# 
# feature_types <- list(
#   active = c("H3K27ac",
#              "Nej",
#              "H3K4me1",
#              "H3K4me3" ,
#              "H4K16ac" ,
#              "H2AV"),
# 
#   Pol_II = c("Rpb3",
#              "PolII.pS2",
#              "H3K36me3",
#              "SSRP1",
#              "Spt16"),
# 
#   histones = c("H3", "H1"),
# 
#   polycomb = c("Pho" ,
#                "Ez" ,
#                "Pc" ,
#                "Psc",
#                "Ph" ,
#                "dRing" ,
#                "H3K27me3"),
# 
#   heterochromatin = c("H3K9me3",
#                       "HP1a"),
# 
#   TF = c("M1BP",
#          "GAF",
#          "BEAF32"),
# 
#   insulator = c("CTCF" ,
#                 "CP190",
#                 "Su_hw",
#                 "Mod_mdg4")
# )
# 
# zld_chip_classes <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")
# grh_chip_classes <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")
# twi_chip_classes <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")
# 
# zld_ame_results_fn <- "results/ChIP_tissue_classes/motifs/zld_AME_results.tsv"
# zld_motif_factor_RNAseq_fn <- "results/ChIP_tissue_classes/motifs/zld_motif_factor_RNAseq.tsv"
# grh_ame_results_fn <- "results/ChIP_tissue_classes/motifs/grh_AME_results.tsv"
# grh_motif_factor_RNAseq_fn <- "results/ChIP_tissue_classes/motifs/grh_motif_factor_RNAseq.tsv"
# twi_ame_results_fn <- "results/ChIP_tissue_classes/motifs/twi_AME_results.tsv"
# twi_motif_factor_RNAseq_fn <- "results/ChIP_tissue_classes/motifs/twi_motif_factor_RNAseq.tsv"

# get input files from snakemake
CR_spikeIn_counts_fn <-snakemake@input[["CR_spikeIn_counts_fn"]]
taz_blot_aH3K27me3_image <- snakemake@input[["taz_blot_aH3K27me3_image"]]
taz_blot_aTub_image <- snakemake@input[["taz_blot_aTub_image"]]

bw <- c(
  H3K27ac = snakemake@input[["H3K27ac_bw"]],
  Nej = snakemake@input[["Nej_bw"]],
  H3K4me1 = snakemake@input[["H3K4me1_bw"]],
  H3K4me3 = snakemake@input[["H3K4me3_bw"]],
  H4K16ac = snakemake@input[["H4K16ac_bw"]],
  H2AV = snakemake@input[["H2AV_bw"]],
  
  Rpb3 = snakemake@input[["Rpb3_bw"]],
  `PolII-pS2` = snakemake@input[["PolII_pS2_bw"]],
  H3K36me3 = snakemake@input[["H3K36me3_bw"]],
  SSRP1 = snakemake@input[["SSRP1_bw"]],
  Spt16 = snakemake@input[["Spt16_bw"]],
  
  H3 = snakemake@input[["H3_bw"]],
  H1 = snakemake@input[["H1_bw"]],
  
  Pho = snakemake@input[["Pho_bw"]],
  Ez = snakemake@input[["Ez_bw"]],
  Pc = snakemake@input[["Pc_bw"]],
  Psc = snakemake@input[["Psc_bw"]],
  Ph = snakemake@input[["Ph_bw"]],
  dRing =snakemake@input[["dRing_bw"]],
  H3K27me3 = snakemake@input[["H3K27me3_bw"]],
  
  H3K9me3 = snakemake@input[["H3K9me3_bw"]],
  HP1a = snakemake@input[["HP1a_bw"]],
  
  M1BP = snakemake@input[["M1BP_bw"]],
  GAF = snakemake@input[["GAF_bw"]],
  BEAF32 = snakemake@input[["BEAF32_bw"]],
  
  CTCF =snakemake@input[["CTCF_bw"]],
  CP190 = snakemake@input[["CP190_bw"]],
  Su_hw = snakemake@input[["Su_hw_bw"]],
  Mod_mdg4 = snakemake@input[["Mod_mdg4_bw"]]
  
)

feature_types <- list(
  active = c("H3K27ac",
             "Nej",
             "H3K4me1",
             "H3K4me3" ,
             "H4K16ac" ,
             "H2AV"),
  
  Pol_II = c("Rpb3",
             "PolII.pS2",
             "H3K36me3",
             "SSRP1",
             "Spt16"),
  
  histones = c("H3", "H1"),
  
  polycomb = c("Pho" ,
               "Ez" ,
               "Pc" ,
               "Psc",
               "Ph" ,
               "dRing" ,
               "H3K27me3"),
  
  heterochromatin = c("H3K9me3",
                      "HP1a"),
  
  TF = c("M1BP",
         "GAF",
         "BEAF32"),
  
  insulator = c("CTCF" ,
                "CP190",
                "Su_hw",
                "Mod_mdg4")
)

zld_chip_classes <- read_tsv(snakemake@input[["zld_chip_classes"]]) 
grh_chip_classes <- read_tsv(snakemake@input[["grh_chip_classes"]]) 
twi_chip_classes <- read_tsv(snakemake@input[["twi_chip_classes"]]) 

zld_ame_results_fn <- snakemake@input[["zld_ame_results_fn"]]
zld_motif_factor_RNAseq_fn <- snakemake@input[["zld_motif_factor_RNAseq_fn"]]
grh_ame_results_fn <- snakemake@input[["grh_ame_results_fn"]]
grh_motif_factor_RNAseq_fn <- snakemake@input[["grh_motif_factor_RNAseq_fn"]]
twi_ame_results_fn <- snakemake@input[["twi_ame_results_fn"]]
twi_motif_factor_RNAseq_fn <- snakemake@input[["twi_motif_factor_RNAseq_fn"]]

# create blank layout for plot ===============================================
fig_width <-  18
fig_height <- 19

# open pdf
# pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
pdf("manuscript/figures/extended_data_fig8.pdf", useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

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


# reference points for positioning figure components
x_offset_class_label <- 0.25
x_offset_browser_label <- 1
x_offset_browser <- 1.5

# set genome browser height
gb_height <- 0.3

# set heatmap parameters
hm_upstream <-  500
hm_downstream <-  500


# panel A ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 0.5

plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# import images of western blot of H3K27me3 and tubulin in DMSO or control cells
# read in tiff of western blot
H3K27me3_blot_image <- readImage(taz_blot_aH3K27me3_image)

# rotate image
H3K27me3_blot_image <- H3K27me3_blot_image |>
  rotate(91, bg.col = "white") |>
  flop()

# crop image
H3K27me3_blot_image <- H3K27me3_blot_image[1536:2188,880:1102]

# anti-Tubulin
tub_blot_image <- readImage(taz_blot_aTub_image)

# rotate image
tub_blot_image <- tub_blot_image |>
  rotate(91, bg.col = "white") |>
  flop()

# crop image
tub_blot_image <- tub_blot_image[3280:3932,1184:1406]

# adjust brightness and contrast
tub_blot_image <- (tub_blot_image - 0.3)

# get aspect ratio
H3K27me3_blot_dim <- dim(H3K27me3_blot_image)
H3K27me3_blot_aspect_ratio <- H3K27me3_blot_dim[2] / H3K27me3_blot_dim[1]

tub_blot_dim <- dim(tub_blot_image)
tub_blot_aspect_ratio <- tub_blot_dim[2] / tub_blot_dim[1]


# check that images were cropped to the same width
if (H3K27me3_blot_dim[1] != tub_blot_dim[1]) {stop("K27me3 and anti-tub blot are different widths")}


# place blot on page
plot_width <- 4

plotRaster(
  H3K27me3_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 0.5,
  width = plot_width,
  height =  plot_width * H3K27me3_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")

)

plotRaster(
  tub_blot_image,
  x = ref_x + 1.5,
  y = ref_y + 2.25,
  width = plot_width,
  height =  plot_width * tub_blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")

)

plotSegments(
  x0 = ref_x + 1.8, y0 = ref_y + 0.25, x1 = ref_x + 3.4, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 3.6, y0 = ref_y + 0.25, x1 = ref_x + 5.1, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

# panel label
plotText(
  label = "DMSO", params = large_text_params, fontface = "bold",
  x = ref_x + 2.7, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "taz", params = large_text_params, fontface = "bold",
  x = ref_x + 4.25, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "H3K27me3", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "tubulin", params = large_text_params, fontface = "bold",
  x = ref_x + 1.25, y = ref_y + 2.8, just = c("right","center"), default.units = "cm"
)

# panel B ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 5

plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in barcode counts
CR_spikeIn_counts <- read_tsv(CR_spikeIn_counts_fn)

AB_specificity <- CR_spikeIn_counts |>
  group_by(sample_name, target) |>
  summarise(barcode_count = sum(count))

barcode_totals <- CR_spikeIn_counts |>
  group_by(sample_name) |>
  summarise(total_count = sum(count))

AB_specificity <- AB_specificity |>
  left_join(barcode_totals) |>
  mutate(barcode_percent = barcode_count / total_count * 100)


sample_order <- c(
  "S2-Zld_DMSO_aH3K27me3_rep1",
  "S2-Zld_DMSO_aH3K27me3_rep2",
  "S2-Grh_DMSO_aH3K27me3_rep1",
  "S2-Grh_DMSO_aH3K27me3_rep2",
  "S2-HA-Twi_DMSO_aH3K27me3_rep1",
  "S2-HA-Twi_DMSO_aH3K27me3_rep2",

  "S2-Zld_Taz_aH3K27me3_rep1",
  "S2-Zld_Taz_aH3K27me3_rep2",
  "S2-Grh_Taz_aH3K27me3_rep1",
  "S2-Grh_Taz_aH3K27me3_rep2",
  "S2-HA-Twi_Taz_aH3K27me3_rep1",
  "S2-HA-Twi_Taz_aH3K27me3_rep2",

  "S2-Zld_DMSO_IgG_rep1",
  "S2-Zld_DMSO_IgG_rep2",
  "S2-Grh_DMSO_IgG_rep1",
  "S2-Grh_DMSO_IgG_rep2",
  "S2-HA-Twi_DMSO_IgG_rep1",
  "S2-HA-Twi_DMSO_IgG_rep2"
)

sample_labels <- sample_order |>
  tibble() |>
  separate(sample_order, into = c("cell_line", "treatment", "antibody", "replicate"), sep = "_") |>
  unite("label", cell_line, replicate) |>
  pull(label)

target_order <- c(
  "Unmodified",
  "H3K4me1",
  "H3K4me2",
  "H3K4me3",
  "H3K9me1",
  "H3K9me2",
  "H3K9me3",
  "H3K27me1",
  "H3K27me2",
  "H3K27me3",
  "H3K36me1",
  "H3K36me2",
  "H3K36me3",
  "H4K20me1",
  "H4K20me2",
  "H4K20me3"
)


# generate plot
b_plot <- AB_specificity |>
  mutate(sample_name = factor(sample_name, levels = sample_order)) |>
  mutate(sample_name = fct_rev(sample_name)) |>
  mutate(target = factor(target, levels = target_order)) |>
  ggplot(aes(x = target, y = sample_name, fill = barcode_percent)) + geom_tile() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    text = element_text(size = small_text_params$fontsize),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10)
  ) +
  scale_y_discrete(labels = sample_labels)

# place heatmap on page
plotGG(
  plot = b_plot,
  x = (ref_x + 1),
  y = (ref_y),
  width = 4.7,
  height = 4,
  just = c("left", "top"),
  default.units = "cm"
)

# add sample labels to heatmap
plotSegments(
  x0 = (ref_x + 1), y0 = (ref_y + 0.25), x1 = (ref_x + 1), y1 = (ref_y + 0.95),
  default.units = "cm"
)

plotSegments(
  x0 = (ref_x + 1), y0 = (ref_y + 1.05), x1 = (ref_x + 1), y1 = (ref_y + 1.85),
  default.units = "cm"
)

plotSegments(
  x0 = (ref_x + 0.25), y0 = (ref_y + 0.25), x1 = (ref_x + 0.25), y1 = (ref_y + 1.85),
  default.units = "cm"
)


plotSegments(
  x0 = (ref_x + 0.25), y0 = (ref_y + 1.9), x1 = (ref_x + 0.25), y1 = (ref_y + 2.7),
  default.units = "cm"
)

plotText(
  label = "DMSO",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x + 0.8,
  y = ref_y + 0.6,
  default.units = "cm"
)

plotText(
  label = "taz",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x + 0.8,
  y = ref_y + 1.4,
  default.units = "cm"
)

plotText(
  label = "anti-H3K27me3",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x,
  y = ref_y + 1,
  default.units = "cm"
)

plotText(
  label = "IgG",
  fontsize = small_text_params$fontsize,
  rot = 90,
  x = ref_x,
  y = ref_y + 2.2,
  default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 0.5

plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

regions <- zld_chip_classes |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


chrom_mark_averages <-bw |>
  map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |>
  map(rowMeans) |>
  as.data.frame()

p_data <- zld_chip_classes |>
  bind_cols(chrom_mark_averages)

c_plot <- p_data |>
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |>
  mutate(class = str_replace(class, "class ", "")) |>
  dplyr::select(class, H3K27ac:Mod_mdg4) |>
  pivot_longer(H3K27ac:Mod_mdg4, names_to = "chromatin_feature", values_to = "normalized_signal") |>
  group_by(class, chromatin_feature) |> summarise(normalized_signal = mean(normalized_signal)) |>

  mutate(feature_type = case_when(
    chromatin_feature %in% feature_types$active ~ "active",
    chromatin_feature %in% feature_types$Pol_II ~ "Pol_II",
    chromatin_feature %in% feature_types$histones ~ "histones",
    chromatin_feature %in% feature_types$polycomb ~ "polycomb",
    chromatin_feature %in% feature_types$heterochromatin ~ "heterochromatin",
    chromatin_feature %in% feature_types$TF ~ "TF",
    chromatin_feature %in% feature_types$insulator ~ "insulator",
  )) |>
  mutate(chromatin_feature = factor(chromatin_feature, levels = unlist(feature_types, use.names = FALSE))) |>
  ungroup() |>
  mutate(class = fct_rev(as.factor(class))) |>
  ggplot(aes(x=chromatin_feature, y = class, fill = normalized_signal)) +
  # facet_grid(feature_type ~ ., scales = "free") +
  geom_tile() +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,0)
  ) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(fill = paste0("normalized", "\n", "signal"))


plotGG(
  plot = c_plot,
  x = (ref_x + 0.2), y = (ref_y - 0.25),
  width = 10.6, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 3

plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


regions <- grh_chip_classes |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


chrom_mark_averages <-bw |>
  map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |>
  map(rowMeans) |>
  as.data.frame()

p_data <- grh_chip_classes |>
  bind_cols(chrom_mark_averages)

d_plot <- p_data |>
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |>
  mutate(class = str_replace(class, "class ", "")) |>
  dplyr::select(class, H3K27ac:Mod_mdg4) |>
  pivot_longer(H3K27ac:Mod_mdg4, names_to = "chromatin_feature", values_to = "normalized_signal") |>
  filter(!is.na(normalized_signal)) |>
  group_by(class, chromatin_feature) |> summarise(normalized_signal = mean(normalized_signal)) |>

  mutate(feature_type = case_when(
    chromatin_feature %in% feature_types$active ~ "active",
    chromatin_feature %in% feature_types$Pol_II ~ "Pol_II",
    chromatin_feature %in% feature_types$histones ~ "histones",
    chromatin_feature %in% feature_types$polycomb ~ "polycomb",
    chromatin_feature %in% feature_types$heterochromatin ~ "heterochromatin",
    chromatin_feature %in% feature_types$TF ~ "TF",
    chromatin_feature %in% feature_types$insulator ~ "insulator",
  )) |>
  mutate(chromatin_feature = factor(chromatin_feature, levels = unlist(feature_types, use.names = FALSE))) |>
  ungroup() |>
  mutate(class = fct_rev(as.factor(class))) |>
  ggplot(aes(x=chromatin_feature, y = class, fill = normalized_signal)) +
  # facet_grid(feature_type ~ ., scales = "free") +
  geom_tile() +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "right",
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,0)
  ) +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  labs(fill = paste0("normalized", "\n", "signal"))



plotGG(
  plot = d_plot,
  x = (ref_x + 0.2), y = (ref_y - 0.25),
  width = 10.6, height = 2.5, just = c("left", "top"),
  default.units = "cm"
)



# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 5.5

plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

regions <- twi_chip_classes |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


chrom_mark_averages <-bw |>
  map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |>
  map(rowMeans) |>
  as.data.frame()

p_data <- twi_chip_classes |>
  bind_cols(chrom_mark_averages)

e_plot <- p_data |>
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |>
  mutate(class = str_replace(class, "class ", "")) |>
  dplyr::select(class, H3K27ac:Mod_mdg4) |>
  pivot_longer(H3K27ac:Mod_mdg4, names_to = "chromatin_feature", values_to = "normalized_signal") |>
  group_by(class, chromatin_feature) |> summarise(normalized_signal = mean(normalized_signal)) |>

  mutate(feature_type = case_when(
    chromatin_feature %in% feature_types$active ~ "active",
    chromatin_feature %in% feature_types$Pol_II ~ "Pol_II",
    chromatin_feature %in% feature_types$histones ~ "histones",
    chromatin_feature %in% feature_types$polycomb ~ "polycomb",
    chromatin_feature %in% feature_types$heterochromatin ~ "heterochromatin",
    chromatin_feature %in% feature_types$TF ~ "TF",
    chromatin_feature %in% feature_types$insulator ~ "insulator",
  )) |>
  mutate(chromatin_feature = factor(chromatin_feature, levels = unlist(feature_types, use.names = FALSE))) |>
  ungroup() |>
  mutate(class = fct_rev(as.factor(class))) |>
  ggplot(aes(x=chromatin_feature, y = class, fill = normalized_signal)) +
  # facet_grid(feature_type ~ ., scales = "free") +
  geom_tile() +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.text.x = element_text(angle = 90, hjust=1),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "right",
    legend.justification = c(0.9),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,0)
  ) +
  scale_fill_distiller(palette = "GnBu", direction = 1) +
  labs(fill = paste0("normalized", "\n", "signal"))


plotGG(
  plot = e_plot,
  x = (ref_x + 0.2), y = (ref_y - 0.25),
  width = 10.6, height = 3.23, just = c("left", "top"),
  default.units = "cm"
)

# add labels to bottom of heatmap
plotSegments(
  x0 = ref_x + 0.6, y0 = ref_y + 3.25, x1 = ref_x + 4.2, y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 4.3, y0 = ref_y + 3.25, x1 = ref_x + 4.8, y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 4.9, y0 = ref_y + 3.25, x1 =ref_x + 7.1 , y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 7.2, y0 = ref_y + 3.25, x1 = ref_x + 7.8, y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 7.9, y0 = ref_y + 3.25, x1 = ref_x + 8.8, y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 8.9, y0 = ref_y + 3.25, x1 = ref_x + 10, y1 = ref_y + 3.25,
  default.units = "cm",
  lwd = 1
)

plotText(
  label = "active",
  params = large_text_params,
  x = ref_x + 2.45,
  y = ref_y + 3.5,
  just = c("center", "bottom"),
  default.units = "cm"
)

plotText(
  label = "histones",
  params = large_text_params,
  x = ref_x + 4.6,
  y = ref_y + 3.5,
  just = c("center", "bottom"),
  default.units = "cm"
)

plotText(
  label = "polycomb",
  params = large_text_params,
  x = ref_x + 5.8,
  y = ref_y + 3.5,
  just = c("center", "bottom"),
  default.units = "cm"
)

plotText(
  label = "heterochromatin",
  params = large_text_params,
  x = ref_x + 7.4,
  y = ref_y + 3.5,
  just = c("center", "bottom"),
  default.units = "cm"
)

plotText(
  label = "TFs",
  params = large_text_params,
  x = ref_x + 8.3,
  y = ref_y + 3.6,
  just = c("center", "top"),
  default.units = "cm"
)

plotText(
  label = "insulators",
  params = large_text_params,
  x = ref_x + 9.45,
  y = ref_y + 3.5,
  just = c("center", "bottom"),
  default.units = "cm"
)

# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 9.5

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

zld_ame_results <- zld_ame_results_fn |> 
  read_tsv() |> 
filter(rank < 25) |>
  arrange(`S2-WT_noCuSO4`) |> 
  mutate(gene_symbol = fct_inorder(gene_symbol)) |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |> 
  mutate(class = str_replace(class, "class ", ""))

# visualize enrichment results
f_motif_plot <- zld_ame_results |> 
  plot_ame_heatmap(group = class, id = gene_symbol, value = "normalize") +
  coord_flip() +
  scale_x_discrete(limits = unique(zld_ame_results$gene_symbol)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,1), breaks = c(0,1)) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom",
    legend.justification = 0.9
  )

# place plot on plotGardener page
plotGG(
  plot = f_motif_plot,
  x = (ref_x), y = (ref_y),
  width = 3.5, height = 9, just = c("left", "top"),
  default.units = "cm"
)

zld_motif_factor_RNAseq <- zld_motif_factor_RNAseq_fn |> 
  read_tsv()

# plot heatmap alongside motif plot
hm_data <- zld_motif_factor_RNAseq |> 
  filter(gene_id %in% zld_ame_results$gene_id) |> 
  mutate(gene_symbol = factor(gene_symbol, levels = unique(zld_ame_results$gene_symbol))) |> 
  select(-gene_id, -`embryo_cycle-14C`) |> 
  pivot_longer(2:4, names_to = "sample", values_to = "RPKM") |> 
  mutate(RPKM = log2(RPKM + 0.01))

f_RNAseq_heatmap <- hm_data |> 
  mutate(sample = replace(sample, sample == "S2-WT_noCuSO4", "S2-WT")) |>
  mutate(sample = replace(sample, sample == "embryo_cycle-14D", "embryo")) |>
  mutate(sample = replace(sample, sample == "neuroblast", "NSC")) |>
  mutate(sample = factor(sample, levels = c("S2-WT","embryo", "NSC"))) |>
  dplyr::rename(tissue = sample) |> 
  ggplot(aes(x = tissue, y = gene_symbol, fill = RPKM)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1, ) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom"
  )

# place plot on plotGardener page
plotGG(
  plot = f_RNAseq_heatmap,
  x = (ref_x + 3.5), y = (ref_y),
  width = 2, height = 9, just = c("left", "top"),
  default.units = "cm"
) 

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 9.5

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

grh_ame_results <- grh_ame_results_fn |> 
  read_tsv() |> 
  filter(rank < 25) |>
  arrange(`S2-WT_noCuSO4`) |> 
  mutate(gene_symbol = fct_inorder(gene_symbol)) |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |> 
  mutate(class = str_replace(class, "class ", ""))

# visualize enrichment results
g_motif_plot <- grh_ame_results |> 
  plot_ame_heatmap(group = class, id = gene_symbol, value = "normalize") +
  coord_flip() +
  scale_x_discrete(limits = unique(grh_ame_results$gene_symbol)) +
  scale_fill_distiller(palette = "Oranges", direction = 1, limits = c(0,1), breaks = c(0,1)) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom",
    legend.justification = 0.9
  )

# place plot on plotGardener page
plotGG(
  plot = g_motif_plot,
  x = (ref_x), y = (ref_y),
  width = 3.5, height = 9, just = c("left", "top"),
  default.units = "cm"
)


# import RNA-seq data
grh_motif_factor_RNAseq <- grh_motif_factor_RNAseq_fn |> 
  read_tsv()

# plot RNA-seq heatmap alongside motifs
hm_data <- grh_motif_factor_RNAseq |> 
  filter(gene_id %in% grh_ame_results$gene_id) |> 
  mutate(gene_symbol = factor(gene_symbol, levels = unique(grh_ame_results$gene_symbol))) |> 
  select(-gene_id) |> 
  pivot_longer(2:4, names_to = "sample", values_to = "RPKM") |> 
  mutate(RPKM = log2(RPKM + 0.01))

g_RNAseq_heatmap <- hm_data |> 
  mutate(sample = replace(sample, sample == "S2-WT_noCuSO4", "S2-WT")) |>
  mutate(sample = replace(sample, sample == "embryo_15-16H", "embryo")) |>
  mutate(sample = replace(sample, sample == "wing-disc_L3", "wing disc")) |>
  mutate(sample = factor(sample, levels = c("S2-WT","embryo", "wing disc"))) |>
  dplyr::rename(tissue = sample) |> 
  ggplot(aes(x = tissue, y = gene_symbol, fill = RPKM)) +
  geom_tile() +
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom"
  )

# place plot on plotGardener page
plotGG(
  plot = g_RNAseq_heatmap,
  x = (ref_x + 3.5), y = (ref_y),
  width = 2, height = 9, just = c("left", "top"),
  default.units = "cm"
) 


# panel H ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 9.5

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


twi_ame_results <- twi_ame_results_fn |> 
  read_tsv() |> 
  filter(rank < 25) |>
  arrange(`S2-WT_noCuSO4`) |> 
  mutate(gene_symbol = fct_inorder(gene_symbol)) |> 
  mutate(class = replace(class, class == "repressed_H3K27me3", "class IV")) |>
  mutate(class = replace(class, class == "repressed_H3K9me3", "class V")) |>
  mutate(class = replace(class, class == "repressed_other", "class VI")) |> 
  mutate(class = str_replace(class, "class ", ""))

# visualize enrichment results
h_motif_plot <- twi_ame_results |> 
  plot_ame_heatmap(group = class, id = gene_symbol, value = "normalize") +
  coord_flip() +
  scale_x_discrete(limits = unique(twi_ame_results$gene_symbol)) +
  scale_fill_distiller(palette = "GnBu", direction = 1, limits = c(0,1), breaks = c(0,1)) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom",
    legend.justification = 0.9
  )


# place plot on plotGardener page
plotGG(
  plot = h_motif_plot,
  x = (ref_x), y = (ref_y),
  width = 3.5, height = 9, just = c("left", "top"),
  default.units = "cm"
)

# import RNAseq data
twi_motif_factor_RNAseq <- twi_motif_factor_RNAseq_fn |> 
  read_tsv()

# plot RNA-seq heatmap alongside motifs
hm_data <- twi_motif_factor_RNAseq |> 
  filter(gene_id %in% twi_ame_results$gene_id) |> 
  mutate(gene_symbol = factor(gene_symbol, levels = unique(twi_ame_results$gene_symbol))) |> 
  select(-gene_id, -`embryo_cycle-14C`) |> 
  pivot_longer(2:3, names_to = "sample", values_to = "RPKM") |> 
  mutate(RPKM = log2(RPKM + 0.01))

h_RNAseq_heatmap <- hm_data |> 
  mutate(sample = replace(sample, sample == "S2-WT_noCuSO4", "S2-WT")) |>
  mutate(sample = replace(sample, sample == "embryo_cycle-14D", "embryo")) |>
  mutate(sample = factor(sample, levels = c("S2-WT","embryo"))) |>
  dplyr::rename(tissue = sample) |> 
  ggplot(aes(x = tissue, y = gene_symbol, fill = RPKM)) +
  geom_tile() +
  scale_fill_distiller(palette = "GnBu", direction = 1) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.key.size = unit(2, 'mm'),
    legend.position = "bottom"
  )

# place plot on plotGardener page
plotGG(
  plot = h_RNAseq_heatmap,
  x = (ref_x + 3.5), y = (ref_y),
  width = 2, height = 9, just = c("left", "top"),
  default.units = "cm"
) 



# close graphics device ========================================================
dev.off()

