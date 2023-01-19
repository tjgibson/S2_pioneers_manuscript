# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
# Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
# Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]
# Zld_WT_RNAseq_bw <- snakemake@input[["Zld_WT_RNAseq_bw"]]
# Zld_Zld_RNAseq_bw <- snakemake@input[["Zld_Zld_RNAseq_bw"]]
# 
# Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
# Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
# Grh_Grh_ATAC_bw <- snakemake@input[["Grh_Grh_ATAC_bw"]]
# Grh_WT_RNAseq_bw <- snakemake@input[["Grh_WT_RNAseq_bw"]]
# Grh_Grh_RNAseq_bw <- snakemake@input[["Grh_Grh_RNAseq_bw"]]
# 
# zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes"]]
# grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes"]]

Zld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
Zld_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw"
Zld_Zld_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw"

Zld_nc14_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/individual/embryo-nc14_aZld.bw"

Grh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
Grh_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw"
Grh_Grh_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw"

Grh_embryo_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw"

Twi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
Twi_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw"
Twi_Twi_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Twi_40uM_small.bw"

Twi_embryo_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw"

H3K27me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/individual/GSE151983_S2_aH3K27me3_IP.bw"
H3K9me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw"

CR_spikeIn_counts_fn <-"CUTandRUN/results/scaling_factors/epiCypher_barcode_counts.tsv"


# open graphics device =========================================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig8.pdf")
# create blank layout for plot =================================================
pageCreate(width = 18, height = 18, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
  grh_color <- "#F98400"
    twi_color <- "#00A08A"
      
    H3K27me3_color <- "gray40"      
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
      
  

# panel M ==================================================================
# panel label
ref_x <- 0.5
ref_y <- 16.5

plotText(
  label = "m", params = panel_label_params, fontface = "bold",
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
m_plot <- AB_specificity |>
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
  plot = m_plot,
  x = (ref_x + 1),
  y = (ref_y),
  width = 4.5,
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
      
      # close graphics device ========================================================
      dev.off()
      
      