# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(rtracklayer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# define input files explicitly for testing
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
# zld_chip_classes <- read_tsv("results/ChIP_peak_classes/zld_ChIP_classes.tsv")
# grh_chip_classes <- read_tsv("results/ChIP_peak_classes/grh_ChIP_classes.tsv")
# twi_chip_classes <- read_tsv("results/ChIP_peak_classes/twi_ChIP_classes.tsv")
# 
# Zld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# Zld_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw"
# Zld_Zld_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw"
# 
# Grh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# Grh_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw"
# Grh_Grh_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw"
# 
# Twi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# Twi_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw"
# Twi_Twi_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Twi_40uM_small.bw"
# 
# H3K27me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/individual/GSE151983_S2_aH3K27me3_IP.bw"
# H3K9me3_ChIP_bw <- "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw"

# get input files from snakemake
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

Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]
Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
Zld_Zld_ATAC_bw <- snakemake@input[["Zld_Zld_ATAC_bw"]]

Grh_ChIP_bw <- snakemake@input[["Grh_ChIP_bw"]]
Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
Grh_Grh_ATAC_bw <-snakemake@input[["Grh_Grh_ATAC_bw"]]

Twi_ChIP_bw <- snakemake@input[["Twi_ChIP_bw"]]
Twi_WT_ATAC_bw <- snakemake@input[["Twi_WT_ATAC_bw"]]
Twi_Twi_ATAC_bw <- snakemake@input[["Twi_Twi_ATAC_bw"]]

H3K27me3_ChIP_bw <- snakemake@input[["H3K27me3_ChIP_bw"]]
H3K9me3_ChIP_bw <- snakemake@input[["H3K9me3_ChIP_bw"]]

# create blank layout for plot ===============================================
fig_width <-  18
fig_height <- 12.5

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig4.pdf", useDingbats = FALSE)

# generate plotGardener page
pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)


# general figure settings ======================================================
# text parameters
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
  zld_heatmap_colors <- brewer.pal(9, "Blues")
  grh_color <- "#F98400"
    grh_heatmap_colors <- brewer.pal(9, "Oranges")
    
    twi_color <- "#00A08A"
      # twi_heatmap_colors <- c("white","#00A08A", "#154734")
    twi_heatmap_colors <- brewer.pal(9, "GnBu")
    
    H3K27me3_color <- "gray40"      
      
    # reference points for positioning figure components
    x_offset_class_label <- 0.25
    x_offset_browser_label <- 1
    x_offset_browser <- 1.5
    
    # set genome browser height
    gb_height <- 0.3
    
    
    
    # panel A ======================================================================
    # reference points for positioning figure components
    ref_x <- 0.5
    ref_y <- 0.5
    
    # panel label
    plotText(
      label = "a",
      params = panel_label_params,
      fontface = "bold",
      x = ref_x,
      y = ref_y,
      just = "bottom",
      default.units = "cm"
    )
    
    
    regions <- zld_chip_classes |>
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    
    chrom_mark_averages <-bw |> 
      map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |> 
      map(rowMeans) |> 
      as.data.frame()
    
    p_data <- zld_chip_classes |> 
      bind_cols(chrom_mark_averages)
    
    a_plot <- p_data |> 
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
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            axis.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.key.size = unit(2, 'mm'),
            legend.position = "top",
            legend.justification = c(0.9),
            legend.margin=margin(0,0,4,0),
            legend.box.margin=margin(-10,-10,-10,-10)
      ) +
      scale_fill_distiller(palette = "Blues", direction = 1)
    
    
    plotGG(
      plot = a_plot,
      x = (ref_x), y = ref_y + 0.5,
      width = 10, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 3, x1 = ref_x + 3.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 3, x1 = ref_x + 4.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 3, x1 =ref_x + 6.5 , y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 3, x1 = ref_x + 7.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 3, x1 = ref_x + 8.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 3, x1 = ref_x + 9.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "histones",
      params = large_text_params,
      x = ref_x + 4.2,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 3.3,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    # panel B ======================================================================
    # panel label
    ref_x <- 11
    ref_y <- 0.5
    
    plotText(
      label = "b", params = panel_label_params, fontface = "bold",
      x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    plotText(
      label = "class III", params = large_text_params, fontface = "bold",
      x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    # add broad H3K27me3 region
    plotText(
      label = paste("H3K27me3"),
      params = small_text_params,
      fontface = "bold",
      x = (ref_x + x_offset_browser_label),
      y = (ref_y + 0.25),
      just = c("right", "center"),
      default.units = "cm"
    )
    
    region <- pgParams(
      chrom = "chr3R",
      chromstart = 6647050, chromend = 6823266,
      assembly = "dm6"
    )
    
    `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                  params = region
    )
    
    ChIP_range <- signal_range(c(H3K27me3_ChIP$score))
    
    s1 <- plotSignal(
      data = `H3K27me3_ChIP`, params = region,
      x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s1, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    # add zoomed in view of ChIP and ATAC
    zoomRegion <- pgParams(
      chrom = "chr3R",
      chromstart = 6751637, chromend = 6760890,
      assembly = "dm6"
    )
    
    annoZoomLines(
      plot = s1, params = zoomRegion,
      y0 = (s1$y + s1$height / 2), x1 = c(s1$x, s1$x + s1$width), y1 = (ref_y + 0.9), extend = c(gb_height, 1.5),
      default.units = "cm",
      lty = 2
    )
    
    # plotText(
    #   label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
    #   x = (ref_x + x_offset_browser_label), y = (ref_y + 0.75), just = c("right","center"), default.units = "cm"
    # )
    # 
    
    plotText(
      label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
    )
    
    
    plotText(
      label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
    )
    
    plotText(
      label = "ZLD ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
    )
    
    
    
    
    
    `Zld_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                             params = zoomRegion
    )
    
    
    `S2-WT_ATAC` <- readBigwig(file = Zld_WT_ATAC_bw, 
                               params = zoomRegion
    )
    
    `S2-Zld_ATAC` <- readBigwig(file = Zld_Zld_ATAC_bw, 
                                params = zoomRegion
    )
    
    ChIP_range <- signal_range(c(Zld_ChIP$score))
    ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Zld_ATAC`$score))
    
    
    s2 <- plotSignal(
      data = `Zld_ChIP`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s2, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s3 <- plotSignal(
      data = `S2-WT_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s3, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s4 <- plotSignal(
      data = `S2-Zld_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s4, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    plotGenes(
      params = zoomRegion, assembly = "dm6",
      x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
      fill = c("gray50", "gray50"),
      fontcolor = c("gray50", "gray50"),
      just = c("left", "top"),
      default.units = "cm"
    )
    
    annoHighlight(
      plot = s2,
      chrom = "chr3R",
      chromstart = 6755358,
      chromend = 6755908,
      y = (ref_y + 1), height = 1.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    
    # panel C ======================================================================
    # reference points for positioning figure components
    ref_x <- 0.5
    ref_y <- 4.5
    
    # panel label
    plotText(
      label = "c",
      params = panel_label_params,
      fontface = "bold",
      x = ref_x,
      y = ref_y,
      just = "bottom",
      default.units = "cm"
    )
    
    
    regions <- grh_chip_classes |>
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    
    chrom_mark_averages <-bw |> 
      map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |> 
      map(rowMeans) |> 
      as.data.frame()
    
    p_data <- grh_chip_classes |> 
      bind_cols(chrom_mark_averages)
    
    c_plot <- p_data |> 
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
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            axis.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.key.size = unit(2, 'mm'),
            legend.position = "top",
            legend.justification = c(0.9),
            legend.margin=margin(0,0,4,0),
            legend.box.margin=margin(-10,-10,-10,-10)
      ) +
      scale_fill_distiller(palette = "Oranges", direction = 1)
    
    
    plotGG(
      plot = c_plot,
      x = (ref_x), y = ref_y + 0.5,
      width = 10, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 3, x1 = ref_x + 3.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 3, x1 = ref_x + 4.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 3, x1 =ref_x + 6.5 , y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 3, x1 = ref_x + 7.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 3, x1 = ref_x + 8.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 3, x1 = ref_x + 9.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "histones",
      params = large_text_params,
      x = ref_x + 4.2,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 3.3,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    # panel D ======================================================================
    # panel label
    ref_x <- 11
    ref_y <- 4.5
    
    plotText(
      label = "d", params = panel_label_params, fontface = "bold",
      x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    plotText(
      label = "class III", params = large_text_params, fontface = "bold",
      x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    # add broad H3K27me3 region
    plotText(
      label = paste("H3K27me3"),
      params = small_text_params,
      fontface = "bold",
      x = (ref_x + x_offset_browser_label),
      y = (ref_y + 0.25),
      just = c("right", "center"),
      default.units = "cm"
    )
    
    region <- pgParams(
      chrom = "chr2L",
      chromstart = 14042440, chromend = 14236112,
      assembly = "dm6"
    )
    
    `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                  params = region
    )
    
    ChIP_range <- signal_range(c(H3K27me3_ChIP$score))
    
    s1 <- plotSignal(
      data = `H3K27me3_ChIP`, params = region,
      x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s1, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    # add zoomed in view of ChIP and ATAC
    zoomRegion <- pgParams(
      chrom = "chr2L",
      chromstart = 14100727, chromend = 14113159,
      assembly = "dm6"
    )
    
    annoZoomLines(
      plot = s1, params = zoomRegion,
      y0 = (s1$y + s1$height / 2), x1 = c(s1$x, s1$x + s1$width), y1 = (ref_y + 0.9), extend = c(gb_height, 1.5),
      default.units = "cm",
      lty = 2
    )
    
    
    plotText(
      label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
    )
    
    
    plotText(
      label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
    )
    
    plotText(
      label = "Grh ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
    )
    
    
    
    
    
    `Grh_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                             params = zoomRegion
    )
    
    
    `S2-WT_ATAC` <- readBigwig(file = Grh_WT_ATAC_bw, 
                               params = zoomRegion
    )
    
    `S2-Grh_ATAC` <- readBigwig(file = Grh_Grh_ATAC_bw, 
                                params = zoomRegion
    )
    
    ChIP_range <- signal_range(c(Grh_ChIP$score))
    ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Grh_ATAC`$score))
    
    
    s2 <- plotSignal(
      data = `Grh_ChIP`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s2, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s3 <- plotSignal(
      data = `S2-WT_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s3, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s4 <- plotSignal(
      data = `S2-Grh_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s4, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    plotGenes(
      params = zoomRegion, assembly = "dm6",
      x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
      fill = c("gray50", "gray50"),
      fontcolor = c("gray50", "gray50"),
      just = c("left", "top"),
      default.units = "cm"
    )
    
    
    # panel E ======================================================================
    # reference points for positioning figure components
    ref_x <- 0.5
    ref_y <- 8.5
    
    # panel label
    plotText(
      label = "e",
      params = panel_label_params,
      fontface = "bold",
      x = ref_x,
      y = ref_y,
      just = "bottom",
      default.units = "cm"
    )
    
    
    regions <- twi_chip_classes |>
      makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    
    chrom_mark_averages <-bw |> 
      map(coverage_matrix, regions = regions, upstream = 500, downstream = 500, region_width = 1000) |> 
      map(rowMeans) |> 
      as.data.frame()
    
    p_data <- twi_chip_classes |> 
      bind_cols(chrom_mark_averages)
    
    b_plot <- p_data |> 
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
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            axis.title = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.key.size = unit(2, 'mm'),
            legend.position = "top",
            legend.justification = c(0.9),
            legend.margin=margin(0,0,4,0),
            legend.box.margin=margin(-10,-10,-10,-10)
      ) +
      scale_fill_distiller(palette = "GnBu", direction = 1)
    
    
    plotGG(
      plot = b_plot,
      x = (ref_x), y = ref_y + 0.5,
      width = 10, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 3, x1 = ref_x + 3.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 3, x1 = ref_x + 4.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 3, x1 =ref_x + 6.5 , y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 3, x1 = ref_x + 7.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 3, x1 = ref_x + 8.5, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 3, x1 = ref_x + 9.75, y1 = ref_y + 3,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "histones",
      params = large_text_params,
      x = ref_x + 4.2,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 3.3,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 3.25,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
  
    
    # panel F ======================================================================
    # panel label
    ref_x <- 11
    ref_y <- 8.5
    
    plotText(
      label = "f", params = panel_label_params, fontface = "bold",
      x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    plotText(
      label = "class III", params = large_text_params, fontface = "bold",
      x = ref_x + 3, y = ref_y, just = "bottom", default.units = "cm"
    )
    
    # add broad H3K27me3 region
    plotText(
      label = paste("H3K27me3"),
      params = small_text_params,
      fontface = "bold",
      x = (ref_x + x_offset_browser_label),
      y = (ref_y + 0.25),
      just = c("right", "center"),
      default.units = "cm"
    )
    
    region <- pgParams(
      chrom = "chr2R",
      chromstart = 24259395, chromend = 24314647,
      assembly = "dm6"
    )
    
    `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                  params = region
    )
    
    ChIP_range <- signal_range(c(H3K27me3_ChIP$score))
    
    
    s1 <- plotSignal(
      data = `H3K27me3_ChIP`, params = region,
      x = (ref_x + x_offset_browser), y = (ref_y + 0.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s1, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    # add zoomed in view of ChIP and ATAC
    zoomRegion <- pgParams(
      chrom = "chr2R",
      chromstart = 24277315, chromend = 24284445,
      assembly = "dm6"
    )
    
    annoZoomLines(
      plot = s1, params = zoomRegion,
      y0 = (s1$y + s1$height / 2), x1 = c(s1$x, s1$x + s1$width), y1 = (ref_y + 0.9), extend = c(gb_height, 1.5),
      default.units = "cm",
      lty = 2
    )
    
    
    plotText(
      label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.25), just = c("right","center"), default.units = "cm"
    )
    
    
    plotText(
      label = "S2 WT ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
    )
    
    plotText(
      label = "Twi ATAC", params = small_text_params, fontface = "bold",
      x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
    )
    
    
    
    
    
    `Twi_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                             params = zoomRegion
    )
    
    
    `S2-WT_ATAC` <- readBigwig(file = Twi_WT_ATAC_bw, 
                               params = zoomRegion
    )
    
    `S2-Twi_ATAC` <- readBigwig(file = Twi_Twi_ATAC_bw, 
                                params = zoomRegion
    )
    
    ChIP_range <- signal_range(c(Twi_ChIP$score))
    ATAC_range <- signal_range(c(`S2-WT_ATAC`$score, `S2-Twi_ATAC`$score))
    
    
    s2 <- plotSignal(
      data = `Twi_ChIP`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
      range = ChIP_range
    )
    
    annoYaxis(
      plot = s2, at = round(ChIP_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s3 <- plotSignal(
      data = `S2-WT_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s3, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    
    
    s4 <- plotSignal(
      data = `S2-Twi_ATAC`, params = zoomRegion,
      x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
      just = c("left", "center"), default.units = "cm",
      linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
      range = ATAC_range
    )
    
    annoYaxis(
      plot = s4, at = round(ATAC_range, 1),
      axisLine = TRUE, fontsize = 5, lwd = 0.5, 
    )
    
    plotGenes(
      params = zoomRegion, assembly = "dm6",
      x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
      fill = c("gray50", "gray50"),
      fontcolor = c("gray50", "gray50"),
      just = c("left", "top"),
      default.units = "cm"
    )
    
    
    
    
    # close graphics device ========================================================
    dev.off()
    
    
    
    
    