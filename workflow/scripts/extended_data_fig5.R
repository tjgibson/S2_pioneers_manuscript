# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(rtracklayer)

source("workflow/scripts/plot_heatmap.R")

# define input files ===========================================================
# Zld_ChIP_bw <- snakemake@input[["Zld_ChIP_bw"]]

bw <- c(
  H3K27ac = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw",
  Nej = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw",
  H3K4me1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw",
  H3K4me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw",
  H4K16ac = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE37865_aH4K16ac.bw",
  H2AV = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw",
  
  Rpb3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Rpb3_IP.bw",
  `PolII-pS2` = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_PolII_phosphoSer2_IP.bw",
  H3K36me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3K36me3_IP.bw",
  SSRP1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_SSRP1_IP.bw",
  Spt16 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Spt16_IP.bw",
  
  H3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3_IP.bw",
  H1 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE127227_aH1.bw",
  
  Pho = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE84502_aPho_IP.bw",
  Ez = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_Ez_IP.bw",
  Pc = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Pc_IP.bw",
  Psc = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Psc_IP.bw",
  Ph = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Ph_IP.bw",
  dRing = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_dRing_IP.bw",
  H3K27me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE151983_S2_aH3K27me3_IP.bw",
  
  H3K9me3 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw",
  HP1a = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE56101_aHP1a.bw",
  
  M1BP = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_M1BP_IP.bw",
  GAF = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_GAGA_IP.bw",
  BEAF32 =
    "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE52887_aBEAF32_IP.bw",
  
  CTCF = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCTCF.bw",
  CP190 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCP190_IP.bw",
  Su_hw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aSu-Hw_IP.bw",
  Mod_mdg4 = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aMod-mdg4.bw"
  
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

zld_chip_classes <- read_tsv("results/ChIP_peak_classes/zld_ChIP_classes.tsv")
grh_chip_classes <- read_tsv("results/ChIP_peak_classes/grh_ChIP_classes.tsv")
twi_chip_classes <- read_tsv("results/ChIP_peak_classes/twi_ChIP_classes.tsv")


# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig5.pdf", useDingbats = FALSE)
pageCreate(
  width = 18,
  height = 12.5,
  default.units = "cm",
  showGuides = FALSE
)

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
      # twi_heatmap_colors <- c("white","#00A08A", "#154734")
    twi_heatmap_colors <- brewer.pal(9, "GnBu")
    
  
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
            legend.key.size = unit(2, 'mm')
            ) +
      scale_fill_distiller(palette = "Blues", direction = 1)
    
    
    plotGG(
      plot = a_plot,
      x = (ref_x), y = ref_y,
      width = 12, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 2.5, x1 = ref_x + 3.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 2.5, x1 = ref_x + 4.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 2.5, x1 =ref_x + 6.5 , y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 2.5, x1 = ref_x + 7.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 2.5, x1 = ref_x + 8.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 2.5, x1 = ref_x + 9.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )

    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 2.8,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    # panel B ======================================================================
    # reference points for positioning figure components
    ref_x <- 0.5
    ref_y <- 3.5
    
    # panel label
    plotText(
      label = "b",
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
            legend.key.size = unit(2, 'mm')
      ) +
      scale_fill_distiller(palette = "Oranges", direction = 1)
    
    
    plotGG(
      plot = b_plot,
      x = (ref_x), y = ref_y,
      width = 12, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 2.5, x1 = ref_x + 3.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 2.5, x1 = ref_x + 4.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 2.5, x1 =ref_x + 6.5 , y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 2.5, x1 = ref_x + 7.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 2.5, x1 = ref_x + 8.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 2.5, x1 = ref_x + 9.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 2.8,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    # panel C ======================================================================
    # reference points for positioning figure components
    ref_x <- 0.5
    ref_y <- 6.5
    
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
            legend.key.size = unit(2, 'mm')
      ) +
      scale_fill_distiller(palette = "GnBu", direction = 1)
    
    
    plotGG(
      plot = b_plot,
      x = (ref_x), y = ref_y,
      width = 12, height = 2.5, just = c("left", "top"),
      default.units = "cm"
    )
    
    # add labels to bottom of heatmap
    plotSegments(
      x0 = ref_x, y0 = ref_y + 2.5, x1 = ref_x + 3.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4, y0 = ref_y + 2.5, x1 = ref_x + 4.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 4.75, y0 = ref_y + 2.5, x1 =ref_x + 6.5 , y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 6.75, y0 = ref_y + 2.5, x1 = ref_x + 7.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 7.75, y0 = ref_y + 2.5, x1 = ref_x + 8.5, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotSegments(
      x0 = ref_x + 8.75, y0 = ref_y + 2.5, x1 = ref_x + 9.75, y1 = ref_y + 2.5,
      default.units = "cm",
      lwd = 1
    )
    
    plotText(
      label = "active",
      params = large_text_params,
      x = ref_x + 1.875,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "polycomb",
      params = large_text_params,
      x = ref_x + 5.5,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "heterochromatin",
      params = large_text_params,
      x = ref_x + 7.1,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    plotText(
      label = "TFs",
      params = large_text_params,
      x = ref_x + 8.1,
      y = ref_y + 2.8,
      just = c("center", "top"),
      default.units = "cm"
    )
    
    plotText(
      label = "insulators",
      params = large_text_params,
      x = ref_x + 9.25,
      y = ref_y + 2.75,
      just = c("center", "bottom"),
      default.units = "cm"
    )
    
    
    # close graphics device ========================================================
    dev.off()
    
    
    
    
    