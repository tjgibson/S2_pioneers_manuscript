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
# pdf("manuscript/figures/extended_data_fig8.pdf")
# create blank layout for plot =================================================
pageCreate(width = 18, height = 21, default.units = "cm", showGuides = TRUE)

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
      
      # panel A ======================================================================
      # panel label
      ref_x <- 0.5
      ref_y <- 0.5
      
      plotText(
        label = "a", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      # panel B ======================================================================
      # panel label
      ref_x <- 6.5
      ref_y <- 0.5
      
      plotText(
        label = "b", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      # panel C ======================================================================
      # panel label
      ref_x <- 12.5
      ref_y <- 0.5
      
      plotText(
        label = "c", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      
      
      # panel D ==================================================================
      # panel label
      ref_x <- 0.5
      ref_y <- 4.5
      
      plotText(
        label = "d", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chr3R",
        chromstart = 30824716, chromend = 30888891,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr3R",
        chromstart = 30846356, chromend = 30858282,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                                    params = zoomRegion
      )
      
      `Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Zld_ChIP$score, Zld_nc14_ChIP$score))
      
      s2 <- plotSignal(
        data = `Zld_nc14_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Zld_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      # panel E ==================================================================
      # panel label
      ref_x <- 6.5
      ref_y <- 4.5
      
      plotText(
        label = "e", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chr3L",
        chromstart = 11796001, chromend = 12064216,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr3L",
        chromstart = 11919547, chromend = 11940115,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                                    params = zoomRegion
      )
      
      `Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Zld_ChIP$score, Zld_nc14_ChIP$score))
      
      s2 <- plotSignal(
        data = `Zld_nc14_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Zld_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      # panel F ==================================================================
      # panel label
      ref_x <- 12.5
      ref_y <- 4.5
      
      plotText(
        label = "f", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chrX",
        chromstart = 393315, chromend = 402092,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- c(-1.2,8)
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chrX",
        chromstart = 393449, chromend = 400413,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Zld embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Zld S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Zld_nc14_ChIP` <- readBigwig(file = Zld_nc14_ChIP_bw, 
                                    params = zoomRegion
      )
      
      `Zld_S2_ChIP` <- readBigwig(file = Zld_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Zld_ChIP$score, Zld_nc14_ChIP$score))
      
      s2 <- plotSignal(
        data = `Zld_nc14_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Zld_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = zld_color, fill = zld_color, baseline = FALSE, baseline.lwd = 0, baseline.color = zld_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      # panel G ==================================================================
      # panel label
      ref_x <- 0.5
      ref_y <- 8.5
      
      plotText(
        label = "g", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chr3L",
        chromstart = 3232675, chromend = 3324202,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr3L",
        chromstart = 3254519, chromend = 3272949,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Grh_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Grh_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      # panel H ==================================================================
      # panel label
      ref_x <- 6.5
      ref_y <- 8.5
      
      plotText(
        label = "h", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chr2L",
        chromstart = 4443392, chromend = 4669912,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr2L",
        chromstart = 4562557, chromend = 4626155,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Grh_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Grh_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      # panel I ==================================================================
      # panel label
      ref_x <- 12.5
      ref_y <- 8.5
      
      plotText(
        label = "i", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      # region candidates:
      # chr3R:22,748,753-22,807,827
      # chrX:13,216,271-13,269,214
      
      region <- pgParams(
        chrom = "chrX",
        chromstart = 13181251, chromend = 13295616,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- c(-1.1,8)
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chrX",
        chromstart = 13218840, chromend = 13238073,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Grh embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Grh S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Grh_embryo_ChIP` <- readBigwig(file = Grh_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Grh_S2_ChIP` <- readBigwig(file = Grh_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Grh_S2_ChIP$score, Grh_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Grh_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Grh_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = grh_color, fill = grh_color, baseline = FALSE, baseline.lwd = 0, baseline.color = grh_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      # panel J ==================================================================
      # panel label
      ref_x <- 0.5
      ref_y <- 12.5
      
      plotText(
        label = "j", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      region <- pgParams(
        chrom = "chr2R",
        chromstart = 12753782, chromend = 12796150,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr2R",
        chromstart = 12773275, chromend = 12791478,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Twi_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Twi_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      # panel K ==================================================================
      # panel label
      ref_x <- 6.5
      ref_y <- 12.5
      
      plotText(
        label = "k", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      # candidate regions
      # chrX:18,315,069-18,367,911
      # chrX:17891901-18095637
      # chr3L:13,423,435-14,067,450
      
      region <- pgParams(
        chrom = "chr3L",
        chromstart = 13423435, chromend = 14067450,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- signal_range(c(H3K27me3_ChIP$score, H3K9me3_ChIP$score))
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr3L",
        chromstart = 13658008, chromend = 13683890,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Twi_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Twi_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
      
      # panel L ==================================================================
      # panel label
      ref_x <- 12.5
      ref_y <- 12.5
      
      plotText(
        label = "l", params = panel_label_params, fontface = "bold",
        x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
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
      
      plotText(
        label = paste("H3K9me3"),
        params = small_text_params,
        fontface = "bold",
        x = (ref_x + x_offset_browser_label),
        y = (ref_y + 0.75),
        just = c("right", "center"),
        default.units = "cm"
      )
      
      
      region <- pgParams(
        chrom = "chr3L",
        chromstart = 10650450, chromend = 10684646,
        assembly = "dm6"
      )
      
      `H3K27me3_ChIP` <- readBigwig(file = H3K27me3_ChIP_bw, 
                                    params = region
      )
      
      `H3K9me3_ChIP` <- readBigwig(file = H3K9me3_ChIP_bw, 
                                   params = region
      )
      
      ChIP_range <- c(-1.1,8)
      
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
      
      s2 <- plotSignal(
        data = `H3K9me3_ChIP`, params = region,
        x = (ref_x + x_offset_browser), y = (ref_y + 0.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = H3K27me3_color, fill = H3K27me3_color, baseline = FALSE, baseline.lwd = 0, baseline.color = H3K27me3_color,
        range = c(-0.5, 10)
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      # add zoomed in view of ChIP
      zoomRegion <- pgParams(
        chrom = "chr3L",
        chromstart = 10667170, chromend = 10675448,
        assembly = "dm6"
      )
      
      annoZoomLines(
        plot = s2, params = zoomRegion,
        y0 = (s2$y + s2$height / 2), x1 = c(s2$x, s2$x + s2$width), y1 = (ref_y + 0.9), extend = c(1, 1),
        default.units = "cm",
        lty = 2
      )
      
      plotText(
        label = "Twi embryo ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 1.75), just = c("right","center"), default.units = "cm"
      )
      
      
      plotText(
        label = "Twi S2 ChIP", params = small_text_params, fontface = "bold",
        x = (ref_x + x_offset_browser_label), y = (ref_y + 2.25), just = c("right","center"), default.units = "cm"
      )
      
      
      
      `Twi_embryo_ChIP` <- readBigwig(file = Twi_embryo_ChIP_bw, 
                                      params = zoomRegion
      )
      
      `Twi_S2_ChIP` <- readBigwig(file = Twi_ChIP_bw, 
                                  params = zoomRegion
      )
      
      
      
      ChIP_range <- signal_range(c(Twi_S2_ChIP$score, Twi_embryo_ChIP$score))
      
      s2 <- plotSignal(
        data = `Twi_embryo_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 1.75), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ChIP_range
      )
      
      annoYaxis(
        plot = s2, at = round(ChIP_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      s3 <- plotSignal(
        data = `Twi_S2_ChIP`, params = zoomRegion,
        x = (ref_x + x_offset_browser), y = (ref_y + 2.25), width = 3.75, height = gb_height,
        just = c("left", "center"), default.units = "cm",
        linecolor = twi_color, fill = twi_color, baseline = FALSE, baseline.lwd = 0, baseline.color = twi_color,
        range = ATAC_range
      )
      
      annoYaxis(
        plot = s3, at = round(ATAC_range, 1),
        axisLine = TRUE, fontsize = 5, lwd = 0.5, 
      )
      
      
      
      plotGenes(
        params = zoomRegion, assembly = "dm6",
        x = (ref_x + x_offset_browser), y = (ref_y + 2.5), height = 1, width = 3.75, fontsize = small_text_params$fontsize,
        just = c("left", "top"),
        default.units = "cm"
      )
      
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
      
      plotSegments()
      
      # close graphics device ========================================================
      dev.off()
      
      