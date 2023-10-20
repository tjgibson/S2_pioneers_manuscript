# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)
library(RColorBrewer)
library(EBImage)

source("workflow/scripts/plot_heatmap.R")

# define input files ===========================================================
# define input files explicitly for interactve testing
# RPKM_table_fn <- "RNAseq/results/count_tables/S2-Grh_RNAseq_RPKM.tsv"
# 
# zld_blot_image <- "data/immunoblot_raw_images/2018-10-17_Zld_induction/2018-1018-144408_pub.tif"
# grh_blot_image <- "data/immunoblot_raw_images/2018-11-08_Grh_induction/2018-1108-132834_pub.tif"
# 
# S2_Zld_aZld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw"
# S2_WT_aZld_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aZld_IP.bw"
# S2_Zld_IgG_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aIgG_IP.bw"
# Zld_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw"
# Grh_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw"
# Twi_WT_ATAC_bw <- "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw"
# 
# zld_ChIP_classes_fn <- "results/ChIP_peak_classes/zld_ChIP_classes.tsv"
# 
# S2_Grh_aGrh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw"
# S2_WT_aGrh_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aGrh_IP.bw"
# S2_Grh_IgG_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aIgG_IP.bw"
# 
# grh_ChIP_classes_fn <- "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
# 
# S2_Twi_aTwi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw"
# S2_WT_aTwi_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aTwi_IP.bw"
# S2_Twi_IgG_ChIP_bw <- "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aIgG_IP.bw"
# 
# twi_ChIP_classes_fn <- "results/ChIP_peak_classes/twi_ChIP_classes.tsv"

# get input from snakemake
RPKM_table_fn <- snakemake@input[["RPKM_table_fn"]]

zld_blot_image <- snakemake@input[["zld_blot_image"]]
grh_blot_image <- snakemake@input[["grh_blot_image"]]

S2_Zld_aZld_ChIP_bw <- snakemake@input[["S2_Zld_aZld_ChIP_bw"]]
S2_WT_aZld_ChIP_bw <- snakemake@input[["S2_WT_aZld_ChIP_bw"]]
S2_Zld_IgG_ChIP_bw <- snakemake@input[["S2_Zld_IgG_ChIP_bw"]]
Zld_WT_ATAC_bw <- snakemake@input[["Zld_WT_ATAC_bw"]]
Grh_WT_ATAC_bw <- snakemake@input[["Grh_WT_ATAC_bw"]]
Twi_WT_ATAC_bw <- snakemake@input[["Twi_WT_ATAC_bw"]]

zld_ChIP_classes_fn <- snakemake@input[["zld_ChIP_classes_fn"]]

S2_Grh_aGrh_ChIP_bw <- snakemake@input[["S2_Grh_aGrh_ChIP_bw"]]
S2_WT_aGrh_ChIP_bw <- snakemake@input[["S2_WT_aGrh_ChIP_bw"]]
S2_Grh_IgG_ChIP_bw <- snakemake@input[["S2_Grh_IgG_ChIP_bw"]]

grh_ChIP_classes_fn <- snakemake@input[["grh_ChIP_classes_fn"]]

S2_Twi_aTwi_ChIP_bw <- snakemake@input[["S2_Twi_aTwi_ChIP_bw"]]
S2_WT_aTwi_ChIP_bw <- snakemake@input[["S2_WT_aTwi_ChIP_bw"]]
S2_Twi_IgG_ChIP_bw <- snakemake@input[["S2_Twi_IgG_ChIP_bw"]]

twi_ChIP_classes_fn <- snakemake@input[["twi_ChIP_classes_fn"]]


# # create blank layout for plot ===============================================
fig_width <-  18
fig_height <- 18.5

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig1.pdf", useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

pageCreate(width = fig_width, height = fig_height, default.units = "cm", showGuides = FALSE)

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
hm_upstream <-  1000
hm_downstream <-  1000


# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# get model for EPS file
PostScriptTrace(
  "manuscript/figures/models/system_schematic.eps",
  outfilename = "manuscript/figures/models/system_schematic.eps.xml",
  charpath = FALSE,
  charpos = FALSE
)
schematic <- readPicture("manuscript/figures/models/system_schematic.eps.xml")
schematic_gtree <- grid.grabExpr(grid.picture(schematic))


# plot model on page
plotGG(schematic_gtree, x = 0.1, y = 0.5, width = 12, height = 4.8, default.units = "cm")


# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 12
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in sample table 
sample_table <- read_tsv("config/RNAseq_units.tsv") |> 
  dplyr::select(sample_name, sample_group)

# read in RPKM data
zld_RPKM_table <- read_tsv("RNAseq/results/count_tables/S2-Zld_RNAseq_RPKM.tsv") |> 
  pivot_longer(3:10, names_to = "sample_name",values_to =  "RPKM") |> 
  mutate(RPKM = RPKM + 0.01) |> 
  left_join(sample_table, by = "sample_name")

grh_RPKM_table <- read_tsv("RNAseq/results/count_tables/S2-Grh_RNAseq_RPKM.tsv") |> 
  pivot_longer(3:10, names_to = "sample_name",values_to =  "RPKM") |> 
  mutate(RPKM = RPKM + 0.01) |> 
  left_join(sample_table, by = "sample_name")

# generate plots
limits <- c(-7, 14)

# histogram of gene expression levels in WT S2 cells
b_p1 <- zld_RPKM_table |>
  filter(sample_group == "S2-WT_noCuSO4") |>
  ggplot(aes(x = log2(RPKM), y=..density..)) + 
  geom_histogram(bins = 50, color="black", fill="gray", lwd = 0.1) +
  xlim(limits) +
  ylim(c(0,0.15)) +
  scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(axis.title.x = element_blank()) +
  geom_vline(xintercept = log2(1), lty = 2, lwd = 0.1)

# boxplot of Zld rpkm values
b_p2 <- zld_RPKM_table |>
  filter(gene_symbol == "zld") |>
  mutate(condition = fct_relevel(sample_group, c("S2-WT_noCuSO4", "S2-WT_1000uM", "S2-Zld_noCuSO4", "S2-Zld_1000mM"))) |>
  mutate(condition = fct_rev(condition)) |>
  mutate(
    condition = recode(condition,
                          "S2-WT_noCuSO4" ="S2-WT uninduced", 
                          "S2-WT_1000uM" = "S2-WT induced", 
                          "S2-Zld_noCuSO4" = "S2-Zld uninduced", 
                          "S2-Zld_1000mM" = "S2-Zld induced"
    )
  ) |> 
  ggplot(aes(x=condition, y=log2(RPKM)),) + 
  geom_boxplot(fill = "deepskyblue", lwd = 0.05) + 
  coord_flip() +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(axis.title.x = element_blank()) +
  ylim(limits) +
  xlab("Zld") +
  geom_hline(yintercept = log2(1), lty = 2, lwd = 0.1)

# boxplot of Grh rpkm levels
b_p3 <- grh_RPKM_table |>
  filter(gene_symbol == "grh") |>
  mutate(condition = fct_relevel(sample_group, c("S2-WT_noCuSO4", "S2-WT_100uM", "S2-Grh_noCuSO4", "S2-Grh_100mM"))) |>
  mutate(condition = fct_rev(condition)) |>
  mutate(
    condition = recode(condition,
                          "S2-WT_noCuSO4" ="S2-WT uninduced", 
                          "S2-WT_100uM" = "S2-WT induced", 
                          "S2-Grh_noCuSO4" = "S2-Grh uninduced", 
                          "S2-Grh_100mM" = "S2-Grh induced"
    )
  ) |> 
  ggplot(aes(x=condition, y=log2(RPKM)),) + 
  geom_boxplot(fill = "darkorange", lwd = 0.05) + 
  coord_flip() +
  theme_classic(base_size = small_text_params$fontsize) +
  ylim(limits) +
  xlab("Grh") +
  geom_hline(yintercept = log2(1), lty = 2, lwd = 0.1) +
  ylab(bquote(log[2]("RPKM")))

# combine plots into composite plot
b_plot <- rbind(ggplotGrob(b_p1), ggplotGrob(b_p2), ggplotGrob(b_p3), size = "last")

# place plots on plotGardener page
plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 5.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5.5

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot in grayscale
blot_image <- readImage(zld_blot_image) |> 
  channel("gray")

# crop image
blot_image <- blot_image[511:1154,421:580]

# get blot aspect ratio
blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 1.5,
  width = plot_width,
  height = plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# add labels to western blot
plotSegments(
  x0 = ref_x + 1.1, y0 = ref_y + 0.25, x1 = ref_x + 3.25, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 3.5, y0 = ref_y + 0.25, x1 = ref_x + 5.6, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 6, y0 = ref_y + 0.25, x1 = ref_x + 7.75, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotText(
  label = "S2-Zld line A", params = large_text_params, fontface = "bold",
  x = ref_x + 2.15, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "S2-Zld line B", params = large_text_params, fontface = "bold",
  x = ref_x + 4.5, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "2-3H embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 6.9, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 0.75, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "n embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, rot = 45,
  x = ref_x + 1.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "750", params = large_text_params, rot = 45,
  x = ref_x + 1.9, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1000", params = large_text_params, rot = 45,
  x = ref_x + 2.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1250", params = large_text_params, rot = 45,
  x = ref_x + 2.9, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1500", params = large_text_params, rot = 45,
  x = ref_x + 3.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, rot = 45,
  x = ref_x + 3.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "750", params = large_text_params, rot = 45,
  x = ref_x + 4.3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1000", params = large_text_params, rot = 45,
  x = ref_x + 4.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1250", params = large_text_params, rot = 45,
  x = ref_x + 5.3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "1500", params = large_text_params, rot = 45,
  x = ref_x + 5.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "5", params = large_text_params, rot = 45,
  x = ref_x + 6.2, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "5", params = large_text_params, rot = 45,
  x = ref_x + 6.7, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)


plotText(
  label = "10", params = large_text_params, rot = 45,
  x = ref_x + 7.3, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "10", params = large_text_params, rot = 45,
  x = ref_x + 7.8, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

# add arrowheads to blot
plotPolygon(x = c(1, 1.25, 1), y = c(7.2,7.3,7.4), default.units = "cm", fill = "black")

plotPolygon(x = c(1, 1.25, 1), y = c(8.35,8.45,8.55), default.units = "cm", fill = "grey")

# # add band intensity values under blot
# blot_quantification <- "data/immunoblot_raw_images/2018-10-17_Zld_induction/quantification.tsv" |> 
#   read_tsv()
# 
# plotText(
#   label = pull(blot_quantification[1,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 1.25, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[2,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 1.7, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[3,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 2.1, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[4,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 2.6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[5,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 3.1, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[6,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 3.6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[7,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 4.1, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[8,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 4.6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[9,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 5.1, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[10,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 5.5, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# 
# plotText(
#   label = pull(blot_quantification[11,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[12,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 6.6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# 
# plotText(
#   label = pull(blot_quantification[13,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 7.1, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[14,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 7.6, y = ref_y + 3.4, just = c("center"), default.units = "cm"
# )


# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 9.5
ref_y <- 5.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# read in tiff of western blot in grayscale
blot_image <- readImage(grh_blot_image) |> 
  channel("gray")

# crop image
blot_image <- blot_image[486:1101,432:625]

# get blot aspect ratio
blot_dim <- dim(blot_image)
blot_aspect_ratio <- blot_dim[2] / blot_dim[1]



# place blot on page
plot_width <- 7

plotRaster(
  blot_image,
  x = ref_x + 1,
  y = ref_y + 1.5,
  width = plot_width,
  height = plot_width * blot_aspect_ratio,
  default.units = "cm",
  just = c("left, top")
  
)

# panel label
plotSegments(
  x0 = ref_x + 1.1, y0 = ref_y + 0.25, x1 = ref_x + 3.25, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 3.5, y0 = ref_y + 0.25, x1 = ref_x + 5.6, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)

plotSegments(
  x0 = ref_x + 6, y0 = ref_y + 0.25, x1 = ref_x + 7.75, y1 = ref_y + 0.25,
  default.units = "cm",
  lwd = 1
)


plotText(
  label = "S2-Grh line A", params = large_text_params, fontface = "bold",
  x = ref_x + 2.15, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "S2-Grh line B", params = large_text_params, fontface = "bold",
  x = ref_x + 4.5, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "2-3H embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 6.9, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "[CuSO4]", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 0.75, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "n embryos", params = large_text_params, fontface = "bold",
  x = ref_x + 1, y = ref_y + 1.25, just = c("right","center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, rot = 45,
  x = ref_x + 1.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "50", params = large_text_params, rot = 45,
  x = ref_x + 1.9, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "100", params = large_text_params, rot = 45,
  x = ref_x + 2.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "200", params = large_text_params, rot = 45,
  x = ref_x + 2.9, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "400", params = large_text_params, rot = 45,
  x = ref_x + 3.4, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "0", params = large_text_params, rot = 45,
  x = ref_x + 3.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "50", params = large_text_params, rot = 45,
  x = ref_x + 4.3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "100", params = large_text_params, rot = 45,
  x = ref_x + 4.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "200", params = large_text_params, rot = 45,
  x = ref_x + 5.3, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)

plotText(
  label = "400", params = large_text_params, rot = 45,
  x = ref_x + 5.8, y = ref_y + 0.75, just = c("center"), default.units = "cm"
)


plotText(
  label = "5", params = large_text_params, rot = 45,
  x = ref_x + 6.2, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "5", params = large_text_params, rot = 45,
  x = ref_x + 6.7, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)


plotText(
  label = "10", params = large_text_params, rot = 45,
  x = ref_x + 7.3, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

plotText(
  label = "10", params = large_text_params, rot = 45,
  x = ref_x + 7.8, y = ref_y + 1.25, just = c("center"), default.units = "cm"
)

# add arrowheads to blot
plotPolygon(x = c(10, 10.25, 10), y = c(7.3,7.4,7.5), default.units = "cm", fill = "black")

plotPolygon(x = c(10, 10.25, 10), y = c(8.9,9,9.1), default.units = "cm", fill = "grey")

# # add band intensity values under blot
# blot_quantification <- "data/immunoblot_raw_images/2018-11-08_Grh_induction/quantification.tsv" |> 
#   read_tsv()
# 
# plotText(
#   label = pull(blot_quantification[1,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 1.2, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[2,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 1.7, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[3,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 2.1, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[4,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 2.6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[5,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 3.1, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[6,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 3.6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[7,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 4.1, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[8,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 4.6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[9,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 5.1, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[10,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 5.5, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# 
# plotText(
#   label = pull(blot_quantification[11,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[12,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 6.6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# 
# plotText(
#   label = pull(blot_quantification[13,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 7.1, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )
# 
# plotText(
#   label = pull(blot_quantification[14,"Raw volume"]), params = large_text_params, rot = 45,
#   x = ref_x + 7.6, y = ref_y + 3.9, just = c("center"), default.units = "cm"
# )



# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 10

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define bigWig files to use for heatmap
bw <- c(
  S2_WT_ATAC = Zld_WT_ATAC_bw,
  S2_ZLD_ChIP =  S2_Zld_aZld_ChIP_bw,
  S2_WT_aZld_ChIP = S2_WT_aZld_ChIP_bw,
  S2_Zld_IgG_ChIP = S2_Zld_IgG_ChIP_bw
)

# define regions to use for heatmap
regions <- zld_ChIP_classes_fn |>
  read_tsv() |> 
  mutate(class = str_to_upper(class)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# generate metaplots
plot_range <- c(-0.5,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)


metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_3 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Blues") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

# generate heatmap
hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = zld_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  scale_group = c(1,2,2,2),
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

e_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap and metaplots on page
plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.0542), y = (ref_y + 0.25),
  width = (1.1 + 0.1958), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 1.55), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_3,
  x = (ref_x + 2.85), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_4,
  x = (ref_x + 4.15), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

# add legend
plot_colors <- c(
  I = "#DEEBF7",
  II = "#9ECAE1",
  III = "#3182BD"
)


plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = (ref_x + 4.12), y = ref_y + 0.3, width = 0.5, height = 1,
  just = c("center", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

plotGG(
  plot = e_hm,
  x = (ref_x + 0.25), y = (ref_y + 1.4),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_1_heatmap_body_3_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = paste0("S2-WT", "\n", "ATAC"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 0.8), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-Zld", "\n", "anti-Zld"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-WT", "\n", "anti-Zld"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "IgG", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 1.55), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 3.65), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 5.9), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 6.35), just = c("center"), default.units = "cm"
)

# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 6.5
ref_y <- 10

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define bigWig files to use for heatmap
bw <- c(
  S2_WT_ATAC = Grh_WT_ATAC_bw,
  S2_Grh_ChIP =  S2_Grh_aGrh_ChIP_bw,
  S2_WT_aGrh_ChIP = S2_WT_aGrh_ChIP_bw,
  S2_Grh_IgG_ChIP = S2_Grh_IgG_ChIP_bw
)

# define regions for heatmap
regions <- grh_ChIP_classes_fn |>
  read_tsv() |> 
  mutate(class = str_to_upper(class)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# generate metaplots
plot_range <- c(-0.5,5)

metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_3 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "Oranges") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

# generate heatmap
hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = grh_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  scale_group = c(1,2,2,2),
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

f_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap and metaplots on page
plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.0542), y = (ref_y + 0.25),
  width = (1.1 + 0.1958), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 1.55), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_3,
  x = (ref_x + 2.85), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_4,
  x = (ref_x + 4.15), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

# add legend
plot_colors <- c(
  `I` = "#FDBE85",
  `II` = "#FD8D3C",
  `III` = "#D94701"
)


plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = (ref_x + 4.12), y = ref_y + 0.3, width = 0.5, height = 1,
  just = c("center", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

plotGG(
  plot = f_hm,
  x = (ref_x + 0.25), y = (ref_y + 1.4),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_5_heatmap_body_3_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
plotText(
  label = paste0("S2-WT", "\n", "ATAC"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 0.8), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-Grh", "\n", "anti-Grh"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-WT", "\n", "anti-Grh"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "IgG", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 1.55), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 3.15), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 5.5), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 6.25), just = c("center"), default.units = "cm"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 12.5
ref_y <- 10

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# define bigWigs to use for heatmap
bw <- c(
  S2_WT_ATAC = Twi_WT_ATAC_bw,
  S2_Twi_ChIP =  S2_Twi_aTwi_ChIP_bw,
  S2_WT_aTwi_ChIP = S2_WT_aTwi_ChIP_bw,
  S2_Twi_IgG_ChIP = S2_Twi_IgG_ChIP_bw
)

# define regions to use for heatmap
regions <- twi_ChIP_classes_fn |>
  read_tsv() |> 
  mutate(class = str_to_upper(class)) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



# generate metaplots
plot_range <- c(-0.5,8)


metaplot_1 <- plot_average(bw[1], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "GnBu") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)



metaplot_2 <- plot_average(bw[2], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "GnBu") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_3 <- plot_average(bw[3], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "GnBu") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

metaplot_4 <- plot_average(bw[4], regions = regions, row_split = regions$class, line_width = 0.2) +
  scale_color_brewer(palette = "GnBu") +
  theme(text = element_text(size = 5),
        line = element_line(size = 0.1),
        axis.title.y = element_blank(),
        legend.key.size = unit(2, 'mm'),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length.y = unit(0, "cm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  ylim(plot_range)

# generate heatmap
hm <- plot_heatmap_minimal(
  bw, regions, 
  upstream = hm_upstream, downstream = hm_downstream, 
  colors  = twi_heatmap_colors, 
  row_split = regions$class, 
  # order_by_samples = 1, 
  scale_group = c(1,2,2,2),
  use_raster = TRUE, raster_by_magick = TRUE, raster_magick_filter = "Triangle",
  border = "black",
  border_gp = gpar(lwd = 0.5),
  show_heatmap_legend = FALSE,
  return_heatmap_list = TRUE
)

g_hm <- grid.grabExpr(draw(hm, show_heatmap_legend = FALSE, padding = unit(c(0, 0, 0, 0), "mm")))

# place heatmap and metaplots on page
plotGG(
  plot = metaplot_1,
  x = (ref_x + 0.0542), y = (ref_y + 0.25),
  width = (1.1 + 0.1958), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_2,
  x = (ref_x + 1.55), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_3,
  x = (ref_x + 2.85), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

plotGG(
  plot = metaplot_4,
  x = (ref_x + 4.15), y = (ref_y + 0.25),
  width = (1.1), height = 1.1, just = c("left", "top"),
  default.units = "cm"
)

# add legend
plot_colors <- c(
  `I` = "#E0F3DB",
  `II` = "#A8DDB5",
  `III` = "#43A2CA"
)


plotLegend(
  legend = names(plot_colors),
  fill = plot_colors,
  border = FALSE,
  x = (ref_x + 4.12), y = ref_y + 0.3, width = 0.5, height = 1,
  just = c("center", "top"),
  default.units = "cm",
  fontsize = small_text_params$fontsize,
  lty = 1,
  orientation = "v"
)

plotGG(
  plot = g_hm,
  x = (ref_x + 0.25), y = (ref_y + 1.4),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# add axes to heatmaps
seekViewport(name = "matrix_9_heatmap_body_3_1")
grid.xaxis(at = c(0, 1), label = c(paste0("-", hm_upstream / 1000, "KB"), paste0("+",hm_downstream / 1000, "KB")), gp = gpar(lwd = 0.5, fontsize = small_text_params$fontsize))
seekViewport(name = "page")

# heatmap labels
# heatmap labels
plotText(
  label = paste0("S2-WT", "\n", "ATAC"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 0.8), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-Twi", "\n", "anti-Twi"), params = small_text_params, fontface = "bold", 
  x = (ref_x + 2.1), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = paste0("S2-WT", "\n", "anti-Twi"), params = small_text_params, fontface = "bold",
  x = (ref_x + 3.4), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "IgG", params = small_text_params, fontface = "bold",
  x = (ref_x + 4.7), y = (ref_y), just = c("center"), default.units = "cm"
)

plotText(
  label = "class:", params = small_text_params, fontface = "bold",
  x = (ref_x - 0.1), y = (ref_y + 1.55), just = c("center"), default.units = "cm"
)

plotText(
  label = "I", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 3.15), just = c("center"), default.units = "cm"
)

plotText(
  label = "II", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 5.65), just = c("center"), default.units = "cm"
)

plotText(
  label = "III", params = small_text_params, fontface = "bold",
  x = (ref_x), y = (ref_y + 6.35), just = c("center"), default.units = "cm"
)


# close graphics device ========================================================
dev.off()
