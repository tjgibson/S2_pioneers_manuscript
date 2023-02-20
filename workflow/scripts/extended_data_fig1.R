# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)
library(RColorBrewer)
library(EBImage)


# define input files ===========================================================
# for testing script
# RPKM_table_fn <- "RNAseq/results/count_tables/S2-Grh_RNAseq_RPKM.tsv"
# 
# zld_blot_image <- "data/immunoblot_raw_images/2018-10-17_Zld_induction/2018-1018-144408_pub.tif"
# grh_blot_image <- "data/immunoblot_raw_images/2018-11-08_Grh_induction/2018-1108-132834_pub.tif"

# get input from snakemake
RPKM_table_fn <- snakemake@input[["RPKM_table_fn"]]

zld_blot_image <- snakemake@input[["zld_blot_image"]]
grh_blot_image <- snakemake@input[["grh_blot_image"]]


# # create blank layout for plot ===============================================
pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/extended_data_fig1.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 18.5, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors




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
  ggplot(aes(x=condition, y=log2(RPKM)),) + 
  geom_boxplot(fill = "darkorange", lwd = 0.05) + 
  coord_flip() +
  theme_classic(base_size = small_text_params$fontsize) +
  ylim(limits) +
  xlab("Grh") +
  geom_hline(yintercept = log2(1), lty = 2, lwd = 0.1)

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
  label = "S2 Zld line A", params = large_text_params, fontface = "bold",
  x = ref_x + 2.15, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "S2 Zld line B", params = large_text_params, fontface = "bold",
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
  label = "S2 Grh line A", params = large_text_params, fontface = "bold",
  x = ref_x + 2.15, y = ref_y, just = "top", default.units = "cm"
)

plotText(
  label = "S2 Grh line B", params = large_text_params, fontface = "bold",
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

# close graphics device ========================================================
dev.off()




