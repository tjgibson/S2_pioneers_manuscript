# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(grImport)
library(grid)
library(RColorBrewer)


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

RPKM_table_fn <- "RNAseq/results/count_tables/S2-Grh_RNAseq_RPKM.tsv"

# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig1.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 18.5, default.units = "cm", showGuides = TRUE)

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
# grid.picture(schematic)
# schematic_grob <- pictureGrob(schematic, exp = 0)
schematic_gtree <- grid.grabExpr(grid.picture(schematic))


# plot model on page
plotGG(schematic_gtree, x = 0.1, y = 0.5, width = 12, height = 4.8, default.units = "cm")

# # panel B ======================================================================
# # reference points for positioning figure components
# ref_x <- 7
# ref_y <- 0.5
# 
# # panel label
# plotText(
#   label = "b", params = panel_label_params, fontface = "bold",
#   x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
# )
# 
# # read in RPKM data
# WT_RPKM_table <- RPKM_table_fn |> 
#   read_tsv() |> 
#   dplyr::select(gene_id, gene_symbol, starts_with("S2-WT_noCuSO4")) |> 
#   pivot_longer(3:4, names_to = "sample", values_to = "RPKM")
# 
# # define S2 marker genes 
# s2_markers <- c("FBgn0263772", "FBgn0263772", "FBgn0038294")
# kc_markers <- c("FBgn0259175", "FBgn0261987")
# markers <- c(s2_markers, kc_markers)
# 
# b_plot <- WT_RPKM_table |> 
#   filter(gene_id %in% markers) |>
#   mutate(marker_type = case_when(gene_id %in% s2_markers ~ "S2 markers",
#                                  gene_id %in% kc_markers ~ "Kc markers")) |> 
#   
#   
#   ggplot(aes(x = gene_symbol, y = RPKM)) + 
#   geom_boxplot() +
#   coord_flip() +
#   theme_classic(base_size = small_text_params$fontsize) +
#   theme(axis.title.y = element_blank())
# 
# plotGG(
#   plot = b_plot,
#   x = (ref_x + 0.8), y = (ref_y),
#   width = 3.5, height = 2.5, just = c("left", "top"),
#   default.units = "cm"
# )
# 
# # add additional label to plot
# plotSegments(
#   x0 = (ref_x + 0.75), y0 = (ref_y + 0.1), x1 = (ref_x + 0.75), y1 = (ref_y + 0.9),
#   default.units = "cm",
#   lwd = 1
# )
# 
# plotSegments(
#   x0 = (ref_x + 0.75), y0 = (ref_y + 1.1), x1 = (ref_x + 0.75), y1 = (ref_y + 1.9),
#   default.units = "cm",
#   lwd = 1
# )
# 
# plotText(
#   label = paste0("Kc cell", "\n", "markers"), 
#   x = (ref_x + 0.25), y = (ref_y + 0.5),
#   default.units = "cm", 
#   fontsize = small_text_params$fontsize
# )
# 
# plotText(
#   label = paste0("S2 cell", "\n", "markers"), 
#   x = (ref_x + 0.25), y = (ref_y + 1.5),
#   default.units = "cm", 
#   fontsize = small_text_params$fontsize
# )

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

# placeholder for western blot
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 8, height = 3, default.units = "cm", just = c("top","left"))

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 9.5
ref_y <- 5.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# placeholder for western blot
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 8, height = 3, default.units = "cm", just = c("top","left"))

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 9.5
ref_y <- 5.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# placeholder for western blot
plotRect(x = (ref_x + 0.25), y = (ref_y + 0.25), width = 8, height = 3, default.units = "cm", just = c("top","left"))

# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 9

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# close graphics device ========================================================
dev.off()




