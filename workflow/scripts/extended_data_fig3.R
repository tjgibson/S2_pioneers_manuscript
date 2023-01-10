# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)


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

# zld_ATAC_results_fn <- "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results_annotated.tsv"
# grh_ATAC_results_fn  <- "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results_annotated.tsv"


zld_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results_annotated.tsv"
grh_RNAseq_results_fn  <- "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results_annotated.tsv"

# # create blank layout for plot ===============================================
# pdf(snakemake@output[[1]], useDingbats = FALSE)
pdf("manuscript/figures/extended_data_fig3.pdf", useDingbats = FALSE)
pageCreate(width = 18, height = 12, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
# colors
zld_color <- "#5BBCD6"
grh_color <- "#F98400"
  



# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5


# read in RNAseq_results
zld_RNAseq_results <- read_tsv(zld_RNAseq_results_fn)


# define categories for coloring points
zld_RNAseq_results <- zld_RNAseq_results |> 
  mutate(diff_class = case_when(
    is_diff & has_ChIP_peak ~ "Zld bound",
    is_diff & !has_ChIP_peak ~ "not Zld bound",
    !is_diff ~ "ns"
  ))

# volcano plot
a_plot <- zld_RNAseq_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Zld bound", "Zld bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("grey", "black", zld_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

# place plots on plotGardener page
plotGG(
  plot = a_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

test_genes <- zld_RNAseq_results |> 
  dplyr::filter(is_diff, log2FoldChange > 0) |> 
  pull(gene_id)


ego <- enrichGO(gene          = test_genes,
                # universe      = background_genes,
                OrgDb         = org.Dm.eg.db,
                ont           = "BP",
                keyType = "FLYBASE",
                qvalueCutoff = 0.05,
                readable      = TRUE)


# generate GO enrichment plot
b_plot <- barplot(ego, showCategory=10) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 6, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5


# read in RNAseq_results
grh_RNAseq_results <- read_tsv(grh_RNAseq_results_fn)


# define categories for coloring points
grh_RNAseq_results <- grh_RNAseq_results |> 
  mutate(diff_class = case_when(
    is_diff & has_ChIP_peak ~ "Grh bound",
    is_diff & !has_ChIP_peak ~ "not Grh bound",
    !is_diff ~ "ns"
  ))

# volcano plot
c_plot <- grh_RNAseq_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Grh bound", "Grh bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1) + 
  scale_color_manual(values = c("grey", "black", grh_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

# place plots on plotGardener page
plotGG(
  plot = c_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 4, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

test_genes <- grh_RNAseq_results |> 
  dplyr::filter(is_diff, log2FoldChange > 0) |> 
  pull(gene_id)


ego <- enrichGO(gene          = test_genes,
                # universe      = background_genes,
                OrgDb         = org.Dm.eg.db,
                ont           = "BP",
                keyType = "FLYBASE",
                qvalueCutoff = 0.05,
                readable      = TRUE)


# generate GO enrichment plot
b_plot <- barplot(ego, showCategory=10) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 6, height = 4, just = c("left", "top"),
  default.units = "cm"
)



# close graphics device ========================================================
dev.off()




