# setup ========================================================================
library(plotgardener)
library(tidyverse)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(rtracklayer)
library(ChIPseeker)
library(ggpubr)
library(rstatix)

# define input files ===========================================================
# define input files explicitly for testing
# zld_RNAseq_results_fn <- "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results_annotated.tsv"
# grh_RNAseq_results_fn  <- "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results_annotated.tsv"
# 
# zld_ATAC_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Zld_ATACseq_S2-Zld-FL-vs-S2-WT_results.tsv"
# grh_ATAC_results_fn <- "ATACseq/results/DEseq2_results_filtered/S2-Grh_ATACseq_S2-Grh-FL-vs-S2-WT_results.tsv"
# 
# zld_ChIP_peaks_fn <- "ChIPseq/results/peaks/final/S2-Zld_aZld_IP.narrowPeak"
# grh_ChIP_peaks_fn <- "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP.narrowPeak"

# get input files from snakemake
zld_RNAseq_results_fn <- snakemake@input[["zld_RNAseq_results_fn"]]
grh_RNAseq_results_fn  <- snakemake@input[["grh_RNAseq_results_fn"]]

zld_ATAC_results_fn <- snakemake@input[["zld_ATAC_results_fn"]]
grh_ATAC_results_fn <- snakemake@input[["grh_ATAC_results_fn"]]

zld_ChIP_peaks_fn <- snakemake@input[["zld_ChIP_peaks_fn"]]
grh_ChIP_peaks_fn <- snakemake@input[["grh_ChIP_peaks_fn"]]


# create blank layout for plot ===============================================
fig_width <-  16
fig_height <- 18

# open pdf
pdf(snakemake@output[[1]], useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)
# pdf("manuscript/figures/extended_data_fig2.pdf", useDingbats = FALSE, width = fig_width / 2.54,height = fig_height / 2.54)

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



# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5


# read in ATAC-seq results
zld_ATAC_results <- read_tsv(zld_ATAC_results_fn)

# annotate ATAC peaks as overlapping a ChIP-seq peak
zld_ChIP_peaks <- import(zld_ChIP_peaks_fn)

zld_ATAC_gr <- zld_ATAC_results |> 
  makeGRangesFromDataFrame()

zld_ATAC_results$has_ChIP_peak <- FALSE
overlaps <- findOverlaps(zld_ATAC_gr, zld_ChIP_peaks)@from
zld_ATAC_results$has_ChIP_peak[overlaps] <- TRUE

# define categories for coloring points
zld_ATAC_results <- zld_ATAC_results |> 
  mutate(diff_class = case_when(
    is_diff & has_ChIP_peak ~ "Zld bound",
    is_diff & !has_ChIP_peak ~ "not Zld bound",
    !is_diff ~ "ns"
  ))

# volcano plot
a_plot <- zld_ATAC_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Zld bound", "Zld bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1, alpha = 0.5, shape = 16) + 
  scale_color_manual(values = c("grey", "black", zld_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size=5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  guides(colour = guide_legend(override.aes = list(size=0.6))) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))

# place plots on plotGardener page
plotGG(
  plot = a_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "ATAC-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
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
b_plot <- zld_RNAseq_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Zld bound", "Zld bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1, alpha = 0.5, shape = 16) + 
  scale_color_manual(values = c("grey", "black", zld_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size=5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  guides(colour = guide_legend(override.aes = list(size=0.6))) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))

# place plots on plotGardener page
plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "RNA-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 10.5
ref_y <- 0.5


# assign peaks to nearest tss
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
peaks <- makeGRangesFromDataFrame(zld_ATAC_results)

promoter_interval <- c(-500, 100)
peakAnno <- annotatePeak(peaks, tssRegion=promoter_interval,
                         TxDb=txdb)

anno_df <- as.data.frame(peakAnno@anno) |>
  dplyr::select(seqnames, start, end, geneId, distanceToTSS) |>
  dplyr::rename(
    peak_chrom = seqnames,
    peak_start = start,
    peak_end = end,
    gene_id = geneId
  )

p_data <- zld_ATAC_results |> 
  left_join(anno_df, by = c("peak_chrom","peak_start", "peak_end")) |>
  mutate(ATAC_direction_of_change = case_when(
    is_diff & (log2FoldChange > 0) ~ "increased",
    is_diff & (log2FoldChange < 0) ~ "decreased",
    !is_diff ~ "ns"
  )) |>
  mutate(ATAC_direction_of_change = fct_relevel(ATAC_direction_of_change, c("decreased", "ns", "increased"))) |>
  dplyr::select(ATAC_direction_of_change, gene_id, distanceToTSS) |>
  left_join(dplyr::select(zld_RNAseq_results, gene_id, log2FoldChange), by = "gene_id")

# perform statistical comparison with wilcoxon test
stat_test <- compare_means(
  log2FoldChange ~ ATAC_direction_of_change,
  ref.group = "ns",
  p.adjust.method = "bonferroni",
  method='wilcox.test',
  data = p_data
) |>
  mutate(y.position = c(10,11))

# generate plot
c_plot <- p_data |>
  ggplot(aes(x = ATAC_direction_of_change, y = log2FoldChange)) +
  geom_violin(fill = zld_color, lwd = 0.1) + 
  geom_boxplot(width = 0.05,color = "gray50", outlier.shape = NA, outlier.size = 0.01, lwd = 0.1) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  xlab("change in ATAC") +
  ylab(bquote("RNA-seq "~log[2]("fold change"))) +
  theme_classic(base_size = small_text_params$fontsize) +
  stat_pvalue_manual(stat_test, label = "p.adj", size = small_text_params$fontsize * 0.35)

plotGG(
  plot = c_plot,
  x = (ref_x), y = (ref_y),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)



# panel D ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 6


# read in ATAC-seq results
grh_ATAC_results <- read_tsv(grh_ATAC_results_fn)

# annotate ATAC peaks as overlapping a ChIP-seq peak
grh_ChIP_peaks <- import(grh_ChIP_peaks_fn)

grh_ATAC_gr <- grh_ATAC_results |> 
  makeGRangesFromDataFrame()

grh_ATAC_results$has_ChIP_peak <- FALSE
overlaps <- findOverlaps(grh_ATAC_gr, grh_ChIP_peaks)@from
grh_ATAC_results$has_ChIP_peak[overlaps] <- TRUE

# define categories for coloring points
grh_ATAC_results <- grh_ATAC_results |> 
  mutate(diff_class = case_when(
    is_diff & has_ChIP_peak ~ "Grh bound",
    is_diff & !has_ChIP_peak ~ "not Grh bound",
    !is_diff ~ "ns"
  ))

# volcano plot
d_plot <- grh_ATAC_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Grh bound", "Grh bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1, alpha = 0.5, shape = 16) + 
  scale_color_manual(values = c("grey", "black", grh_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size=5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  guides(colour = guide_legend(override.aes = list(size=0.6))) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))

# place plots on plotGardener page
plotGG(
  plot = d_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "ATAC-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 5.5
ref_y <- 6


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
e_plot <- grh_RNAseq_results |>
  mutate(diff_class = fct_relevel(diff_class, c("ns", "not Grh bound", "Grh bound"))) |>
  # mutate(diff_class = fct_rev(diff_class)) |>
  ggplot(aes(x=log2FoldChange, y=-log10(padj), color = diff_class)) + 
  geom_point(size = 0.1, alpha = 0.5, shape = 16) + 
  scale_color_manual(values = c("grey", "black", grh_color)) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size=5),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0)) +
  guides(colour = guide_legend(override.aes = list(size=0.6))) +
  xlab(bquote(log[2]("fold change"))) +
  ylab(bquote(-log[10]("adj. p-value")))


# place plots on plotGardener page
plotGG(
  plot = e_plot,
  x = (ref_x), y = (ref_y),
  width = 4.5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

plotText(
  label = "RNA-seq", params = large_text_params, fontface = "bold",
  x = ref_x + 2.5, y = ref_y, just = "bottom", default.units = "cm"
)

# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 10.5
ref_y <- 6


# assign peaks to nearest tss
peaks <- makeGRangesFromDataFrame(grh_ATAC_results)

promoter_interval <- c(-500, 100)
peakAnno <- annotatePeak(peaks, tssRegion=promoter_interval,
                         TxDb=txdb)

anno_df <- as.data.frame(peakAnno@anno) |>
  dplyr::select(seqnames, start, end, geneId, distanceToTSS) |>
  dplyr::rename(
    peak_chrom = seqnames,
    peak_start = start,
    peak_end = end,
    gene_id = geneId
  )

p_data <- grh_ATAC_results |> 
  left_join(anno_df, by = c("peak_chrom","peak_start", "peak_end")) |>
  mutate(ATAC_direction_of_change = case_when(
    is_diff & (log2FoldChange > 0) ~ "increased",
    is_diff & (log2FoldChange < 0) ~ "decreased",
    !is_diff ~ "ns"
  )) |>
  mutate(ATAC_direction_of_change = fct_relevel(ATAC_direction_of_change, c("decreased", "ns", "increased"))) |>
  dplyr::select(ATAC_direction_of_change, gene_id, distanceToTSS) |>
  left_join(dplyr::select(grh_RNAseq_results, gene_id, log2FoldChange), by = "gene_id")

# perform statistical comparison with wilcoxon test
stat_test <- compare_means(
  log2FoldChange ~ ATAC_direction_of_change,
  ref.group = "ns",
  p.adjust.method = "bonferroni",
  method='wilcox.test',
  data = p_data
) |>
  mutate(y.position = c(13,14))


f_plot <- p_data |>
  ggplot(aes(x = ATAC_direction_of_change, y = log2FoldChange)) +
  geom_violin(fill = grh_color, lwd = 0.1) + 
  geom_boxplot(width = 0.05,color = "gray50", outlier.shape = NA, outlier.size = 0.01, lwd = 0.1) +
  geom_hline(yintercept = 0, color = "gray", lty = 2) +
  xlab("change in ATAC") +
  ylab(bquote("RNA-seq "~log[2]("fold change"))) +
  theme_classic(base_size = small_text_params$fontsize) +
  stat_pvalue_manual(stat_test, label = "p.adj", size = small_text_params$fontsize * 0.35)

plotGG(
  plot = f_plot,
  x = (ref_x), y = (ref_y),
  width = 5, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 11.5

# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
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
g_plot <- barplot(ego, showCategory=10) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = g_plot,
  x = (ref_x), y = (ref_y),
  width = 7, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# panel H ======================================================================
# reference points for positioning figure components
ref_x <- 8
ref_y <- 11.5

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
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
h_plot <- barplot(ego, showCategory=10) +
  theme_minimal(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'))

plotGG(
  plot = h_plot,
  x = (ref_x), y = (ref_y),
  width = 7, height = 5, just = c("left", "top"),
  default.units = "cm"
)

# close graphics device ========================================================
dev.off()