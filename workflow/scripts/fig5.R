# setup ========================================================================
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Dm.eg.db))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))
# library(grImport)
suppressPackageStartupMessages(library(grid))
library(RColorBrewer)

source("workflow/scripts/plot_heatmap.R")
source("workflow/scripts/utils.R")

# define input files ===========================================================
# zld_tissue_occupancy_fn <- "results/ChIP_tissue_classes/zld_tissue_classes.tsv"
# grh_tissue_occupancy_fn <- "results/ChIP_tissue_classes/grh_tissue_classes.tsv"
# twi_tissue_occupancy_fn <- "results/ChIP_tissue_classes/twi_tissue_classes.tsv"
# zld_titration_classes_fn <- "results/ChIP_titration_classes/zld_titration_classes.tsv"
# grh_titration_classes_fn <- "results/ChIP_titration_classes/grh_titration_classes.tsv"
# twi_titration_classes_fn <- "results/ChIP_titration_classes/twi_titration_classes.tsv"
# zld_titration_ChIP_rpkm_fn <- "results/ChIP_tissue_classes/zld_classes_titration_ChIP_rpkm.tsv"
# zld_titration_ATAC_rpkm_fn <- "results/ChIP_tissue_classes/zld_classes_titration_ATAC_rpkm.tsv"
# grh_titration_ChIP_rpkm_fn <- "results/ChIP_tissue_classes/grh_classes_titration_ChIP_rpkm.tsv"
# grh_titration_ATAC_rpkm_fn <- "results/ChIP_tissue_classes/grh_classes_titration_ATAC_rpkm.tsv"
# twi_titration_ChIP_rpkm_fn <- "results/ChIP_tissue_classes/twi_classes_titration_ChIP_rpkm.tsv"
# twi_titration_ATAC_rpkm_fn <- "results/ChIP_tissue_classes/twi_classes_titration_ATAC_rpkm.tsv"


zld_titration_ChIP_rpkm_fn <- snakemake@input[["zld_titration_ChIP_rpkm_fn"]]
zld_titration_ATAC_rpkm_fn <- snakemake@input[["zld_titration_ATAC_rpkm_fn"]]

grh_titration_ChIP_rpkm_fn <- snakemake@input[["grh_titration_ChIP_rpkm_fn"]]
grh_titration_ATAC_rpkm_fn <- snakemake@input[["grh_titration_ATAC_rpkm_fn"]]

twi_titration_ChIP_rpkm_fn <- snakemake@input[["twi_titration_ChIP_rpkm_fn"]]
twi_titration_ATAC_rpkm_fn <- snakemake@input[["twi_titration_ATAC_rpkm_fn"]]

zld_tissue_occupancy_fn <- snakemake@input[["zld_tissue_occupancy"]]
grh_tissue_occupancy_fn <- snakemake@input[["grh_tissue_occupancy"]]
twi_tissue_occupancy_fn <- snakemake@input[["twi_tissue_occupancy"]]

zld_titration_classes_fn <- snakemake@input[["zld_titration_classes_fn"]]
grh_titration_classes_fn <- snakemake@input[["grh_titration_classes_fn"]]
twi_titration_classes_fn <- snakemake@input[["twi_titration_classes_fn"]]

## create blank layout for plot =================================================
pdf(snakemake@output[[1]], useDingbats = FALSE)
# pdf("manuscript/figures/fig5.pdf")
pageCreate(width = 18, height = 16, default.units = "cm", showGuides = FALSE)

# general figure settings ======================================================
# text parameters for Nature Genetics
panel_label_params <- pgParams(fontsize = 8)
large_text_params <- pgParams(fontsize = 7)
small_text_params <- pgParams(fontsize = 5)

# colors
zld_color <- "#5BBCD6"
grh_color <- "#F98400"
twi_color <- "#00A08A"


# panel A ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 0.5

# panel label
plotText(
  label = "a", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


a_plot <- zld_titration_ChIP_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` = `S2-Zld-0uM_aZld_IP`,
    `500` = `S2-Zld-500uM_aZld_IP`,
    `1000` = `S2-Zld-1000uM_aZld_IP`,
    `1500` = `S2-Zld-1500uM_aZld_IP`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "500", "1000", "1500"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ChIP-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = a_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)


# panel B ======================================================================
# reference points for positioning figure components
ref_x <- 7.5
ref_y <- 0.5

# panel label
plotText(
  label = "b", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


b_plot <- zld_titration_ATAC_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` = `titration_S2-Zld_0uM`,
    `500` = `titration_S2-Zld_500uM`,
    `1000` = `titration_S2-Zld_1000uM`,
    `1500` = `titration_S2-Zld_1500uM`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "500", "1000", "1500"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = zld_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ATAC-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = b_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)

# panel C ======================================================================
# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 0.5


# generate plot
zld_class_plot <- zld_titration_classes_fn |> 
  read_tsv() |> 
  filter(class != "i") |> 
    ggplot(aes(x = CuSO4, fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm'),
        legend.position = "right") +
  scale_fill_manual(values = brewer.pal(3, "Blues")[2:3]) +
  ylab("n peaks") 

# place chart on page
plotGG(
  plot = zld_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel label
plotText(
  label = "c", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

# panel D ======================================================================
library(eulerr)

# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 2.5

# panel label
plotText(
  label = "d", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "top", default.units = "cm"
)

# define input files
input_files <- c(
  zld_0uM = "ChIPseq/results/peaks/filtered/S2-Zld-0uM_aZld_IP.narrowPeak",
  zld_500uM = "ChIPseq/results/peaks/filtered/S2-Zld-500uM_aZld_IP.narrowPeak",
  zld_1000uM = "ChIPseq/results/peaks/filtered/S2-Zld-1000uM_aZld_IP.narrowPeak",
  zld_1500uM = "ChIPseq/results/peaks/filtered/S2-Zld-1500uM_aZld_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


zld_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = zld_euler_plot,
  x = (ref_x + 0.5), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel E ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 5.5

# panel label
plotText(
  label = "e", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

e_plot <- grh_titration_ChIP_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` = `S2-Grh-0uM_aGrh_IP`,
    `25` = `S2-Grh-25uM_aGrh_IP`,
    `100` = `S2-Grh-100uM_aGrh_IP`,
    `400` = `S2-Grh-400uM_aGrh_IP`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "25", "100", "400"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = grh_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ChIP-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = e_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)


# panel F ======================================================================
# reference points for positioning figure components
ref_x <- 7.5
ref_y <- 5.5

# panel label
plotText(
  label = "f", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


e_plot <- grh_titration_ATAC_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` =  `titration_S2-Grh_0uM`,
    `25` = `titration_S2-Grh_25uM`,
    `100` = `titration_S2-Grh_100uM`,
    `400` = `titration_S2-Grh_400uM`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "25", "100", "400"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = grh_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ATAC-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = e_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)

# panel G ======================================================================
# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 5.5


# generate plot
grh_class_plot <- grh_titration_classes_fn |> 
  read_tsv() |> 
  filter(class != "i") |> 
  ggplot(aes(x = as.factor(CuSO4), fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_manual(values = brewer.pal(3, "Oranges")[2:3]) +
  ylab("n peaks")

# place chart on page
# place chart on page
plotGG(
  plot = grh_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)


# panel label
plotText(
  label = "g", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# panel H ======================================================================

# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 7.5

# panel label
plotText(
  label = "h", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "top", default.units = "cm"
)

# define input files
input_files <- c(
  grh_0uM = "ChIPseq/results/peaks/filtered/S2-Grh-0uM_aGrh_IP.narrowPeak",
  grh_25uM = "ChIPseq/results/peaks/filtered/S2-Grh-25uM_aGrh_IP.narrowPeak",
  grh_100uM = "ChIPseq/results/peaks/filtered/S2-Grh-100uM_aGrh_IP.narrowPeak",
  grh_400uM = "ChIPseq/results/peaks/filtered/S2-Grh-400uM_aGrh_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


grh_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = grh_euler_plot,
  x = (ref_x + 0.5), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

# panel I ======================================================================
# reference points for positioning figure components
ref_x <- 0.5
ref_y <- 10.5

# panel label
plotText(
  label = "i", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)

i_plot <- twi_titration_ChIP_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` = `S2-HA-Twi-0uM_aHA_IP`,
    `10` = `S2-HA-Twi-10uM_aHA_IP`,
    `40` = `S2-HA-Twi-40uM_aHA_IP`,
    `160` = `S2-HA-Twi-160uM_aHA_IP`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "10", "40", "160"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = twi_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ChIP-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = i_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)

# panel J ======================================================================
# reference points for positioning figure components
ref_x <- 7.5
ref_y <- 10.5

# panel label
plotText(
  label = "j", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


j_plot <- twi_titration_ATAC_rpkm_fn |> 
  read_tsv() |> 
  filter(class %in% c("class I", "class II", "class III")) |> 
  dplyr::rename(
    `0` =  `titration_S2-HA-Twi_0uM`,
    `10` = `titration_S2-HA-Twi_10uM`,
    `40` = `titration_S2-HA-Twi_40uM`,
    `160` = `titration_S2-HA-Twi_160uM`
  ) |> 
  pivot_longer(-c("seqnames","start","end","class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("0", "10", "40", "160"))) |>
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = twi_color, outlier.size = 0.01, lwd = 0.1) +
  xlab(bquote("[" ~ CuSO[4] ~ "]")) +
  ylab(bquote("ATAC-seq" ~ log[2](RPKM))) +
  theme_bw(base_size = small_text_params$fontsize) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)


plotGG(
  plot = j_plot,
  x = (ref_x), y = (ref_y),
  width = 6.5, height = 4.5, just = c("left", "top"),
  default.units = "cm"
)

# panel K ======================================================================
# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 10.5


# generate plot
twi_class_plot <- twi_titration_classes_fn |> 
  read_tsv() |> 
  filter(class != "i") |> 
  ggplot(aes(x = as.factor(CuSO4), fill = class)) + 
  geom_bar(color = "black", linewidth = 0.1) +
  theme_classic(base_size = small_text_params$fontsize) +
  theme(legend.key.size = unit(2, 'mm')) +
  scale_fill_manual(values = brewer.pal(3, "GnBu")[2:3]) +
  ylab("n peaks")

# place chart on page
plotGG(
  plot = twi_class_plot,
  x = (ref_x), y = ref_y,
  width = 3, height = 2, just = c("left", "top"),
  default.units = "cm"
)


# panel label
plotText(
  label = "k", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "bottom", default.units = "cm"
)


# panel L ======================================================================

# reference points for positioning figure components
ref_x <- 14.5
ref_y <- 12.5

# panel label
plotText(
  label = "l", params = panel_label_params, fontface = "bold",
  x = ref_x, y = ref_y, just = "top", default.units = "cm"
)

# define input files
input_files <- c(
  twi_0uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-0uM_aHA_IP.narrowPeak",
  twi_10uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-10uM_aHA_IP.narrowPeak",
  twi_40uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-40uM_aHA_IP.narrowPeak",
  twi_160uM = "ChIPseq/results/peaks/filtered/S2-HA-Twi-160uM_aHA_IP.narrowPeak"
  
)

peak_overlaps <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table()


twi_euler_plot <- peak_overlaps %>% 
  dplyr::select(8,9) %>% 
  as.matrix() %>% 
  euler() %>% 
  plot(quantities = list(fontsize = small_text_params$fontsize), labels = list(fontsize = small_text_params$fontsize))

plotGG(
  plot = twi_euler_plot,
  x = (ref_x + 0.5), y = (ref_y + 0.25),
  width = 2, height = 2, just = c("left", "top"),
  default.units = "cm"
)

dev.off()
