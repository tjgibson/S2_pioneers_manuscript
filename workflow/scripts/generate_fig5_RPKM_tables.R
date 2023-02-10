library(tidyverse)
library(Rsubread)
library(GenomicRanges)

source("workflow/scripts/utils.R")

# Zld ChIP =====================================================================
# import unit names
units <- read_tsv("config/ChIP_units.tsv")
units <- units[c(37:84),]

units <- units |> 
  filter(call_peaks) |> 
  filter(grepl("Zld", sample_name))

# import regions of interest
zld_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")

saf <- zld_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(zld_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ChIPseq/results/aligned_reads/filtered/", units$sample_name, ".bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9].bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(zld_tissue_occupancy)

# prepare output table
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  select(-peak_id) |> 
  write_tsv("results/ChIP_tissue_classes/zld_classes_titration_ChIP_rpkm.tsv")

# Zld ATAC =====================================================================
# import unit names
units <- read_tsv("config/ATACseq_units.tsv")
units <- units[25:48,]

units <- units |> 
  filter(grepl("Zld", sample_name))

# import regions of interest
zld_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/zld_tissue_classes.tsv")

saf <- zld_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(zld_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ATACseq/results/aligned_reads/split_fragments/", units$sample_name, "_small.bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE, isPairedEnd = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9]_small.bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(zld_tissue_occupancy)

# prepare output table
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  select(-peak_id) |> 
  write_tsv("results/ChIP_tissue_classes/zld_classes_titration_ATAC_rpkm.tsv")



# Grh ==========================================================================
# import unit names
units <- read_tsv("config/ChIP_units.tsv")
units <- units[37:84,]

units <- units |> 
  filter(call_peaks) |> 
  filter(grepl("Grh", sample_name))

# import regions of interest
grh_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")

saf <- grh_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(grh_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ChIPseq/results/aligned_reads/filtered/", units$sample_name, ".bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9].bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(grh_tissue_occupancy)

# prepare output table
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  select(-peak_id) |> 
  write_tsv("results/ChIP_tissue_classes/grh_classes_titration_ChIP_rpkm.tsv")


# Grh ATAC =====================================================================
# import unit names
units <- read_tsv("config/ATACseq_units.tsv")
units <- units[25:48,]

units <- units |> 
  filter(grepl("Grh", sample_name))

# import regions of interest
grh_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/grh_tissue_classes.tsv")

saf <- grh_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(grh_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ATACseq/results/aligned_reads/split_fragments/", units$sample_name, "_small.bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE, isPairedEnd = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9]_small.bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(grh_tissue_occupancy)

# prepare output table
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  select(-peak_id) |> 
  write_tsv("results/ChIP_tissue_classes/grh_classes_titration_ATAC_rpkm.tsv")

# Twi ==========================================================================
# import unit names
units <- read_tsv("config/ChIP_units.tsv")
units <- units[37:84,]

units <- units |> 
  filter(call_peaks) |> 
  filter(grepl("Twi", sample_name))

# import regions of interest
twi_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")

saf <- twi_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(twi_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ChIPseq/results/aligned_reads/filtered/", units$sample_name, ".bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9].bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(twi_tissue_occupancy)

# prepare output table
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  select(-peak_id) |> 
  write_tsv("results/ChIP_tissue_classes/grh_classes_titration_ATAC_rpkm.tsv")

# Twi ATAC =====================================================================
# import unit names
units <- read_tsv("config/ATACseq_units.tsv")
units <- units[25:48,]

units <- units |> 
  filter(grepl("Twi", sample_name))

# import regions of interest
twi_tissue_occupancy <- read_tsv("results/ChIP_tissue_classes/twi_tissue_classes.tsv")

saf <- twi_tissue_occupancy |> 
  add_column(GeneID = paste0("peak_", seq(nrow(twi_tissue_occupancy)))) |> 
  add_column("strand" = "+") |> 
  select(GeneID, seqnames, start, end, strand) |> 
  dplyr::rename(Chr = seqnames, Start = start, End = end)

# get list of BAM files
bam_files <- paste0("ATACseq/results/aligned_reads/split_fragments/", units$sample_name, "_small.bam")

# get read count_table for ChIP-seq data
read_counts <- featureCounts(bam_files, annot.ext = saf, allowMultiOverlap = TRUE, isPairedEnd = TRUE)

# combine replicate counts
annotation <- read_counts$annotation
read_counts <- read_counts$counts |> 
  as.data.frame() |> 
  rownames_to_column("peak_id") 

colnames(read_counts) <- str_replace(colnames(read_counts), "_rep[0-9]_small.bam", "")

read_counts <- read_counts |> 
  pivot_longer(-c("peak_id"), names_to = "sample_name", values_to = "read_count") |> 
  group_by(peak_id, sample_name) |> 
  summarise(read_count = sum(read_count)) |> 
  pivot_wider(names_from = sample_name, values_from = read_count)

# add length information
peak_widths <- annotation |> 
  select(GeneID, Length) |> 
  dplyr::rename(peak_id = "GeneID", peak_width = "Length")

read_counts <- read_counts |> 
  left_join(peak_widths) |> 
  column_to_rownames("peak_id")

# RPKM normalize counts
rpkm_counts <- rpkm(select(read_counts, -c("peak_width")), widths = read_counts$peak_width)

# add class information
peak_class <- annotation |> 
  select(GeneID,Chr, Start, End) |> 
  dplyr::rename(peak_id = "GeneID", seqnames = "Chr", start = "Start", end = "End") |> 
  left_join(twi_tissue_occupancy) |> 
  select(peak_id, class)

# plot boxplot
rpkm_counts |> 
  rownames_to_column("peak_id") |> 
  left_join(peak_class) |> 
  pivot_longer(-c("peak_id", "class"), names_to = "sample_name", values_to = "RPKM") |> 
  mutate(sample_name = factor(sample_name, levels = c("titration_S2-HA-Twi_0uM","titration_S2-HA-Twi_10uM", "titration_S2-HA-Twi_40uM", "titration_S2-HA-Twi_160uM"))) |> 
  ggplot(aes(x = sample_name, y = log2(RPKM))) +
  geom_boxplot(fill = twi_color) +
  theme_bw(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(~class)

