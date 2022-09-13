# setup ========================================================================
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# import input files ===========================================================
chip_classes <- snakemake@input[[1]] %>% 
  read_tsv()

nonspecific_sites <- snakemake@input[[2]] %>% 
  rtracklayer::import()

# filter ChIP classes to exclude nonspecific class I sites =====================
chip_classes_gr <- chip_classes %>% 
  makeGRangesFromDataFrame()

overlaps <- findOverlaps(chip_classes_gr, nonspecific_sites)

nonspecific_sites <- chip_classes[overlaps@from,]

chip_classes_filtered <- chip_classes %>% 
  filter(!((peak_id %in% nonspecific_sites$peak_id) & (class == "i")))

# write output data ============================================================
chip_classes_filtered %>% 
  write_tsv(snakemake@output[[1]])

chip_classes_filtered %>% 
  rtracklayer::export(snakemake@output[[2]])
