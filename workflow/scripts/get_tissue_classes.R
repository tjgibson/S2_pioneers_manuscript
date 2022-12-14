# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))

source("workflow/scripts/utils.R")


# define input files ===========================================================
# S2 ChIP peaks and classes
S2_bound_regions <- c(
  all_peaks = snakemake@input[["all_peaks"]],
  class_I_peaks = snakemake@input[["class_I_peaks"]],
  class_II_peaks = snakemake@input[["class_II_peaks"]],
  class_III_peaks = snakemake@input[["class_III_peaks"]]
)

# other tissue ChIP peaks
other_tissue_bound_regions <- c(
  snakemake@input[["other_tissue_peaks"]]
)
names(other_tissue_bound_regions) <- str_replace(basename(other_tissue_bound_regions), ".bed", "")

# Regions with repressive histone modifications
repressed_regions <- c(
  H3K27me3_peaks = snakemake@input[["H3K27me3_peaks"]],
  H3K9me3_peaks = snakemake@input[["H3K9me3_peaks"]]
)

# nonspecific sites to exclude
ns_sites <- c(
  ns_sites = snakemake@input[["nonspecific_sites"]]
)

# combine all inputs
input_files <- c(
  S2_bound_regions,
  other_tissue_bound_regions,
  repressed_regions,
  ns_sites
)


# Define potential binding sites ===============================================
overlap_table <- input_files |>
  map(rtracklayer::import) |> 
  GRangesList()  |> 
  peak_overlap_table()

# annotate peaks ==============================================================
other_tissue_names <- names(other_tissue_bound_regions)

overlap_table <- overlap_table |>
  rowwise() |> 
  mutate(n_tissues_bound = sum(c_across(all_of(c("all_peaks",other_tissue_names))))) |>
  ungroup() |> 
  mutate(class = case_when(
    (!all_peaks & n_tissues_bound >= 1 & H3K27me3_peaks) ~ "repressed_H3K27me3",
    (!all_peaks & n_tissues_bound >= 1 & H3K9me3_peaks) ~ "repressed_H3K9me3",
    (!all_peaks & n_tissues_bound >= 1 & !H3K27me3_peaks & !H3K9me3_peaks) ~ "repressed_other",
    (class_I_peaks) ~ "class I",
    (class_II_peaks) ~ "class II",
    (class_III_peaks) ~ "class III"
  )) |>
  filter(!is.na(class)) |> 
  filter(!ns_sites)

# # # generate random background regions ==================================================
keep_chroms <- readLines(snakemake@input[["keep_chroms"]])
# 
# tile_seq_info <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
# tile_seq_info <- tile_seq_info[seqnames(tile_seq_info)[seqnames(tile_seq_info) %in% keep_chroms]]
# 
# genome_tiles <- tileGenome(tile_seq_info, tilewidth=500,
#                            cut.last.tile.in.chrom=TRUE) |> 
#   subsetByOverlaps(rtracklayer::import(ns_sites), invert = TRUE)
# 
# 
# background_sites <- genome_tiles[sample(seq(genome_tiles), size = sum(overlap_table$all_peaks))] |> 
#   as.data.frame() |> 
#   add_column(class = "background")
# 



## output all annotated sites ==================================================
# get output sites
out_sites <- overlap_table |> 
  # bind_rows(background_sites) |> 
  dplyr::select(seqnames, start, end, class)

# write output sites to file
out_sites  |>  
  filter(seqnames %in% keep_chroms) |> 
  write_tsv(snakemake@output[[1]])
