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
  class_III_peaks = snakemake@input[["class_III_peaks"]],
  motif_instances = snakemake@input[["motifs"]]
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
  map(rtracklayer::import) 

overlap_table$motif_instances <- resize(overlap_table$motif_instances, width = 500, fix = "center")

overlap_table <- overlap_table |>
  GRangesList() |>
  peak_overlap_table()

# annotate motifs ==============================================================
other_tissue_names <- names(other_tissue_bound_regions)

overlap_table <- overlap_table |>
  mutate(n_tissues_bound = rowSums(select(., all_peaks, all_of(other_tissue_names)))) |>
  mutate(class = case_when(
    (motif_instances & n_tissues_bound == 0) ~ "motif_never_bound",
    (motif_instances & !all_peaks & n_tissues_bound >= 1 & H3K27me3_peaks) ~ "repressed_H3K27me3",
    (motif_instances & !all_peaks & n_tissues_bound >= 1 & H3K9me3_peaks) ~ "repressed_H3K9me3",
    (motif_instances & !all_peaks & n_tissues_bound >= 1 & !H3K27me3_peaks & !H3K9me3_peaks) ~ "repressed_other",
    (class_I_peaks) ~ "class I",
    (class_II_peaks) ~ "class II",
    (class_III_peaks) ~ "class III"
  )) |>
  filter(!is.na(class)) |> 
  filter(!ns_sites)

# # generate background regions ==================================================
keep_chroms <- readLines(snakemake@input[["keep_chroms"]])

tile_seq_info <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
tile_seq_info <- tile_seq_info[seqnames(tile_seq_info)[seqnames(tile_seq_info) %in% keep_chroms]]

genome_tiles <- tileGenome(tile_seq_info, tilewidth=500,
                           cut.last.tile.in.chrom=TRUE) |> 
  subsetByOverlaps(rtracklayer::import(ns_sites), invert = TRUE)


background_sites <- genome_tiles[sample(seq(genome_tiles), size = sum(overlap_table$all_peaks))] |> 
  as.data.frame() |> 
  add_column(class = "background")




## output all annotated sites ==================================================
# get output sites
out_sites <- overlap_table |> 
  bind_rows(background_sites) |> 
  select(seqnames, start, end, class)

# filter out motifs within 500bp of chromosome end 
# this will prevent errors when trying to get 500bp of sequence around motifs
out_sites_gr <- out_sites |> 
  makeGRangesFromDataFrame() |> 
  resize(width = snakemake@params[["bichrom_window_size"]], fix = "center")

seqlevels(out_sites_gr) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
seqinfo(out_sites_gr) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)

# function to get out of bound ranges (function taken from GenomicRanges source code)
get_out_of_bound_index <- function(x)
{
  if (length(x) == 0L)
    return(integer(0))
  x_seqnames_id <- as.integer(seqnames(x))
  x_seqlengths <- unname(seqlengths(x))
  seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
  seqlength_is_na <- is.na(x_seqlengths)
  seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
  which(seqlevel_has_bounds[x_seqnames_id] &
          (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

OOB_sites <- get_out_of_bound_index(out_sites_gr)

message(paste("filtering out",length(OOB_sites), "sites that are within", snakemake@params[["bichrom_window_size"]], "bp of chromosome ends"))

if (length(OOB_sites) > 0) {
  out_sites <- out_sites[-OOB_sites,]
}


# write output sites to file
out_sites  |>  
  filter(seqnames %in% keep_chroms) |> 
  distinct() |> 
  write_tsv(snakemake@output[[1]])
