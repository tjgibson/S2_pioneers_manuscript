# setup -----------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

source(snakemake@params[["r_source"]])

# read peaks -------------------------------------------------------------------
peak_fn <- snakemake@input[[1]]
peaks <- rtracklayer::import(peak_fn)


# extend peak summits ----------------------------------------------------------
extended_summits <- peaks %>% 
  extract_summits(extend_width = as.integer(snakemake@params[["extend_width"]]))

# filter out artifactual regions -----------------------------------------------
# get name of TF
factor <- basename(peak_fn) %>% 
  str_split(pattern = "_")
factor <- factor[[1]][2] %>% 
  str_remove("a") %>% 
  str_to_lower()

# import artifacts file
artifact_fn <- paste0("resources/", factor, "_regions_to_exclude.bed")

artifacts_gr <- rtracklayer::import(artifact_fn)

filtered_peaks <- extended_summits %>% 
  subsetByOverlaps(artifacts_gr, invert = TRUE)

# export filtered peaks to file  -----------------------------------------------
filtered_peaks %>% 
  as.data.frame() %>% 
  select(1:3, 6:7, 5, 8:10) %>% 
  mutate(strand = ".") %>% 
  mutate(peak = -1) %>% 
  write_tsv(snakemake@output[[1]], col_names = FALSE)

rtracklayer::export(filtered_peaks, snakemake@output[[2]])