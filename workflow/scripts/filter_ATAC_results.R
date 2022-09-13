# setup -----------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

# read DESeq2 results ----------------------------------------------------------
res_fn <- snakemake@input[[1]]
results <- read_tsv(res_fn)


# read in artifactual regions --------------------------------------------------
# get name of TF
factor <- basename(res_fn) %>% 
  str_split(pattern = "_")
factor <- factor[[1]][1] %>% 
  str_split(pattern = "-")

factor <- factor[[1]][2] %>% 
  str_to_lower()

# import artifacts file
artifact_fn <- paste0("config/", factor, "_artifacts.bed")

artifacts_gr <- rtracklayer::import(artifact_fn)

# filter results table to exclude artifactual regions --------------------------
results_gr <- makeGRangesFromDataFrame(results)

overlaps <- findOverlaps(results_gr, artifacts_gr) %>% 
  queryHits()

filtered_results <- results[-overlaps,] 

# export filtered peaks to file  -----------------------------------------------
filtered_results %>% 
  
  write_tsv(snakemake@output[[1]])
