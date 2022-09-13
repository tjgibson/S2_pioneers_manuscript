# setup ========================================================================
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(Biostrings))

source(snakemake@params[["r_source"]])


# define input files ===========================================================
input_files <- snakemake@input %>% 
  unlist()

# Define nonspecific binding sites =============================================
nonspecific_sites <- input_files %>%
  map(rtracklayer::import) %>%
  GRangesList() %>%
  peak_overlap_table() %>% 
  mutate(n_samples = rowSums(select(.,6:8))) %>%
  filter(n_samples == 3) %>% 
  makeGRangesFromDataFrame()

# write nonspecific sites to file ==============================================
nonspecific_sites$peak_id <- paste0("nonspecific_classI_region_", 1:length(nonspecific_sites))
seqlevels(nonspecific_sites) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
seqinfo(nonspecific_sites) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
nonspecific_sites <- trim(nonspecific_sites)

seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, nonspecific_sites)
names(seqs) <- nonspecific_sites$peak_id

# write peaks and sequences to files
writeXStringSet(seqs, snakemake@output[[1]])

nonspecific_sites$name <- nonspecific_sites$peak_id
rtracklayer::export(nonspecific_sites, snakemake@output[[2]])
