# setup ========================================================================
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(rtracklayer))

source(snakemake@params[["r_source"]])

# get peak sequences ===========================================================
# get input peaks
peak_file <- snakemake@input[[1]]
sites_for_seqs <- rtracklayer::import(peak_file)

# extend peak summits ----------------------------------------------------------
peaks <- sites_for_seqs %>% 
  extract_summits(extend_width = 201L)

# get underlying peak sequences ------------------------------------------------
seqlevels(peaks) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
seqinfo(peaks) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
peaks <- trim(peaks)

peak_seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, peaks)
names(peak_seqs) <- peaks$name



# export to bed file ===========================================================
writeXStringSet(peak_seqs, snakemake@output[[1]])

