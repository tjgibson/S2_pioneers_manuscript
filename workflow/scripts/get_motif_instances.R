# setup ========================================================================
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(rtracklayer))

# get motif instances ==========================================================
# get input motif
motif <- DNAString(snakemake@params[["motif"]])

# check if motif is palindrome
is_palindrome <- ifelse(motif == reverseComplement(motif), TRUE, FALSE)

# check if motif is palindrome 
if (is_palindrome) {
  motif_instances <- GRanges(vmatchPattern(motif, BSgenome.Dmelanogaster.UCSC.dm6, fixed = "subject"))
  motif_instances <- unique(motif_instances)
  names(motif_instances) <- paste0("motif_", 1:length(motif_instances))
} else {
  for_matches <- GRanges(vmatchPattern(motif, BSgenome.Dmelanogaster.UCSC.dm6, fixed = "subject"))
  rev_matches <- GRanges(vmatchPattern(reverseComplement(motif), BSgenome.Dmelanogaster.UCSC.dm6, fixed = "subject"))
  motif_instances <- c(for_matches, rev_matches)
  motif_instances <- unique(motif_instances)
  names(motif_instances) <- paste0("motif_", 1:length(motif_instances))
  
}

# export to bed file ===========================================================
motif_instances %>% 
  rtracklayer::export(snakemake@output[[1]])
