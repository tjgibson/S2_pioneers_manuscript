# setup ------------------------------------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(Biostrings))

# get input file names ---------------------------------------------------------
tissue_classes <- snakemake@input[[1]] %>% 
  read_tsv()

# export each class to bed file ------------------------------------------------
classes <- 
  c(
    "class I",           
    "class II",       
    "class III",
    "repressed_H3K27me3",
    "repressed_H3K9me3",
    "repressed_other"
    )

for (i in seq(classes)) {
tissue_classes %>% 
  filter(class ==  classes[i]) %>% 
  makeGRangesFromDataFrame() %>% 
  rtracklayer::export(snakemake@output[[i]])
}



# get underlying sequences and write to fasta file -----------------------------
for (i in seq(classes)) {
  peaks <- tissue_classes %>% 
    filter(class == classes[i]) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  seqlevels(peaks) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
  seqinfo(peaks) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
  peaks <- trim(peaks)
  
  peak_seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, peaks)
  names(peak_seqs) <-  paste0("peak_", 1:length(peak_seqs))
  
  writeXStringSet(peak_seqs, snakemake@output[[i + 6]])
}
