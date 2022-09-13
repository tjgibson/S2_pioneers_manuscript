# setup ------------------------------------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(Biostrings))

# get input file names ---------------------------------------------------------
ChIP_classes <- snakemake@input[[1]] %>% 
  read_tsv()

# export each class to bed file ------------------------------------------------
ChIP_classes %>% 
  filter(class == "i") %>% 
  makeGRangesFromDataFrame() %>% 
  rtracklayer::export(snakemake@output[[1]])

ChIP_classes %>% 
  filter(class == "ii") %>% 
  makeGRangesFromDataFrame() %>% 
  rtracklayer::export(snakemake@output[[2]])

ChIP_classes %>% 
  filter(class == "iii") %>% 
  makeGRangesFromDataFrame() %>% 
  rtracklayer::export(snakemake@output[[3]])

ChIP_classes %>% 
  filter(class != "i") %>% 
  makeGRangesFromDataFrame() %>% 
  rtracklayer::export(snakemake@output[[4]])

# get underlying sequences and write to fasta file -----------------------------
classes <- c("i", "ii", "iii")

for (i in seq(classes)) {
peaks <- ChIP_classes %>% 
  filter(class == classes[i]) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(peaks) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
seqinfo(peaks) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
peaks <- trim(peaks)

peak_seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, peaks)
names(peak_seqs) <- peaks$peak_id

writeXStringSet(peak_seqs, snakemake@output[[i + 4]])
}

peaks <- ChIP_classes %>% 
  filter(class != "i") %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqlevels(peaks) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
seqinfo(peaks) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
peaks <- trim(peaks)

peak_seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, peaks)
names(peak_seqs) <- peaks$peak_id

writeXStringSet(peak_seqs, snakemake@output[[8]])