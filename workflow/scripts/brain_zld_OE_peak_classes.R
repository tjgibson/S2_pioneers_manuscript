# setup ========================================================================
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(here)
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(Biostrings))

source(here("../useful_functions/peak_processing.R"))

# define input files ===========================================================
# define input files explicitly for testing script
brat_aZld_peaks_fn <-  "TF_CUTandRUN/results/peaks/filtered/brain_brat_aZld.narrowPeak"
OE_aZld_peaks_fn <-  "TF_CUTandRUN/results/peaks/filtered/brain_UAS-HA-Zld_aZld.narrowPeak"

# get input files from snakemake
snakemake@input[[1]]
snakemake@input[[2]]

# define peak classes ==========================================================
NB_peaks <- c(
  brat_aZld = brat_aZld_peaks_fn,
  OE_aZld = OE_aZld_peaks_fn) |> 
  map(rtracklayer::import) |>
  GRangesList() |>
  peak_overlap_table() |> 
  mutate(
    class = case_when(
      brat_aZld ~ "endogenous",
      !brat_aZld ~ "novel"
    )
  )

# export novel binding sites ===================================================
NB_peaks |> 
  filter(class == "endogenous") |> 
  makeGRangesFromDataFrame() |> 
  rtracklayer::export(snakemake@output[[1]])

NB_peaks |> 
  filter(class == "novel") |> 
  makeGRangesFromDataFrame() |> 
  rtracklayer::export(snakemake@output[[2]])

# export sequences for novel binding sites =====================================
classes <- c("endogenous", "novel")

for (i in seq(classes)) {
  peaks <- NB_peaks %>% 
    filter(class == classes[i]) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  seqlevels(peaks) <- seqlevels(BSgenome.Dmelanogaster.UCSC.dm6)
  seqinfo(peaks) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)
  peaks <- trim(peaks)
  
  peak_seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, peaks)
  names(peak_seqs) <- paste0("peak_", 1:length(peaks))
  
  out_name <- snakemake@output[[(2 + i)]]
  writeXStringSet(peak_seqs, out_name)
}

# export table of binding sites ================================================
NB_peaks |> 
  select(seqnames, start, end, class) |> 
  dplyr::rename(peak_chr = seqnames, peak_start = start, peak_end = end, class = class) |> 
  write_tsv(snakemake@output[[5]])
