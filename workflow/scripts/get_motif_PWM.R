# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))

# get input from snakemake =====================================================
meme_file <- snakemake@input[["meme_file"]]
motif_of_interest <- snakemake@params[["meme_motif_number"]]

# extract PWM for motif of interest ============================================
meme <- read_meme(meme_file)
pwm <- meme[[motif_of_interest]]@motif

# write output =================================================================
pwm %>% 
  as.data.frame() %>% 
  write.table(snakemake@output[[1]], sep = "\t", quote = FALSE)