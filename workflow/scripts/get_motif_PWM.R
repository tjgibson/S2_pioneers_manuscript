# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))

# get input from snakemake =====================================================
meme_file <- snakemake@input[["meme_file"]]
motif_of_interest <- snakemake@params[["meme_motif_number"]]

# extract PWM for motif of interest ============================================
meme <- read_meme(meme_file) %>% 
  convert_type("PWM", pseudocount = 0.1)
pwm <- meme[[motif_of_interest]]
pwm@pseudocount <- 0.1
# write output =================================================================
pwm %>% 
  as.data.frame() %>% 
  write.table(snakemake@output[[1]], sep = "\t", quote = FALSE)