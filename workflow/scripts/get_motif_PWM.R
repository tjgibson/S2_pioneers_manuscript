# load packages
library(tidyverse)
library(universalmotif)
library(Biostrings)

# get input from snakemake
meme_file <- snakemake@input[["meme_file"]]
motif_of_interest <- snakemake@input[["meme_motif_number"]]

# get motif of interest
meme <- read_meme(meme_file)
pwm <- twi_meme[[motif_of_interest]]@motif

# write output
write_tsv(pwm, snakemake@output[[1]])