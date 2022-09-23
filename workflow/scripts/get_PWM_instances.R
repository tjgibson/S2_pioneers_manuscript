# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))

# get input from snakemake =====================================================
pwm_file <- snakemake@input[["PWM"]]


# score motif occurrences in genome ============================================
# read in PWM file
message("reading PWM from file")
pwm <- pwm_file %>% 
  read.table() %>% 
  as.matrix()

# get motif instances
message("getting motif instances")
motif_instances <- matchPWM(pwm, BSgenome.Dmelanogaster.UCSC.dm6, min.score = snakemake@params[["threshold"]], with.score = TRUE)

# filter out duplicates caused by palindromic motifs
message("filtering out duplicate motifs")
n_motifs <- length(motif_instances)
motif_instances <- motif_instances %>% 
  as.data.frame() %>% 
  distinct(seqnames, start, end, .keep_all = TRUE) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

seqinfo(motif_instances) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)

message(paste("filtered out", n_motifs - length(motif_instances), "duplicated motifs"))

# transform motif score into percentage
motif_instances <- motif_instances %>% 
  plyranges::mutate(score = score / maxScore(pwm) * 100)


# export motif occurences ======================================================
message("writing output files")
motif_instances %>% 
  as.data.frame() %>% 
  write_tsv(snakemake@output[["tsv"]])

motif_instances %>% 
  export(snakemake@output[["bed"]])


motif_instances %>% 
  plyranges::filter(strand == "+") %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export(snakemake@output[["bw_plus"]])

motif_instances %>% 
  plyranges::filter(strand == "-") %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export(snakemake@output[["bw_minus"]])

motif_instances %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export(snakemake@output[["bw_all"]])