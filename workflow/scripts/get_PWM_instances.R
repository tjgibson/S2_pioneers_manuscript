# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(memes))


# get input from snakemake =====================================================
pwm_file <- snakemake@input[["PWM"]]


# score motif occurrences in genome ============================================
# read in PWM file
message("reading PWM from file")
pwm <- pwm_file |> 
  readRDS()

# get motif instances
message("getting motif instances")
# motif_instances <- matchPWM(pwm, BSgenome.Dmelanogaster.UCSC.dm6, min.score = snakemake@params[["threshold"]], with.score = TRUE)
motif_instances <- 
  get_sequence(GRangesForBSGenome("dm6", chrom = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")), BSgenome.Dmelanogaster.UCSC.dm6) |> 
  runFimo(pwm, thresh = snakemake@params[["threshold"]])


seqinfo(motif_instances) <- seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)

# export motif occurences ======================================================
message("writing output files")
motif_instances |> 
  as.data.frame() |> 
  write_tsv(snakemake@output[["tsv"]])

motif_instances |> 
  export(snakemake@output[["bed"]])


motif_instances |> 
  plyranges::mutate(score = -log10(pvalue)) |>  
  plyranges::reduce_ranges(score = mean(score)) |> 
  export(snakemake@output[["bw_motif_pvalue"]])

motif_instances |> 
  plyranges::reduce_ranges(score = mean(score)) |> 
  export(snakemake@output[["bw__motif_score"]])