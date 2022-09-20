# load packages
library(tidyverse)
library(universalmotif)
library(Biostrings)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# test reading meme output
twi_meme <- read_meme("published_ChIPseq/results/motifs/embryo-1-3H_aTwi/meme.txt")
twi_pwm <- twi_meme[[8]]@motif

grh_meme <- read_meme("published_ChIPseq/results/motifs/embryo-15-16H_aGrh/meme.txt")
grh_pwm <- grh_meme[[3]]@motif

zld_meme <- read_meme("published_ChIPseq/results/motifs/embryo-nc14_aZld/meme.txt")
zld_pwm <- zld_meme[[3]]@motif

# score motif occurences in genome
twi_motif_instances <- matchPWM(twi_pwm, BSgenome.Dmelanogaster.UCSC.dm6, min.score = "80%", with.score = TRUE)
grh_motif_instances <- matchPWM(grh_pwm, BSgenome.Dmelanogaster.UCSC.dm6, min.score = "80%", with.score = TRUE)
zld_motif_instances <- matchPWM(zld_pwm, BSgenome.Dmelanogaster.UCSC.dm6, min.score = "80%", with.score = TRUE)



hist( mcols(twi_motif_instances)$score / maxScore(twi_pwm) )
hist( mcols(grh_motif_instances)$score / maxScore(grh_pwm) )
hist( mcols(zld_motif_instances)$score / maxScore(zld_pwm) )

# export motif occurences
twi_motif_instances %>% 
  as.data.frame() %>% 
  write_tsv("results/motif_instances/twi_motifs_80percent.tsv")

zld_motif_instances %>% 
  as.data.frame() %>% 
  write_tsv("results/motif_instances/zld_motifs_80percent.tsv")

grh_motif_instances %>% 
  as.data.frame() %>% 
  write_tsv("results/motif_instances/grh_motifs_80percent.tsv")


twi_motif_instances %>% 
  plyranges::filter(strand == "+") %>% 
  plyranges::mutate(score = score / maxScore(twi_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export("results/motif_instances/twi_motifs_80percent_plus.bw")

twi_motif_instances %>% 
  plyranges::filter(strand == "-") %>% 
  plyranges::mutate(score = score / maxScore(twi_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export("results/motif_instances/twi_motifs_80percent_minus.bw")

twi_motif_instances %>% 
  plyranges::filter(strand == "-") %>% 
  plyranges::mutate(score = score / maxScore(twi_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>% 
  export("results/motif_instances/twi_motifs_80percent_all.bw")

zld_motif_instances %>% 
  plyranges::filter(strand == "+") %>% 
  plyranges::mutate(score = score / maxScore(zld_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>%
  export("results/motif_instances/zld_motifs_80percent_plus.bw")

zld_motif_instances %>% 
  plyranges::filter(strand == "-") %>% 
  plyranges::mutate(score = score / maxScore(zld_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>%
  export("results/motif_instances/zld_motifs_80percent_minus.bw")



grh_motif_instances %>% 
  plyranges::filter(strand == "+") %>% 
  plyranges::mutate(score = score / maxScore(grh_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>%
  export("results/motif_instances/grh_motifs_80percent_plus.bw")

grh_motif_instances %>% 
  plyranges::filter(strand == "-") %>% 
  plyranges::mutate(score = score / maxScore(grh_pwm) * 100 - 80) %>% 
  plyranges::reduce_ranges(score = mean(score)) %>%
  export("results/motif_instances/grh_motifs_80percent_minus.bw")

