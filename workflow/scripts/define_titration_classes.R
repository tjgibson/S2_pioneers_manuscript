# setup ------------------------------------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))


source("workflow/scripts/utils.R")

## Zld 500 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Zld-500uM_aZld_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-500uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
zld_500_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)



## Zld 1000 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Zld-1000uM_aZld_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-1000uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
zld_1000_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)


## Zld 1500 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Zld-1500uM_aZld_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-1500uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
zld_1500_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)

## combine results for Zld =====================================================
zld_500_ChIP_classes$CuSO4 <- 500
zld_1000_ChIP_classes$CuSO4 <- 1000
zld_1500_ChIP_classes$CuSO4 <- 1500


zld_comb_results <- bind_rows(zld_500_ChIP_classes, zld_1000_ChIP_classes, zld_1500_ChIP_classes)

# # plot of class distribution for each experiment
# zld_comb_results |> 
#   ggplot(aes(x = CuSO4, fill = class)) + geom_bar() +
#   scale_fill_brewer(palette = "Blues")

## Grh 25 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Grh-25uM_aGrh_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-25uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
grh_25_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)

## Grh 100 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Grh-100uM_aGrh_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-100uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
grh_100_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)


## Grh 400 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-Grh-400uM_aGrh_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-400uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
grh_400_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)

## combine results for Grh =====================================================
grh_25_ChIP_classes$CuSO4 <- 25
grh_100_ChIP_classes$CuSO4 <- 100
grh_400_ChIP_classes$CuSO4 <- 400


grh_comb_results <- bind_rows(grh_25_ChIP_classes, grh_100_ChIP_classes, grh_400_ChIP_classes)

# # plot of class distribution for each experiment
# grh_comb_results |> 
#   ggplot(aes(x = as.factor(CuSO4), fill = class)) + geom_bar() +
#   scale_fill_brewer(palette = "Oranges")


## Twi 10 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-HA-Twi-10uM_aHA_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-10uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
twi_10_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)

## Twi 40 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-HA-Twi-40uM_aHA_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-40uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
twi_40_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)


## Twi 10 uM ==================================================================
# get input file names ---------------------------------------------------------
ChIP_peaks <- "ChIPseq/results/peaks/filtered/S2-HA-Twi-160uM_aHA_IP.narrowPeak"
wt_ATAC_peaks <- "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed"
ATAC_results <- "ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-160uM-vs-0uM_results.tsv"

# annotate ChIP classes --------------------------------------------------------
# import ChIP peaks
peaks_gr <- rtracklayer::import(ChIP_peaks)

peaks_df <- peaks_gr %>% 
  as.data.frame() %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end, peak_id = name, MACS2_enrichment = signalValue, MACS2_pvalue = pValue, MACS2_qValue = qValue) %>% 
  select(-c(width, strand,score,peak))

# read in WT ATAC peaks
wt_atac <- rtracklayer::import(wt_ATAC_peaks) %>% 
  resize(width = 201, fix = "center")

# read in diff ATAC results
diff_atac_df <- read_tsv(ATAC_results)
diff_atac_gr <- makeGRangesFromDataFrame(diff_atac_df, keep.extra.columns = TRUE)

atac_increased <- diff_atac_df %>% 
  filter(is_diff, log2FoldChange > 0) %>% 
  makeGRangesFromDataFrame()

# generate overlap table and annotate ChIP classes
twi_160_ChIP_classes <- GRangesList(
  ChIP_peak = peaks_gr,
  wt_atac = wt_atac,
  atac_increased = atac_increased
) %>% 
  
  peak_overlap_table() %>% 
  filter(ChIP_peak) %>% 
  mutate(class = case_when(
    ChIP_peak & wt_atac ~ "i",
    ChIP_peak & !wt_atac & !atac_increased ~ "ii",
    ChIP_peak & !wt_atac & atac_increased ~ "iii"
  )) %>% 
  dplyr::rename(peak_chrom = seqnames, peak_start = start, peak_end = end) %>% 
  select(peak_chrom, peak_start, peak_end, class)

## combine results for Twi =====================================================
twi_10_ChIP_classes$CuSO4 <- 10
twi_40_ChIP_classes$CuSO4 <- 40
twi_160_ChIP_classes$CuSO4 <- 160


twi_comb_results <- bind_rows(twi_10_ChIP_classes, twi_40_ChIP_classes, twi_160_ChIP_classes)

# # plot of class distribution for each experiment
# twi_comb_results |>
#   ggplot(aes(x = as.factor(CuSO4), fill = class)) + geom_bar() +
#   scale_fill_brewer(palette = "GnBu")
# 
# 

# export output files ==========================================================
zld_comb_results |> 
  write_tsv("results/ChIP_titration_classes/zld_titration_classes.tsv")

grh_comb_results |> 
  write_tsv("results/ChIP_titration_classes/grh_titration_classes.tsv")

twi_comb_results |> 
  write_tsv("results/ChIP_titration_classes/twi_titration_classes.tsv")
