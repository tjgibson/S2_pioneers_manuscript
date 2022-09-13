# setup ------------------------------------------------------------------------
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))

source(snakemake@params[["r_source"]])

# get input file names ---------------------------------------------------------
ChIP_peaks <- snakemake@input[["ChIP_peaks"]]
wt_ATAC_peaks <- snakemake@input[["WT_ATAC_peaks"]]
ATAC_results <- snakemake@input[["ATAC_results"]]
gtf <- snakemake@input[["gtf_genome_annotation"]]
RNAseq_results <- snakemake@input[["RNAseq_results"]]
ChIP_features <- snakemake@input[["ChIP_feature_annotation"]]

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
ChIP_classes <- GRangesList(
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

# annotate ChIP peaks based on ChIP classes
peaks_annotated <- peaks_df %>% 
  left_join(ChIP_classes)

# add additional information to peak table -------------------------------------
# add diff ATAC data
overlaps <- findOverlaps(peaks_gr, diff_atac_gr, minoverlap = 1L)
peaks_annotated$atac_is_diff <- as.logical(NA)
peaks_annotated$atac_padj <- as.logical(NA)
peaks_annotated$atac_FC <- as.logical(NA)

peaks_annotated$atac_is_diff[overlaps@from] <- diff_atac_gr$is_diff[overlaps@to] 
peaks_annotated$atac_padj[overlaps@from] <- diff_atac_gr$padj[overlaps@to] 
peaks_annotated$atac_FC[overlaps@from] <- diff_atac_gr$log2FoldChange[overlaps@to] 

# annotate with nearest gene
gtf <- rtracklayer::import(gtf) %>% 
  mcols() %>%
  as.data.frame() %>%
  filter(type == "gene") %>%
  dplyr::select(gene_id, gene_symbol)

peak_feature_annotation <- read_tsv(ChIP_features) %>% 
  mutate(feature = case_when(annotation == "Promoter" ~ "promoter", annotation != "Promoter" ~ "distal")) %>% 
  rename(peak_id = name, gene_id = geneId) %>% 
  dplyr::select(peak_id, feature, gene_id) %>% 
  left_join(gtf, by = "gene_id")

peaks_annotated <- peaks_annotated %>%
  left_join(peak_feature_annotation, by = "peak_id")

# read in DE genes from RNAseq
DE_genes <- read_tsv(RNAseq_results) %>% 
  dplyr::select(gene_id, is_diff, log2FoldChange, padj) %>% 
  rename(RNAseq_is_diff = is_diff, RNAseq_FC = log2FoldChange, RNAseq_padj = padj)

# add diff expression data to peaks
peaks_annotated <- peaks_annotated %>% 
  left_join(DE_genes, by = "gene_id")

# export table with annotated ChIP classes -------------------------------------
peaks_annotated %>% 
  write_tsv(snakemake@output[[1]])