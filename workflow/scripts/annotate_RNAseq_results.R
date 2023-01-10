# setup ========================================================================
library(tidyverse)

# define input files ===========================================================
RNAseq_results_fn <- snakemake@input[[1]]
annotated_ChIP_peaks_fn <-  snakemake@input[[2]]


# annotate RNAseq genes based on binding in ChIP-seq data ======================
RNAseq_results <- RNAseq_results_fn |> 
  read_tsv()

annotated_ChIP_peaks <- annotated_ChIP_peaks_fn |> 
  read_tsv()

RNAseq_results$has_ChIP_peak <- RNAseq_results$gene_id %in% annotated_ChIP_peaks$geneId

# export results ===============================================================
RNAseq_results |> 
  write_tsv(snakemake@output[[1]])