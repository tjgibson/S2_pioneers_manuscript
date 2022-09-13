# setup ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(TxDb.Dmelanogaster.UCSC.dm6.ensGene))

# import peaks -----------------------------------------------------------------
peak_fn <- snakemake@input[[1]]
peaks <- rtracklayer::import(peak_fn)

# assign peaks to nearest tss --------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene


promoter_interval <- c(as.integer(snakemake@params[["promoter_upstream"]]), as.integer(snakemake@params[["promoter_downstream"]]))
peakAnno <- annotatePeak(peaks, tssRegion=promoter_interval,
                         TxDb=txdb)

# export table with peak annotations -------------------------------------------
anno_df <- as.data.frame(peakAnno@anno)
write.table(anno_df, snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)
