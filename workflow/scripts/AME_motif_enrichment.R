# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
library(pheatmap)
library(GenomicRanges)
library(memes)

# define input files ===========================================================
zld_tissue_classes_fn <- "results/ChIP_tissue_classes/zld_tissue_classes.tsv"
grh_tissue_classes_fn <- "results/ChIP_tissue_classes/grh_tissue_classes.tsv"
twi_tissue_classes_fn <- "results/ChIP_tissue_classes/twi_tissue_classes.tsv"

embryo_nc14_RNAseq_fn <- "RNAseq/results/count_tables/embryo_RNAseq_RPKM.tsv"
NB_RNAseq_fn <- "RNAseq/results/count_tables/neuroblast_RNA-seq_RPKM.tsv"

`embryo_15-16H_RNAseq_fn` <- "RNAseq/results/count_tables/embryo-15-16H_RNAseq_RPKM.tsv"
wing_disc_RNAseq_fn <- "RNAseq/results/count_tables/wing-disc_RNAseq_RPKM.tsv"

flybase_synonym_table_fn <- "resources/fb_synonyms.tsv.gz"


# import and clean motif database ==============================================
# import Drosophila motif database
meme_db_path <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE)
# meme_db_path <- "/Users/tylergibson/meme/motif_databases/FLY/combined_drosophila_database.meme"

meme_db <- read_meme(meme_db_path) |> 
  to_df() 

# import gene synonym table
fb_synonym_table <- read_tsv(flybase_synonym_table_fn, comment = "#", 
                             col_names = c("primary_FBid",  
                                           "organism_abbreviation",
                                           "current_symbol",
                                           "current_fullname",
                                           "fullname_synonym",
                                           "symbol_synonym")) |> 
  filter(startsWith(primary_FBid, "FBgn")) |> 
  filter(organism_abbreviation == "Dmel") |> 
  separate_longer_delim(symbol_synonym, delim = "|") |> 
  select(primary_FBid, current_symbol, symbol_synonym) |> 
  dplyr::rename(gene_id = primary_FBid)

fb_synonym_table <- fb_synonym_table |> 
  mutate(gene_symbol = current_symbol)

# get list of current symbols
current_gene_symbols <- fb_synonym_table |> 
  select(gene_id, current_symbol, gene_symbol) |> 
  distinct()

# add gene_id to motif table
meme_db <- meme_db |> 
  left_join(current_gene_symbols, by = join_by(altname == current_symbol))

# fill in gene_ids for symbols that don't match the current gene symbol
meme_db_w_gene_id <- meme_db |> 
  filter(!is.na(gene_id))

meme_db_wo_gene_id <- meme_db |> 
  filter(is.na(gene_id)) |> 
  select(-gene_id, -gene_symbol) |> 
  left_join(fb_synonym_table, by = join_by(altname == symbol_synonym)) |> 
  select(-current_symbol) |> 
  distinct(name, .keep_all = TRUE)

meme_db <- bind_rows(meme_db_w_gene_id, meme_db_wo_gene_id)

# ==============================================================================
# =================== prepare Zld motif and RNA-seq data =======================
# ==============================================================================
class_names <- c(
  "class_I",
  "class_II",
  "class_III",
  "class_IV",
  "class_V",
  "class_VI"
)

tissue_classes <- zld_tissue_classes_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


embryo_RNAseq <- embryo_nc14_RNAseq_fn |> 
  read_tsv()

NB_RNAseq <- NB_RNAseq_fn |> 
  read_tsv()

RNAseq_RPKM <- NB_RNAseq |> 
  select(-c(2:4)) |> 
  left_join(embryo_RNAseq, by = "gene_id") |> 
  select(gene_id, gene_symbol, everything())

# get sequences underlying peaks ===============================================
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

# Get sequences in a 100bp window around the peak summit
tissue_classes_100bp <- tissue_classes |>
  # this ensures we take 50bp up & downstream of the summit for a total width of 100bp
  plyranges::anchor_center() |> 
  plyranges::mutate(width = 100)

# get sequences for regions
peak_seqs <- tissue_classes_100bp %>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$class) |>
  # look up the DNA sequence of each peak within each group
  get_sequence(dm.genome)

# filter for motif factors expressed in S2 cells or embryo =====================
expressed_genes <- RNAseq_RPKM |> 
  filter(if_any(3:ncol(RNAseq_RPKM), \(x) x > 1))

meme_db_expressed <- meme_db |> 
  filter(gene_id %in% expressed_genes$gene_id) |> 
  to_list()

# test for motif enrichment using AME ==========================================
zld_ame_results_db <- peak_seqs |> 
  runAme(database = meme_db_expressed) |>  
  bind_rows(.id = "class") |> 
  dplyr::group_by(class, motif_alt_id) |> 
  left_join(select(meme_db, name, gene_id, gene_symbol), by = join_by(motif_id == name))

# get RNA-seq RPKM data for motif factors =======================================
# get gene_ids for factors associated with each motif
motif_factors <- zld_ame_results_db |> 
  # filter(rank < 25) |>
  pull(gene_id)

# take the average RPKM across replicates in each tissue
zld_motif_factor_RNAseq <- RNAseq_RPKM |> 
  filter(gene_id %in% motif_factors) |> 
  pivot_longer(3:ncol(RNAseq_RPKM), names_to = "sample", values_to = "RPKM") |> 
  separate(sample, into = c("sample", "replicate"), sep="_(?=[^_]+$)") |> 
  select(-replicate) |> 
  mutate(sample = str_replace(sample, "_female", "")) |> 
  mutate(sample = str_replace(sample, "_male", "")) |> 
  group_by(gene_id, gene_symbol, sample) |> 
  summarise(mean_RPKM = mean(RPKM)) |> 
  pivot_wider(names_from = "sample", values_from = "mean_RPKM")

# zld_gene_order <- zld_motif_factor_RNAseq |> 
#   arrange(`S2-WT_noCuSO4`) |> 
#   pull(gene_symbol) |> 
#   unique()

# write output files ===========================================================
tmp <- zld_motif_factor_RNAseq |> 
  ungroup() |> 
  select(gene_id, `S2-WT_noCuSO4`)

zld_ame_results_db |> 
  left_join(tmp, by = "gene_id") |> 
  write_tsv("results/ChIP_tissue_classes/motifs/zld_AME_results.tsv")

zld_motif_factor_RNAseq |> 
  write_tsv("results/ChIP_tissue_classes/motifs/zld_motif_factor_RNAseq.tsv")

# ==============================================================================
# =================== prepare Grh motif and RNA-seq data =======================
# ==============================================================================

# define input files ===========================================================

tissue_classes <- grh_tissue_classes_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


embryo_RNAseq <- `embryo_15-16H_RNAseq_fn` |> 
  read_tsv()

wing_disc_RNAseq <- wing_disc_RNAseq_fn |> 
  read_tsv()

RNAseq_RPKM <- wing_disc_RNAseq |> 
  select(-c(2:4)) |> 
  left_join(embryo_RNAseq, by = "gene_id") |> 
  select(gene_id, gene_symbol, everything())

# get sequences underlying peaks ===============================================
# Get sequences in a 100bp window around the peak summit
tissue_classes_100bp <- tissue_classes |>
  # this ensures we take 50bp up & downstream of the summit for a total width of 100bp
  plyranges::anchor_center() |> 
  plyranges::mutate(width = 100)

# get sequences for regions
peak_seqs <- tissue_classes_100bp %>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$class) |>
  # look up the DNA sequence of each peak within each group
  get_sequence(dm.genome)


# filter for motif factors expressed in S2 cells or embryo =====================
expressed_genes <- RNAseq_RPKM |> 
  filter(if_any(3:ncol(RNAseq_RPKM), \(x) x > 1))

meme_db_expressed <- meme_db |> 
  filter(gene_id %in% expressed_genes$gene_id) |> 
  to_list()

# test for motif enrichment using AME ==========================================
grh_ame_results_db <- peak_seqs |> 
  runAme(database = meme_db_expressed) |>  
  bind_rows(.id = "class") |> 
  dplyr::group_by(class, motif_alt_id) |> 
  left_join(select(meme_db, name, gene_id, gene_symbol), by = join_by(motif_id == name))

# get RNA-seq RPKM data for motif factors =======================================
# get gene_ids for factors associated with each motif
motif_factors <- grh_ame_results_db |> 
  # filter(rank < 25) |>
  pull(gene_id)

# take the average RPKM across replicates in each tissue
grh_motif_factor_RNAseq <- RNAseq_RPKM |> 
  filter(gene_id %in% motif_factors) |> 
  pivot_longer(3:ncol(RNAseq_RPKM), names_to = "sample", values_to = "RPKM") |> 
  separate(sample, into = c("sample", "replicate"), sep="_(?=[^_]+$)") |> 
  select(-replicate) |> 
  group_by(gene_id, gene_symbol, sample) |> 
  summarise(mean_RPKM = mean(RPKM)) |> 
  pivot_wider(names_from = "sample", values_from = "mean_RPKM")

# grh_gene_order <- grh_motif_factor_RNAseq |> 
#   arrange(`S2-WT_noCuSO4`) |> 
#   pull(gene_symbol) |> 
#   unique()

# write output files ===========================================================
tmp <- grh_motif_factor_RNAseq |> 
  ungroup() |> 
  select(gene_id, `S2-WT_noCuSO4`)

grh_ame_results_db |> 
  left_join(tmp, by = "gene_id") |> 
  write_tsv("results/ChIP_tissue_classes/motifs/grh_AME_results.tsv")

grh_motif_factor_RNAseq |> 
  write_tsv("results/ChIP_tissue_classes/motifs/grh_motif_factor_RNAseq.tsv")


# ==============================================================================
# =================== prepare Twi motif and RNA-seq data =======================
# ==============================================================================

# define input files ===========================================================

tissue_classes <- twi_tissue_classes_fn |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

RNAseq_RPKM <- embryo_nc14_RNAseq_fn |> 
  read_tsv()

# get sequences underlying peaks ===============================================
# Get sequences in a 100bp window around the peak summit
tissue_classes_100bp <- tissue_classes |>
  # this ensures we take 50bp up & downstream of the summit for a total width of 100bp
  plyranges::anchor_center() |> 
  plyranges::mutate(width = 100)

# get sequences for regions
peak_seqs <- tissue_classes_100bp %>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$class) |>
  # look up the DNA sequence of each peak within each group
  get_sequence(dm.genome)

# filter for motif factors expressed in S2 cells or embryo =====================
expressed_genes <- RNAseq_RPKM |> 
  filter(if_any(3:ncol(RNAseq_RPKM), \(x) x > 1))

meme_db_expressed <- meme_db |> 
  filter(gene_id %in% expressed_genes$gene_id) |> 
  to_list()

# test for motif enrichment using AME ==========================================
twi_ame_results_db <- peak_seqs |> 
  runAme(database = meme_db_expressed) |>  
  bind_rows(.id = "class") |> 
  dplyr::group_by(class, motif_alt_id) |> 
  left_join(select(meme_db, name, gene_id, gene_symbol), by = join_by(motif_id == name))

# get RNA-seq RPKM data for motif factors =======================================
# get gene_ids for factors associated with each motif
motif_factors <- twi_ame_results_db |> 
  # filter(rank < 25) |>
  pull(gene_id)

# take the average RPKM across replicates in each tissue
twi_motif_factor_RNAseq <- RNAseq_RPKM |> 
  filter(gene_id %in% motif_factors) |> 
  pivot_longer(3:ncol(RNAseq_RPKM), names_to = "sample", values_to = "RPKM") |> 
  separate(sample, into = c("sample", "replicate"), sep="_(?=[^_]+$)") |> 
  select(-replicate) |> 
  mutate(sample = str_replace(sample, "_female", "")) |> 
  mutate(sample = str_replace(sample, "_male", "")) |> 
  group_by(gene_id, gene_symbol, sample) |> 
  summarise(mean_RPKM = mean(RPKM)) |> 
  pivot_wider(names_from = "sample", values_from = "mean_RPKM")

# twi_gene_order <- twi_motif_factor_RNAseq |> 
#   arrange(`S2-WT_noCuSO4`) |> 
#   pull(gene_symbol) |> 
#   unique()

# write output files ===========================================================
tmp <- twi_motif_factor_RNAseq |> 
  ungroup() |> 
  select(gene_id, `S2-WT_noCuSO4`)

twi_ame_results_db |> 
  left_join(tmp, by = "gene_id") |> 
  write_tsv("results/ChIP_tissue_classes/motifs/twi_AME_results.tsv")
twi_motif_factor_RNAseq |> 
  write_tsv("results/ChIP_tissue_classes/motifs/twi_motif_factor_RNAseq.tsv")
