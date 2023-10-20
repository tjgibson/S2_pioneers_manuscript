# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
library(pheatmap)
library(GenomicRanges)
library(memes)
library(RColorBrewer)

# define input files ===========================================================
meme_files <- c(
  "results/Zld_NB_OE/motifs/zld_NB_endogenous/combined.meme",
  "results/Zld_NB_OE/motifs/zld_NB_novel/combined.meme"

)

tissue_classes <- "results/Zld_NB_OE/peak_classes.tsv" |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# import motfs for all classes =================================================
class_names <- c(
  "endogenous",
  "OE"
)

all_motifs <- meme_files |> 
  map(read_meme) 
names(all_motifs) <- class_names


# combine all motifs into data frame ===========================================
all_motifs <- all_motifs |> 
  map(convert_type, "PWM", pseudocount = 0.1) |> 
  map(to_df) |> 
  bind_rows(.id = "class") |> 
  mutate(name = paste0(class, "_", name))

# compute motif similarity =====================================================
motif_similarity <- all_motifs %>%
  # Convert to universalmotif format 
  to_list() %>%
  # Compute the pearson correlation for each motif with all other motifs
  universalmotif::compare_motifs(method = "PCC")

# group motifs based on motif similarity =======================================
threshold <- 0.9
similarity_binary <- motif_similarity > threshold
d <- dist(similarity_binary, method = "binary")
hc <- hclust(d)
motif_groups <- cutree(hc, h = threshold)
motif_groups |> unique() |> length()

all_motifs$group <- motif_groups

# make heatmap of motif similarity =============================================
# This is for adding the colored annotation blocks indicating group membership
# to the heatmap
anno.df <- all_motifs %>% 
  dplyr::select(name, class) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("name")

# perform hierarchical clustering and add clusters to anno-df
anno.df$group <- as.character(motif_groups)

# Plot the correlation matrix along with the annotations
motif_similarity |> 
  pheatmap::pheatmap(
    annotation_col = anno.df,
    show_colnames = FALSE) 


# filter motifs based on similarity groups =====================================
unique_motifs <- all_motifs |> 
  dplyr::group_by(group)  |> 
  # within each binding type, select the TF hit with the lowest adjusted p-value
  dplyr::slice(which.min(eval))


# write motifs to table
unique_motifs |> 
  saveRDS("results/Zld_NB_OE/motifs/meme_chip_motifs.rds")

# test for enrichment of these motifs across tissue classes ====================
dm.genome <- BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6

# Get sequences in a 100bp window around the peak summit
tissue_classes_100bp <- tissue_classes %>%
  # this ensures we take 50bp up & downstream of the summit for a total width of 100bp
  plyranges::anchor_center() %>% 
  plyranges::mutate(width = 100)

# get sequences for regions
peak_seqs <- tissue_classes_100bp %>%
  # Get a list of chip peaks belonging to each set
  split(mcols(.)$class) %>%
  # look up the DNA sequence of each peak within each group
  get_sequence(dm.genome)

# run AME to calculate motif enrichment
ame_results <- peak_seqs |> 
  runAme(database = to_list(unique_motifs))

# output AME results to table
ame_results |>
  bind_rows(.id = "class") |>
  write_tsv("results/Zld_NB_OE/motifs/meme_chip_motifs_AME_enrichment.tsv")

