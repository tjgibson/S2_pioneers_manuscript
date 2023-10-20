# setup ========================================================================
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(Biostrings))
library(pheatmap)
library(GenomicRanges)
library(memes)

# define input files ===========================================================
meme_files <- c(
  # "results/ChIP_tissue_classes/motifs/grh_class_I/combined.meme",
  # "results/ChIP_tissue_classes/motifs/grh_class_II/combined.meme",
  # "results/ChIP_tissue_classes/motifs/grh_class_III/combined.meme",
  "results/ChIP_tissue_classes/motifs/grh_class_IV/combined.meme",
  "results/ChIP_tissue_classes/motifs/grh_class_V/combined.meme",
  "results/ChIP_tissue_classes/motifs/grh_class_VI/combined.meme"
)

class_names <- c(
  # "class_I",
  # "class_II",
  # "class_III",
  "class_IV",
  "class_V",
  "class_VI"
)

tissue_classes <- "results/ChIP_tissue_classes/grh_tissue_classes.tsv" |> 
  read_tsv() |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# import motfs for all classes =================================================
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


# visualize tree of all motifs to assess motif grouping ========================
# library(cowplot)
# pdf("exploratory_analysis/grh_motifs_tree.pdf")
# for (i in unique(all_motifs$group)) {
# # for (i in 1) {
#   ## Get the tree: make sure it's a horizontal type layout
#   motifs <- all_motifs |> 
#     filter(group == i) |> 
#     to_list()
#   
#   print(paste("group =", i))
#   print(paste("n motifs =", length(motifs)))
#   
#   if (length(motifs) == 1) {
#     print(view_motifs(motifs))
#     print("plotting single motif")
#   } else {
#     tree <- motif_tree(motifs, layout = "rectangular", linecol = "none")
#     ## Now, make sure we order our list of motifs to match the order of tips:
#     mot.names <- sapply(motifs, function(x) x["name"])
#     names(motifs) <- mot.names
#     new.order <- tree$data$label[tree$data$isTip]
#     new.order <- rev(new.order[order(tree$data$y[tree$data$isTip])])
#     motifs <- motifs[new.order]
#     ## Plot the two together (finessing of margins and positions may be required):
#     print("plotting motif tree")
#     print(
#       plot_grid(nrow = 1, rel_widths = c(1, -0.15, 1),
#                 tree + xlab(""), NULL,
#                 view_motifs(motifs, names.pos = "right") +
#                   ylab(element_blank()) +
#                   theme(
#                     axis.line.y = element_blank(),
#                     axis.ticks.y = element_blank(),
#                     axis.text.y = element_blank(),
#                     axis.text = element_text(colour = "white")
#                   )
#       )
#     )
#   }
# }
# 
# dev.off()

# filter motifs based on similarity groups =====================================
unique_motifs <- all_motifs |> 
  dplyr::group_by(group)  |> 
  # within each binding type, select the TF hit with the lowest adjusted p-value
  dplyr::slice(which.min(eval))

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

pdf("exploratory_analysis/grh_denovo_motifs_heatmap.pdf", height = 12, width = 8)
# visualize enrichment results
ame_results |>
  bind_rows(.id = "class") |>
  plot_ame_heatmap(group = class, scale_max = 50, id = motif_alt_id) +
  coord_flip()
dev.off()

# get instances of these motifs across classes
peak_seqs <- 
  tissue_classes |>
  get_sequence(dm.genome)

motif_counts <- "results/ChIP_tissue_classes/twi_tissue_classes.tsv" |> 
  read_tsv()

for (i in seq(nrow(unique_motifs)))  {
  message <- paste(i, "out of", nrow(unique_motifs), "started")
  print(message)
  motif <- unique_motifs[i,] |> to_list()
  motif_instances <- runFimo(peak_seqs, motif)
  if (!is.null(motif_instances)) {
    count <- GenomicRanges::countOverlaps(tissue_classes, motif_instances)
  } else {
    count <- 0
  }
  
  motif_name <- unique_motifs[i, "consensus"] |> pull()
  motif_counts[,motif_name] <- count
  
}

motif_consensus <- unique_motifs |> 
  ungroup() |> 
  dplyr::select(name, consensus, altname) |> 
  dplyr::rename(motif_name = "name")

motif_counts |> 
  pivot_longer(5:ncol(motif_counts), names_to = "motif_name", values_to = "n_motifs") |> 
  mutate(has_motif = n_motifs > 0) |> 
  group_by(class, motif_name) |> 
  summarise(percent_with_motif = round(mean(has_motif) * 100, 2)) |> 
  left_join(motif_consensus, by = "motif_name") |>  
  ggplot(aes(x=class, y = consensus, fill = percent_with_motif)) + 
  geom_tile(color = "black") +
  scale_fill_distiller(palette = "GnBu", direction = 1, limits=c(0,100), breaks = c(0,100), name = NULL) +
  # scale_x_discrete(labels = c("background", "class I", "class II", "class III")) +
  geom_text(aes(label = percent_with_motif)) +
  theme(
    axis.text.x = element_text(angle=45, vjust=1, hjust=1),
    axis.title = element_blank()
  )

# test enrichment of known motifs from flyfactorsurvey =========================
meme_db_path <- system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE)

control_seqs <- "results/ChIP_tissue_classes/twi_tissue_classes.tsv" |> 
  read_tsv() |> 
  dplyr::filter(class %in% c("class II", "class III")) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |> 
  get_sequence(dm.genome)

ame_results_db <- peak_seqs |>
  runAme(database = meme_db_path, control = control_seqs)

# visualize enrichment results
ame_results_db |>
  bind_rows(.id = "class") |>
  filter(rank < 25) |>
  plot_ame_heatmap(group = class, scale_max = 100)
