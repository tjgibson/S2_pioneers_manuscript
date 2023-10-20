# setup ========================================================================
library(tidyverse)

# ==============================================================================
# for ZLD
# ==============================================================================

pdf("data/immunoblot_raw_images/2018-10-17_Zld_induction/quantification_plot.pdf", useDingbats = FALSE)
# defin input files ============================================================
raw_band_intensities <- "data/immunoblot_raw_images/2018-10-17_Zld_induction/quantification.tsv"

# analysis of band intensities =================================================
# get normalized intensity
gel_quantification <- raw_band_intensities |> 
  read_tsv(col_names = c("lane", "band", "Rf", "raw_volume", "cal_volume", "MW"), skip = 1) |> 
  select(lane, band, raw_volume) |> 
  pivot_wider(names_from = band, names_prefix = "band_", values_from = raw_volume) |> 
  add_column(sample =
               c(
                 "Zld_lineA_0", 
                 "Zld_lineA_750", 
                 "Zld_lineA_1000", 
                 "Zld_lineA_1250",
                 "Zld_lineA_1500",
                 
                 "Zld_lineB_0", 
                 "Zld_lineB_750", 
                 "Zld_lineB_1000", 
                 "Zld_lineB_1250",
                 "Zld_lineB_1500",
                 
                 "5 2-3H embryos A",
                 "5 2-3H embryos B",
                 "10 2-3H embryos A",
                 "10 2-3H embryos B"
               )
  ) |>
  mutate(normalized_intensity = band_1 / band_2) 


# get estimated WT level
wt_level <- gel_quantification[11:14, "normalized_intensity"] |> 
  pull() |> 
  mean()

# plot band intensities
gel_quantification |> 
  mutate(normalized_intensity = normalized_intensity / wt_level) |> 
  ggplot(aes(x = fct_inorder(sample), y = normalized_intensity)) + geom_bar(stat = "identity") +
  geom_hline(yintercept = wt_level, lty = 2, color = "grey") +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,), axis.title.x = element_blank()) 

dev.off()
# ==============================================================================
# for Grh
# ==============================================================================

pdf("data/immunoblot_raw_images/2018-11-08_Grh_induction/quantification_plot.pdf", useDingbats = FALSE)
# defin input files ============================================================
raw_band_intensities <- "data/immunoblot_raw_images/2018-11-08_Grh_induction/quantification.tsv"

# analysis of band intensities =================================================
# get normalized intensity
gel_quantification <- raw_band_intensities |> 
  read_tsv(col_names = c("lane", "band", "Rf", "raw_volume", "cal_volume", "MW"), skip = 1) |> 
  select(lane, band, raw_volume) |> 
  pivot_wider(names_from = band, names_prefix = "band_", values_from = raw_volume) |> 
  add_column(sample =
               c(
                 "Grh_lineA_0", 
                 "Grh_lineA_50", 
                 "Grh_lineA_100", 
                 "Grh_lineA_200",
                 "Grh_lineA_400",
                 
                 "Grh_lineB_0", 
                 "Grh_lineB_50", 
                 "Grh_lineB_100", 
                 "Grh_lineB_200",
                 "Grh_lineB_400",
                 
                 "5 2-3H embryos A",
                 "5 2-3H embryos B",
                 "10 2-3H embryos A",
                 "10 2-3H embryos B"
               )
  ) |>
  mutate(normalized_intensity = band_1 / band_2) 


# get estimated WT level
wt_level <- gel_quantification[11:14, "normalized_intensity"] |> 
  pull() |> 
  mean()

# plot band intensities
gel_quantification |> 
  mutate(normalized_intensity = normalized_intensity / wt_level) |> 
  ggplot(aes(x = fct_inorder(sample), y = normalized_intensity)) + geom_bar(stat = "identity") +
  geom_hline(yintercept = wt_level, lty = 2, color = "grey") +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,), axis.title.x = element_blank()) 


dev.off()

# ==============================================================================
# for Twi
# ==============================================================================

pdf("data/immunoblot_raw_images/2021-9-13_Twi-vs-embryo/quantification_plot.pdf", useDingbats = FALSE)
# defin input files ============================================================
raw_band_intensities <- "data/immunoblot_raw_images/2021-9-13_Twi-vs-embryo/quantification.tsv"

# analysis of band intensities =================================================
# get normalized intensity
gel_quantification <- raw_band_intensities |> 
  read_tsv(col_names = c("lane", "band", "Rf", "raw_volume", "cal_volume", "MW"), skip = 1) |> 
  select(lane, band, raw_volume) |> 
  pivot_wider(names_from = band, names_prefix = "band_", values_from = raw_volume) |> 
  add_column(sample =
               c(
                 "Twi_lineA_0", 
                 "Twi_lineA_5", 
                 "Twi_lineA_10", 
                 "Twi_lineA_20",
                 "Twi_lineA_40",
                 
                 "HA-Twi_lineA_0", 
                 "HA-Twi_lineA_5", 
                 "HA-Twi_lineA_10", 
                 "HA-Twi_lineA_20",
                 "HA-Twi_lineA_40",
                 
                 "5 3-4H embryos",
                 "10 3-4H embryos",
                 "20 3-4H embryos"
               )
  ) |>
  mutate(normalized_intensity = band_2 / band_1) 


# get estimated WT level
wt_level <- gel_quantification[11:13, "normalized_intensity"] |> 
  pull() |> 
  mean()

# plot band intensities
gel_quantification |> 
  mutate(normalized_intensity = normalized_intensity / wt_level) |> 
  ggplot(aes(x = fct_inorder(sample), y = normalized_intensity)) + geom_bar(stat = "identity") +
  geom_hline(yintercept = wt_level, lty = 2, color = "grey") +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 90, hjust=1,), axis.title.x = element_blank()) 


dev.off()
