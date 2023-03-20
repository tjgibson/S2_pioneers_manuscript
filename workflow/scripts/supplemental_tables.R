# setup ========================================================================
library(tidyverse)
library(openxlsx)

# supplemental table 1 =========================================================
zld_ChIP_classes <- read_tsv("results/ChIP_peak_classes/zld_ChIP_classes.tsv") |> 
 mutate(class = str_to_upper(class))
grh_ChIP_classes <- read_tsv("results/ChIP_peak_classes/grh_ChIP_classes.tsv") |> 
  mutate(class = str_to_upper(class))
twi_ChIP_classes <- read_tsv("results/ChIP_peak_classes/twi_ChIP_classes.tsv") |> 
  mutate(class = str_to_upper(class))

supp_table_1 <- list(
  Zld_classes = zld_ChIP_classes,
  Grh_classes = grh_ChIP_classes,
  Twi_classes = twi_ChIP_classes
)

openxlsx::write.xlsx(supp_table_1, file = "manuscript/supplemental_tables/supplemental_table_1.xlsx") 


# supplemental table 2 =========================================================
published_ChIP_data <- read_tsv("config/published_ChIP_units_all.tsv") |> 
  distinct(sample_group, .keep_all = TRUE) |> 
  select(sample_group, GSE_accession)

openxlsx::write.xlsx(published_ChIP_data, file = "manuscript/supplemental_tables/supplemental_table_2.xlsx") 
