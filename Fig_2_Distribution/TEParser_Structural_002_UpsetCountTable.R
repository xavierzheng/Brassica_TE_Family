#------------------------------
# AIM
#	List the TE family in each upset-group 
#
#=============================

# setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine")

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
library(UpSetR)
library(ggupset)


setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)


df <- fread(
  file = "/Data/Fig_2_Distribution/StructuralTE_TotalCount_FamilyLevel.txt", 
  header = T, sep = "\t"
)

glimpse(df)

df_list <- df %>%
  summarise(
    .by = c(SPECIES, FULL_NAME, TE_CLASS, TE_FAM), 
    TOTAL = n()
  ) %>%
  summarise(
    .by = c(TE_CLASS, TE_FAM), 
    LIST_ACCESS = list(FULL_NAME)
  ) %>%
  glimpse()

df_key <- df_list %>%
  summarise(
    .by = LIST_ACCESS, 
    PREVEL = n()
  ) %>%
  arrange(
    desc(PREVEL)
  ) %>%
  slice_head(n = 20) %>%
  mutate(
    CAT = 1:20
  ) %>%
  glimpse()

df_use <- left_join(
  x = df_key, 
  y = df_list, 
  by = "LIST_ACCESS"
) %>%
  unnest(
    cols = c(LIST_ACCESS)
  ) %>%
  rename(
    FULL_NAME = LIST_ACCESS
  ) %>%
  select(-PREVEL) %>%
  left_join(
    x = .,
    y = df, 
    by = c("FULL_NAME", "TE_CLASS", "TE_FAM")
  ) %>%
  mutate(
    FULL_NAME = factor(
      FULL_NAME, 
      levels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2"), 
      labels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2")
    ), 
    TE_CLASS = base::factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")
    ),
    SPECIES = factor(
      SPECIES,
      levels = c("B. rapa", "B. oleracea"),
      labels = c("B. rapa", "B. oleracea")
    )
  )  

df_use %>%
  mutate(
    COUNT_LEN = paste0(COUNT, ";", SUM_LEN)
  ) %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    CAT + TE_CLASS + TE_FAM + Share_or_Unique ~ FULL_NAME, 
    value.var = "COUNT", fill = 0
  ) %>%
  fwrite(
    file = "/Data/Fig_2_Distribution/StructuralTE_UpsetCount_FamilyLevel_Excel.txt", 
    col.names = T, sep = "\t"
  )

  