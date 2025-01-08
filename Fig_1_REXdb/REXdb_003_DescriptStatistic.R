#=============================
# AIM
#	Estimate how many family can be assigned to REXdb
#
#============================


setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/Compare_REXdb")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggalluvial))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

INPUT <- "unknown"

df <- fread(
  file = paste0("Summary_REXdb_besthit_", INPUT, ".txt"), 
  header = T, sep = "\t", nThread = 12
)

#str(df)

df_report2 <- df[
  order(-LENGTH, -pident, evalue -bitscore), head(.SD, 1) ,by = .(qseqid)
][
  ,KEY_ORDER := ifelse(Q_ORDER==REX_ORDER, Q_ORDER, "DIFF")
][
  ,KEY_SUPFAM := ifelse(Q_SUPFAM==REX_SUPFAM, Q_SUPFAM, "DIFF")
]

#glimpse(df_report2)

fam_list <- sort(unique(df_report2$REX_FAM))

name_list <- sort(unique(df$Q_FAM))

df_draw <- df_report2 %>%
  group_by(
    Q_FAM, REX_FAM
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  mutate(
    Q_FAM = factor(
      Q_FAM, 
      levels = name_list, 
      labels = name_list
    ), 
    REX_FAM = factor(
      REX_FAM, 
      levels = fam_list, 
      labels = fam_list
    )
  ) 
  
df_allu <- df_draw %>%
  group_by(
    Q_FAM
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  mutate(
    GROUP = ifelse(
      COUNT == 1, 
      "Single", "Multiple"
    )
  ) %>%
  select(
    -COUNT
  ) %>%
  left_join(
    x =  df_draw, 
    y = ., 
    by = "Q_FAM"
  ) %>%
  arrange(
    Q_FAM, -COUNT
  ) %>%
  group_by(
    Q_FAM
  ) %>%
  mutate(
    PERC = COUNT/sum(COUNT)
  ) %>%
  ungroup() %>%
  group_by(
    Q_FAM
  ) %>%
  mutate(
    ANY = any(PERC > 0.95), 
    GROUP2 = ifelse(
      GROUP=="Multiple" && ANY==TRUE, 
      "Near_single", GROUP
    )
  ) %>%
  ungroup() %>%
  mutate(
    GROUP2 = factor(
      GROUP2,
      levels = c("Single", "Near_single", "Multiple"),
      labels = c("Single", "Near_single", "Multiple")
    )
  ) %>%
  select(
    -PERC, -ANY, -GROUP
  )

#glimpse(df_allu)

df_allu %>%
  group_by(
    Q_FAM
  ) %>%
  arrange(
    desc(COUNT)
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(
    GROUP2
  ) %>%
  summarise(
    COUNT2 = n()
  )


# Copia
# A tibble: 3 × 2
# GROUP2      COUNT2
# <fct>        <int>
# 1 Single          49
# 2 Near_single      8
# 3 Multiple         7


# Gypsy
# A tibble: 3 × 2
# GROUP2      COUNT2
# <fct>        <int>
# 1 Single          21
# 2 Near_single      5
# 3 Multiple         9


# unknown
# A tibble: 2 × 2
# GROUP2   COUNT2
# <fct>     <int>
# 1 Single        8
# 2 Multiple      6