#-------------------
# AIM
#   Analysis the TE similarity difference and association between nucleotide identity and count prevalence
#
#-------------------

# setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

PREFIX <- "Relax"

df_stat <- fread(
  file = paste0("/Data/Fig_3_Nucleotide_Identity/AgeComparison_Top10_statistic_", PREFIX, ".txt"), 
  header = T, sep = "\t"
)


print("# Any correlation betweeen TE count and similarity ? --------------")

df_count <- fread(
  file = "/Data/Fig_3_Nucleotide_Identity/Compare_AA_CC_Count.txt",
  header = T, sep = "\t"
)

glimpse(df_count)

glimpse(df_stat)

df_stat %>%
  filter(
    str_detect(Prevelance, "^Unique", negate = T)
  ) %>%
  mutate(
    Prevelance = as.numeric(Prevelance), 
    SIG_FDR = ifelse(
      nchar(SIG_FDR) < 1,
      NA, SIG_FDR
    ),
    LAB4 = ifelse(
      is.na(SIG_FDR),
      "No_different",
      ifelse(
        Prevelance > 0,
        "AA_higher", 
        "CC_higher"
      )
    )
  ) %>%
  dplyr::left_join(
    x = ., 
    y = df_count %>% dplyr::rename(S_FAM = TE_NAME) %>% dplyr::select(S_FAM, CAT),
    by = "S_FAM"
  ) -> df_sim_count
glimpse(df_sim_count)

  
table_sim_count <- table(df_sim_count$LAB4, df_sim_count$CAT)
table_sim_count

##### using Relax parameter:
#               AA > CC CC > AA No different
# AA_higher          8       5           10
# CC_higher          1      23            2
# No_different       6      14           18

set.seed(1234)
chisq.test(table_sim_count, simulate.p.value = T, B = 1000)

##### using Relax parameter:
# Pearson's Chi-squared test with simulated p-value (based on 1000 replicates)
# 
# data:  table_sim_count
# X-squared = 27.719, df = NA, p-value = 0.000999