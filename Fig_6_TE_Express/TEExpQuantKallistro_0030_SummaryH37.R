#-------------------------
# AIM
#   could I used a better way to estimate the expression difference of TE family level 
#
#------------------------

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(ggsci))
library(broom)

setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

# Read Raw result from kallisto (No included)

name_list <- fread(
  file = "/Data/Fig_6_TE_Express/name.list", 
  header = F, sep = "\t"
) %>%
  select(V1) %>%
  pull

temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0("kallisto_RF_", NAM, "/abundance.tsv"), 
    header = T, sep = "\t"
  )
  
  df_temp$SAMPLE <- NAM
  
  return(df_temp)
  
})

df <- bind_rows(temp_list) %>%
  mutate(
    TE_ID = str_extract(target_id, "ID=[A-Za-z0-9_]{1,}") %>% str_remove(., "ID="),
    TE_NAME = str_extract(target_id, "Name=[A-Za-z0-9_]{1,}") %>% str_remove(., "Name="), 
    TE_TEMP = str_extract(target_id, "Classification=[A-Za-z0-9_/:\\-]{1,}"),  #%>% str_split_fixed(., "::", 2)[,1], 
    TE_CLASS = str_split_fixed(TE_TEMP, "::", 2)[,1] %>% str_remove(., "Classification="), 
    TE_TEMP2 = str_split_fixed(TE_TEMP, "::", 2)[,2], 
    TE_CHROM = str_split_fixed(TE_TEMP2, ":", 2)[,1],
    TE_TEMP3 = str_split_fixed(TE_TEMP2, ":", 2)[,2], #%>% str_split(., "\\-", n=2)[1]
    TE_BED_START = str_split_fixed(TE_TEMP3, "\\-", 2)[,1] %>% as.numeric(), 
    TE_BED_END = str_split_fixed(TE_TEMP3, "\\-", 2)[,2] %>% as.numeric()
  ) %>%
  select(
    -contains("TEMP")
  ) 

glimpse(df)
fwrite(
  df,
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_GeneTE_H37.txt", 
  col.names = T, sep = "\t"
)


## Start From this part !!!!!!!!!!!!!

df <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_GeneTE_H37.txt.gz", 
  header = T, sep = '\t'
)


df_total <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/summary.kallistro.map.H37.txt.gz", 
  header = T, sep = "\t"
) %>%
  mutate(
    SAMPLE = str_remove_all(READ1, "/nfs/project_ssd/project3/pxzhe/TE_exp/00_ReadILMN/") %>%
      str_remove_all(., "\\.R1.qc.fq.gz")
  ) %>%
  select(SAMPLE, TOTAL_READ, ALIGNED_READ)

glimpse(df_total)

glimpse(df)


print("# estimate FPKM and TPM by hand ----")
df_exp <- df %>%
  left_join(
    x = .,
    y = df_total, 
    by = "SAMPLE"
  ) %>%
  mutate(
    length = as.numeric(length), # 不要讓數值以整數形態存在，會遇到2^31上限
    FPKM = (est_counts*1e+9)/(ALIGNED_READ*length), # estimate FPKM
    TEMP = (est_counts*1e+3)/length # 這是為了TPM計算使用
  ) %>%
  group_by(
    SAMPLE
  ) %>%
  mutate(
    TPM = (TEMP/sum(TEMP))*1e+6
  ) %>%
  ungroup() %>%
  select(-TEMP)

glimpse(df_exp)

fwrite(
  df_exp, 
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_GeneTE_TPM_FPKM_H37.txt", 
  col.names = T, sep = "\t"
)

df_exp <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_GeneTE_TPM_FPKM_H37.txt.gz", 
  header = T, sep = "\t"
)

glimpse(df_exp)

print("# focus on TE expression ===================")
df_exp_TEFAM <- df_exp %>%
  filter(
    nchar(TE_ID)>0
  ) %>%
  mutate(
    SPECIES = ifelse(
      str_detect(TE_CHROM, "^D"), 
      "B_rapa", "B_oleracea"
    )
  ) %>%
  summarise(
    .by = c(SAMPLE, TE_CLASS, TE_NAME), 
    TE_COUNT = n(),
    TEFAM_TPM_MAX = max(TPM, na.rm = T),
    TEFAM_TPM_RANGE = (max(TPM, na.rm = T) - min(TPM, na.rm = T)),
    TEFAM_TPM_SUM = sum(TPM, na.rm = T),
    TEFAM_TPM_MED = median(TPM, na.rm = T), 
    TEFAM_TPM_AVE = mean(TPM, na.rm = T),
    TEFAM_TPM_25 = quantile(TPM, probs = 0.25, na.rm = T),
    TEFAM_TPM_75 = quantile(TPM, probs = 0.75, na.rm = T),
    TEFAM_SUM_LEN = sum(length, na.rm = T)
  ) %>%
  mutate(
    NORMALIZE_TPM = (TEFAM_TPM_SUM / TEFAM_SUM_LEN) * 1e+6
  )

glimpse(df_exp_TEFAM)

fwrite(
  x = df_exp_TEFAM,
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_H37.txt",
  col.names = T, 
  sep = "\t"
)

print("# Test New Normalize TPM parameter, their theshold -------------------")
summary(df_exp_TEFAM$NORMALIZE_TPM)

df_exp_TEFAM %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  arrange(
    desc(NORMALIZE_TPM)
  ) %>%
  slice_head(n=50)

print("### Test per superfamily quantile distribution -----------")
df_thres <- df_exp_TEFAM %>%
  summarise(
    .by = TE_CLASS, 
    T05 = quantile(NORMALIZE_TPM, probs = 0.05), 
    T50 = quantile(NORMALIZE_TPM, probs = 0.5), 
    T95 = quantile(NORMALIZE_TPM, probs = 0.95)
  )
df_thres

# TE_CLASS T05       T50        T95
# 1    LTR/Copia   0 0.4790757  10.598066
# 2      TIR/DTA   0 0.0000000  60.045336
# 3      TIR/DTC   0 0.3583670  75.155586
# 4      TIR/DTH   0 0.0000000 266.563090
# 5      TIR/DTM   0 0.0000000 219.085989
# 6      TIR/DTT   0 0.0000000 443.920815
# 7    LTR/Gypsy   0 0.2239069   5.794848
# 8 DNA/Helitron   0 2.5237794  62.215851
# 9  LTR/unknown   0 0.0000000  56.389983

fwrite(
  x = df_thres, 
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_threshold_H37.txt",
  col.names = T, sep = "\t"
)
