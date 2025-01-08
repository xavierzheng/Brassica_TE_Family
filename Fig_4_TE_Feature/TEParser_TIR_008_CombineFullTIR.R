#=====================================
# AIM
#	Combine full length of GC and methylation motif in TIR
#
#------------------------------------

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))

setDTthreads(threads = 8)
options(dplyr.summarise.inform=FALSE)

print("# read GC of TIR full length =============================")
list_TIR <- c("DTA", "DTC", "DTH", "DTM", "DTT", "Helitron")

temp_list <- lapply(list_TIR, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0("FullLength_TIR_", NAM, ".gc"), 
    header = T, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = str_split_fixed(FASTA_HEAD, ";", 6)[,1],
      CHROM = str_split_fixed(FASTA_HEAD, ";", 6)[,2],
      TE_ID = str_split_fixed(FASTA_HEAD, ";", 6)[,3],
      TE_NAME = str_split_fixed(FASTA_HEAD, ";", 6)[,4],
      TE_CLASS = str_split_fixed(FASTA_HEAD, ";", 6)[,5],
      TE_CLASS = str_remove_all(TE_CLASS, "\\(\\.\\)")
    ) %>%
    select(
      -FASTA_HEAD
    )
  
  return(df_temp)
})

df_TIR_GC <- bind_rows(temp_list)

glimpse(df_TIR_GC)
tail(df_TIR_GC)

fwrite(
  df_TIR_GC, 
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_GC.txt", 
  col.names = T, sep = "\t"
)

print("# read methylation motif, CG, CHH, and CHG ==========================")
rm(temp_list)
temp_list <- lapply(list_TIR, FUN = function(NAM){
  
  temp_1 <- lapply(c("CG", "CHH", "CHG"), FUN = function(MET){
    
    df_temp <- fread(
      file = paste0("TIR_", NAM,".methyl_", MET,".bed"), 
      header = F, sep = "\t"
    )
    return(df_temp)
  })
  
  df_ret <- bind_rows(temp_1)
  return(df_ret)
})

df_TIR_meth <- bind_rows(temp_list)


print("# take a long time --------------")

df_TIR_meth_temp <- df_TIR_meth %>%
  rename(
    FASTA_HEAD = V1, 
    METHYL = V4
  ) %>%
  mutate(
    METHYL_BASE = nchar(METHYL), 
    PLANT_ID = str_split_fixed(FASTA_HEAD, ";", 6)[,1],
    CHROM = str_split_fixed(FASTA_HEAD, ";", 6)[,2],
    TE_ID = str_split_fixed(FASTA_HEAD, ";", 6)[,3],
    TE_NAME = str_split_fixed(FASTA_HEAD, ";", 6)[,4],
    TE_CLASS = str_split_fixed(FASTA_HEAD, ";", 6)[,5],
    TE_CLASS = str_remove_all(TE_CLASS, "\\(\\.\\)")
  ) %>%
  group_by(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, METHYL
  ) %>%
  summarise(
    ACCUM_METHYL_LEN = sum(METHYL_BASE, na.rm = T)
  ) %>%
  ungroup() %>%
  left_join(
    x = .,
    y = df_TIR_GC, 
    by = c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS")
  ) 


df_temp1 <- df_TIR_meth_temp %>%
  mutate(
    TE_CLASS = str_remove_all(TE_CLASS, "\\([+-]\\)"), 
    ACCUM_METHYL_PERC = (ACCUM_METHYL_LEN/LEN)*100
  ) %>%
  glimpse()

fwrite(
  df_temp1, 
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_methyl_raw.txt", 
  col.names = T, sep = "\t"
)

print("# change to friendly table ==================")
df_temp_CG <- df_temp1 %>%
  filter(
    METHYL %in% "CG"
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, ACCUM_METHYL_PERC
  ) %>%
  rename(
    ACCUM_METHYL_CG_PERC = ACCUM_METHYL_PERC
  )

df_temp_CHH <- df_temp1 %>%
  filter(
    METHYL %in% "CHH"
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, ACCUM_METHYL_PERC
  ) %>%
  rename(
    ACCUM_METHYL_CHH_PERC = ACCUM_METHYL_PERC
  )

df_temp_CHG <- df_temp1 %>%
  filter(
    METHYL %in% "CHG"
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, ACCUM_METHYL_PERC
  ) %>%
  rename(
    ACCUM_METHYL_CHG_PERC = ACCUM_METHYL_PERC
  )

df_out_meth <- df_temp_CG %>%
  left_join(
    x = .,
    y = df_temp_CHG, 
    by = c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS")
  ) %>%
  left_join(
    x = .,
    y = df_temp_CHH, 
    by = c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS")
  )

fwrite(
  df_out_meth, 
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_methyl.txt", 
  col.names = T, sep = "\t"
)  
