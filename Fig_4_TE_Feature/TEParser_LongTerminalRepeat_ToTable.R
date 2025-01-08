#=====================
# Aim
#	Combine LTR long terminal repeat GC and methylation site into a table
#	Key should be a universal name
#
#=====================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

name_list <- fread(
  file = "name.list", 
  header = F, sep = "\t"
) %>% select(V1) %>% pull

print("# Add length===========================")
rm(temp_list)

temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0(NAM, ".LongTerminalRepeat.length"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = NAM
    )
  return(df_temp)
  
})

df_len <- bind_rows(temp_list)

df_len <- df_len %>%
  mutate(
    TE_ID_LTR = str_split_fixed(V1, ";", 3)[,1] %>% str_remove_all(., "^ID="),
    TE_NAME = str_split_fixed(V1, ";", 3)[,2] %>% str_remove_all(., "^Name="),
    TEMP = str_split_fixed(V1, ";", 3)[,3]
  ) %>%
  mutate(
    TE_ID = str_replace(TE_ID_LTR, "[lr]LTR", replacement = "RR"),
    TE_CLASS = str_split_fixed(TEMP, ":", 4)[,1] %>% str_remove_all(., "Classification="), 
    CHROM = str_split_fixed(TEMP, ":", 4)[,3], 
    LEN = V2
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, LEN
  ) %>%
  glimpse()

df_len_uniq <- df_len %>%
  group_by(
    PLANT_ID, CHROM, TE_NAME, TE_CLASS, TE_ID
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() 


print("# read GC ======================")

rm(temp_list)
temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0(NAM, ".LongTerminalRepeat.gc"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = NAM
    )
  
  return(df_temp)
  
})

df_gc <- bind_rows(temp_list)

df_gc <- df_gc %>%
  mutate(
    TE_ID_LTR = str_split_fixed(V1, ";", 3)[,1] %>% str_remove_all(., "^ID="),
    TE_NAME = str_split_fixed(V1, ";", 3)[,2] %>% str_remove_all(., "^Name="),
    TEMP = str_split_fixed(V1, ";", 3)[,3]
  ) %>%
  mutate(
    TE_ID = str_replace(TE_ID_LTR, "[lr]LTR", replacement = "RR"),
    TE_CLASS = str_split_fixed(TEMP, ":", 4)[,1] %>% str_remove_all(., "Classification="), 
    CHROM = str_split_fixed(TEMP, ":", 4)[,3], 
    GC = V2
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, GC
  ) %>%
  left_join(
    x = .,
    y = df_len_uniq, 
    by = c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS")
  ) %>%
  glimpse()

fwrite(
  df_gc, 
  file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_GC.txt", 
  col.names = T, row.names = F, sep = "\t"
)

df_gc <- fread(
  file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_GC.txt", 
  header = T, sep = "\t"
)

# glimpse(df_gc)

print("## family median and IQR============================")
df_gc %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  summarise(
    MEDIAN_GC = median(GC), 
    Q25_GC = quantile(GC, probs = 0.25), 
    Q75_GC = quantile(GC, probs = 0.75), 
    No_of_StructualTE = n()/2
  ) %>%
  ungroup() %>%
  arrange(
    MEDIAN_GC
  ) %>%
  fwrite(
    file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_GC_PerFamily.txt", 
    col.names = T, row.names = F, sep = "\t", quote = F
  )


print("# read methylation file =====================================")
print("## Only keep the no. of row")

rm(temp_list)

temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_CG <- fread(
    file = paste0(NAM, ".LongTerminalRepeat.methyl_CG.bed"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = NAM
    )
  
  df_CHH <- fread(
    file = paste0(NAM, ".LongTerminalRepeat.methyl_CHH.bed"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = NAM
    )
  
  df_CHG <- fread(
    file = paste0(NAM, ".LongTerminalRepeat.methyl_CHG.bed"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      PLANT_ID = NAM
    )
  
  df_temp <- rbind(df_CG, df_CHH, df_CHG)
  return(df_temp)
  
})

df_methyl <- bind_rows(temp_list)

df_methyl <- df_methyl %>%
  mutate(
    TE_ID_LTR = str_split_fixed(V1, ";", 3)[,1] %>% str_remove_all(., "^ID="),
    TE_NAME = str_split_fixed(V1, ";", 3)[,2] %>% str_remove_all(., "^Name="),
    TEMP = str_split_fixed(V1, ";", 3)[,3]
  ) %>%
  mutate(
    TE_ID = str_replace(TE_ID_LTR, "[lr]LTR", replacement = "RR"),
    TE_CLASS = str_split_fixed(TEMP, ":", 4)[,1] %>% str_remove_all(., "Classification="), 
    CHROM = str_split_fixed(TEMP, ":", 4)[,3], 
    METHYL_LEN = nchar(V4)
  ) %>%
  rename(
    BED_START = V2, 
    BED_END = V3, 
    METHYL = V4, 
    STRAND = V6
  ) %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, BED_START, BED_END, METHYL, STRAND, METHYL_LEN
  ) 

# glimpse(df_methyl)

df_methyl_len <- df_methyl %>%
  group_by(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, METHYL
  ) %>%
  summarise(
    ACCUM_METHYL_LEN = sum(METHYL_LEN, na.rm = T)
  ) %>%
  ungroup() %>%
  left_join(
    x = ., 
    y = df_len_uniq, 
    by = c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS")
  ) %>%
  rename(
    ACCUM_LTR_LEN = LEN
  ) %>%
  ## we consider forward and reverse strand for methylation site
  ## so length should be double
  mutate(
    ACCUM_METHYL_PERC = (ACCUM_METHYL_LEN/(ACCUM_LTR_LEN*2))*100
  ) %>%
  arrange(
    TE_NAME, TE_ID, METHYL
  ) %>%
  glimpse()

View(df_methyl_len)

fwrite(
  df_methyl_len,
  file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_MethylSite_PerInsert.txt",
  col.name = T, row.name = F, sep = "\t", quote = F
)

