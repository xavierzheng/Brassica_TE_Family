#======================================
# AIM
#	From TE_INFO_AA_CC.txt to generate TE frequency heatmap
#
#======================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(pheatmap, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

df <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_AA_CC.txt.gz", 
  header = T, sep = "\t"
) %>%
  filter(
    # remove NA
    !is.na(TE_NAME)
  ) %>%
  filter(
    # remove strange TE family, they should be NA
    nchar(TE_NAME) > 1
  )

df_plantid <- fread(
  file = "/Data/Fig_4_TE_Feature/Plant_IDNAME.txt", 
  header = T, sep = "\t"
)

df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt", 
  header = T, sep = "\t"
) %>%
  rename(
    TE_NAME = TE_Family, 
    TE_CLASS = TE_superfamily,
    ShareUniq = Share_or_Unique
  )

# fwrite(
#   df,
#   file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt",
#   col.names = T, row.names = F, sep = "\t", quote = F
# )

CLAS <- "Total" 

df <- df %>%
  left_join(
    x = .,
    y = df_plantid,
    by = "PLANT_ID"
  )

lapply(c("homology", "structural"), FUN = function(NAM){
  
  df_temp <- df %>%
    filter(
      TE_METHOD %in% NAM
    ) 
  
  # How many TE?========================="
  df_temp %>%
    group_by(
      PLANT_NAME, TE_CLASS
    ) %>%
    summarise(
      COUNT = n()
    ) %>%
    ungroup() %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_CLASS ~ PLANT_NAME, value.var = "COUNT"
    ) %>%
    fwrite(
      file = paste0("TE", NAM,"_", CLAS,"_No.txt"), # WithinFam
      col.names = T, row.names = F, sep = "\t", quote = F
    )
  
  # what is the lenght of TE?=========================
  df_temp %>%
    mutate(
      TE_LEN = GFF_END - GFF_START + 1
    ) %>%
    group_by(
      PLANT_NAME, TE_CLASS
    ) %>%
    summarise(
      SUM_TE_LEN = sum(TE_LEN)
    ) %>%
    ungroup() %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_CLASS ~ PLANT_NAME, value.var = "SUM_TE_LEN"
    ) %>%
    fwrite(
      file = paste0("TE", NAM,"_", CLAS,"_Len.txt"), # WithinFam
      col.names = T, row.names = F, sep = "\t", quote = F
    )
  
})





print("# prepare data frame for heatmap=======================================")


#===============================================================
#===============================================================
#===============================================================

#### HERE IS STRUCTURAL ONLY
lapply(c("DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"), FUN = function(NAM){
  
  OUT_PREFIX <- str_replace_all(NAM, "/", "_")
  
  df_draw_temp <- df %>%
    filter(
      TE_METHOD %in% "structural"
    ) %>%
    filter(
      TE_CLASS %in% NAM
    ) %>%
    group_by(
      PLANT_NAME, TE_NAME
    ) %>%
    summarise(
      COUNT = n()
    ) %>%
    ungroup() %>%
    as.data.table() %>%
    data.table::dcast(
      TE_NAME ~ PLANT_NAME, value.var = "COUNT", fill = 0
    ) %>%
    mutate(
      ROWSUM = rowSums(across(where(is.numeric)))
    ) %>%
    # filter(
    #   ROWSUM > 20
    # ) %>%
    select(
      -ROWSUM
    ) %>%
    select(
      TE_NAME, Z1, Chiifu, CXA, MIZ, OIB, PCA, TUA, 
      RCBO, JZSv2, OX_heart, HDEM, Korso
    )
  
  
  # z scale by me. If any raw value is 0, I forced z-scaled results as -3
  df_draw_temp_Z <- df_draw_temp %>%
    as.data.table() %>%
    data.table::melt.data.table(
      id.vars = "TE_NAME"
    ) %>%
    group_by(
      TE_NAME
    ) %>%
    mutate(
      SD = sd(value, na.rm = T)
    ) %>%
    ungroup() %>%
    filter(
      !(SD %in% 0)
    ) %>%
    filter(
      is.finite(SD)
    ) %>%
    group_by(
      TE_NAME
    ) %>%
    mutate(
      ZSCALE = (value -mean(value))/sd(value), 
      variable = factor(
        variable, 
        levels = c("Z1", "Chiifu", "CXA", "MIZ", "OIB", "PCA", "TUA", "RCBO", "JZSv2", "OX_heart", "HDEM", "Korso"), 
        labels = c("Z1", "Chiifu", "CXA", "MIZ", "OIB", "PCA", "TUA", "RCBO", "JZSv2", "OX_heart", "HDEM", "Korso")
      )
    ) %>%
    mutate(
      ZSCALE = ifelse(
        value == 0,
        -3, 
        ZSCALE
      )
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_NAME ~ variable, value.var = "ZSCALE"
    )
  
  fwrite(
    df_draw_temp_Z, 
    file = paste0("/Data/Fig_4_TE_Feature/Heatmap_TE_Structural/heatmap_TEStructural_AddShareUniq_Zscale_", OUT_PREFIX, ".txt"), 
    col.names = T, row.names = F, sep = "\t", quote = F
  )
  
  
})


