#==========================
# AIM
#	classify TE and TE relationship
#	  1. intact TE contains another TE
#	  2. intact TE within another TE
#	  3. TE and TE have partial overlap
#	  4. Complex event: TE has more than 1 even list above
#	  5. TE doest not overlap with TE
#
#===========================


suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

name_list <- fread(
  file = "name.list", 
  header = F, sep = "\t"
) %>% dplyr::select(V1) %>% pull


temp_list <- lapply(name_list, FUN = function(NAM){
  
  df <- fread(
    file = paste0("IntactTE_TotalTE_relation.", NAM, ".tab"), 
    header = F, sep = "\t"
  ) %>%
    select(V1, V4, V5, V7, V9, V10, V13, V14, V16, V18, V19)
  
  colnames(df) <- c("TE_CHROM", "TE_GFF_START", "TE_GFF_END", "TE_STRAND", "TEMP1", "OTHERTE_CHROM", "OTHERTE_GFF_START", "OTHERTE_GFF_END", "OTHERTE_STRAND", "TEMP2", "LEN_OVERLAP")
  
  #glimpse(df)
  
  # modify format
  df_format <- df %>%
    dplyr::mutate(
      TE_ID = str_extract(TEMP1, pattern = "ID=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "ID="),
      TE_NAME = str_extract(TEMP1, pattern = "Name=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "Name="),
      TE_CLASS = str_extract(TEMP1, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>% 
        str_remove(., pattern = "Classification="),
      OTHERTE_ID = str_extract(TEMP2, pattern = "ID=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "ID="),
      OTHERTE_NAME = str_extract(TEMP2, pattern = "Name=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "Name="),
      OTHERTE_CLASS = str_extract(TEMP2, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>% 
        str_remove(., pattern = "Classification="),
    ) %>%
    dplyr::select(-TEMP1, -TEMP2) %>%
    dplyr::mutate(
      TE_LEN = TE_GFF_END - TE_GFF_START +1,
      OTHERTE_LEN = OTHERTE_GFF_END - OTHERTE_GFF_START +1
    )
  
  #glimpse(df_format)  
  
  # add tag
  df_info <- df_format %>%
    dplyr::mutate(
      RELATION_OTHERTE_TE = ifelse(
        !(OTHERTE_LEN==LEN_OVERLAP | TE_LEN==LEN_OVERLAP), 
        "CLASS3_partial", 
        ifelse(
          LEN_OVERLAP==OTHERTE_LEN & OTHERTE_GFF_START>= TE_GFF_START & OTHERTE_GFF_END <=TE_GFF_END, 
          "CLASS1_TEcontainOTHERTE", 
          ifelse(
            LEN_OVERLAP==TE_LEN & TE_GFF_START>= OTHERTE_GFF_START & OTHERTE_GFF_END>=TE_GFF_END, 
            "CLASS2_TEwithinOTHERTE", 
            "OTHER"
          )
        )
      )
    ) %>%
    dplyr::select(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END, RELATION_OTHERTE_TE, 
      OTHERTE_ID, OTHERTE_NAME, OTHERTE_CLASS
    ) 
  
  df_temp1 <- df_info %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # make sure if TE has complex relationship to OTHERTE, which means that they have more than 1 target OTHERTE, 
    # anothor tag should be added
    dplyr::mutate(
      COUNT = n()
    ) %>%
    dplyr::mutate(
      RELATION_OTHERTE_TE = ifelse(
        COUNT > 1, 
        "CLASS4_ComplexEvent", 
        RELATION_OTHERTE_TE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-COUNT) %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # change format, put all OTHERTE ID to a vector, which seperated by __;
    dplyr::mutate(
      temp_col = paste0("V", 1:n())
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_ID+TE_NAME+TE_CLASS+TE_CHROM+TE_GFF_START+TE_GFF_END+RELATION_OTHERTE_TE~temp_col, 
      value.var = "OTHERTE_ID"
    ) %>%
    tidyr::unite(
      OTHERTE_ID, 
      starts_with("V"), 
      sep = "__;",
      na.rm = T
    ) 
  
  df_temp2 <- df_info %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # make sure if TE has complex relationship to OTHERTE, which means that they have more than 1 target OTHERTE, 
    # anothor tag should be added
    dplyr::mutate(
      COUNT = n()
    ) %>%
    dplyr::mutate(
      RELATION_OTHERTE_TE = ifelse(
        COUNT > 1, 
        "CLASS4_ComplexEvent", 
        RELATION_OTHERTE_TE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-COUNT) %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # change format, put all OTHERTE ID to a vector, which seperated by __;
    dplyr::mutate(
      temp_col = paste0("V", 1:n())
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_ID+TE_NAME+TE_CLASS+TE_CHROM+TE_GFF_START+TE_GFF_END+RELATION_OTHERTE_TE~temp_col, 
      value.var = "OTHERTE_NAME"
    ) %>%
    tidyr::unite(
      OTHERTE_NAME, 
      starts_with("V"), 
      sep = "__;",
      na.rm = T
    ) 
  
  df_temp3 <- df_info %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # make sure if TE has complex relationship to OTHERTE, which means that they have more than 1 target OTHERTE, 
    # anothor tag should be added
    dplyr::mutate(
      COUNT = n()
    ) %>%
    dplyr::mutate(
      RELATION_OTHERTE_TE = ifelse(
        COUNT > 1, 
        "CLASS4_ComplexEvent", 
        RELATION_OTHERTE_TE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-COUNT) %>%
    dplyr::group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # change format, put all OTHERTE ID to a vector, which seperated by __;
    dplyr::mutate(
      temp_col = paste0("V", 1:n())
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_ID+TE_NAME+TE_CLASS+TE_CHROM+TE_GFF_START+TE_GFF_END+RELATION_OTHERTE_TE~temp_col, 
      value.var = "OTHERTE_CLASS"
    ) %>%
    tidyr::unite(
      OTHERTE_CLASS, 
      starts_with("V"), 
      sep = "__;",
      na.rm = T
    ) 
  
  df_temp <- left_join(
    x = df_temp1, 
    y = df_temp2, 
    by = c("TE_ID", "TE_NAME", "TE_CLASS", "TE_CHROM", "TE_GFF_START", "TE_GFF_END", "RELATION_OTHERTE_TE")
  ) %>%
    left_join(
      x = . ,
      y = df_temp3, 
      by = c("TE_ID", "TE_NAME", "TE_CLASS", "TE_CHROM", "TE_GFF_START", "TE_GFF_END", "RELATION_OTHERTE_TE")
    )
  
  
  # read total TE=========
  
  df_singleLine <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/", NAM, ".TE.InChr.IntactSingleLine.sort.gff3"), 
    header = F, sep = "\t"
  ) %>%
    dplyr::mutate(
      TE_ID = str_extract(V9, pattern = "ID=[A-Za-z0-9_]{1,}") %>% str_remove(., pattern = "ID="),
      TE_NAME = str_extract(V9, pattern = "Name=[A-Za-z0-9_]{1,}") %>% str_remove(., pattern = "Name="),
      TE_CLASS = str_extract(V9, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>% str_remove(., pattern = "Classification=") 
    ) %>%
    dplyr::rename(
      TE_CHROM = V1,
      TE_GFF_START = V4, 
      TE_GFF_END = V5
    ) %>%
    dplyr::select(
      starts_with("TE")
    )
  
  # left joint to find the intact TE which was not overlapped with OTHERTE
  left_join(
    x = df_singleLine, 
    y = df_temp,
    by = c("TE_ID", "TE_NAME", "TE_CLASS", "TE_CHROM", "TE_GFF_START", "TE_GFF_END")
  ) %>%
    dplyr::mutate(
      RELATION_OTHERTE_TE = ifelse(
        is.na(RELATION_OTHERTE_TE),
        "CLASS5_NotOverlapWithOTHERTE", 
        RELATION_OTHERTE_TE
      ), 
      PLANT_ID = NAM
    ) -> df_return
  
  return(df_return)
  
})

df_out <- bind_rows(temp_list)

glimpse(df_out)

df_out <- df_out %>%
  dplyr::select(
    PLANT_ID, TE_CHROM, TE_GFF_START, TE_GFF_END, TE_ID, TE_NAME, TE_CLASS, RELATION_OTHERTE_TE, 
    OTHERTE_ID, OTHERTE_NAME, OTHERTE_CLASS
  ) 

fwrite(
  df_out,
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Relationship.txt", 
  col.names = T, row.names = F, sep = "\t", quote = T
)
