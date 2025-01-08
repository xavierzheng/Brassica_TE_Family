#==========================
# AIM
#	classify TE gene relationship
#	  1. TE contains complete gene
#	  2. TE within gene
#	  3. TE and gene have partial overlap
#	  4. Complex event: TE has more than 1 even list above
#	  5. TE doest not overlap with gene
#
#===========================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

name_list <- fread(
  file = "name.list", 
  header = F, sep = "\t"
) %>% select(V1) %>% pull


temp_list <- lapply(name_list, FUN = function(NAM){
  
  df <- fread(
    file = paste0("Gene_IntactTE_relation.", NAM, ".tab"), 
    header = F, sep = "\t"
  ) %>%
    select(V1, V4, V5, V7, V9, V10, V13, V14, V16, V18, V19)
  
  colnames(df) <- c("TE_CHROM", "TE_GFF_START", "TE_GFF_END", "TE_STRAND", "TEMP", "GENE_CHROM", "GENE_GFF_START", "GENE_GFF_END", "GENE_STRAND", "GENE_ID", "LEN_OVERLAP")
  
  #glimpse(df)
  
  # modify format
  df_format <- df %>%
    mutate(
      TE_ID = str_extract(TEMP, pattern = "ID=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "ID="),
      TE_NAME = str_extract(TEMP, pattern = "Name=[A-Za-z0-9_]{1,}") %>% 
        str_remove(., pattern = "Name="),
      TE_CLASS = str_extract(TEMP, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>% 
        str_remove(., pattern = "Classification="),
      GENE_ID = str_remove(GENE_ID, pattern = "ID=")
    ) %>%
    select(-TEMP) %>%
    mutate(
      TE_LEN = TE_GFF_END - TE_GFF_START +1,
      GENE_LEN = GENE_GFF_END - GENE_GFF_START +1
    )
  
  #glimpse(df_format)  
  
  # add tag
  df_info <- df_format %>%
    mutate(
      RELATION_GENE_TE = ifelse(
        !(GENE_LEN==LEN_OVERLAP | TE_LEN==LEN_OVERLAP), 
        "CLASS3_partial", 
        ifelse(
          LEN_OVERLAP==GENE_LEN & GENE_GFF_START>= TE_GFF_START & GENE_GFF_END<=TE_GFF_END, 
          "CLASS1_TEcontainGENE", 
          ifelse(
            LEN_OVERLAP==TE_LEN & TE_GFF_START>= GENE_GFF_START & GENE_GFF_END>=TE_GFF_END, 
            "CLASS2_TEwithinGENE", 
            "OTHER"
          )
        )
      )
    ) %>%
    select(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END, RELATION_GENE_TE, GENE_ID
    ) 
  
  df_temp <- df_info %>%
    group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # make sure if TE has complex relationship to gene, which means that they have more than 1 target gene, 
    # anothor tag should be added
    mutate(
      COUNT = n()
    ) %>%
    mutate(
      RELATION_GENE_TE = ifelse(
        COUNT > 1, 
        "CLASS4_ComplexEvent", 
        RELATION_GENE_TE
      )
    ) %>%
    ungroup() %>%
    select(-COUNT) %>%
    group_by(
      TE_ID, TE_NAME, TE_CLASS, TE_CHROM, TE_GFF_START, TE_GFF_END
    ) %>%
    # change format, put all gene ID to a vector, which seperated by __;
    mutate(
      temp_col = paste0("V", 1:n())
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      TE_ID+TE_NAME+TE_CLASS+TE_CHROM+TE_GFF_START+TE_GFF_END+RELATION_GENE_TE~temp_col, 
      value.var = "GENE_ID"
    ) %>%
    tidyr::unite(
      GENE_ID, 
      starts_with("V"), 
      sep = "__;",
      na.rm = T
    ) 
  
  
  # read total TE=========
  
  df_singleLine <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/", NAM, ".TE.InChr.IntactSingleLine.sort.gff3"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      TE_ID = str_extract(V9, pattern = "ID=[A-Za-z0-9_]{1,}") %>% str_remove(., pattern = "ID="),
      TE_NAME = str_extract(V9, pattern = "Name=[A-Za-z0-9_]{1,}") %>% str_remove(., pattern = "Name="),
      TE_CLASS = str_extract(V9, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>% str_remove(., pattern = "Classification=")
    ) %>%
    rename(
      TE_CHROM = V1,
      TE_GFF_START = V4, 
      TE_GFF_END = V5
    ) %>%
    select(
      starts_with("TE")
    )
  
  # left joint to find the intact TE which was not overlapped with gene
  left_join(
    x = df_singleLine, 
    y = df_temp,
    by = c("TE_ID", "TE_NAME", "TE_CLASS", "TE_CHROM", "TE_GFF_START", "TE_GFF_END")
  ) %>%
    mutate(
      RELATION_GENE_TE = ifelse(
        is.na(RELATION_GENE_TE),
        "CLASS5_NotOverlapWithGene", 
        RELATION_GENE_TE
      ), 
      PLANT_ID = NAM
    ) -> df_return
    
  return(df_return)
  
})

df_out <- bind_rows(temp_list)

df_out <- df_out %>%
  select(
    PLANT_ID, TE_CHROM, TE_GFF_START, TE_GFF_END, TE_ID, TE_NAME, TE_CLASS, RELATION_GENE_TE, GENE_ID
  ) 

fwrite(
  df_out,
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt", 
  col.names = T, row.names = F, sep = "\t", quote = T
)
  


# #==================================================
# print("# ANY characteristic?==============================")
# df_out <- fread(
#   file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt", 
#   header = T, sep = "\t"
# )

# df_out %>%
#   glimpse()

# table(df_out$TE_CLASS, df_out$RELATION_GENE_TE)

# chisq.test(
#   x = df_out$TE_CLASS, 
#   y = df_out$RELATION_GENE_TE
# )

# # Pearson's Chi-squared test
# # 
# # data:  df_out$TE_CLASS and df_out$RELATION_GENE_TE
# # X-squared = 16491, df = 32, p-value < 2.2e-16

# df_out_copia <- df_out %>%
#   filter(
#     TE_CLASS %in% "LTR/Copia"
#   )

# table(df_out_copia$TE_NAME, df_out_copia$RELATION_GENE_TE)
# chisq.test(
#   x = df_out_copia$TE_NAME, 
#   y = df_out_copia$RELATION_GENE_TE
# )

# # Pearson's Chi-squared test
# # 
# # data:  df_out_copia$TE_NAME and df_out_copia$RELATION_GENE_TE
# # X-squared = 3721.3, df = 316, p-value < 2.2e-16


# df_out_gypsy <- df_out %>%
#   filter(
#     TE_CLASS %in% "LTR/Gypsy"
#   )
# table(df_out_gypsy$TE_NAME, df_out_gypsy$RELATION_GENE_TE)
# chisq.test(
#   x = df_out_gypsy$TE_NAME, 
#   y = df_out_gypsy$RELATION_GENE_TE
# )

# # Pearson's Chi-squared test
# # 
# # data:  df_out_gypsy$TE_NAME and df_out_gypsy$RELATION_GENE_TE
# # X-squared = 2896.3, df = 180, p-value < 2.2e-16