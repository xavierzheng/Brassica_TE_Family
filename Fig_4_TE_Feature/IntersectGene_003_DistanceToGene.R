#========================================
# AIM
# 	Estimate the distance between intack TE and closest gene
#	Only in class 5
#
#========================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

name_list <- fread(
  file = "name.list", 
  header = F, sep = "\t"
) %>% select(V1) %>% pull

# import summary =======================================================
df_info <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt", 
  header = T, sep = "\t"
)

#glimpse(df_info)

df_info_class5 <- df_info %>%
  filter(
    RELATION_GENE_TE %in% "CLASS5_NotOverlapWithGene"
  )


# per-plant ============================================================

temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_gene <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/GeneAnnotate_", NAM,".InChr.sort.gff.gz"), 
    header = F, sep = "\t"
  ) %>%
    filter(
      V3 %in% "gene"
    ) %>%
    mutate(
      GENE_ID = str_extract(V9, pattern = "ID=[A-Za-z0-9._]{1,}") %>%
        str_remove(., pattern = "ID="), 
      FEATURE = "GENE"
    ) %>%
    select(
      V1, V4, V5, GENE_ID, FEATURE
    ) %>%
    rename(
      CHROM = V1, 
      GFF_START = V4, 
      GFF_END = V5, 
      FEATURE_ID = GENE_ID
    ) 
  
  df_singleLine <- df_info_class5 %>%
    filter(
      PLANT_ID %in% NAM
    ) %>%
    rename(
      CHROM = TE_CHROM, 
      GFF_START = TE_GFF_START, 
      GFF_END = TE_GFF_END, 
      FEATURE_ID = TE_ID
    ) %>%
    mutate(
      FEATURE = "TE"
    ) %>%
    select(
      CHROM, GFF_START, GFF_END, FEATURE_ID, TE_NAME, TE_CLASS, FEATURE
    )
  
  df_dist_TE_GENE <- df_singleLine %>%
    select(
      CHROM, GFF_START, GFF_END, FEATURE, FEATURE_ID
    ) %>%
    rbind(., df_gene) %>%
    arrange(
      CHROM, GFF_START
    ) %>%
    mutate(
      PREVIOUS_GENE_ID = ifelse(
        FEATURE == "TE" & CHROM==lag(CHROM, n = 1) & lag(FEATURE, n = 1)=="GENE",
        lag(FEATURE_ID, n = 1), 
        ifelse(
          FEATURE == "TE" & CHROM==lag(CHROM, n = 2) & lag(FEATURE, n = 2)=="GENE",
          lag(FEATURE_ID, n = 2),
          ifelse(
            FEATURE == "TE" & CHROM==lag(CHROM, n = 3) & lag(FEATURE, n = 3)=="GENE",
            lag(FEATURE_ID, n = 3),
            ifelse(
              FEATURE == "TE" & CHROM==lag(CHROM, n = 4) & lag(FEATURE, n = 4)=="GENE",
              lag(FEATURE_ID, n = 4),
              ifelse(
                FEATURE == "TE" & CHROM==lag(CHROM, n = 5) & lag(FEATURE, n = 5)=="GENE",
                lag(FEATURE_ID, n = 5),
                NA
              )
            )
          )
        )
      ), 
      PREVIOUS_GENE_DISTANCE = ifelse(
        FEATURE == "TE" & CHROM==lag(CHROM, n = 1) & lag(FEATURE, n = 1)=="GENE",
        GFF_START - lag(GFF_END, n = 1) +1, 
        ifelse(
          FEATURE == "TE" & CHROM==lag(CHROM, n = 2) & lag(FEATURE, n = 2)=="GENE",
          GFF_START - lag(GFF_END, n = 2) +1,
          ifelse(
            FEATURE == "TE" & CHROM==lag(CHROM, n = 3) & lag(FEATURE, n = 3)=="GENE",
            GFF_START - lag(GFF_END, n = 3) +1,
            ifelse(
              FEATURE == "TE" & CHROM==lag(CHROM, n = 4) & lag(FEATURE, n = 4)=="GENE",
              GFF_START - lag(GFF_END, n = 4) +1,
              ifelse(
                FEATURE == "TE" & CHROM==lag(CHROM, n = 5) & lag(FEATURE, n = 5)=="GENE",
                GFF_START - lag(GFF_END, n = 5) +1,
                NA
              )
            )
          )
        )
      ), 
      NEXT_GENE_ID = ifelse(
        FEATURE == "TE" & CHROM==lead(CHROM, n = 1) & lead(FEATURE, n = 1)=="GENE",
        lead(FEATURE_ID, n = 1),
        ifelse(
          FEATURE == "TE" & CHROM==lead(CHROM, n = 2) & lead(FEATURE, n = 2)=="GENE",
          lead(FEATURE_ID, n = 2),
          ifelse(
            FEATURE == "TE" & CHROM==lead(CHROM, n = 3) & lead(FEATURE, n = 3)=="GENE",
            lead(FEATURE_ID, n = 3),
            ifelse(
              FEATURE == "TE" & CHROM==lead(CHROM, n = 4) & lead(FEATURE, n = 4)=="GENE",
              lead(FEATURE_ID, n = 4),
              ifelse(
                FEATURE == "TE" & CHROM==lead(CHROM, n = 5) & lead(FEATURE, n = 5)=="GENE",
                lead(FEATURE_ID, n = 5),
                NA
              )
            )
          )
        )
      ),
      NEXT_GENE_DISTANCE = ifelse(
        FEATURE == "TE" & CHROM==lead(CHROM, n = 1) & lead(FEATURE, n = 1)=="GENE",
        lead(GFF_START, n = 1) - GFF_END +1,
        ifelse(
          FEATURE == "TE" & CHROM==lead(CHROM, n = 2) & lead(FEATURE, n = 2)=="GENE",
          lead(GFF_START, n = 2) - GFF_END +1,
          ifelse(
            FEATURE == "TE" & CHROM==lead(CHROM, n = 3) & lead(FEATURE, n = 3)=="GENE",
            lead(GFF_START, n = 3) - GFF_END +1,
            ifelse(
              FEATURE == "TE" & CHROM==lead(CHROM, n = 4) & lead(FEATURE, n = 4)=="GENE",
              lead(GFF_START, n = 4) - GFF_END +1,
              ifelse(
                FEATURE == "TE" & CHROM==lead(CHROM, n = 5) & lead(FEATURE, n = 5)=="GENE",
                lead(GFF_START, n = 5) - GFF_END +1,
                NA
              )
            )
          )
        )
      )
    ) %>%
    filter(
      FEATURE %in% "TE"
    )
  
  df_dist_TE_GENE  %>%
    mutate(
      CLOSEST_GENE_ID = ifelse(
        PREVIOUS_GENE_DISTANCE > NEXT_GENE_DISTANCE & !is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
        NEXT_GENE_ID, 
        ifelse(
          PREVIOUS_GENE_DISTANCE < NEXT_GENE_DISTANCE & !is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
          PREVIOUS_GENE_ID, 
          ifelse(
            is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
            NEXT_GENE_ID, 
            ifelse(
              is.na(NEXT_GENE_DISTANCE) & !is.na(PREVIOUS_GENE_DISTANCE), 
              PREVIOUS_GENE_ID, 
              NA
            )
          )
        )
      ),
      CLOSEST_GENE_DISTANCE = ifelse(
        PREVIOUS_GENE_DISTANCE > NEXT_GENE_DISTANCE & !is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
        NEXT_GENE_DISTANCE, 
        ifelse(
          PREVIOUS_GENE_DISTANCE < NEXT_GENE_DISTANCE & !is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
          PREVIOUS_GENE_DISTANCE, 
          ifelse(
            is.na(PREVIOUS_GENE_DISTANCE) & !is.na(NEXT_GENE_DISTANCE), 
            NEXT_GENE_DISTANCE, 
            ifelse(
              is.na(NEXT_GENE_DISTANCE) & !is.na(PREVIOUS_GENE_DISTANCE), 
              PREVIOUS_GENE_DISTANCE, 
              NA
            )
          )
        )
      )
    ) %>%
    left_join(
      x = .,
      y = df_singleLine,
      by = c("CHROM", "GFF_START", "GFF_END", "FEATURE", "FEATURE_ID")
    ) %>%
    rename(
      TE_CHROM = CHROM, 
      TE_GFF_START = GFF_START, 
      TE_GFF_END = GFF_END,
      TE_ID = FEATURE_ID
    ) %>%
    mutate(
      PLANT_ID = NAM
    ) %>%
    select(
      PLANT_ID, TE_CHROM, TE_GFF_START, TE_GFF_END, TE_ID, TE_NAME, TE_CLASS, 
      PREVIOUS_GENE_ID, PREVIOUS_GENE_DISTANCE, NEXT_GENE_ID, NEXT_GENE_DISTANCE, 
      CLOSEST_GENE_ID, CLOSEST_GENE_DISTANCE
    ) -> df_ret
  
  return(df_ret)
  
})


df_save <- bind_rows(temp_list)

fwrite(
  df_save, 
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Distance.txt",
  col.names = T, row.names = F, sep = "\t", quote = T
)

