#========================================
# AIM
# 	Estimate the distance between intack TE and closest TE
#	  Only in class 5
#
#========================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

name_list <- fread(
  file = "name.list", 
  header = F, sep = "\t"
) %>% dplyr::select(V1) %>% pull

# import summary =======================================================
df_info <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Relationship.txt", 
  header = T, sep = "\t"
)

#glimpse(df_info)

df_info_class5 <- df_info %>%
  filter(
    RELATION_OTHERTE_TE %in% "CLASS5_NotOverlapWithOTHERTE"
  )

# per-plant ======================================================

temp_list <- lapply(name_list, FUN = function(NAM){
  
  # import homo and intact TE
  df_OTHERTE_complete <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/", NAM, ".TE.InChr.IntactSingleLineAndHomoTE.sort.gff3"), 
    header = F, sep = "\t"
  ) %>%
    mutate(
      OTHERTE_ID = str_extract(V9, pattern = "ID=[A-Za-z0-9_]{1,}") %>%
        str_remove(., pattern = "ID="),
      OTHERTE_NAME = str_extract(V9, pattern = "Name=[A-Za-z0-9_]{1,}") %>%
        str_remove(., pattern = "Name="),
      OTHERTE_CLASS = str_extract(V9, pattern = "Classification=[A-Za-z0-9_/]{1,}") %>%
        str_remove(., pattern = "Classification="), 
      FEATURE = "OTHERTE"
    ) %>%
    rename(
      CHROM = V1, 
      GFF_START = V4, 
      GFF_END = V5, 
      FEATURE_ID = OTHERTE_ID
    ) 
  df_OTHERTE <- df_OTHERTE_complete %>%
    select(
      CHROM, GFF_START, GFF_END, FEATURE_ID, FEATURE
    )
  
  # class 5 in from Summary file
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
  
  # remove the same record ======================
  df_OTHERTE <- df_OTHERTE %>%
    filter(
      !(FEATURE_ID %in% df_singleLine$FEATURE_ID)
    ) 
  
  df_dist_TE_GENE <- df_singleLine %>%
    select(
      CHROM, GFF_START, GFF_END, FEATURE, FEATURE_ID
    ) %>%
    rbind(., df_OTHERTE) %>%
    arrange(
      CHROM, GFF_START
    ) %>%
    mutate(
      PREVIOUS_OTHERTE_ID = ifelse(
        FEATURE == "TE" & CHROM==lag(CHROM, n = 1),
        lag(FEATURE_ID, n = 1),
        ifelse(
          FEATURE == "TE" & CHROM==lag(CHROM, n = 2),
          lag(FEATURE_ID, n = 2),
          ifelse(
            FEATURE == "TE" & CHROM==lag(CHROM, n = 3),
            lag(FEATURE_ID, n = 3),
            ifelse(
              FEATURE == "TE" & CHROM==lag(CHROM, n = 4),
              lag(FEATURE_ID, n = 4),
              ifelse(
                FEATURE == "TE" & CHROM==lag(CHROM, n = 5),
                lag(FEATURE_ID, n = 5),
                NA
              )
            )
          )
        )
      ),
      PREVIOUS_OTHERTE_DISTANCE = ifelse(
        FEATURE == "TE" & CHROM==lag(CHROM, n = 1),
        GFF_START - lag(GFF_END, n = 1) +1,
        ifelse(
          FEATURE == "TE" & CHROM==lag(CHROM, n = 2),
          GFF_START - lag(GFF_END, n = 2) +1,
          ifelse(
            FEATURE == "TE" & CHROM==lag(CHROM, n = 3),
            GFF_START - lag(GFF_END, n = 3) +1,
            ifelse(
              FEATURE == "TE" & CHROM==lag(CHROM, n = 4),
              GFF_START - lag(GFF_END, n = 4) +1,
              ifelse(
                FEATURE == "TE" & CHROM==lag(CHROM, n = 5),
                GFF_START - lag(GFF_END, n = 5) +1,
                NA
              )
            )
          )
        )
      ),
      NEXT_OTHERTE_ID = ifelse(
        FEATURE == "TE" & CHROM==lead(CHROM, n = 1),
        lead(FEATURE_ID, n = 1),
        ifelse(
          FEATURE == "TE" & CHROM==lead(CHROM, n = 2),
          lead(FEATURE_ID, n = 2),
          ifelse(
            FEATURE == "TE" & CHROM==lead(CHROM, n = 3),
            lead(FEATURE_ID, n = 3),
            ifelse(
              FEATURE == "TE" & CHROM==lead(CHROM, n = 4),
              lead(FEATURE_ID, n = 4),
              ifelse(
                FEATURE == "TE" & CHROM==lead(CHROM, n = 5),
                lead(FEATURE_ID, n = 5),
                NA
              )
            )
          )
        )
      ),
      NEXT_OTHERTE_DISTANCE = ifelse(
        FEATURE == "TE" & CHROM==lead(CHROM, n = 1),
        lead(GFF_START, n = 1) - GFF_END +1,
        ifelse(
          FEATURE == "TE" & CHROM==lead(CHROM, n = 2),
          lead(GFF_START, n = 2) - GFF_END +1,
          ifelse(
            FEATURE == "TE" & CHROM==lead(CHROM, n = 3),
            lead(GFF_START, n = 3) - GFF_END +1,
            ifelse(
              FEATURE == "TE" & CHROM==lead(CHROM, n = 4),
              lead(GFF_START, n = 4) - GFF_END +1,
              ifelse(
                FEATURE == "TE" & CHROM==lead(CHROM, n = 5),
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
  
  df_dist_TE_GENE %>%
    mutate(
      ## IF distance is the same, I report the NEXT feature
      CLOSEST_OTHERTE_ID = ifelse(
        PREVIOUS_OTHERTE_DISTANCE >= NEXT_OTHERTE_DISTANCE & !is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
        NEXT_OTHERTE_ID,
        ifelse(
          PREVIOUS_OTHERTE_DISTANCE < NEXT_OTHERTE_DISTANCE & !is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
          PREVIOUS_OTHERTE_ID,
          ifelse(
            is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
            NEXT_OTHERTE_ID,
            ifelse(
              is.na(NEXT_OTHERTE_DISTANCE) & !is.na(PREVIOUS_OTHERTE_DISTANCE),
              PREVIOUS_OTHERTE_ID,
              NA
            )
          )
        )
      ),
      CLOSEST_OTHERTE_DISTANCE = ifelse(
        PREVIOUS_OTHERTE_DISTANCE >= NEXT_OTHERTE_DISTANCE & !is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
        NEXT_OTHERTE_DISTANCE,
        ifelse(
          PREVIOUS_OTHERTE_DISTANCE < NEXT_OTHERTE_DISTANCE & !is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
          PREVIOUS_OTHERTE_DISTANCE,
          ifelse(
            is.na(PREVIOUS_OTHERTE_DISTANCE) & !is.na(NEXT_OTHERTE_DISTANCE),
            NEXT_OTHERTE_DISTANCE,
            ifelse(
              is.na(NEXT_OTHERTE_DISTANCE) & !is.na(PREVIOUS_OTHERTE_DISTANCE),
              PREVIOUS_OTHERTE_DISTANCE,
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
      PREVIOUS_OTHERTE_ID, PREVIOUS_OTHERTE_DISTANCE, NEXT_OTHERTE_ID, NEXT_OTHERTE_DISTANCE,
      CLOSEST_OTHERTE_ID, CLOSEST_OTHERTE_DISTANCE
    ) -> df_ret
  
})

df_save <- bind_rows(temp_list)

# glimpse(df_save)

fwrite(
  df_save, 
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Distance.txt", 
  col.names = T, sep = "\t", quote = T, row.names = F
)
