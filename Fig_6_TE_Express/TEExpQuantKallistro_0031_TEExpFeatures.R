#====================================
# AIM
#	Using party::cforest to identify importance features to predict TE family expression
#
#====================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
#library(party)

setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)


print("# HINT: dplyr::where is masked by party::where")

print("# import ID table ========================================")
name_list <- fread(
  file = "/Data/Fig_4_TE_Feature/name.list",
  header = F, sep = "\t"
) %>% select(V1) %>% pull

df_plant_name <- fread(
  file = "/Data/Fig_4_TE_Feature/Plant_IDNAME.txt",
  header = T, sep = "\t"
)

df_plant_species <- fread(
  file = "/Data/Fig_4_TE_Feature/Plant_Species.txt",
  header = T, sep = "\t"
)

print("# import prevalence in Z score format =====================")
cat_list <- c("DNA_Helitron", "LTR_Copia", "LTR_Gypsy", "LTR_unknown", "TIR_DTA", "TIR_DTC", "TIR_DTH", "TIR_DTM", "TIR_DTT")

temp_list <- lapply(cat_list, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/Heatmap_TE_Structural/heatmap_TEStructural_AddShareUniq_Zscale_", NAM, ".txt"),
    header = T, sep = "\t"
  )
  
  OUT <- str_replace(NAM, "_", "/")
  
  df_ret <- df_temp %>%
    mutate(
      TE_CLASS = OUT
    )
  
  return(df_ret)
  
})

df_z_score <- bind_rows(temp_list)
rm(temp_list)

df_temp_z <- df_z_score %>%
  as.data.table() %>%
  melt.data.table(
    id.vars = c("TE_CLASS", "TE_NAME"),
    variable.name = "PLANT_NAME",
    variable.factor = F, 
    value.name = "PREVALENCE"
  ) %>%
  left_join(
    x = .,
    y = df_plant_species,
    by = "PLANT_NAME"
  ) %>%
  glimpse()

# the original name of this variable is df_temp_z_all
df_fam_prevelance_z <- df_temp_z %>%
  summarise(
    .by = c(TE_CLASS, TE_NAME), 
    PREVALENCE_AVE_ALL = mean(PREVALENCE, na.rm = T),
    PREVALENCE_SD_ALL = sd(PREVALENCE, na.rm = T)
  ) %>%
  glimpse()


rm(temp_list)
rm(df_temp_z)
rm(df_temp_z_all)
rm(df_temp_z_spec)

print("# import share and uniq info=======================")
df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt",
  header = T, sep = "\t"
)
colnames(df_share) <- c("TE_CLASS", "TE_NAME", "SHARE_UNIQUE")

print("# import LTR intergation time  ======================")
df_age <- fread(
  file = "/Data/Supplement_Data/Insert_Age_Brassica_AA_CC_TEIDChange.txt",
  header = T, sep = "\t"
)

df_fam_age <- df_age[
  , .(mean(INSERT_AGE, na.rm = T), sd(INSERT_AGE, na.rm = T)), by = .(TE_CLASS, TE_NAME)
]
colnames(df_fam_age) <- c("TE_CLASS", "TE_NAME", "LTR_AGE_AVE", "LTR_AGE_SD")


print("# import count and length ============================")
df_count <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt.gz", 
  header = T, sep = "\t"
)

df_fam_count_len <- df_count[
  TE_METHOD == "structural"
][
  , "LEN" := GFF_END - GFF_START + 1
][
  , .(.N, mean(LEN, na.rm = T), sd(LEN, na.rm = T)), by = .(TE_CLASS, TE_NAME)
]
colnames(df_fam_count_len) <- c("TE_CLASS", "TE_NAME", "COUNT", "LEN_AVE", "LEN_SD")

glimpse(df_fam_count_len)

print("# import LTR GC content ===============================")
df_gc_LTR <- fread(
  file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_GC.txt.gz",
  header = T, sep = "\t"
)

print("# import LTR methylation data ==========================")
df_meth_LTR <- fread(
  file = "/Data/Fig_4_TE_Feature/Combine_LongTerminalRepeat_MethylSite_PerInsert.txt.gz",
  header = T, sep = "\t"
)

print("# import type II GC content, full length ===================================")
df_gc_TIR <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_GC.txt.gz",
  header = T, sep = "\t"
)

print("# Combine GC Content =============================")
df_gc <- df_gc_TIR %>%
  select(
    PLANT_ID, CHROM, TE_ID, TE_NAME, TE_CLASS, GC, LEN
  ) %>%
  rbind(., df_gc_LTR) %>%
  as.data.table()

glimpse(df_gc)
df_fam_gc <- df_gc[
  , .(mean(GC*100, na.rm = T), sd(GC*100, na.rm = T)), by = .(TE_CLASS, TE_NAME)
]
colnames(df_fam_gc) <- c("TE_CLASS", "TE_NAME", "GC_AVE", "GC_SD")
glimpse(df_fam_gc)

df_fam_LTR_LEN <- df_gc %>%
  filter(
    str_detect(TE_CLASS, "^LTR")
  ) %>%
  summarise(
    .by = c(TE_CLASS, TE_NAME),
    LTR_LEN_AVE = mean(LEN, na.rm = T),
    LTR_LEN_SD = sd(LEN, na.rm = T)
  ) 



print("# import ClassII methylation data, full length ==============================")
df_meth_TIR <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_methyl_raw.txt.gz",
  header = T, sep = "\t"
) %>%
  select(-GC)

print("# Combine Methylation info =====================================")
colnames(df_meth_LTR) <- c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS", "METHYL", "ACCUM_METHYL_LEN", "ACCUM_FLANK_LEN", "ACCUM_METHYL_PERC")
colnames(df_meth_TIR) <- c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS", "METHYL", "ACCUM_METHYL_LEN", "ACCUM_FLANK_LEN", "ACCUM_METHYL_PERC")

df_meth <- rbind(df_meth_LTR, df_meth_TIR) %>% as.data.table()

glimpse(df_meth)

df_meth_temp <- df_meth[
  ,.(mean(ACCUM_METHYL_PERC, na.rm = T), sd(ACCUM_METHYL_PERC, na.rm = T)), by = .(TE_CLASS, TE_NAME, METHYL)
] 


temp_list <- lapply(c("CG", "CHG", "CHH"), FUN = function(NAM){
  
  df_ret <- df_meth_temp %>%
    filter(
      METHYL %in% NAM
    ) %>%
    select(
      -METHYL
    )
  
  colnames(df_ret) <- c("TE_CLASS", "TE_NAME", 
                        paste0("METHYL_", NAM, "_AVE"), 
                        paste0("METHYL_", NAM, "_SD"))
  return(df_ret)
})

df_fam_methyl <- left_join(
  x = temp_list[[1]],
  y = temp_list[[2]],
  by = c("TE_CLASS", "TE_NAME")
) %>%
  left_join(
    x = .,
    y = temp_list[[3]], 
    by = c("TE_CLASS", "TE_NAME")
  )

rm(temp_list)
rm(df_meth_temp)

print("# import the relaionship to gene ==============================")
df_gene_relation <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt.gz",
  header = T, sep = "\t"
)

df_temp <- df_gene_relation[
  , .N, by = .(TE_CLASS, TE_NAME, RELATION_GENE_TE)
][
  , "PERC" := (N/sum(N) * 100), by = .(TE_CLASS, TE_NAME)
]

df_fam_TE_GENE <- df_temp %>%
  dcast.data.table(
    TE_CLASS + TE_NAME ~ RELATION_GENE_TE, value.var = "PERC", fill = 0
  ) %>%
  rename_with(
    .fn = ~paste0("GENE_", .x), 
    .cols = starts_with("CLASS")
  )

rm(df_temp)


print("# import the distnace to gene in class 5 ========================")
df_gene_dist <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Distance.txt.gz",
  header = T,  sep = "\t"
)

df_fam_gene_dist <- df_gene_dist[
  ,.(mean(CLOSEST_GENE_DISTANCE, na.rm = T), 
    sd(CLOSEST_GENE_DISTANCE, na.rm = T)),
  , by = .(TE_CLASS, TE_NAME)
]

colnames(df_fam_gene_dist) <- c("TE_CLASS", "TE_NAME", 
                                "CLOSEST_GENE_DISTANCE_AVE",
                                "CLOSEST_GENE_DISTANCE_SD")

df_fam_gene_dist <- df_fam_gene_dist %>%
  mutate(
    CLOSEST_GENE_DISTANCE_AVE = ifelse(
      is.na(CLOSEST_GENE_DISTANCE_AVE), 
      0, CLOSEST_GENE_DISTANCE_AVE
    ),
    CLOSEST_GENE_DISTANCE_SD = ifelse(
      is.na(CLOSEST_GENE_DISTANCE_SD),
      0, CLOSEST_GENE_DISTANCE_SD
    )
  ) 


print("# import the relationship to TE ================================")
df_TE_relation <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Relationship.txt.gz",
  header = T,  sep = "\t"
)

df_temp <- df_TE_relation[
  , .N, by = .(TE_CLASS, TE_NAME, RELATION_OTHERTE_TE)
][
  , "PERC" := N/sum(N)*100 , by = .(TE_CLASS, TE_NAME)
]

df_fam_TE_TE <- df_temp %>%
  dcast.data.table(
    TE_CLASS + TE_NAME ~ RELATION_OTHERTE_TE, value.var = "PERC", fill = 0
  ) %>%
  rename_with(
    .fn = ~ paste0("TE_", .x), 
    .cols = starts_with("CLASS")
  )

rm(df_temp)
  
print("# import the distance to TE in class 5 =========================")
df_TE_dist <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Distance.txt.gz",
  header = T,  sep = "\t"
)

df_fam_TE_dist <- df_TE_dist[
  , .(mean(CLOSEST_OTHERTE_DISTANCE, na.rm = T),
      sd(CLOSEST_OTHERTE_DISTANCE, na.rm = T))
  , by = .(TE_CLASS, TE_NAME)
]

colnames(df_fam_TE_dist) <- c("TE_CLASS", "TE_NAME", "CLOSEST_OTHERTE_DISTANCE_AVE",
                              "CLOSEST_OTHERTE_DISTANCE_SD")

df_fam_TE_dist <- df_fam_TE_dist %>%
  mutate(
    CLOSEST_OTHERTE_DISTANCE_AVE = ifelse(
      is.na(CLOSEST_OTHERTE_DISTANCE_AVE),
      0, CLOSEST_OTHERTE_DISTANCE_AVE
    ),
    CLOSEST_OTHERTE_DISTANCE_SD = ifelse(
      is.na(CLOSEST_OTHERTE_DISTANCE_SD),
      0, CLOSEST_OTHERTE_DISTANCE_SD
    )
  )


print("# import outside 2000 bp GC content =================================")
df_out2k_GC <- fread(
  file = "/Data/Fig_4_TE_Feature/Both_sidesTE_2K_GC.txt.gz",
  header = T, sep = "\t"
)

df_fam_out2k_GC <- df_out2k_GC[
  , .(mean(GC_content*100, na.rm = T), sd(GC_content*100, na.rm = T))
  , by = .(TE_CLASS, TE_NAME)
]

colnames(df_fam_out2k_GC) <- c("TE_CLASS", "TE_NAME", "OUTER_GC_AVE", 
                               "OUTER_GC_SD")

df_fam_out2k_GC <- df_fam_out2k_GC %>%
  mutate(
    OUTER_GC_AVE = ifelse(
      is.na(OUTER_GC_AVE),
      0, OUTER_GC_AVE
    ),
    OUTER_GC_SD = ifelse(
      is.na(OUTER_GC_SD),
      0, OUTER_GC_SD
    )
  )

# Summarize to family specific ===================================================================
print("# Summarize to family specific ======================================")

df_fam_feature <- left_join(
  x = df_fam_count_len, 
  y = df_fam_age,
  by = c("TE_CLASS", "TE_NAME")
) %>%
  left_join(
    x = .,
    y = df_fam_LTR_LEN,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_prevelance_z,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_gc,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_methyl,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_TE_GENE,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_gene_dist,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_TE_TE,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_TE_dist,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  left_join(
    x = .,
    y = df_fam_out2k_GC,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  glimpse()

print("# make sure there is no NA in data frame, especially I use sd as feature")
df_fam_feature[is.na(df_fam_feature)] <- 0

fwrite(
  df_fam_feature,
  file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix.txt",
  col.names = T, sep = "\t"
)

df_fam_feature <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix.txt.gz",
  header = T, sep = "\t"
)


print("# read expression data =============================")
df_exp_TEFAM <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_H37.txt.gz",
  header = T, sep = "\t"
)

df_thres <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_threshold_H37.txt.gz", 
  header = T, sep = "\t"
)

glimpse(df_exp_TEFAM)
glimpse(df_thres)

print("# Using the highest expression level from 5 conditions to represent TE family expression ======")
df_fam_exp_rank <- df_exp_TEFAM %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  arrange(
    .by_group = T,
    desc(NORMALIZE_TPM)
  ) %>%
  slice_head(
    n = 1
  ) %>%
  select(
    SAMPLE, TE_CLASS, TE_NAME, NORMALIZE_TPM
  ) %>%
  group_by(
    TE_CLASS
  ) %>%
  arrange(
    NORMALIZE_TPM, 
    .by_group = T
  ) %>%
  mutate(
    NORMALIZE_EXP = (NORMALIZE_TPM - min(NORMALIZE_TPM))/(max(NORMALIZE_TPM) - min(NORMALIZE_TPM)), 
    NORMALIZE_EXP_RANK = rank(NORMALIZE_TPM, ties.method = "min")
  ) %>%
  ungroup() %>%
  select(
    -SAMPLE
  )

glimpse(df_fam_exp_rank)

print("# Using YES or NO (through threshold) to represent expression ============")
df_fam_exp <- left_join(
  x = df_exp_TEFAM,
  y = df_thres,
  by = "TE_CLASS"
) %>%
  mutate(
    TE_FAM_EXP = ifelse(
      NORMALIZE_TPM >= T95,
      "YES", "NO"
    ),
    TE_FAM_EXP = factor(
      TE_FAM_EXP, 
      levels = c("YES", "NO"),
      labels = c("YES", "NO")
    )
  ) %>%
  summarise(
    .by = c(TE_CLASS, TE_NAME, TE_FAM_EXP),
    COUNT = n()
  ) %>%
  slice_head(
    by = c(TE_CLASS, TE_NAME),
    n = 1
  ) %>%
  select(
    TE_CLASS, TE_NAME, TE_FAM_EXP
  ) 

glimpse(df_fam_exp)


print("# put expression and features data together =====================")

df_for_tree <- left_join(
  x = df_fam_exp,
  y = df_fam_feature,
  by = c("TE_CLASS", "TE_NAME")
) %>%
  glimpse()

fwrite(
  df_for_tree,
  file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix_AddExp.txt",
  col.names = T, sep = "\t"
)

print("# saving rank result ----------------------------")
left_join(
  x = df_fam_exp_rank,
  y = df_fam_feature,
  by = c("TE_CLASS", "TE_NAME")
) %>%
  fwrite(
    file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix_AddExpRank.txt",
    col.names = T, sep = "\t"
  )
