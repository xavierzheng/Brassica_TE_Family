#------------------------------------
# AIM
#	Focus on the difference between AA and CC
#
#===================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(corrplot, quietly = F))
suppressPackageStartupMessages(library(broom))

setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

print("# import ID table ========================================")
name_list <- fread(
  file = "/Data/Fig_4_TE_Feature/name.list", 
  header = F, sep = "\t"
) %>% select(V1) %>% pull

df_plant_name <- fread(
  file = "/Data/Fig_4_TE_Feature/Plant_IDNAME.txt", 
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

print("# import share and uniq info=======================")
df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt",
  header = T, sep = "\t"
)
colnames(df_share) <- c("TE_CLASS", "TE_NAME", "Share_or_Unique")

print("# import LTR intergation time  ======================")
df_age <- fread(
  file = "/Data/Supplement_Data/Insert_Age_Brassica_AA_CC_TEIDChange.txt", 
  header = T, sep = "\t"
)

print("# import count and length ============================")
df_count <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt.gz",
  header = T, sep = "\t"
)

print("# import structural TE count and length ========================")
df_count_intact <- fread(
  file = "/Data/Fig_2_Distribution/StructuralTE_TotalCount_FamilyLevel.txt", 
  header = T, sep = "\t"
)

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
  rbind(., df_gc_LTR)

print("# import ClassII methylation data, full length ==============================")
df_meth_TIR <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_TypeII_FullLength_methyl_raw.txt.gz",
  header = T, sep = "\t"
) %>%
  select(-GC)


print("# Combine Methylation info =====================================")
colnames(df_meth_LTR) <- c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS", "METHYL", "ACCUM_METHYL_LEN", "ACCUM_FLANK_LEN", "ACCUM_METHYL_PERC")
colnames(df_meth_TIR) <- c("PLANT_ID", "CHROM", "TE_ID", "TE_NAME", "TE_CLASS", "METHYL", "ACCUM_METHYL_LEN", "ACCUM_FLANK_LEN", "ACCUM_METHYL_PERC")

df_meth <- rbind(df_meth_LTR, df_meth_TIR)

print("# import the relaionship to gene ==============================")
df_gene_relation <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt.gz", 
  header = T, sep = "\t"
)

print("# import the distnace to gene in class 5 ========================")
df_gene_dist <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Distance.txt.gz", 
  header = T,  sep = "\t"
)

print("# import the relationship to TE ================================")
df_TE_relation <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Relationship.txt.gz", 
  header = T,  sep = "\t"
)

print("# import the distance to TE in class 5 =========================")
df_TE_dist <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_IntactTE_TE_Distance.txt.gz", 
  header = T,  sep = "\t"
)

print("# import outside 2000 bp GC content =================================")
df_out2k_GC <- fread(
  file = "/Data/Fig_4_TE_Feature/Both_sidesTE_2K_GC.txt.gz", 
  header = T, sep = "\t"
)

print("# upset table for synteny TE ==================================")
df_synteny_AA_TEfamily <- read.table(
  file = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_AA.txt.gz", 
  header = T, sep = "\t", row.names = 1
)

df_synteny_CC_TEfamily <- read.table(
  file = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_CC.txt.gz", 
  header = T, sep = "\t", row.names = 1
)

df_synteny_AA_TEclass <- read.table(
  file = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_AA.txt.gz", 
  header = T, sep = "\t", row.names = 1
)

df_synteny_CC_TEclass <- read.table(
  file = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_CC.txt.gz", 
  header = T, sep = "\t", row.names = 1
)

df_species <- fread(
  file = "/Data/Fig_4_TE_Feature/Plant_Species.txt", 
  header = T, sep = "\t"
)


# define function ========================================================================
print("# Define function ================================")
estimate_diff <- function(DF, COL, OUTPUT_COL){
  df_temp_raw <- DF %>%
    left_join(
      x = ., 
      y = df_species, 
      by = "PLANT_ID"
    ) %>%
    left_join(
      x = .,
      y = df_share, 
      by = c("TE_CLASS", "TE_NAME")
    ) %>%
    mutate(
      PLANT_SPECIES = factor(
        PLANT_SPECIES, 
        levels = c("B_rapa", "B_oleracea"), 
        labels = c("B_rapa", "B_oleracea")
      )
    )
  
  # for shared 
  df_temp <- df_temp_raw %>%
    filter(
      Share_or_Unique %in% "Share_in_AA_CC"
    ) 
  
  df_key <- df_temp %>%
    filter(
      !is.na(get(COL))
    ) %>%
    group_by(
      TE_CLASS, TE_NAME, PLANT_SPECIES
    ) %>%
    summarise(
      COUNT = n()
    ) %>%
    ungroup() %>%
    group_by(
      TE_CLASS, TE_NAME
    ) %>%
    summarise(
      COUNT = n()
    ) %>%
    ungroup() %>%
    filter(
      COUNT >=2
    ) %>%
    select(-COUNT) %>%
    mutate(
      SELECT = "KEY"
    )
  
  df_temp_pval <- df_temp %>%
    left_join(
      x = ., 
      y = df_key, 
      by = c("TE_CLASS", "TE_NAME")
    ) %>%
    filter(
      SELECT %in% "KEY"
    ) %>%
    nest_by(
      TE_CLASS, TE_NAME
    ) %>%
    summarise(
      PVAL = wilcox.test(get(COL) ~ PLANT_SPECIES, data = data, exact = F)$p.value
    ) %>%
    ungroup() 
  
  df_temp %>%
    group_by(
      TE_CLASS, TE_NAME, PLANT_SPECIES
    ) %>%
    summarise(
      MED = median(get(COL), na.rm = T)
    ) %>%
    ungroup() %>%
    reshape2::dcast(
      TE_CLASS + TE_NAME ~ PLANT_SPECIES, 
      value.var = "MED"
    ) %>%
    mutate(
      DIFF_AA_CC = B_rapa - B_oleracea
    ) %>%
    left_join(
      x = .,
      y = df_temp_pval, 
      by = c("TE_CLASS", "TE_NAME")
    ) %>%
    mutate(
      FDR = p.adjust(PVAL, method = "fdr")
    ) %>%
    mutate(
      CAT = ifelse(
        FDR < 0.01 & DIFF_AA_CC > 0, 
        "AA > CC", 
        ifelse(
          FDR < 0.01 & DIFF_AA_CC < 0, 
          "CC > AA", 
          "No different"
        )
      ), 
      COMPARE = OUTPUT_COL
    ) -> df_ret
  
  return(df_ret)
}


plot_diff <- function(DF, COL, Y_AXIS){
  
  # label
  df_temp_raw <- DF %>%
    left_join(
      x = ., 
      y = df_species, 
      by = "PLANT_ID"
    ) %>%
    left_join(
      x = .,
      y = df_share, 
      by = c("TE_CLASS", "TE_NAME")
    ) %>%
    mutate(
      PLANT_SPECIES = factor(
        PLANT_SPECIES, 
        levels = c("B_rapa", "B_oleracea"), 
        labels = c("B_rapa", "B_oleracea")
      )
    )
  
  # keep shared
  df_temp <- df_temp_raw %>%
    filter(
      Share_or_Unique %in% "Share_in_AA_CC"
    ) 
  
  # only use Copia0002, Copia0011, Copia0006, the most prevelant Copia family
  df_key <- data.frame(
    TE_CLASS = rep("LTR/Copia", 4), 
    TE_NAME = c("Copia0002", "Copia0011", "Copia0006", "Copia0035"), 
    SELECT = rep("KEY", 4)
  )
  
  # drawing
  df_temp %>%
    left_join(
      x = ., 
      y = df_key, 
      by = c("TE_CLASS", "TE_NAME")
    ) %>%
    filter(
      SELECT %in% "KEY"
    ) %>%
    mutate(
      TE_NAME = factor(
        TE_NAME, 
        levels = c("Copia0002", "Copia0011", "Copia0006", "Copia0035"), 
        labels = c("Copia0002", "Copia0011", "Copia0006", "Copia0035")
      )
    ) %>%
    ggplot(
      aes(
        x = TE_NAME, 
        y = get(COL), 
        fill = PLANT_SPECIES
      )
    )+
    geom_boxplot(
      outlier.shape = NA,
      linewidth = 0.1
    )+
    scale_fill_manual(
      breaks = c("B_rapa", "B_oleracea"), 
      values = c("#998EC3", "#F1A340")
    )+
    theme_cowplot()+
    theme(
      axis.title.x = element_blank(), 
      axis.title.y = element_text(
        size = 6
      ),
      axis.text.x = element_text(
        size = 5, 
        angle = 90, 
        hjust = 1, 
        vjust = 0.5
      ), 
      axis.text.y = element_text(
        size = 5
      ), 
      legend.position = "none", 
      plot.margin = margin(t = 1, r = 0, b = 1, l = 2, unit = "pt")
    )+
    labs(
      y = Y_AXIS
    ) -> plot_p
  
  return(plot_p)
  
}

# For real data =================================================================
print("whole length------------------------------------------------------------")
glimpse(df_count)
rm(df_count_temp)
df_count_temp <- df_count %>%
  mutate(
    TE_LEN = GFF_END - GFF_START +1
  ) %>%
  filter(
    TE_METHOD %in% "structural"
  )
glimpse(df_count_temp)

df_result_whole_length <- estimate_diff(df_count_temp, "TE_LEN", "Whole_length")

plot_whole_length <- plot_diff(df_count_temp, "TE_LEN", "Whole_length")

print("# Count ----------------------------------------------------------------")
glimpse(df_count)
rm(df_count_temp)
df_count_temp <- df_count %>%
  filter(
    TE_METHOD %in% "structural"
  ) %>%
  summarise(
    .by = c(PLANT_ID, TE_CLASS, TE_NAME),
    NO = n()
  ) %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    TE_CLASS + TE_NAME ~ PLANT_ID, value.var = "NO", fill = 0
  ) %>%
  data.table::melt.data.table(
    id.vars = c("TE_CLASS", "TE_NAME"), 
    variable.name = "PLANT_ID", 
    variable.factor = F,
    value.name = "NO"
  ) 
glimpse(df_count_temp)

print("##### ONLY TRY spearman correlation --------------------------------------------------")
df_count_temp2 <- df_count_temp %>%
  left_join(
    x = .,
    y = df_species,
    by = "PLANT_ID"
  ) %>%
  left_join(
    x = .,
    y = df_share, 
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  filter(
    Share_or_Unique %in% "Share_in_AA_CC"
  ) 

df_count_temp_p <- df_count_temp2 %>%
  mutate(
    CAT = ifelse(
      PLANT_SPECIES %in% "B_rapa", 
      0, 1
    )
  ) %>%
  nest_by(
    TE_CLASS, TE_NAME
  ) %>%
  summarise(
    PVAL = cor.test( ~ CAT + NO, data = data, method = "spearman", exact = F)$p.value
  ) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(PVAL, method = "fdr")
  ) 

df_result_count <- df_count_temp2 %>%
  summarise(
    .by = c(PLANT_SPECIES, TE_CLASS, TE_NAME), 
    MED = median(NO, na.rm = T)
  ) %>%
  mutate(
    PLANT_SPECIES = factor(
      PLANT_SPECIES,
      levels = c("B_rapa", "B_oleracea"), 
      labels = c("B_rapa", "B_oleracea")
    )
  ) %>%
  reshape2::dcast(
    TE_CLASS + TE_NAME ~ PLANT_SPECIES, value.var = "MED", fill = 0
  ) %>%
  mutate(
    DIFF_AA_CC = B_rapa - B_oleracea
  ) %>%
  left_join(
    x = ., 
    y = df_count_temp_p, 
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    CAT = ifelse(
      FDR < 0.01 & DIFF_AA_CC > 0, 
      "AA > CC", 
      ifelse(
        FDR < 0.01 & DIFF_AA_CC < 0, 
        "CC > AA", 
        "No different"
      )
    ), 
    COMPARE = "COUNT"
  ) 

glimpse(df_result_count)
print("##### ABOVE is strange, but I used @@@@@@@@@@@@@@@@@@@@ ----------------------------")

print("# LTR length ============================================================")  
glimpse(df_gc)
rm(df_temp)

df_gc_temp <- df_gc %>%
  filter(
    str_detect(TE_CLASS, "^LTR")
  )
df_result_LTR_length <- estimate_diff(df_gc_temp, "LEN", "LTR_length")
plot_LTR_length <- plot_diff(df_gc_temp, "LEN", "LTR_length")

print("GC is different =========================================================")
glimpse(df_gc)
#rm(df_temp)
df_result_GC <- estimate_diff(df_gc, "GC", "GC_content")
plot_GC <- plot_diff(df_gc, "GC", "GC_content")


print("# Insertion age is different===================================================")
glimpse(df_age)
df_age_temp <- df_age %>%
  mutate(
    FAKE = str_split_fixed(CHROM, "_", 2)[,1],
    PLANT_ID = ifelse(
      str_detect(FAKE, "^D"), 
      paste0(FAKE, ".AA"), 
      paste0(FAKE, ".CC")
    )
  ) %>%
  select(-FAKE) 

df_result_LTR_AGE <- estimate_diff(df_age_temp, "INSERT_AGE", "LTR_INSERT_AGE")
plot_LTR_AGE <- plot_diff(df_age_temp, "INSERT_AGE", "LTR_INSERT_AGE")

print("# methyl motif ===========================================================")
glimpse(df_meth)
rm(df_meth_temp)
df_meth_temp <- df_meth %>%
  filter(
    METHYL %in% "CG"
  )
df_result_meth_CG <- estimate_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CG_motif")
plot_meth_CG <- plot_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CG_motif")


rm(df_meth_temp)
df_meth_temp <- df_meth %>%
  filter(
    METHYL %in% "CHG"
  )
df_result_meth_CHG <- estimate_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CHG_motif")
plot_meth_CHG <- plot_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CHG_motif")

rm(df_meth_temp)
df_meth_temp <- df_meth %>%
  filter(
    METHYL %in% "CHH"
  )
df_result_meth_CHH <- estimate_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CHH_motif")
plot_meth_CHH <- plot_diff(df_meth_temp, "ACCUM_METHYL_PERC", "CHH_motif")

print("# TE-gene relationship ============================================================")
glimpse(df_gene_relation)

rm(df_temp)
df_temp <- df_gene_relation %>%
  left_join(
    x = .,
    y = df_species, 
    by = "PLANT_ID"
  ) %>%
  left_join(
    x = ., 
    y = df_share,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  filter(
    Share_or_Unique %in% "Share_in_AA_CC"
  ) %>%
  filter(
    !is.na(RELATION_GENE_TE)
  ) %>%
  group_by(
    TE_CLASS, TE_NAME, RELATION_GENE_TE, PLANT_SPECIES
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  reshape2::dcast(
    TE_CLASS + TE_NAME + PLANT_SPECIES ~ RELATION_GENE_TE, 
    value.var = "COUNT", fill = 0
  ) 

temp_list <- lapply(unique(df_temp$TE_NAME), FUN = function(NAM){
  
  df_chi <- df_temp %>%
    filter(
      TE_NAME %in% NAM
    ) %>%
    select(
      PLANT_SPECIES, starts_with("CLASS")
    ) %>%
    column_to_rownames(
      var = "PLANT_SPECIES"
    )
  set.seed(1234)
  df_P <- chisq.test(df_chi, simulate.p.value = T, B = 2000)$p.value
  df_ret <- data.frame(
    TE_NAME = NAM, 
    PVAL = df_P
  )
  return(df_ret)
  
})
  
bind_rows(temp_list) %>%
  mutate(
    FDR = p.adjust(PVAL, method = "fdr"), 
    CAT = ifelse(
      FDR < 0.01, 
      "Differnt", 
      "Similar"
    ), 
    COMPARE = "TE_GENE_RELATION"
  ) %>%
  left_join(
    x = df_temp,
    y = ., 
    by = "TE_NAME"
  ) -> df_result_TE_gene
glimpse(df_result_TE_gene)


print("# TE-gene distance ==================================================================")
glimpse(df_gene_dist)
df_result_gene_dist <- estimate_diff(df_gene_dist, "CLOSEST_GENE_DISTANCE", "CLOSEST_GENE_DISTANCE")
plot_gene_dist <- plot_diff(df_gene_dist, "CLOSEST_GENE_DISTANCE", "CLOSEST_GENE_DISTANCE")


print(" TE-TE relationship ===================================================================")
glimpse(df_TE_relation)

rm(df_temp)
df_temp <- df_TE_relation %>%
  left_join(
    x = .,
    y = df_species, 
    by = "PLANT_ID"
  ) %>%
  left_join(
    x = ., 
    y = df_share,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  filter(
    Share_or_Unique %in% "Share_in_AA_CC"
  ) %>%
  filter(
    !is.na(RELATION_OTHERTE_TE)
  ) %>%
  group_by(
    TE_CLASS, TE_NAME, RELATION_OTHERTE_TE, PLANT_SPECIES
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  reshape2::dcast(
    TE_CLASS + TE_NAME + PLANT_SPECIES ~ RELATION_OTHERTE_TE, 
    value.var = "COUNT", fill = 0
  ) 

temp_list <- lapply(unique(df_temp$TE_NAME), FUN = function(NAM){
  
  df_chi <- df_temp %>%
    filter(
      TE_NAME %in% NAM
    ) %>%
    select(
      PLANT_SPECIES, starts_with("CLASS")
    ) %>%
    column_to_rownames(
      var = "PLANT_SPECIES"
    )
  set.seed(1234)
  df_P <- chisq.test(df_chi, simulate.p.value = T, B = 2000)$p.value
  df_ret <- data.frame(
    TE_NAME = NAM, 
    PVAL = df_P
  )
  return(df_ret)
  
})

bind_rows(temp_list) %>%
  mutate(
    FDR = p.adjust(PVAL, method = "fdr"), 
    CAT = ifelse(
      FDR < 0.01, 
      "Differnt", 
      "Similar"
    ),
    COMPARE = "TE_TE_RELATION"
  ) %>%
  left_join(
    x = df_temp,
    y = ., 
    by = "TE_NAME"
  ) -> df_result_TE_TE
glimpse(df_result_TE_TE)


print("# TE-TE distance ========================================================================")
glimpse(df_TE_dist)
df_result_TE_dist <- estimate_diff(df_TE_dist, "CLOSEST_OTHERTE_DISTANCE", "CLOSEST_OTHERTE_DISTANCE")
plot_TE_dist <- plot_diff(df_TE_dist, "CLOSEST_OTHERTE_DISTANCE", "CLOSEST_OTHERTE_DISTANCE")

print("# Outside 2kb, GC content ==============================================================")
glimpse(df_out2k_GC)
rm(df_temp)
df_temp <- df_out2k_GC %>%
  mutate(
    PLANT_SPECIES = factor(
      PLANT_Species, 
      levels = c("B. rapa", "B. oleracea"), 
      labels = c("B_rapa", "B_oleracea")
    )
  ) %>%
  left_join(
    x = .,
    y = df_share, 
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  filter(
    Share_or_Unique %in% "Share_in_AA_CC", 
    type %in% "Both sides"
  ) %>%
  # filter(
  #   direction %in% "right"
  # ) %>%
  # select(
  #   TE_CLASS, TE_NAME, TE_ID, PLANT_ID, PLANT_SPECIES, GC_content
  # ) %>%
  group_by(
    TE_CLASS, TE_NAME, TE_ID, PLANT_ID, PLANT_SPECIES
  ) %>%
  summarise(
    AVE = mean(GC_content, na.rm = T)
  ) %>%
  ungroup() %>%
  select(
    -PLANT_SPECIES
  ) %>%
  glimpse()
#glimpse(df_temp)
df_result_OUTGC_2kb <- estimate_diff(df_temp, "AVE", "Outside_2kb_GC")
#df_result_OUTGC_2kb <- estimate_diff(df_temp, "GC_content", "Outside_right_2kb_GC")

plot_OUTGC_2kb <- plot_diff(df_temp, "AVE", "Outside_2kb_GC")


glimpse(df_result_OUTGC_2kb)


print("# Save out =================================================================================")
fwrite(
  df_result_GC, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_GCContent.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_gene_dist, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_GeneDistance.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_LTR_AGE, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_LTRAge.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_LTR_length, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_LTRLength.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_meth_CG, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCG.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_meth_CHG, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCHG.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_meth_CHH, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCHH.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_TE_dist, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_TEDistance.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_TE_gene, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_RelationGeneTE.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_TE_TE, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_RelationTETE.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_whole_length, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_WholeLength.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_OUTGC_2kb, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_Outside_2kb_GC.txt",
  col.names = T, sep = "\t"
)

fwrite(
  df_result_count, 
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_Count.txt",
  col.names = T, sep = "\t"
)


#--------------------------------------------------------------------------------------------
# draw example plot -------------------------------------------------------------------------
print("# draw example plot, numeric features --------")

plot_LTR_length <- plot_LTR_length+
  ggplot2::coord_cartesian(
    ylim = c(NA, 400)
  )

cowplot::plot_grid(
  plotlist = list(plot_LTR_AGE, plot_LTR_length, plot_GC, plot_meth_CG), 
  align = "h", ncol = 4
)

ggsave2(
  filename = "/Fig_5_Assymetric_Evolution/ForExample_numeric.pdf", 
  width = 3.5, 
  height = 1.4
)

print("# catelog features --------")
df_gene_relation %>%
  filter(
    str_detect(TE_NAME, "(Copia0002)|(Copia0011)|(Copia0006)|(Copia0035)")
  ) %>%
  left_join(
    x = .,
    y = df_species, 
    by = "PLANT_ID"
  ) %>%
  summarise(
    .by = c(PLANT_SPECIES, TE_NAME, RELATION_GENE_TE), 
    COUNT = n()
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = c("Copia0002", "Copia0011", "Copia0006", "Copia0035"), 
      labels = c("Copia0002", "Copia0011", "Copia0006", "Copia0035")
    ), 
    RELATION_GENE_TE = factor(
      RELATION_GENE_TE,
      levels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial",
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene"),
      labels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial",
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene")
    ), 
    PLANT_SPECIES = factor(
      PLANT_SPECIES, 
      levels = c("B_rapa", "B_oleracea"), 
      labels = c("B. rapa", "B. oleracea")
    )
  ) %>%
  mutate(
    .by = c(PLANT_SPECIES, TE_NAME), 
    PREC = COUNT/sum(COUNT)*100
  ) %>%
  ggplot(
    aes(
      x = PLANT_SPECIES, 
      y = PREC, 
      fill = RELATION_GENE_TE
    )
  )+
  geom_col(
    position = "stack", 
    color = "black",
    linewidth = 0.1
  )+
  scale_fill_manual(
    breaks = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial",
               "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene", NA),
    values = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba", "grey")
  )+
  coord_cartesian(
    expand = F
  )+
  facet_grid(
    cols = vars(TE_NAME)
  )+
  theme_cowplot()+
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(
      size = 6
    ), 
    axis.text.x = element_text(
      size = 6, 
      angle = 90, 
      hjust = 1, 
      vjust = 0.5
    ),
    axis.text.y = element_text(
      size = 6
    ),
    strip.background = element_blank(), 
    strip.text = element_text(
      size = 6
    ),
    legend.position = "none",
    panel.spacing = unit(5, units = "mm")
  )+
  labs(
    y = "%"
  )

ggsave2(
  filename = "/Fig_5_Assymetric_Evolution/ForExample_category.pdf", 
  width = 3.5, 
  height = 1.6
)
