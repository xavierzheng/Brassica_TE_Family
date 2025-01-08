#-------------------------------------
# AIM
#	Summary the result of AA-CC different
#
#------------------------------------

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(corrplot, quietly = F))
suppressPackageStartupMessages(library(broom))

setDTthreads(threads = 8)
options(dplyr.summarise.inform=FALSE)

df_result_GC <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_GCContent.txt",
  header = T, sep = "\t"
)

df_result_gene_dist <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_GeneDistance.txt",
  header = T, sep = "\t"
)

df_result_LTR_AGE <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_LTRAge.txt",
  header = T, sep = "\t"
)

df_result_LTR_length <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_LTRLength.txt",
  header = T, sep = "\t"
)

df_result_meth_CG <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCG.txt",
  header = T, sep = "\t"
)

df_result_meth_CHG <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCHG.txt",
  header = T, sep = "\t"
)

df_result_meth_CHH <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_MethCHH.txt",
  header = T, sep = "\t"
)

df_result_TE_dist <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_TEDistance.txt",
  header = T, sep = "\t"
)

df_result_TE_gene <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_RelationGeneTE.txt",
  header = T, sep = "\t"
)

df_result_TE_TE <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_RelationTETE.txt",
  header = T, sep = "\t"
)

df_result_whole_length <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_WholeLength.txt",
  header = T, sep = "\t"
)

df_result_OUTGC_2kb <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_Outside_2kb_GC.txt",
  header = T, sep = "\t"
)

print("# stranger method, spearman correlation ---------------------")
df_result_count <- fread(
  file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_Count.txt",
  header = T, sep = "\t"
)

print("# Joint-------------------------------------------")
df_result <- rbind(
  df_result_whole_length, df_result_LTR_length, df_result_LTR_AGE, df_result_GC,  df_result_meth_CG, 
  df_result_meth_CHG, df_result_meth_CHH, df_result_gene_dist, df_result_TE_dist, df_result_OUTGC_2kb, 
  df_result_count # add this
)

df_temp1 <- df_result_TE_gene %>%
  select(
    TE_CLASS, TE_NAME, FDR, CAT, COMPARE
  ) %>%
  unique()

df_temp2 <- df_result_TE_TE %>%
  select(
    TE_CLASS, TE_NAME, FDR, CAT, COMPARE
  ) %>%
  unique()
df_result_relation <- rbind(df_temp1, df_temp2)

print("# Select --------------------------------------------")
glimpse(df_result)
unique(df_result$CAT)

rm(df_out)
df_out <- df_result %>%
  filter(
    nchar(CAT) >=1
  ) %>%
  group_by(
    TE_CLASS, COMPARE, CAT
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  reshape2::dcast(
    TE_CLASS + COMPARE ~ CAT, value.var = "COUNT", fill = 0
  ) 
glimpse(df_out)
# df_out %>%
#   fwrite(
#     file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_SummaryExcel.txt",
#     col.names = T, sep = "\t"
#   )

rm(df_out)
df_out <- df_result_relation %>%
  unique() %>%
  filter(
    nchar(CAT) >=1
  ) %>%
  group_by(
    TE_CLASS, COMPARE, CAT
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  reshape2::dcast(
    TE_CLASS + COMPARE ~ CAT, value.var = "COUNT", fill = 0
  )
glimpse(df_out)
# df_out %>%
#   fwrite(
#     file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_SummaryExcelRelation.txt",
#     col.names = T, sep = "\t"
#   )

print("# Full list of TE family =============================================")
glimpse(df_result)
unique(df_result$COMPARE)

# [1] "Whole_length"             "LTR_length"               "LTR_INSERT_AGE"           "GC_content"              
# [5] "CG_motif"                 "CHG_motif"                "CHH_motif"                "CLOSEST_GENE_DISTANCE"   
# [9] "CLOSEST_OTHERTE_DISTANCE" "Outside_2kb_GC"           "COUNT"    

 
df_out_numeric <- df_result %>%
  filter(
    nchar(CAT) >=1
  ) %>%
  mutate(
    Features = factor(
      COMPARE, 
      levels = c("COUNT", "Whole_length" ,"LTR_length" ,"LTR_INSERT_AGE" ,"GC_content", "CG_motif" ,"CHG_motif", "CHH_motif", "CLOSEST_GENE_DISTANCE", "CLOSEST_OTHERTE_DISTANCE" ,"Outside_2kb_GC"), 
      labels = c("Count", "Whole length", "LTR length", "LTR inserted age", "GC content", "CG motif", "CHG motif", "CHH motif", " Closed distance to gene", "Closest distance to TE", "Outer GC (%)")
    )
  ) %>%
  select(
    -COMPARE
  )

glimpse(df_out_numeric)

# df_out_numeric %>%
#   fwrite(
#     file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_SupplementNumeric.txt",
#     col.names = T, sep = "\t", quote = T
#   )

glimpse(df_result_relation)

df_out_catelog <- df_result_relation %>%
  filter(
    nchar(CAT) >=1
  ) %>%
  mutate(
    Features = factor(
      COMPARE, 
      levels = c("TE_GENE_RELATION", "TE_TE_RELATION"), 
      labels = c("TE-gene relation", "TE-TE relation")
    )
  ) %>%
  select(
    -COMPARE
  ) 

glimpse(df_out_catelog)

# df_out_catelog %>%
#   fwrite(
#     file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_SupplementCatelog.txt",
#     col.names = T, sep = "\t", quote = T
#   )

#====================================================
print("# family-wise siginificant features ==============================")
glimpse(df_out_numeric)
glimpse(df_out_catelog)

rm(df_temp_numeric)
df_temp_numeric <- df_out_numeric %>%
  filter(
    !(CAT %in% "No different")
  ) %>%
  select(
    TE_CLASS, TE_NAME, Features
  )
rm(df_temp_cat)
df_temp_cat <- df_out_catelog %>%
  filter(
    !(CAT %in% "Similar")
  ) %>%
  select(
    TE_CLASS, TE_NAME, Features
  )

df_temp_out <- rbind(df_temp_numeric, df_temp_cat) %>%
  mutate(
    FAKE = paste0("XX_", 1:n())
  ) %>%
  reshape2::dcast(
    TE_CLASS + TE_NAME ~ FAKE, value.var = "Features"
  ) %>%
  tidyr::unite(
    FEATURES, starts_with("XX_"), sep = ";", na.rm = T
  ) 

glimpse(df_temp_out)

# df_temp_out %>%
#   fwrite(
#     file = "/Data/Fig_5_Assymetric_Evolution/Compare_Info_AA_CC/Compare_AA_CC_SupplementCombineFeatures.txt",
#     col.names = T, sep = "\t", quote = T
#   )


#===================================
print("# share family = family with diverse features ?? ")

df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt",
  header = T, sep = "\t"
)
colnames(df_share) <- c("TE_CLASS", "TE_NAME", "Share_or_Unique")
glimpse(df_share)

df_check <- df_share %>%
  left_join(
    x = ., 
    y = df_temp_out, 
    by = c("TE_CLASS", "TE_NAME")
  ) 
glimpse(df_check)
df_stat <- df_check %>%
  mutate(
    SHARENESS = ifelse(
      Share_or_Unique %in% "Share_in_AA_CC", 
      "Share", 
      "Unique"
    ), 
    YesNoFeature = ifelse(
      is.na(FEATURES), 
      "NoFeature",
      "AtLeastOneFeature"
    )
  )  %>%
  glimpse()

table(df_stat$SHARENESS, df_stat$YesNoFeature)

# AtLeastOneFeature NoFeature
# Share                547       381
# Unique                 0       950


glimpse(df_stat)

temp_list <- lapply(unique(df_stat$TE_CLASS), FUN = function(NAM){
  
  df_temp <- df_stat %>%
    filter(
      TE_CLASS %in% NAM
    )
  #print(NAM)
  table(df_temp$SHARENESS, df_temp$YesNoFeature) %>%
    as.data.frame() %>%
    rename(
      SHARENESS = Var1, 
      YesNoFeature = Var2
    ) %>%
    mutate(
      TE_CLASS = NAM,
      SHARENESS = as.character(SHARENESS), 
      YesNoFeature = as.character(YesNoFeature), 
    ) -> df_ret
  return(df_ret)
    
})

df_feature <- bind_rows(temp_list)

rm(df_draw)
df_draw <- df_feature %>%
  mutate(
    TE_CLASS = factor(
      TE_CLASS, 
      levels = rev(c("DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "TIR/DTA",
                 "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")), 
      labels = rev(c("DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "TIR/DTA",
                 "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"))
    ), 
    YesNoFeature = factor(
      YesNoFeature, 
      levels = c("NoFeature", "AtLeastOneFeature"), 
      labels = c("No", "AtLeastOne")
    )
  ) %>%
  filter(
    SHARENESS %in% "Share"
  ) %>%
  group_by(
    TE_CLASS
  ) %>%
  mutate(
    PERC = Freq/(sum(Freq))*100
  ) %>%
  ungroup()

df_draw %>%
  ggplot(
    aes(
      x = PERC, 
      y = TE_CLASS,
      fill = YesNoFeature
    )
  )+
  geom_col(
    position = "stack"
  )+
  geom_vline(
    aes(
      xintercept = 547/(547 + 381)*100
    ), 
    color = "red"
  )+
  scale_fill_manual(
    breaks = c("AtLeastOne", "No"), 
    values = c("#F8B500", "#5C636E")
  )+
  coord_cartesian(
    expand = F
  )+
  theme_classic()+
  theme(
    axis.title = element_text(
      size =8
    ), 
    axis.text = element_text(
      size = 7
    ), 
    axis.ticks.y = element_blank(),
    legend.title = element_text(
      size = 8
    ), 
    legend.text = element_text(
      size = 7
    ), 
    legend.position = "bottom", 
    legend.margin = margin(0, 0, 0, 0, "pt"),
    plot.margin = margin(0, 0.5, 0.1, 0.1, "cm")
  )+
  labs(
    x = "%",
    y = "TE superfamily", 
    fill = "Feature Difference between\nB. rapa and B. oleracea"
  ) -> plot_chi
plot_chi

ggsave2(
  plot_chi, 
  filename = "/Fig_5_Assymetric_Evolution/driving_force_AddCount.pdf", 
  width = 3.3, 
  height = 3.5
)


#======================
print("# statistic analysis, goodness of fit model to background")
glimpse(df_feature)
df_bg <- df_feature %>%
  filter(
    SHARENESS %in% "Share"
  ) %>%
  group_by(
    YesNoFeature
  ) %>%
  summarise(
    TOT = sum(Freq)
  ) %>%
  ungroup()

df_bg

# # A tibble: 2 × 2
# YesNoFeature        TOT
# <chr>             <int>
# 1 AtLeastOneFeature   547
# 2 NoFeature           381

print("# using this for background probability")
# may be better for plot veritcal line
# 547/(547 + 381)
# [1] 0.5894397

df_for_chi <- df_feature %>%
  filter(
    SHARENESS %in% "Share"
  ) %>%
  left_join(
    x = .,
    y = df_bg, 
    by = "YesNoFeature"
  )

glimpse(df_for_chi)

order_TE_CLASS <- c("DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "TIR/DTA",
                    "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")

rm(temp_list)
temp_list <- lapply(order_TE_CLASS, FUN = function(NAM){
  
  df_temp <- df_for_chi %>%
    filter(TE_CLASS %in% NAM) 
  df_P <- chisq.test(
    x = df_temp$Freq,
    p = df_temp$TOT, 
    rescale.p = T
  )$p.value
  
  df_ret <- data.frame(
    TE_CLASS = NAM,
    PVAL = df_P
  )
  return(df_ret)
  
})
df_goodness_of_fit <- bind_rows(temp_list) %>%
  mutate(
    FDR = p.adjust(PVAL, method = "fdr"), 
    CAT = ifelse(
      FDR < 0.001,
      "Different", 
      "Similar"
    )
  )

df_goodness_of_fit 

#----------------------------------------
# add count data (from spearman correlation )

# TE_CLASS         PVAL          FDR       CAT
# 1 DNA/Helitron 8.868127e-01 8.868127e-01   Similar
# 2    LTR/Copia 1.846636e-09 8.309862e-09 Different
# 3    LTR/Gypsy 2.679391e-07 8.038172e-07 Different
# 4  LTR/unknown 2.076608e-05 4.672369e-05 Different
# 5      TIR/DTA 1.767894e-02 3.182209e-02   Similar
# 6      TIR/DTC 5.178130e-01 6.657596e-01   Similar
# 7      TIR/DTH 6.726092e-01 7.566854e-01   Similar
# 8      TIR/DTM 3.062360e-18 2.756124e-17 Different
# 9      TIR/DTT 1.044455e-01 1.566682e-01   Similar




