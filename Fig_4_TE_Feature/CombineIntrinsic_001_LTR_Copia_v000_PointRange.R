#====================================
# AIM
#	Focus on Copia
#   Combine Count, length, Insertion Time, GC content, and Methylation frequency into heatmap
#====================================

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))

setDTthreads(threads = 4)
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

print("# import share and uniq info=======================")
df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt",
  header = T, sep = "\t"
)

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


#===========================================================================
print(" draw heatmap using tile according to Family size ========================================")

list_CLASS <- str_replace_all(cat_list, pattern = "_", replacement = "/")
order_PLANT_NAME <- c("Korso", "HDEM", "OX_heart", "JZSv2", "RCBO", "TUA", "PCA", "OIB", "MIZ", "CXA", "Chiifu", "Z1")


print("# COUNT =============================================")
df_draw_temp_count <- df_count %>%
  filter(
    TE_METHOD %in% "structural"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  group_by(
    TE_CLASS, PLANT_ID ,TE_NAME
  ) %>%
  summarize(
    COUNT = n()
  ) %>%
  ungroup() 

print("# ORDER based on COUNT ===============================")
df_order <- df_draw_temp_count %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  summarise(
    SUM = sum(COUNT)
  ) %>%
  ungroup() %>%
  arrange(
    desc(SUM)
  ) %>%
  select(
    -SUM
  )
# fwrite(
#   df_order,
#   file = "Order_Copia.txt",
#   col.names = T, sep = "\t"
# )

print("# Draw count of structural TE ================================")
glimpse(df_draw_temp_count)
rm(df_temp)
df_temp <- df_draw_temp_count %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(log10(COUNT), na.rm = T), 
    Q25 = quantile(log10(COUNT), probs = 0.25, na.rm = T), 
    Q75 = quantile(log10(COUNT), probs = 0.75, na.rm = T)
  ) %>%
  ungroup() 

ggplot()+
  geom_hline(
    yintercept = quantile(log10(df_temp$COUNT), probs = c(0.5), na.rm = T), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = quantile(log10(COUNT), probs = 0.1, na.rm = T)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = quantile(log10(COUNT), probs = 0.9, na.rm = T),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN,
      ymin = Q25, 
      ymax = Q75
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$COUNT), probs = c(0.1, 0.5, 0.9))),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(-0.1, NA)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "log10\nNo. TE") -> plot_count
plot_count

print("# Draw whole length ===============================")
rm(df_temp)
df_temp <- df_count %>%
  filter(
    TE_METHOD %in% "structural"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_LEN = GFF_END - GFF_START +1, 
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  )
rm(df_draw)
df_draw<- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(TE_LEN, na.rm = T), 
    Q25 = quantile(TE_LEN, probs = 0.25, na.rm = T), 
    Q75 = quantile(TE_LEN, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = quantile(df_temp$TE_LEN/1000, probs = c(0.5), na.rm = T), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf, #quantile(df_temp$TE_LEN/1000, probs = 0.1),
      ymax = quantile(TE_LEN/1000, probs = 0.1, na.rm = T)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = quantile(TE_LEN/1000, probs = 0.9, na.rm = T),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN/1000,
      ymin = Q25/1000, 
      ymax = Q75/1000
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$TE_LEN/1000, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.1), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(3, 7)
    #ylim = quantile(df_temp$TE_LEN/1000, probs = c(0.01, 0.95))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Whole\nlength\n(kb)")-> plot_len

plot_len

print("# Draw Insertion time ======================================================")
glimpse(df_age)

rm(df_temp)
df_temp <- df_age %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  filter(
    !(is.na(INSERT_AGE))
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME,
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  )
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(INSERT_AGE, na.rm = T), 
    Q25 = quantile(INSERT_AGE, probs = 0.25, na.rm = T), 
    Q75 = quantile(INSERT_AGE, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$INSERT_AGE/1e+6, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf, #quantile(df_temp$INSERT_AGE/1e+6, probs = 0.2),
      ymax = round(quantile(INSERT_AGE/1e+6, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(INSERT_AGE/1e+6, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN/1e+6,
      ymin = Q25/1e+6, 
      ymax = Q75/1e+6
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0), add = c(0, 0)), 
    breaks = round(unname(quantile(df_temp$INSERT_AGE/1e+6, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
  )+
  coord_cartesian(
    ylim = c(-0.2, 2.5)
    #ylim = quantile(df_temp$INSERT_AGE/1e+6, probs = c(0, 0.99))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Inserted\nage\n(MYA)")-> plot_age

plot_age

print("# GC content of flanking repeat ===========================================")
glimpse(df_gc)
rm(df_temp)
df_temp <- df_gc %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(GC, na.rm = T), 
    Q25 = quantile(GC, probs = 0.25, na.rm = T), 
    Q75 = quantile(GC, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$GC*100, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(GC*100, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(GC*100, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN*100,
      ymin = Q25*100, 
      ymax = Q75*100
    ), 
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$GC*100, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = quantile(df_temp$GC*100, probs = c(0.01, 0.99), na.rm = T)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "GC\ncontent\n(%)")-> plot_GC

plot_GC

print("# flanking repeat length ===========================================")
glimpse(df_gc)
rm(df_temp)
df_temp <- df_gc %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(LEN, na.rm = T), 
    Q25 = quantile(LEN, probs = 0.25, na.rm = T), 
    Q75 = quantile(LEN, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$LEN, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(LEN, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(LEN, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN,
      ymin = Q25, 
      ymax = Q75
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$LEN, probs = c(0.1, 0.5, 0.9), na.rm = T)),1), #c(250, 500, 750),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = quantile(df_temp$LEN, probs = c(0, 0.995))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "LTR\nLength\n(bp)")-> plot_flank_LEN

plot_flank_LEN

print("# Draw methylation CG =======================================")
rm(df_temp)
df_temp <- df_meth %>%
  filter(
    METHYL %in% "CG"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(ACCUM_METHYL_PERC, na.rm = T), 
    Q25 = quantile(ACCUM_METHYL_PERC, probs = 0.25, na.rm = T), 
    Q75 = quantile(ACCUM_METHYL_PERC, probs = 0.75, na.rm = T)
  ) %>%
  ungroup() 

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(ACCUM_METHYL_PERC, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(ACCUM_METHYL_PERC, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN,
      ymin = Q25, 
      ymax = Q75
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(0, NA)
    #ylim = quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0, 0.99))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Methyl\nCG\n(%)") -> plot_meth_CG

plot_meth_CG

print("# Draw methylation CHH =======================================")  
rm(df_temp)
df_temp <- df_meth %>%
  filter(
    METHYL %in% "CHH"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(ACCUM_METHYL_PERC, na.rm = T), 
    Q25 = quantile(ACCUM_METHYL_PERC, probs = 0.25, na.rm = T), 
    Q75 = quantile(ACCUM_METHYL_PERC, probs = 0.75, na.rm = T)
  ) %>%
  ungroup() 

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(ACCUM_METHYL_PERC, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(ACCUM_METHYL_PERC, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN,
      ymin = Q25, 
      ymax = Q75
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Methyl\nCHH\n(%)") -> plot_meth_CHH

plot_meth_CHH

print("# Draw methylation CHG =======================================")
rm(df_temp)
df_temp <- df_meth %>%
  filter(
    METHYL %in% "CHG"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(ACCUM_METHYL_PERC, na.rm = T), 
    Q25 = quantile(ACCUM_METHYL_PERC, probs = 0.25, na.rm = T), 
    Q75 = quantile(ACCUM_METHYL_PERC, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(ACCUM_METHYL_PERC, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(ACCUM_METHYL_PERC, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN,
      ymin = Q25, 
      ymax = Q75
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.1), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.02, 0.95), na.rm = T)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Methyl\nCHG\n(%)") -> plot_meth_CHG

plot_meth_CHG

print("# Draw relationship to gene ===========================================")
glimpse(df_gene_relation)
rm(df_draw)
df_draw <- df_gene_relation %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  group_by(
    TE_CLASS, TE_NAME, RELATION_GENE_TE
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  mutate(
    SUM_COUNT = sum(COUNT, na.rm = T), 
    PERC_COUNT = COUNT/SUM_COUNT*100
  ) %>%
  ungroup() %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    ), 
    RELATION_GENE_TE = factor(
      RELATION_GENE_TE, 
      levels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial", 
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene"), 
      labels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial", 
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene")
    )
  ) 
# fwrite(
#   df_draw, 
#   file = "Copia_RELATION_GENE_TE.txt", 
#   col.names = T, sep = "\t", row.names = F, quote = F
# )

df_draw %>%
  ggplot(
    aes(
      x = TE_NAME,
      y = PERC_COUNT, 
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
  scale_y_continuous(
    breaks = c(25, 50, 75),
    expand = expansion(mult = c(0, 0), add = c(0, 2))
  )+
  coord_cartesian(
    ylim = c(0, 100)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    legend.position = "none",
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Relation\nto\nGene") -> plot_relation_gene

plot_relation_gene  

print("# Distance to closest gene ==================================")
rm(df_temp)
df_temp <- df_gene_dist %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  filter(
    !is.na(CLOSEST_GENE_DISTANCE)
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  )
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(CLOSEST_GENE_DISTANCE, na.rm = T), 
    Q25 = quantile(CLOSEST_GENE_DISTANCE, probs = 0.25, na.rm = T), 
    Q75 = quantile(CLOSEST_GENE_DISTANCE, probs = 0.75, na.rm = T)
  ) %>%
  ungroup() 

ggplot()+
  geom_hline(
    yintercept = round(quantile(log10(df_temp$CLOSEST_GENE_DISTANCE), probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(log10(CLOSEST_GENE_DISTANCE), probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(log10(CLOSEST_GENE_DISTANCE), probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw %>%
      mutate(
        MEDIAN = ifelse(
          is.na(MEDIAN),
          0.01,
          MEDIAN
        ),
        Q25 = ifelse(
          is.na(Q25),
          0.01,
          Q25
        ),
        Q75 = ifelse(
          is.na(Q75),
          0.01,
          Q75
        )
      ),
    mapping = aes(
      x = TE_NAME, 
      y = log10(MEDIAN),
      ymin = log10(Q25), 
      ymax = log10(Q75)
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$CLOSEST_GENE_DISTANCE), probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(1.8, 4.8)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Distance\nto gene\nlog10(bp)") -> plot_dist_gene

plot_dist_gene

print("# Relation to closest TE =========================================")
glimpse(df_TE_relation)
rm(df_draw)
df_draw <- df_TE_relation %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  group_by(
    TE_CLASS, TE_NAME, RELATION_OTHERTE_TE
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  mutate(
    SUM_COUNT = sum(COUNT, na.rm = T), 
    PERC_COUNT = COUNT/SUM_COUNT*100
  ) %>%
  ungroup() %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    ), 
    RELATION_OTHERTE_TE = factor(
      RELATION_OTHERTE_TE, 
      levels = c("CLASS1_TEcontainOTHERTE", "CLASS2_TEwithinOTHERTE", "CLASS3_partial", 
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithOTHERTE"), 
      labels = c("CLASS1_TEcontainOTHERTE", "CLASS2_TEwithinOTHERTE", "CLASS3_partial", 
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithOTHERTE")
    )
  ) 

fwrite(
  df_draw, 
  file = "Copia_RELATION_OTHERTE_TE.txt", 
  col.names = T, sep = "\t", row.names = F, quote = F
)

df_draw %>%
  ggplot(
    aes(
      x = TE_NAME,
      y = PERC_COUNT, 
      fill = RELATION_OTHERTE_TE
    )
  )+
  geom_col(
    position = "stack",
    color = "black",
    linewidth = 0.1
  )+
  scale_fill_manual(
    breaks = c("CLASS1_TEcontainOTHERTE", "CLASS2_TEwithinOTHERTE", "CLASS3_partial", 
               "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithOTHERTE", NA), 
    values = c("#d7191c", "#fdae61", "#ffffbf", "#abdda4", "#2b83ba", "grey")
  )+
  scale_y_continuous(
    breaks = c(25, 50, 75),
    expand = expansion(mult = c(0, 0), add = c(0, 2))
  )+
  coord_cartesian(
    ylim = c(0, 100)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    legend.position = "none",
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Relation\nto\nTE") -> plot_relation_TE
plot_relation_TE

print("# Draw the distance to closest TE =====================================")
glimpse(df_TE_dist)
rm(df_temp)
df_temp <- df_TE_dist %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  filter(
    !is.na(CLOSEST_OTHERTE_DISTANCE)
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) 
rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(CLOSEST_OTHERTE_DISTANCE, na.rm = T), 
    Q25 = quantile(CLOSEST_OTHERTE_DISTANCE, probs = 0.25, na.rm = T), 
    Q75 = quantile(CLOSEST_OTHERTE_DISTANCE, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(log10(df_temp$CLOSEST_OTHERTE_DISTANCE), probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(log10(CLOSEST_OTHERTE_DISTANCE), probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(log10(CLOSEST_OTHERTE_DISTANCE), probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = log10(MEDIAN),
      ymin = log10(Q25), 
      ymax = log10(Q75)
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$CLOSEST_OTHERTE_DISTANCE), probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(0.2, 4)
    #ylim = quantile(log10(df_temp$CLOSEST_OTHERTE_DISTANCE), probs = c(0, 0.99))
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Distance\nto TE\nlog10(bp)") -> plot_dist_TE

plot_dist_TE

print("# Draw outside 2kb GC content ============================================")
glimpse(df_out2k_GC)
rm(df_temp)
df_temp <- df_out2k_GC %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  filter(
    type %in% "Both sides"
  ) %>%
  group_by(
    TE_ID, TE_NAME, TE_CLASS, PLANT_ID, PLANT_Species
  ) %>%
  summarise(
    GC_OUT = mean(GC_content, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(
    TE_NAME = factor(
      TE_NAME,
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  )

rm(df_draw)
df_draw <- df_temp %>%
  group_by(
    TE_NAME
  ) %>%
  summarise(
    MEDIAN = median(GC_OUT, na.rm = T), 
    Q25 = quantile(GC_OUT, probs = 0.25, na.rm = T), 
    Q75 = quantile(GC_OUT, probs = 0.75, na.rm = T)
  ) %>%
  ungroup()

ggplot()+
  geom_hline(
    yintercept = round(quantile(df_temp$GC_OUT*100, probs = c(0.5), na.rm = T),1), 
    linetype = "solid", 
    color = "#bababa",
    linewidth = 0.2
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,#quantile(df_temp$COUNT, probs = 0.1),
      ymax = round(quantile(GC_OUT*100, probs = 0.1, na.rm = T),1)
    ), 
    fill = "#1ecbe1", 
    alpha = 0.5
  )+
  geom_rect(
    data = df_temp, 
    mapping = aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = round(quantile(GC_OUT*100, probs = 0.9, na.rm = T),1),
      ymax = Inf
    ), 
    fill = "#edd612", 
    alpha = 0.5
  )+
  geom_pointrange(
    data = df_draw, 
    mapping = aes(
      x = TE_NAME, 
      y = MEDIAN*100,
      ymin = Q25*100, 
      ymax = Q75*100
    ),
    color = "#a63603",
    linewidth = 0.2, 
    fatten = 1
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$GC_OUT*100, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0), add = c(0, 0))
  )+
  coord_cartesian(
    #ylim = c(26, 45)
    ylim = quantile(df_temp$GC_OUT*100, probs = c(0, 0.99), na.rm = T)
  )+
  theme_classic()+
  theme(
    axis.line = element_line(
      linewidth = 0.5, 
      color = "black"
    ), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    axis.title.y = element_text(
      size = 7, 
      angle = 90,
      hjust = 0.5, 
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    )
  )+
  labs(y = "Outer\n2 kb\nGC(%)") -> plot_outer2k_GC
plot_outer2k_GC

print("# Draw tile for share or unique infomation ============================")
df_share %>%
  filter(
    TE_superfamily %in% "LTR/Copia"
  ) %>%
  rename(
    TE_CLASS = TE_superfamily,
    TE_NAME = TE_Family
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) %>%
  ggplot(
    aes(
      x = TE_NAME, 
      y = "Unique or share",
      fill = Share_or_Unique
    )
  )+
  geom_tile(
    color = "black",
    linewidth = 0.1
  )+
  scale_fill_manual(
    breaks = c("Share_in_AA_CC", "Unique_in_AA", "Unique_in_CC"),
    values = c("#f0f0f0", "#998EC3", "#F1A340")
  )+
  theme_classic()+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(), 
    axis.title = element_blank(), 
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    ),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    legend.position = "none"
  ) -> plot_share
plot_share

print("# Draw heatmap / tile for their Z-scaled aboundance around AA and CC ========================")
df_z_score %>%
  data.table::melt.data.table(
    id.vars = c("TE_CLASS", "TE_NAME"), 
    variable.name = "PLANT_NAME", 
    value.name = "Z_SCORE"
  ) %>%
  filter(
    TE_CLASS %in% "LTR/Copia"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
  ) %>%
  mutate(
    PLANT_NAME = factor(
      PLANT_NAME, 
      levels = order_PLANT_NAME, 
      labels = order_PLANT_NAME
    ), 
    TE_NAME = factor(
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) %>%
  ggplot(
    aes(
      x = TE_NAME, 
      y = PLANT_NAME,
      fill = Z_SCORE
    )
  )+
  geom_tile()+
  scale_fill_viridis(
    option = "D", 
    discrete = F
  )+
  theme_classic()+
  theme(
    axis.line = element_blank(), 
    axis.ticks = element_blank(),
    axis.title = element_blank(), 
    axis.text.x = element_text(
      size = 6,
      angle = 90, 
      hjust = 1, 
      vjust = 0.5
    ), 
    axis.text.y = element_text(
      size = 6, 
      angle = 0, 
      hjust = 1, 
      vjust = 0.5
    ),
    plot.margin = ggplot2::unit(c(0, 0, 0, 0), units = "in"),
    legend.position = "bottom", 
    legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
    legend.text = element_text(
      size = 6
    ),
    legend.title = element_text(
      size = 6
    ),
    legend.key.height = unit(0.3, units = "cm")
  ) -> plot_heat
plot_heat


# plot together ==============================

COMBINE_PLOT <- plot_outer2k_GC/plot_dist_TE/plot_relation_TE/plot_dist_gene/plot_relation_gene/plot_meth_CG/plot_meth_CHG/plot_meth_CHH/plot_GC/plot_age/plot_flank_LEN/plot_len/plot_count/plot_share/plot_heat + 
  patchwork::plot_layout(heights = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0.5, 5))

#COMBINE_PLOT

cowplot::ggsave2(
  plot = COMBINE_PLOT,
  filename = "/Fig_4_TE_Feature/FullInfo_LTR_Copia.pdf", 
  width = 7.2,
  height = 9.7
)

# #=======================================================================
# #=======================================================================
# #=======================================================================
# #=======================================================================
# print("# Statistic drawing, especially focus on different on B. rapa and B. oleracea ===============")

# df_species <- fread(
#   file = "/nfs/project1/result/Brassica_TE/1000_plot_TotalTE/03_ParserTE/20_CombineIntrinsic/Plant_Species.txt", 
#   header = T, sep = "\t"
# )

# df_share <- fread(
#   file = "/nas/nas6/Brassica_TE/1000_plot_TotalTE/03_ParserTE/04_IntactSingleLineAndHomoTE/TE_family_uniqueness.txt",
#   header = T, sep = "\t"
# )
# colnames(df_share) <- c("TE_CLASS", "TE_NAME", "Share_or_Unique")


# rm(df_temp)
# df_temp <- df_count %>%
#   left_join(
#     x = .,
#     y = df_species, 
#     by = c("PLANT_ID")
#   ) %>%
#   left_join(
#     x = .,
#     y = df_share, 
#     by = c("TE_CLASS", "TE_NAME")
#   ) %>%
#   filter(
#     TE_CLASS %in% "LTR/Copia", 
#     Share_or_Unique %in% "Share_in_AA_CC"
#   ) %>%
#   mutate(
#     LEN = GFF_END - GFF_START + 1, 
#     PLANT_SPECIES = as.factor(PLANT_SPECIES)
#   )

# glimpse(df_temp)
# df_stat <- df_temp %>%
#   select(
#     PLANT_SPECIES, PLANT_ID, TE_NAME, TE_ID, LEN
#   ) %>%
#   mutate(
#     TE_NAME = factor(
#       TE_NAME, 
#       df_order$TE_NAME, 
#       df_order$TE_NAME
#     )
#   ) %>%
#   nest_by(
#     col = TE_NAME
#   ) %>%
#   summarise(
#     PVAL_SPECIES = wilcox.test(LEN~PLANT_SPECIES, data = data, exact = F, conf.int = T)$p.value
#   ) %>%
#   ungroup() %>%
#   mutate(
#     PFDR_SPECIES = p.adjust(PVAL_SPECIES, method = "fdr"), 
#     PFDR_SIG = ifelse(
#       PFDR_SPECIES < 0.001, 
#       "SIG", 
#       "NO"
#     )
#   )
  
# glimpse(df_temp)
# glimpse(df_stat)

# ggplot()+
#   geom_col(
#     data = df_stat %>%
#       filter(
#         PFDR_SIG %in% "SIG"
#       ) %>%
#       rename(
#         TE_NAME = col
#       ),
#     mapping = aes(
#       x = TE_NAME, 
#       y = Inf
#     ), 
#     fill = "#feebe2", 
#     color = "#feebe2"
#   )+
#   geom_pointrange(
#     data = df_temp %>%
#       mutate(
#         TE_NAME = factor(
#           TE_NAME, 
#           df_order$TE_NAME, 
#           df_order$TE_NAME
#         )
#       ) %>%
#       group_by(TE_NAME, PLANT_SPECIES) %>%
#       summarise(
#         MED = median(LEN, na.rm = T),
#         Q25 = quantile(LEN, probs = 0.25, na.rm = T), 
#         Q75 = quantile(LEN, probs = 0.75, na.rm = T)
#       ) %>%
#       ungroup(),
#     mapping = aes(
#       x = TE_NAME, 
#       y = MED/1000, 
#       ymin = Q25/1000, 
#       ymax = Q75/1000, 
#       color = PLANT_SPECIES
#     ),
#     size = 0.3
#   )+
#   scale_color_manual(
#     breaks = c("B_rapa", "B_oleracea"), 
#     values = c("#998EC3", "#F1A340")
#   )+
#   scale_y_continuous(
#     expand = c(0, 0)
#   )+
#   coord_cartesian(
#     ylim = c(0, NA)
#   )+
#   theme_classic()+
#   theme(
#     legend.position = "bottom",
#     legend.title = element_blank(),
#     axis.title.x = element_blank(), 
#     axis.text.x = element_text(
#       size = 7, 
#       angle = 90, 
#       hjust = 0.5, 
#       vjust = 0.5
#     )
#   )+
#   labs(
#     y = "Whole length (kb)"
#   ) -> plot_species_whole_len
# plot_species_whole_len

# #glimpse(df_temp)
# #summary(lm(LEN ~ TE_NAME + PLANT_SPECIES + TE_NAME*PLANT_SPECIES, data = df_temp))
# #model <- nnet::multinom(LEN ~ TE_NAME + PLANT_SPECIES + TE_NAME*PLANT_SPECIES, data = df_temp)
