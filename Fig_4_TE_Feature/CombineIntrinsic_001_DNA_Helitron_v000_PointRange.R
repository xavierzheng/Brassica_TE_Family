#====================================
# AIM
#	Focus on Helitron
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
    TE_CLASS %in% "DNA/Helitron"
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
  select(-SUM)
glimpse(df_order)
# fwrite(
#   df_order, 
#   file = "Order_Helitron.txt", 
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$COUNT), probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.05), add = c(0, 0))
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$TE_LEN/1000, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(NA, NA)
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


print("# GC content of flanking repeat ===========================================")
glimpse(df_gc)
rm(df_temp)
df_temp <- df_gc %>%
  filter(
    TE_CLASS %in% "DNA/Helitron"
  )  %>%
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$GC*100, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.05), add = c(0, 0))
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$LEN, probs = c(0.1, 0.5, 0.9), na.rm = T)),1), #c(250, 500, 750),
    expand = expansion(mult = c(0.05, 0), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = quantile(df_temp$LEN, probs = c(0, 0.995), na.rm = T)
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.05), add = c(0, 0))
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0, 0.05), add = c(0, 0))
  )+
  coord_cartesian(
    ylim = c(NA, NA)
    #ylim = quantile(df_temp$ACCUM_METHYL_PERC, probs = c(0.05, 0.95))
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
    TE_CLASS %in% "DNA/Helitron"
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
#   file = "Helitron_RELATION_GENE_TE.txt", 
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$CLOSEST_GENE_DISTANCE), probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))
  )+
  coord_cartesian(
    #ylim = c(NA, NA)
    ylim = c(1.5, 4.2)
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
    TE_CLASS %in% "DNA/Helitron"
  ) %>%
  left_join(
    x = df_order,
    y = .,
    by = c("TE_CLASS", "TE_NAME")
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
  file = "Helitron_RELATION_OTHERTE_TE.txt", 
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(log10(df_temp$CLOSEST_OTHERTE_DISTANCE), probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0.05, 0.05), add = c(0.03, 0))
  )+
  coord_cartesian(
    #ylim = c(NA, NA)
    ylim = quantile(df_temp$CLOSEST_OTHERTE_DISTANCE/1000, probs = c(0, 0.93), na.rm = T)
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
    TE_CLASS %in% "DNA/Helitron"
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
    fatten = 0.4
  )+
  scale_y_continuous(
    breaks = round(unname(quantile(df_temp$GC_OUT*100, probs = c(0.1, 0.5, 0.9), na.rm = T)),1),
    expand = expansion(mult = c(0.01, 0), add = c(0, 0))
  )+
  coord_cartesian(
    #ylim = c(NA, NA)
    ylim = quantile(df_temp$GC_OUT*100, probs = c(0.02, 0.98))
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
    TE_superfamily %in% "DNA/Helitron"
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
    TE_CLASS %in% "DNA/Helitron"
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
      size = 4,
      angle = 90, 
      hjust = 1, 
      vjust = 0.5
    ), 
    axis.text.y = element_text(
      size = 5, 
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

COMBINE_PLOT <- plot_outer2k_GC/plot_dist_TE/plot_relation_TE/plot_dist_gene/plot_relation_gene/plot_meth_CG/plot_meth_CHG/plot_meth_CHH/plot_GC/plot_len/plot_count/plot_share/plot_heat + 
  patchwork::plot_layout(heights = c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0.5, 5))

#COMBINE_PLOT

cowplot::ggsave2(
  plot = COMBINE_PLOT,
  filename = "/Fig_4_TE_Feature/FullInfo_DNA_Helitron.pdf", 
  width = 18,
  height = 7.2
)









