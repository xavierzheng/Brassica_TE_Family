#===================
# AIM
#   Are the distance and relationship different between B. rapa and B. oleracea
#
#===================


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 讓arrow可以多核心讀寫
arrow::set_cpu_count(12)
arrow::set_io_thread_count(12)
arrow::cpu_count()


df_rel <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Relationship.txt",
  header = T, 
  sep = "\t"
)

# glimpse(df_rel)

df_rel %>%
  mutate(
    SPECIES = ifelse(
      str_detect(PLANT_ID, "^D"),
      "AA", "CC"
    )
  ) %>%
  summarise(
    .by = c(SPECIES, RELATION_GENE_TE),
    NO = n()
  ) %>%
  mutate(
    .by = SPECIES,
    PERC = NO/sum(NO)*100
  ) %>%
  mutate(
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. olerecea")
    ),
    RELATION_GENE_TE = factor(
      RELATION_GENE_TE,
      levels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial",
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene"),
      labels = c("CLASS1_TEcontainGENE", "CLASS2_TEwithinGENE", "CLASS3_partial",
                 "CLASS4_ComplexEvent", "CLASS5_NotOverlapWithGene")
    )
  ) %>%
  ggplot(
    aes(
      x = SPECIES,
      y = PERC,
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
  theme_cowplot()+
  theme(
    axis.title = element_text(
      size = 7
    ),
    axis.text = element_text(
      size = 6
    ),
    legend.title = element_text(
      size = 7
    ),
    legend.text = element_text(
      size = 6
    )
  )+
  labs(
    x = "",
    y = "Percentage",
    fill = ""
  ) -> plot_rel
  


# SPECIES          RELATION_GENE_TE    NO      PERC
# 1       AA CLASS5_NotOverlapWithGene 27770 65.436637
# 2       AA      CLASS1_TEcontainGENE  3415  8.047033
# 3       AA       CLASS4_ComplexEvent  4174  9.835525
# 4       AA            CLASS3_partial  3180  7.493284
# 5       AA       CLASS2_TEwithinGENE  3899  9.187521
# 6       CC CLASS5_NotOverlapWithGene 52171 77.316715
# 7       CC      CLASS1_TEcontainGENE  4148  6.147280
# 8       CC            CLASS3_partial  3499  5.185471
# 9       CC       CLASS2_TEwithinGENE  4892  7.249878
# 10      CC       CLASS4_ComplexEvent  2767  4.100657



# distance ============================================
df <- fread(
  file = "/Data/Fig_4_TE_Feature/Summary_Gene_IntactTE_Distance.txt",
  header = T, sep = "\t"
)

glimpse(df)

df %>%
  mutate(
    SPECIES = ifelse(
      str_detect(PLANT_ID, "^D"),
      "AA", "CC"
    )
  ) %>%
  summarise(
    .by = SPECIES,
    AVE = mean(CLOSEST_GENE_DISTANCE, na.rm = T),
    SD = sd(CLOSEST_GENE_DISTANCE, na.rm = T)
  )

#   SPECIES      AVE       SD
# 1      AA 6220.331 26706.60
# 2      CC 6365.541 25821.52  

df_draw <- df %>%
  mutate(
    SPECIES = ifelse(
      str_detect(PLANT_ID, "^D"),
      "AA", "CC"
    )
  ) %>%
  filter(
    !is.na(CLOSEST_GENE_DISTANCE)
  ) %>%
  mutate(
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. olerecea")
    )
  )

df_draw %>%
  ggplot(
    aes(
      x = SPECIES,
      y = log10(CLOSEST_GENE_DISTANCE)
    )
  )+
  geom_boxplot()+
  theme_cowplot()+
  theme(
    axis.title = element_text(
      size = 7
    ),
    axis.text = element_text(
      size = 6
    )
  )+
  labs(
    x = "",
    y = "log10 Distance to closest gene"
  ) -> plot_dis

plot_grid(
  plotlist = list(plot_rel, plot_dis),
  ncol = 2, rel_widths = c(2,1)
)
ggsave2(
  filename = "GeneTERelation_DiffAACC.pdf",
  width = 7.2,
  height = 3
)

glimpse(df_draw)

wilcox.test(CLOSEST_GENE_DISTANCE ~ SPECIES, data = df_draw)


# Wilcoxon rank sum test with continuity correction
# 
# data:  CLOSEST_GENE_DISTANCE by SPECIES
# W = 598846854, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
