#----------------------------------
# AIM
#	check the prevelance of TE family in AA and CC species using Upset
#
#----------------------------------

# setwd("/nfs/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/04_IntactSingleLineAndHomoTE")

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(pheatmap, quietly = F))
suppressPackageStartupMessages(library(ComplexHeatmap, quietly = F))
suppressPackageStartupMessages(library(UpSetR, quietly = F))

setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

# define function--------------------------------
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements  
  return(data)
}


#-------------------------------------------------
df_plant_species <- fread(
  file = "/Data/Fig_2_Distribution/Plant_Species.txt", 
  header = T, sep = "\t"
)

df_plantid <- fread(
  file = "/Data/Fig_2_Distribution/Plant_IDNAME.txt", 
  header = T, sep = "\t"
)

df <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt.gz", 
  header = T, sep = "\t", nThread = 4
) %>%
  left_join(
    x = ., 
    y = df_plantid, 
    by = "PLANT_ID"
  )

glimpse(df)
unique(df$PLANT_NAME)

plant_order <- c("RCBO", "Korso", "HDEM", "OX_heart", "JZSv2", "PCA", "CXA", "OIB", "MIZ", "TUA", "Chiifu", "Z1")

temp_list <- lapply(plant_order, FUN = function(NAM){
  
  temp_vec <- df %>%
    filter(
      TE_METHOD %in% "structural"
    ) %>%
    filter(
      PLANT_NAME %in% NAM
    ) %>%
    select(
      TE_NAME
    ) %>%
    unique() %>%
    pull()
  
  return(temp_vec)
})

names(temp_list) <- plant_order

cairo_pdf(
  filename = "/Fig_2_Distribution/Upset_TEfamily_v2.pdf", 
  width = 6, 
  height = 5, family = "sans"
)

UpSetR::upset(
  data = fromList(temp_list), 
  sets = plant_order, #給他這個list裡面想拿來畫的資料。尤其要定義順序時候很好用
  nsets = 12, #用這參數決定使用幾個資料集來畫圖，就是上面那7個
  nintersects = 20, #只畫出前40個交集的資料。如果想畫全部的交集，記得寫NA
  keep.order = T, #維持sets的順序
  order.by = "freq", #依照交集的數量大小，從大排到小
  empty.intersections = "on", #就算是交集數量是0也要畫出來
  number.angles = 90, #交集數量barplot的上方數字，要不要轉角度，可是轉角度後很醜，要手動改位置
  text.scale = c(1.5, 1, 1, 1, 1.5, 1.2), #文字要不要放大，可以調整6個區域的文字大小，如下圖顯示
  point.size = 3, #底下的點圖的點要多大
  line.size = 1, #底下的點圖的線要多粗
  mb.ratio = c(0.5, 0.5), # 控制上下圖片比例大小
  mainbar.y.label = "No. TE family",  #控制 y axis名字
  sets.x.label = "Total TE family",  # 控制底下的軸名字，原來叫做Set size
  scale.intersections = "identity", #交集數量可以改成log10或是log2顯示，如果要原始數值，就寫identity
  scale.sets = "identity"
)

dev.off()


print("# generate table ------------------------------------")

print("## PAV-----------------------------------")
df_upset <- fromList(temp_list)
df_upset <- df_upset %>%
  rownames_to_column(
    var = "TE_Family"
  )
glimpse(df_upset)

fwrite(
  df_upset, 
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness_upset.txt", 
  col.names = T,
  row.names = F, 
  sep = "\t", 
  quote = F
)

df_upset <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness_upset.txt", 
  header = T, sep = "\t"
)
glimpse(df_upset)



#========================================================================================
print("# test if share or unique infomation is different -----------------------------")

df_share <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt", 
  header = T, sep = "\t"
)

# df_upset %>%
#   rowwise() %>%
#   mutate(
#     B_rapa = sum(PCA, CXA, OIB, MIZ, TUA, Chiifu, Z1), 
#     B_oleracea = sum(RCBO, Korso, HDEM, OX_heart, JZSv2), 
#   ) %>%
#   ungroup() %>%
#   mutate(
#     Share_or_Unique = ifelse(
#       B_rapa >=1 & B_oleracea >=1, 
#       "Share_in_AA_CC", 
#       ifelse(
#         B_rapa >=1 & B_oleracea ==0, 
#         "Unique_in_AA", 
#         ifelse(
#           B_rapa ==0 & B_oleracea >=1,
#           "Unique_in_CC", 
#           "OTHER"
#         )
#       )
#     )
#   ) %>%
#   left_join(
#     x = .,
#     y = df_share, 
#     by = c("TE_Family")
#   ) %>%
#   mutate(
#     TAG = ifelse(
#       Share_or_Unique.x == Share_or_Unique.y,
#       "SAME", 
#       "DIFF"
#     )
#   ) %>%
#   View()
#========================================================================================

# df_upset %>%
#   filter(
#     Z1 %in% 0, 
#     Chiifu  %in% 0, 
#     TUA %in% 0, 
#     MIZ %in% 0, 
#     OIB %in% 0, 
#     CXA %in% 0, 
#     PCA  %in% 0, 
#     HDEM %in% 0, 
#     Korso %in% 0, 
#     RCBO %in% 0, 
#     JZSv2 %in% 1,
#     OX_heart %in% 1
#   ) %>%
#   arrange(
#     TE_Family
#   ) %>%
#   select(
#     TE_Family
#   ) %>%
#   pull -> LIST_cabbage



print("# Total TE per family ----------------------------")

glimpse(df)
df_count <- df %>%
  group_by(
    TE_CLASS, TE_NAME, PLANT_NAME, TE_METHOD
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    TE_CLASS + TE_NAME + PLANT_NAME ~ TE_METHOD, 
    value.var = "COUNT", 
    fill = 0
  ) %>%
  mutate(
    TOTAL = structural + homology
  ) %>%
  rename(
    STRUCTURAL = structural,
    HOMOLOGY = homology
  ) %>%
  select(
    TE_CLASS, TE_NAME, PLANT_NAME, TOTAL, STRUCTURAL, HOMOLOGY
  ) %>%
  glimpse()
  
fwrite(
  df_count, 
  file = "/Data/Fig_2_Distribution/TE_INFO_count_PerFamilyPerPlant.txt", 
  col.names = T,
  sep = "\t"
)

df_count %>%
  mutate(
    OUT = paste0(TOTAL, ";", STRUCTURAL, ";", HOMOLOGY)
  ) %>%
  select(
    TE_CLASS, TE_NAME, PLANT_NAME, OUT
  ) %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    TE_CLASS + TE_NAME ~ PLANT_NAME, value.var = "OUT", fill = "0;0;0"
  ) %>%
  select(
    TE_CLASS, TE_NAME, Chiifu, Z1, CXA, MIZ, OIB, PCA, TUA, HDEM, Korso, OX_heart, RCBO, JZSv2
  ) %>%
  fwrite(
    file = "/Data/Fig_2_Distribution/TE_INFO_count_PerFamilyPerPlant_ForExcel.txt", 
    col.names = T, sep = "\t"
  )
  

glimpse(df_count)

#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------
df_upset_temp <- df_upset %>%
  data.table::melt.data.table(
    id.vars = "TE_Family", 
    variable.name = "PLANT_NAME", 
    value.name = "PAV"
  ) %>%
  left_join(
    x = ., 
    y = df_plant_species, 
    by = c("PLANT_NAME")
  ) 
  
glimpse(df_upset_temp)

# 260 total conserve ------------------
df_total_share <- df_upset_temp %>%
  group_by(TE_Family) %>%
  mutate(
    TOT = sum(PAV)
  ) %>%
  ungroup() %>%
  filter(
    TOT %in% 12
  ) %>%
  select(TE_Family) %>%
  unique() %>%
  mutate(
    TE_NAME = TE_Family, 
    TAG = "Commonly_found_in_Brassica"
  ) %>%
  glimpse()


# not 260 total conserve ============
df_upset_temp %>%
  group_by(TE_Family) %>%
  mutate(
    TOT = sum(PAV)
  ) %>%
  ungroup() %>%
  filter(
    !(TOT %in% 12)
  ) %>%
  group_by(
    TE_Family, PLANT_SPECIES
  ) %>%
  mutate(
    TOT = sum(PAV)
  ) %>%
  ungroup() %>%
  filter(
    TOT >= 2, 
    PAV %in% 1
  ) %>%
  select(
    TE_Family, PLANT_SPECIES, TOT
  ) %>%
  unique() %>%
  reshape2::dcast(
    TE_Family ~ PLANT_SPECIES, value.var = "TOT", fill = 0
  ) -> df_not_conserve
  
glimpse(df_not_conserve)

df_rapa <- df_not_conserve %>%
  filter(
    B_oleracea %in% 0, 
    B_rapa > 0
  ) %>%
  select(
    TE_Family
  ) %>%
  mutate(
    TE_NAME = TE_Family, 
    TAG = "B_rapa_specific"
  )
  
df_oleracea <- df_not_conserve %>%
  filter(
    B_oleracea > 0, 
    B_rapa %in% 0
  ) %>%
  select(
    TE_Family
  ) %>%
  mutate(
    TE_NAME = TE_Family, 
    TAG = "B_oleracea_specific"
  )

df_other <- df_not_conserve %>%
  filter(
    !(TE_Family %in% df_rapa$TE_Family)
  ) %>%
  filter(
    !(TE_Family %in% df_oleracea$TE_Family)
  ) %>%
  select(
    TE_Family
  ) %>%
  mutate(
    TE_NAME = TE_Family, 
    TAG = "Commonly_found_in_Brassica"
  )
  
glimpse(df_rapa)
glimpse(df_oleracea)
glimpse(df_other)
glimpse(df_total_share)
#dim(df_not_conserve)

rm(df_cat)
#df_cat <- rbind(df_rapa, df_oleracea, df_other)
#df_cat <- rbind(df_rapa, df_oleracea, df_total_share)
df_cat <- rbind(df_rapa, df_oleracea, df_other, df_total_share)
glimpse(df_cat)
glimpse(unique(df_cat))


df_age <- fread(
  file = "/Data/Supplement_Data/Insert_Age_Brassica_AA_CC_TEIDChange.txt", 
  header = T, sep = "\t"
)
glimpse(df_age)

rm(df_age_cat)
df_age_cat <- df_age %>%
  left_join(
    x = ., 
    y = df_cat, 
    by = "TE_NAME"
  ) %>%
  filter(
    !is.na(TAG)
  ) %>%
  filter(
    !is.na(INSERT_AGE)
  ) %>%
  glimpse()

df_age_cat %>%
  filter(
    INSERT_AGE > 0
  ) %>%
  filter(
    str_detect(TE_CLASS, pattern = "^LTR")
  ) %>%
  filter(
    str_detect(TE_CLASS, pattern = "LTR/unknown", negate = T)
  ) %>%
  ggplot(
    aes(
      x = INSERT_AGE/1e+6,
      #y = TAG
    )
  )+
  # geom_histogram(
  #   binwidth = 0.1
  # )+
  geom_density(
    fill = "#cccccc"
  )+
  #ggridges::geom_density_ridges()+
  scale_x_continuous(
    expand = c(0, 0), 
    #limits = c(0, NA)
  )+
  scale_y_continuous(
    expand = c(0, 0),
    #limits = c(0, NA)
  )+
  coord_cartesian(
    xlim = c(0, NA),
    ylim = c(0, NA)
  )+
  facet_grid(
    rows = vars(TAG),
    cols = vars(TE_CLASS),
    scales = "free_y", 
    switch = "y"
  )+
  theme_classic()+
  theme(
    strip.background.x = element_blank(),
    strip.background.y = element_blank(), 
    strip.text.y.left = element_text(
      angle = 0, 
      hjust = 1
    ),
    strip.placement = "outside", 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  )+
  labs(
    x = "MYA"
  ) -> plot_share_age

plot_share_age

cowplot::ggsave2(
  plot = plot_share_age, 
  filename = "/Supplement_Data/SuppFig4_Age_share260_vs_AA_CC_specific.pdf", 
  width = 7.2, 
  height = 4
)
