#===================================
#	Aim:
#	summary TFBS frequency, which TE has higher probability to has TFBS
#
#===================================

setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN/TF_FIMO")

suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(UpSetR))

setDTthreads(threads = 12)
options(tidyverse.quiet = TRUE)
options(dplyr.summarise.inform = FALSE)

print("# filter q value < 1e-5 ==================")

name_list <- read.table(
  file = "name.list", 
  header = F, sep = "\t"
) %>% pull()

temp_list <- lapply(name_list, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0("FIMO_Qfiltered_", NAM,".txt"), 
    header = T, sep = "\t", nThread = 12
  )
  df_temp[,"TE_SUPERFAM" := NAM]
  
  return(df_temp)
})

df <- bind_rows(temp_list)

glimpse(df)

df_info <- df %>%
  mutate(
    TE_FAM = str_extract(sequence_name, pattern = "Name=([A-Za-z0-9]*)") %>%
      str_remove(
        ., "Name="
      ), 
    TE_ID = str_split_fixed(sequence_name, "::", 2)[,2], 
    PLANT = str_split_fixed(CHROM, "_", 2)[,1], 
    SPECIES = ifelse(
      str_detect(PLANT, "^I"), 
      "B_oleracea", "B_rapa"
    )
  ) %>%
  as.data.table() %>%
  glimpse()

glimpse(df_info)


table(df_info$motif_alt_id)

df_plot <- df_info[
  order(TE_SUPERFAM, TE_FAM), .N, keyby = .(TE_SUPERFAM, TE_FAM, motif_alt_id)
]

glimpse(df_plot)
print("# Normalize =================")

print("## make sure this include all family ===========")
df_all_fam <- fread(
  file = "Whole_TE_family.txt", 
  header = F, sep = "\t", col.names = c("TE_SUPERFAM", "TE_FAM", "Total_TE_INSERT")
)
glimpse(df_all_fam)


df_temp <- df_plot %>%
  data.table::dcast.data.table(
    TE_SUPERFAM + TE_FAM ~ motif_alt_id, 
    value.var = "N", fill = 0
  ) %>%
  left_join(
    x = df_all_fam, 
    y = ., 
    by = c("TE_SUPERFAM", "TE_FAM")
  ) %>%
  data.table::melt.data.table(
    id.vars = c("TE_SUPERFAM", "TE_FAM", "Total_TE_INSERT"), 
    variable.name = "motif_alt_id", 
    variable.factor = F, 
    value.name = "N"
  ) %>%
  mutate(
    N = ifelse(
      is.na(N), 
      0, N
    )
  ) 
  
print("# write out a table =======================")
glimpse(df_temp)
df_temp %>%
  data.table::dcast.data.table(
    TE_SUPERFAM + TE_FAM + Total_TE_INSERT ~ motif_alt_id, 
    value.var = "N", fill = 0
  ) %>%
  fwrite(
    file = "TFBS_count_per_TE_family.txt", 
    col.names = T, sep = "\t"
  )

print("# transform the number into category data format ==================")
glimpse(df_temp)
df_plot_category <- df_temp %>%
  mutate(
    COUNT_CAT = ifelse(
      N == 0, 
      "0", ifelse(
        N <= 10, 
        "1-10", ifelse(
          N <= 100, 
          "11-100", ifelse(
            N <= 1000, 
            "101-1000", ">1000"
          )
        )
      )
    )
  ) %>%
  mutate(
    COUNT_CAT = factor(
      COUNT_CAT, 
      levels = c("0", "1-10", "11-100", "101-1000", ">1000"), 
      labels = c("0", "1-10", "11-100", "101-1000", ">1000")
    )
  ) 

df_plot_category %>%
  ggplot(
    aes(
      x = motif_alt_id, 
      y = TE_FAM, 
      fill = COUNT_CAT
    )
  )+
  geom_tile()+
  scale_fill_brewer(
    palette = "RdPu"
  )+
  facet_wrap(
    ~TE_SUPERFAM, scales = "free_y"
  )+
  theme_cowplot()+
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(
      size = 7
    ), 
    axis.line = element_blank(), 
    axis.ticks = element_blank(), 
    axis.text = element_blank(), 
    # axis.text.y = element_text(
    #   size = 5
    # ),
    axis.title = element_text(
      size = 7
    ), 
    legend.title = element_text(
      size = 7, vjust = 0.8
    ), 
    legend.text = element_text(
      size = 6
    ),
    legend.key.width = unit(5, "pt"),
    legend.key.height = unit(3, "pt"),
    legend.position = "bottom", 
    panel.spacing = unit(0, "pt"), 
    panel.border = element_rect(
      color = "black"
    )
  )+
  labs(
    x = "TFBS", 
    y = "TE family", 
    fill = "No. of TFBS in each TE family  "
  ) -> plot_tile

plot_tile

cowplot::ggsave2(
  plot = plot_tile,
  filename = "TFBS_count_per_TE_family.pdf", 
  width = 7.2, 
  height = 4
)


print("# zoom in to identify several example (unknown0001, unknown0002, ...etc) ================")
glimpse(df_plot_category)

print("# choose the top 5 most abundance unknown LTR family --------")
df_temp2 <- df_plot_category %>%
  filter(
    TE_SUPERFAM %in% "unknown"
  ) %>%
  as.data.table() %>%
  dcast.data.table(
    TE_SUPERFAM + TE_FAM + Total_TE_INSERT ~ motif_alt_id, 
    value.var = "N", fill = "0"
  ) %>%
  arrange(
    desc(Total_TE_INSERT)
  ) %>%
  slice_head(
    n = 5 # choose the top 5 most abundance unknown LTR family
  ) %>%
  melt.data.table(
    id.vars = c("TE_SUPERFAM", "TE_FAM", "Total_TE_INSERT"), 
    variable.name = "motif_alt_id", 
    variable.factor = F, 
    value.name = "N"
  ) %>%
  glimpse()


print("# only the TF with more than 200 binding site are chosen for analysis =-----")
df_temp2 %>%
  group_by(
    TE_FAM, motif_alt_id
  ) %>%
  summarise(
    SUM = sum(N)
  ) %>%
  ungroup() %>%
  filter(
    SUM >= 200 # only the TF with more than 200 binding site are chosen for analysis
  ) %>%
  select(
    -SUM
  ) %>%
  left_join(
    x = ., 
    y = df_temp2, 
    by = c("TE_FAM", "motif_alt_id")
  ) %>%
  mutate(
    COUNT_CAT = ifelse(
      N == 0, 
      "0", ifelse(
        N <= 10, 
        "1-10", ifelse(
          N <= 100, 
          "11-100", ifelse(
            N <= 1000, 
            "101-1000", ">1000"
          )
        )
      )
    )
  ) %>%
  mutate(
    COUNT_CAT = factor(
      COUNT_CAT, 
      levels = c("0", "1-10", "11-100", "101-1000", ">1000"), 
      labels = c("0", "1-10", "11-100", "101-1000", ">1000")
    ), 
    TFBS = str_split_fixed(motif_alt_id, "\\.", 3)[,3]
  ) %>%
  ggplot(
    aes(
      x = TFBS, 
      y = TE_FAM, 
      fill = COUNT_CAT, 
      label = N
    )
  )+
  geom_tile()+
  geom_text(
    size = 2.5
  )+
  scale_fill_manual(
    # breaks = c("0", "1-10", "11-100", "101-1000", ">1000"), 
    # values = c("#feebe2", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177")
    breaks = c("101-1000", ">1000"), 
    values = c("#f768a1", "#7a0177")
  )+
  coord_cartesian(
    expand = F
  )+
  theme_cowplot()+
  theme(
    axis.title = element_text(
      size = 7
    ), 
    axis.text.y = element_text(
      size = 6
    ),
    axis.text.x = element_text(
      size = 6, 
      angle = 90, 
      vjust = 0.5, 
      hjust = 1
    ), 
    legend.title = element_text(
      size = 7
    ), 
    legend.text = element_text(
      size = 6
    ), 
    legend.position = "bottom"
  )+
  labs(
    x = "TFBS", 
    y = "TE family", 
    fill = "No. of TFBS in each TE family  "
  ) -> plot_top5

plot_top5

ggsave2(
  plot = plot_top5, 
  filename = "TFBS_count_per_TE_family_LTRunknownTop5.pdf", 
  width = 7.2, 
  height = 3
)


#==================================
# this is used for normalization 
  # group_by(
  #   TE_SUPERFAM, TE_FAM  # group by TE_FAM
  # ) %>%
  # mutate(
  #   Z_COUNT = (N - mean(N))/sd(N)
  # ) %>%
  # # group_by(
  # #   motif_alt_id  # group by TFBS
  # # ) %>%
  # # mutate(
  # #   Z_COUNT = (N - mean(N))/sd(N)
  # # ) %>%
  # ungroup() %>%
  # group_by(
  #   TE_SUPERFAM
  # ) %>%
  # arrange(TE_FAM) %>%
  # ungroup() %>%
  # filter(
  #   TE_SUPERFAM %in% "unknown"
  # ) %>%
  # ggplot(
  #   aes(
  #     x = motif_alt_id, 
  #     y = TE_FAM, 
  #     fill = Z_COUNT
  #   )
  # )+
  # geom_tile()+
  # scale_fill_distiller(
  #   palette = "YlGnBu", direction = -1
  # )+
  # facet_wrap(
  #   ~TE_SUPERFAM, scales = "free_y"
  # )+
  # theme_cowplot()+
  # theme(
  #   strip.background = element_blank(), 
  #   strip.text = element_text(
  #     size = 7
  #   ), 
  #   axis.line = element_blank(), 
  #   axis.ticks = element_blank(), 
  #   axis.text = element_blank(), 
  #   # axis.text.y = element_text(
  #   #   size = 5
  #   # ),
  #   axis.title = element_text(
  #     size = 7
  #   ), 
  #   legend.title = element_text(
  #     size = 7, vjust = 0.7
  #   ), 
  #   legend.text = element_text(
  #     size = 5
  #   ),
  #   legend.key.width = unit(20, "pt"),
  #   legend.key.height = unit(2, "pt"),
  #   legend.position = "bottom", 
  #   panel.spacing = unit(5, "pt"), 
  #   panel.border = element_rect(
  #     color = "black"
  #   )
  # )+
  # labs(
  #   x = "TFBS", 
  #   y = "TE family", 
  #   fill = "Z-normalized count within TE family"
  # )
  # 
