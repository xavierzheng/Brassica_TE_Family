#================================
# AIM
#	Analysis the mean identity of each TE insertion within TE family
#
#===============================

gc()

# setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 另一種畫法, facet, 用這個 ------------------------------------------------------------------

NAME_LIST <- c("Copia", "Gypsy", "unknown", "Helitron", "DTA", "DTC", "DTH", "DTM", "DTT")

print("# 這才是真實的count -----")
df_count <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt.gz", 
  header = T, sep = "\t"
) %>% filter(
  TE_METHOD %in% "structural"
) %>%
  mutate(
    TE_CLASS = as.character(factor(
      TE_CLASS, 
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", 
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"), 
      labels = NAME_LIST
    ))
  )

glimpse(df_count)
unique(df_count$TE_CLASS)

temp_list <- lapply(NAME_LIST, FUN = function(NAM){
  
  df <- fread(
    file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEWithinFamily_Relaxed_", NAM, ".txt"), 
    header = T, sep = "\t"
  )
  
  # chose top 10
  df_count_temp <- df_count %>%
    filter(
      TE_CLASS %in% NAM
    ) %>%
    summarise(
      .by = TE_NAME, 
      COUNT = n()
    ) %>%
    arrange(
      desc(COUNT)
    ) %>%
    slice_head(
      n = 10
    )
  
  # df_count <- df[, .N, by = "S_FAM"][
  #   order(-N), head(.SD, 10)
  # ]
  #str(df_count)
  
  # used all
  # df_count <- df[, .N, by = "S_FAM"][
  #   order(-N), 
  # ]
  
  
  # name_list <- sort(unique(df$S_FAM))
  name_list <- df_count_temp$TE_NAME
  
  keep_list <- c("query_id", "Q_FAM", "S_FAM", "AVE")
  
  df_draw <- df %>%
    filter(
      S_FAM %in% name_list
    ) %>%
    mutate(
      Q_ID = str_extract(query_id, "ID=[A-Za-z0-9_]{1,}"), 
      Q_POS = str_split_fixed(query_id, "::", 2)[,2], 
      Q_CHROM = str_split_fixed(Q_POS, ":", 2)[,1], 
      Q_PLANT = str_split_fixed(Q_CHROM, "_", 2)[,1], 
      Q_SPECIES = ifelse(
        str_detect(Q_PLANT, "^D"), 
        "AA", "CC"
      )
    ) %>%
    mutate(
      S_FAM = factor(
        S_FAM, 
        levels = name_list, 
        labels = name_list
      ), 
      SUPFAM = NAM
    ) 
  
  return(df_draw)
  
})

df_draw <- bind_rows(temp_list)

glimpse(df_draw)

print("# statistic analysis: 1. remove the family that only found in one species")
df_shareness <- fread(
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt", 
  header = T, sep = "\t"
) 
colnames(df_shareness) <- c("SUPFAM", "S_FAM", "Share_or_Unique")
glimpse(df_shareness)  

df_stat1 <- df_draw %>%
  select(
    SUPFAM, Q_SPECIES, S_FAM
  ) %>%
  unique() %>%
  summarise(
    .by = c(SUPFAM ,S_FAM), 
    COUNT = n()
  ) %>%
  filter(
    COUNT == 1
  ) %>%
  mutate(
    SIG = NA, 
    PVAL = NA,
    S_FAM = as.character(S_FAM)
  ) %>%
  left_join(
    x = .,
    y = df_shareness %>% select(-SUPFAM), 
    by = "S_FAM"
  ) %>%
  dplyr::rename(
    Prevelance = Share_or_Unique 
  ) %>%
  select(
    SUPFAM, S_FAM, PVAL, SIG, Prevelance
  )

print("# using wilcox t test ==========================")
df_stat2 <- df_draw %>%
  select(
    SUPFAM, Q_SPECIES, S_FAM, AVE
  ) %>%
  filter(
    !(S_FAM %in% df_stat1$S_FAM) # remove TE family found only in one species
  ) %>%
  mutate(
    Q_SPECIES = factor(
      Q_SPECIES, 
      levels = c("AA", "CC"), 
      labels = c("AA", "CC")
    )
  ) %>%
  nest_by(
    SUPFAM , S_FAM
  ) %>%
  summarise(
    PVAL = wilcox.test(AVE ~ Q_SPECIES, data = data, exact = F)$p.value
  ) %>%
  ungroup() %>%
  mutate(
    SIG = ifelse(
      PVAL < 0.001, 
      "*", NA
    ), 
    S_FAM = as.character(S_FAM)
  ) %>%
  select(
    SUPFAM, S_FAM, PVAL ,SIG
  )

# obtain the describe statistic --------------
df_stat_desc <- df_draw %>%
  select(
    SUPFAM, Q_SPECIES, S_FAM, AVE
  ) %>%
  filter(
    !(S_FAM %in% df_stat1$S_FAM) # remove TE family found only in one species
  ) %>%
  mutate(
    Q_SPECIES = factor(
      Q_SPECIES, 
      levels = c("AA", "CC"), 
      labels = c("AA", "CC")
    )
  ) %>%
  summarise(
    .by = c(SUPFAM, Q_SPECIES, S_FAM), 
    MEAN = mean(AVE)
  ) %>%
  reshape2::dcast(
    SUPFAM + S_FAM ~ Q_SPECIES, value.var = "MEAN", fill = 0
  ) %>%
  mutate(
    S_FAM = as.character(S_FAM), 
    Prevelance = AA - CC
  ) %>%
  left_join(
    x = ., 
    y = df_stat2, 
    by = c("SUPFAM", "S_FAM")
  ) 

# joint all info -------------------------
df_stat <- df_stat_desc %>%
  select(SUPFAM, S_FAM, PVAL, SIG, Prevelance) %>%
  rbind(df_stat1, .) %>%
  mutate(
    S_FAM = factor(
      S_FAM, 
      levels = levels(df_draw$S_FAM), 
      labels = levels(df_draw$S_FAM)
    ), 
    SUPFAM = factor(
      SUPFAM, 
      levels = NAME_LIST, 
      labels = NAME_LIST
    ), 
    FDR = p.adjust(PVAL, method = "fdr"), 
    SIG_FDR = ifelse(
      FDR < 0.001, 
      "*", NA
    )
  )

fwrite(
  df_stat, 
  file = "/Data/Fig_3_Nucleotide_Identity/AgeComparison_Top10_statistic_Relax.txt",
  col.names = T, sep = "\t"
)


df_draw <- df_draw %>%
  mutate(
    SUPFAM = factor(
      SUPFAM, 
      levels = NAME_LIST, 
      labels = NAME_LIST
    )
  ) %>%
  glimpse()


ggplot()+
  geom_jitter(
    data = df_draw, 
    mapping = aes(
      x = S_FAM, 
      y = AVE, 
      color = Q_SPECIES, 
      fill = Q_SPECIES
    ),
    position = position_jitterdodge(
      jitter.width = 0.1, 
      jitter.height = 0, 
      dodge.width = 0.5,
      seed = 123
    ),
    size = 0.1, 
    shape = 16, 
    alpha = 0.1
  )+
  geom_boxplot(
    data = df_draw, 
    mapping = aes(
      x = S_FAM, 
      y = AVE, 
      color = Q_SPECIES,
      fill = Q_SPECIES
    ),
    outlier.shape = NA,
    width = 0.5,
    alpha = 0, 
    position = position_dodge2()
  )+
  geom_text(
    data = df_stat, 
    mapping = aes(
      x = S_FAM, 
      y = 94, 
      label = SIG_FDR
    ),
    color = "red", 
    size = 6
  )+
  scale_fill_manual(
    breaks = c("AA", "CC"), 
    values = c("#998EC3", "#F1A340")
  )+
  scale_color_manual(
    breaks = c("AA", "CC"), 
    values = c("#998EC3", "#F1A340")
  )+
  coord_cartesian(
    ylim = c(10, 105), expand = F
  )+
  facet_wrap(
    ~ SUPFAM, nrow = 3, ncol = 3, scales = "free_x"
  )+
  theme_cowplot()+
  theme(
    axis.title.y = element_text(
      size = 7
    ), 
    axis.title.x = element_blank(),
    axis.text.y = element_text(
      size = 6
    ), 
    axis.text.x = element_text(
      angle = 90, 
      hjust = 1, 
      vjust = 0.5,
      size = 5
    ), 
    legend.text = element_text(
      size = 6
    ),
    legend.title = element_text(
      size = 7
    ),
    legend.position = "bottom", 
    legend.margin = margin(t = 1, r = 2, b = 1, l = 2, unit = "pt"),
    plot.margin = margin(t = 1, r = 2, b = 1, l = 2, unit = "pt"), 
    strip.background = element_blank(), 
    strip.text = element_text(
      size = 7
    )
  )+
  labs(
    color = "", 
    fill = "", 
    y = "Mean identity"
  )

cowplot::ggsave2(
  filename = "/Fig_3_Nucleotide_Identity/boxplot_PerTEWithinFamily_AllTEFamily_Relax_facet.pdf", 
  width = 7.2, 
  height = 5.3
)
