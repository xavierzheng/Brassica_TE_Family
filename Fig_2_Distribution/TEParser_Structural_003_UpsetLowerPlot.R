#-------------------
# AIM
#   Generate the lower plot of upset to describe the superfamily 
#   abundance
#
#-------------------


# setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine")

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(ggsci))


setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

df <- fread(
  file = "/Data/Fig_2_Distribution/StructuralTE_UpsetCount_FamilyLevel_Excel.txt", 
  header = T, sep = "\t"
)

# make sure label is total identical to upset plot
# 16 -> 19
# 17 -> 18
# 18 -> 17
# 19 -> 16


glimpse(df)

df %>%
  arrange(
    CAT
  ) %>%
  mutate(
    CAT = factor(
      CAT,
      levels = as.character(1:20), 
      labels = as.character(1:20)
    )
  ) %>%
  summarise(
    .by = c(CAT, TE_CLASS), 
    COUNT = n()
  ) %>%
  group_by(
    CAT
  ) %>%
  mutate(
    PERC = COUNT/sum(COUNT)*100, # relative to the indicated group
    TE_CLASS = factor(
      TE_CLASS, 
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
    )
  ) %>%
  ggplot(
    aes(
      x = CAT, 
      y = PERC, 
      fill = TE_CLASS
    )
  )+
  geom_col(
    position = "stack"
  )+
  scale_fill_manual(
    breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
               "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
    values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8",
               "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
  )+
  theme_cowplot()+
  coord_cartesian(
    expand = F
  )+
  theme(
    axis.line.x = element_blank(), 
    axis.ticks.length.x = unit(0, "pt"), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.y = element_text(
      size = 7
    ), 
    axis.text.y = element_text(
      size = 6
    ), 
    legend.position = "none"
    # legend.title = element_text(
    #   size = 7
    # ),
    # legend.text = element_text(
    #   size = 6
    # )
  )+
  labs(
    y = "Relative abundance\nTE family (%)", 
    #fill = "TE superfamily"
  )
ggsave(
  filename = "/Fig_2_Distribution/upset_lower_superfamily_abundance.pdf", 
  width = 4, 
  height = 1
)



