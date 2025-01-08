#-----------------------------------------------
# AIM
#   Re-draw upset plot, using smaller region
#
#-----------------------------------------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
library(ggupset)
library(patchwork)

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 讓arrow可以多核心讀寫
arrow::set_cpu_count(12)
arrow::set_io_thread_count(12)
arrow::cpu_count()

# prepare function-------------------------------------------------
print("# prepare plot function ----------------------------------")
upset_superfamily <- function(TAB_LOC, ID_LIST, NAME_LIST){
  
  # read data
  df_TE <- read.table(
    file = TAB_LOC,
    header = T, sep = "\t"
  ) %>%
    tibble::rownames_to_column(
      var = "KEY"
    ) %>%
    mutate(
      PAIR = str_split_fixed(KEY, " ", 2)[,1],
      TE_CLASS = str_split_fixed(KEY, " ", 2)[,2],
      GENE1 = str_split_fixed(PAIR, "-", 2)[,1],
      GENE2 = str_split_fixed(PAIR, "-", 2)[,2],
      GENE1 = str_remove(GENE1, "\\.3\\.1C$"),
      GENE2 = str_remove(GENE2, "\\.3\\.1C$")
    )
  
  # reformat for ggupset
  df_draw <- df_TE %>%
    as.data.table() %>%
    melt.data.table(
      id.vars = c("KEY", "PAIR", "TE_CLASS", "GENE1", "GENE2"),
      variable.name = "PLANT",
      variable.factor = F,
      value.name = "PAV"
    ) %>%
    mutate(
      PAV = ifelse(
        PAV >= 1,
        1, 0
      ),
      PLANT = factor(
        PLANT,
        levels = ID_LIST,
        labels = NAME_LIST
      ),
      PLANT = as.character(PLANT),
      TE_CLASS = factor(
        TE_CLASS,
        levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", 
                   "TIR/DTM", "TIR/DTA", "TIR/DTH", "TIR/DTT"),
        labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", 
                   "TIR/DTM", "DNA/Other", "DNA/Other", "DNA/Other")
      )
    ) %>%
    filter(
      PAV > 0
    ) %>%
    summarise(
      .by = c(KEY, PAIR, TE_CLASS, GENE1, GENE2),
      LIST_PLANT = list(PLANT)
    )
  
  # counting table ----
  df_temp <- df_draw %>%
    summarise(
      .by = LIST_PLANT,
      COUNT = n()
    ) %>%
    arrange(
      desc(COUNT)
    ) %>%
    slice_head(
      n = 10
    )
  
  # drawing --
  ggplot()+
    geom_bar(
      data = df_draw,
      mapping = aes(
        x = LIST_PLANT,
        fill = TE_CLASS
      ),
      color = "black",
      linewidth = 0.5
    )+
    geom_text(
      data = df_temp,
      mapping = aes(
        x = LIST_PLANT,
        y = COUNT,
        label = COUNT
      ),
      vjust=-1,
      size = 6*0.35
    )+
    ggupset::scale_x_upset( # 利用scale_x_upset來產生upset需要的點和線
      order_by = "freq",
      n_sets = length(NAME_LIST),
      sets = NAME_LIST,
      n_intersections = 10
    )+
    scale_fill_manual(
      breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", "TIR/DTM","DNA/Other"),
      #values = c("#a63603", "#fd8d3c", "#fdd0a2", "#54278f", "#9e9ac8", "#dadaeb")
      values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8", "#6baed6","#a1d99b", "white")
    )+
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.25))
    )+
    theme_cowplot()+
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(
        size = 0
      ),
      legend.text = element_text(
        size = 6
      ),
      legend.key.height = unit(10, "pt"),
      legend.key.width = unit(10, "pt"),
      axis.title.y = element_text(
        size = 7
      ),
      axis.text.y = element_text(
        size = 6
      ),
      plot.margin = margin(t = 0, r = 0, b = 2, l = 6, unit = "pt")
    )+
    labs(
      x = "",
      y = "No. syntenic pair\nand TE superfamily",
      fill = ""
    ) -> plot_CLAS
  
  return(plot_CLAS)
  
}


upset_family <- function(TAB_LOC, ID_LIST, NAME_LIST){
  
  # read data ---
  df_TE <- read.table(
    file = TAB_LOC,
    header = T, sep = "\t"
  ) %>%
    tibble::rownames_to_column(
      var = "KEY"
    ) %>%
    mutate(
      PAIR = str_split_fixed(KEY, " ", 2)[,1],
      TE_FAM = str_split_fixed(KEY, " ", 2)[,2],
      GENE1 = str_split_fixed(PAIR, "-", 2)[,1],
      GENE2 = str_split_fixed(PAIR, "-", 2)[,2],
      GENE1 = str_remove(GENE1, "\\.3\\.1C$"),
      GENE2 = str_remove(GENE2, "\\.3\\.1C$")
    )
  
  # reformat for ggupset ------
  df_draw <- df_TE %>%
    as.data.table() %>%
    melt.data.table(
      id.vars = c("KEY", "PAIR", "TE_FAM", "GENE1", "GENE2"),
      variable.name = "PLANT",
      variable.factor = F,
      value.name = "PAV"
    ) %>%
    mutate(
      PAV = ifelse(
        PAV >= 1,
        1, 0
      ),
      PLANT = factor(
        PLANT,
        levels = ID_LIST,
        labels = NAME_LIST
      ),
      PLANT = as.character(PLANT)
    ) %>%
    filter(
      PAV > 0
    ) %>%
    summarise(
      .by = c(KEY, PAIR, TE_FAM, GENE1, GENE2),
      LIST_PLANT = list(PLANT)
    )
  
  # count table -------------
  df_temp <- df_draw %>%
    summarise(
      .by = LIST_PLANT,
      COUNT = n(),
      N_TE_FAM = length(unique(unlist(list(TE_FAM))))
    ) %>%
    mutate(
      DENS_TE_FAM = round(COUNT/N_TE_FAM, digits = 2)
    ) %>%
    arrange(
      desc(COUNT)
    ) %>%
    slice_head(
      n = 10
    )
  
  # drawing ----------
  ggplot()+
    geom_bar(
      data = df_draw,
      mapping = aes(
        x = LIST_PLANT
      ),
      color = "black",
      fill = "#d9d9d9",
      linewidth = 0.5
    )+
    geom_text(
      data = df_temp,
      mapping = aes(
        x = LIST_PLANT,
        y = COUNT,
        label = paste0(COUNT, "\n(",DENS_TE_FAM, ")")
      ),
      vjust = -0.1,
      size = 6*0.35
    )+
    ggupset::scale_x_upset(
      n_intersections = 10,
      sets = NAME_LIST,
      n_sets = length(NAME_LIST)
    )+
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.35))
    )+
    theme_cowplot()+
    theme(
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(
        size = 0
      ),
      axis.title.y = element_text(
        size = 7
      ),
      axis.text.y = element_text(
        size = 6
      ),
      plot.margin = margin(t = 1, r = 0, b = 2, l = 6, unit = "pt")
    )+
    labs(
      x = "",
      y = "No. syntenic pair\nand TE family"
    ) -> plot_ret
  
  return(plot_ret)
  
}


# drawing ----------------------------------------------

print("# drawing -----------------------------------------------")
plot_CLA_AA <- upset_superfamily(
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_AA.txt.gz",
  ID_LIST = c("DA", "DB", "DN", "DP", "DQ", "DR", "DU"),
  NAME_LIST = c("Chiifu", "Z1", "CXA", "MIZ", "OIB", "PCA", "TUA")
)

plot_CLA_CC <- upset_superfamily(
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_CC.txt.gz",
  ID_LIST = c("IC", "IG", "IH", "II", "IJ"),
  NAME_LIST = c("HDEM", "Korso", "OX_heart", "BoA", "JZSv2")
)

plot_FAM_AA <- upset_family(
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_AA.txt.gz",
  ID_LIST = c("DA", "DB", "DN", "DP", "DQ", "DR", "DU"),
  NAME_LIST = c("Chiifu", "Z1", "CXA", "MIZ", "OIB", "PCA", "TUA")
)

plot_FAM_CC <- upset_family(
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_CC.txt.gz",
  ID_LIST = c("IC", "IG", "IH", "II", "IJ"),
  NAME_LIST = c("HDEM", "Korso", "OX_heart", "BoA", "JZSv2")
)

plot_save <- plot_CLA_AA + plot_CLA_CC + plot_FAM_AA + plot_FAM_CC + patchwork::plot_layout(
  ncol = 2, nrow = 2, guides = "collect"
)&theme(legend.position = "right")

plot_save

ggsave2(
  filename = "/Fig_7_Synteny_TE/upset_SyntenyGeneAndTEInsertion.pdf",
  plot = plot_save,
  width = 7.2,
  height = 4
)



# define color code----------------------
# scale_fill_manual(
#   breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
#              "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
#   values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8",
#              "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
# )
  