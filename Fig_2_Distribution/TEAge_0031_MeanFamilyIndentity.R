#------------------------------
# AIM
#	Estimate the mean average identity
#
#------------------------------


# setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN")


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))

setDTthreads(threads = 24)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# VERSION <- as.character(x = args[1])
# INPUT_FILE <- as.character(x = args[2])

# VERSION <- "Relaxed"
# INPUT_FILE <- "DTA"

NAME_LIST <- c("Copia", "Gypsy", "unknown", "Helitron", "DTA", "DTC", "DTH", "DTM", "DTT")

temp_list <- lapply(NAME_LIST, FUN = function(NAM){
  
  VERSION <- "Relaxed"
  
  df <- fread(
    file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEWithinFamily_", VERSION, "_", NAM, ".txt"),
    header = T, sep = "\t"
  )
  
  df_out <- df[
    , .(mean(AVE)), by = .(Q_FAM, S_FAM)
  ][
    order(Q_FAM),.(Q_FAM, V1)
  ]
  
  colnames(df_out) <- c("TE_FAM", "AVE_IDENT")
  
  return(df_out)
})

df_out <- bind_rows(temp_list)

str(df_out)

fwrite(
  df_out,
  file = "/Data/Fig_2_Distribution/PerFamily_Relaxed_AveIdentity.txt",
  col.names = T, sep = "\t"
)

print("# Combine with upset plot ---------------------")
df_upset <- fread(
  file = "/Data/Fig_2_Distribution/StructuralTE_UpsetCount_FamilyLevel_Excel.txt", 
  header = T, sep = "\t"
)

str(df_upset)

df_draw_point <- left_join(
  x = df_upset,
  y = df_out,
  by = "TE_FAM"
) %>%
  mutate(
    TE_CLASS = factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
    )
  )

df_draw_ave <- df_draw_point %>%
  summarise(
    .by = CAT,
    AVE = mean(AVE_IDENT, na.rm = T),
    SD = sd(AVE_IDENT, na.rm = T)
  ) 


ggplot()+
  geom_jitter(
    data = df_draw_point, 
    mapping = aes(
      x = CAT, 
      y = AVE_IDENT, 
      color = TE_CLASS
    ), 
    alpha = 0.5,
    size = 0.2,
    position = position_jitter(height = 0, width = 0.3, seed = 134)
  )+
  geom_errorbar(
    data = df_draw_ave,
    mapping = aes(
      x = CAT, 
      ymin = AVE-SD,
      ymax = AVE+SD
    ),
    color = "black",
    width = 0.5
  )+
  geom_point(
    data = df_draw_ave,
    mapping = aes(
      x = CAT, 
      y = AVE,
    ),
    color = "black"
  )+
  geom_hline(
    yintercept = 80,
    color = "black",
    linewidth = 0.5
  )+
  scale_x_continuous(
    limits = c(0.5, 20.5),
    expand = expansion(mult = 0, add = 0)
  )+
  scale_y_continuous(
    breaks = c(20, 40, 60, 80, 100)
  )+
  scale_color_manual(
    breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
               "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
    values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8",
               "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
  )+
  theme_cowplot()+
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
  )+
  labs(
    y = "Mean identity (%)",
    #fill = "TE superfamily"
  )
ggsave2(
  filename = "/Data/Fig_2_Distribution/Upset_plot_bottom_AveIdent.pdf",
  width = 4, 
  height = 1.6
)  
