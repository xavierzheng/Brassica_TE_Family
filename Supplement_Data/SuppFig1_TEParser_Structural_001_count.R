#============================
# AIM:
#	Focus on the prevelence of No. of structural TE, and length in B. rapa and B. oleraces
#	in superfamily and family level
#
#===========================

setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

df_info <- data.frame(
  SPECIES = c(rep("AA", 7), rep("CC", 5)), 
  ID = c("DA", "DB", "DN", "DP", "DQ", "DR", "DU", "IC", "IG", "IH", "II", "IJ"), 
  FULL_NAME = c("Chiffu", "Z1", "CXA", "MIZ", "OIB", "PCA", "TUA", "HDEM", "Korso", "OX", "BoA", "JZSv2")
) %>%
  mutate(
    PREFIX = paste0(ID, ".", SPECIES)
  )

temp_list <- lapply(df_info$PREFIX, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0(NAM, ".TE.InChr.IntactSingleLine.sort.gff3"), 
    header = F, sep = "\t"
  )
  
  colnames(df_temp) <- c("CHROM", "SOURCE", "TYPE", "GFF_START", "GFF_END", "SCORE", "STRAND", "PHASE", "ATTR")
  
  df_temp$PREFIX <- NAM
  
  df_ret <- df_temp %>%
    mutate(
      TE_ID = str_extract(ATTR, "ID=([A-Za-z0-9_]{1,})") %>% str_remove(., "ID="),
      TE_FAM = str_extract(ATTR, "Name=([A-za-z0-9_]{1,})") %>% str_remove(., "Name="), 
      TE_CLASS = str_extract(ATTR, "Classification=([A-za-z0-9/]{1,})") %>% str_remove(., "Classification="), 
      LEN = GFF_END - GFF_START +1
    )
  
  return(df_ret)
})

df <- bind_rows(temp_list)

df <- left_join(
  x = df, 
  y = df_info, 
  by = "PREFIX"
)

glimpse(df)


print("# estimate chromosome length ---------------------")
temp_chr <- lapply(df_info$PREFIX, FUN = function(NAM){
  
  df_temp <- fread(
    file = paste0("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/00_cleanGFF/", NAM, ".genome.fa.fai"), 
    header = F, sep = "\t"
  ) %>%
    filter(
      str_detect(V1, "^[A-Z]{2,2}_[AC][0-9]")
    ) %>%
    select(
      V1, V2
    )
  
  colnames(df_temp) <- c("CHROM", "LEN")
  df_temp$PREFIX <- NAM
  
  return(df_temp)
  
})

df_chr <- bind_rows(temp_chr)
df_chr <- df_chr %>%
  summarise(
    .by = PREFIX, 
    SUM_CHR_LEN = sum(LEN)
  ) %>%
  left_join(
    x = ., 
    y = df_info, 
    by = "PREFIX"
  )

fwrite(
  df_chr, 
  file = "WholeChromosome_size.txt", 
  col.names = T, sep = "\t"
)

print("# read share-uniqueness ------------")
df_uniq <- fread(
  file = "/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/04_IntactSingleLineAndHomoTE/TE_family_uniqueness.txt", 
  header = T, sep = '\t'
)
colnames(df_uniq) <- c("TE_CLASS", "TE_FAM", "Share_or_Unique")
glimpse(df_uniq)

print("# estimate TE superfamily level ----------")
df_plot <- df %>%
  summarise(
    .by = c(PREFIX, TE_CLASS), 
    COUNT = n(), 
    SUM_LEN = sum(LEN)
  ) %>%
  left_join(
    x = .,
    y = df_info,
    by = "PREFIX"
  ) 

glimpse(df_plot)
fwrite(
  df_plot, 
  file = "StructuralTE_TotalCount_SuperFamilyLevel.txt", 
  col.names = T, sep = "\t"
)

df_excel <- df_plot %>%
  mutate(
    FULL_NAME = factor(
      FULL_NAME, 
      levels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX", "BoA", "JZSv2"), 
      labels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2")
    ),
    TE_CLASS = base::factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")
    ), 
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. oleracea")
    )
  )

df_excel %>%
  reshape2::dcast(
    TE_CLASS ~ FULL_NAME, 
    value.var = "COUNT"
  ) %>%
  fwrite(
    file = "StructuralTE_TotalCount_SuperFamilyLevel_Excel.txt", 
    col.names = T, sep = "\t"
  )

df_excel %>%
  reshape2::dcast(
    TE_CLASS ~ FULL_NAME, 
    value.var = "SUM_LEN"
  ) %>%
  fwrite(
    file = "StructuralTE_TotalLength_SuperFamilyLevel_Excel.txt", 
    col.names = T, sep = "\t"
  )

print("# estimate TE family level ----------------")
df_fam <- df %>%
  summarise(
    .by = c(PREFIX, TE_CLASS ,TE_FAM), 
    COUNT = n(), 
    SUM_LEN = sum(LEN)
  ) %>%
  left_join(
    x = .,
    y = df_info,
    by = "PREFIX"
  ) %>%
  left_join(
    x = .,
    y = df_uniq, 
    by = c("TE_CLASS", "TE_FAM")
  ) %>%
  mutate(
    FULL_NAME = factor(
      FULL_NAME, 
      levels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX", "BoA", "JZSv2"), 
      labels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2")
    ),
    TE_CLASS = base::factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")
    ), 
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. oleracea")
    )
  )
  

glimpse(df_fam)

df_fam %>%
  fwrite(
    file = "StructuralTE_TotalCount_FamilyLevel.txt", 
    col.names = T, sep = "\t"
  )

df_fam %>%
  arrange(
    TE_CLASS, TE_FAM
  ) %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    TE_CLASS + TE_FAM + Share_or_Unique ~ FULL_NAME, 
    value.var = "COUNT", fill = 0
  ) %>%
  fwrite(
    file = "StructuralTE_TotalCount_FamilyLevel_Excel.txt", 
    col.names = T, sep = "\t"
  )

df_fam %>%
  arrange(
    TE_CLASS, TE_FAM
  ) %>%
  as.data.table() %>%
  data.table::dcast.data.table(
    TE_CLASS + TE_FAM + Share_or_Unique ~ FULL_NAME, 
    value.var = "SUM_LEN", fill = 0
  ) %>%
  fwrite(
    file = "StructuralTE_TotalLength_FamilyLevel_Excel.txt", 
    col.names = T, sep = "\t"
  )

#------------------------------------------------------------------------
print("# Draw prevelance -------------------")
df_plot %>%
  mutate(
    FULL_NAME = factor(
      FULL_NAME, 
      levels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX", "BoA", "JZSv2"), 
      labels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2")
    ),
    TE_CLASS = base::factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")
    ), 
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. oleracea")
    )
  ) %>%
  ggplot(
    aes(
      x = COUNT, 
      y = FULL_NAME,
      fill = TE_CLASS
    )
  )+
  geom_col(
    position = "stack"
  )+
  facet_grid(
    rows = vars(SPECIES), 
    scales = "free", 
    space = "free"
  )+
  scale_fill_manual(
    breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",  
               "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"), 
    values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8", 
               "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
  )+
  scale_x_continuous(
    limits = c(0, NA), expand = expansion(mult = c(0, 0.1), add = c(0, 0))
  )+
  theme_cowplot()+
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(
      size = 7,
      face = "italic"
    ),
    axis.title = element_text(
      size = 7
    ),
    axis.text = element_text(
      size = 7
    ),
    legend.title = element_text(
      size = 7
    ), 
    legend.text = element_text(
      size = 6
    ), 
    legend.position = "none"
  )+
  labs(
    x = "Count", 
    y = "", 
    fill = ""
  ) -> plot_C


df_plot %>%
  mutate(
    FULL_NAME = factor(
      FULL_NAME, 
      levels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX", "BoA", "JZSv2"), 
      labels = c("Chiffu", "CXA", "MIZ", "OIB", "PCA", "TUA", "Z1", "HDEM", "Korso", "OX-heart", "BoA", "JZSv2")
    ),
    TE_CLASS = base::factor(
      TE_CLASS,
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
                 "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")
    ), 
    SPECIES = factor(
      SPECIES, 
      levels = c("AA", "CC"),
      labels = c("B. rapa", "B. oleracea")
    )
  ) %>%
  ggplot(
    aes(
      x = SUM_LEN/1e+6, 
      y = FULL_NAME,
      fill = TE_CLASS
    )
  )+
  geom_col(
    position = "stack"
  )+
  facet_grid(
    rows = vars(SPECIES), 
    scales = "free", 
    space = "free"
  )+
  scale_fill_manual(
    breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",  
               "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"), 
    values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8", 
               "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
  )+
  scale_x_continuous(
    limits = c(0, NA), expand = expansion(mult = c(0, 0.1), add = c(0, 0))
  )+
  theme_cowplot()+
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(
      size = 7, 
      face = "italic"
    ),
    axis.title = element_text(
      size = 7
    ),
    axis.text = element_text(
      size = 7
    ),
    legend.title = element_text(
      size = 7
    ), 
    legend.text = element_text(
      size = 6
    )
  )+
  labs(
    x = "Length (Mb)", 
    y = "", 
    fill = "TE Superfamily"
  ) -> plot_L

plot_grid(
  plotlist = list(plot_C, plot_L), 
  align = "hv", rel_widths = c(1,1.2), labels = "AUTO"
)
ggsave2(
  filename = "Prevelene_StructuralTE.pdf", 
  width = 7.2,
  height = 3
)
