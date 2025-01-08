#===============================
# AIM
#	1. combine the REXdb BLASTX result
#	2. Only Focus on TE within the same superfamily
#	3. Check similarity on our family and classification on REXdb
#
#===============================

setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/Compare_REXdb")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggalluvial))


setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

#============================
print("# this code work for Copia and Gypsy and unknown")

NAME <- "Copia"
NAME <- "Gypsy"
NAME <- "unknown"

temp_list <- list.files(
  path = "BLAST_OUT/", 
  pattern = paste0("BLASTX_REXdb_IntackTE_", NAME,".fa_p.*")
)

list2 <- lapply(temp_list, FUN = function(NAM){
  df_temp <- fread(
    cmd = paste0("grep -v '^#' BLAST_OUT/", NAM), 
    header = F, sep = "\t", nThread = 12
  )
  return(df_temp)
})

df <- bind_rows(list2)

rm(list2)

colnames(df) <- c("qseqid", "sseqid", "qlen", "slen", "sstrand", "pident", "LENGTH", "mismatch", "gapopen", "gaps", "qframe", "sframe", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

df <- as.data.table(df)

str(df)
rm(df_proc)

df_proc <- df[
  ,c("Q_ID", "Q_FAM", "Q_TEMP"):=tstrsplit(qseqid, ";", fixed =T)
][
  ,c("Q_TEMP2", "Q_POS"):= tstrsplit(Q_TEMP, "::", fixed = T)
][
  ,c("Q_TEMP3", "Q_TEMP4") := tstrsplit(Q_TEMP2, "=", fixed = T)
][
  ,c("Q_ORDER", "Q_SUPFAM"):= tstrsplit(Q_TEMP4, "/", fixed = T)
][
  , c("S_PROT", "S_XXX", "S_REXIDTEMP"):=tstrsplit(sseqid, "_", fixed = T)
][
  , c("REX_ID", "S_CLAS", "S_ORDER", "S_SUPFAM", "S_FAM"):=tstrsplit(S_REXIDTEMP, ":", fixed = T)
][
  ,!c("Q_TEMP", "Q_TEMP2", "Q_TEMP3", "Q_TEMP4", "S_XXX", "S_REXIDTEMP", "S_CLAS", "S_ORDER", "S_SUPFAM", "S_FAM")
][
  ,c("Q_ID", "Q_FAM"):= .(str_remove_all(Q_ID, "^ID="), str_remove_all(Q_FAM, "^Name="))
]

str(df_proc)


print("# import rexdb =========================================")
df_rex <- read.delim2(
  file = "Viridiplantae_v3.0_ALL_classification_edited",
  header = F, sep = "\t", fill = T, 
  col.names = c("REX_ID", "REX_CLASS", "REX_ORDER", "REX_SUPFAM", 
                "REX_FAM1", "REX_FAM2", "REX_FAM3", "REX_FAM4")
) %>%
  mutate(
    REX_FAM = paste0(REX_FAM1, "/", REX_FAM2, "/", REX_FAM3, "/", REX_FAM4), 
    REX_FAM = str_remove_all(REX_FAM, "/*$")
  )

glimpse(df_rex)

print("# add rexdb information to blast result======================")
df_proc_rex <- left_join(
  x = df_proc, 
  y = df_rex, 
  by = "REX_ID"
)

str(df_proc_rex)

print("# filter target with 80-80-80 rule -----------------------")
df_temp1 <- df_proc_rex %>%
  filter(
    LENGTH/slen > 0.8, # alignment coverage > 80%
    pident > 80, # identity > 80%
    qend-qstart > 80 # alignment nt more than 80nt
  )

df_temp1 <- as.data.table(df_temp1)

glimpse(df_temp1)


print("# select the best hit per query and subject --------------")
df_report <- df_temp1[
  order(-LENGTH, -pident, evalue -bitscore), head(.SD, 1) ,by = .(qseqid, sseqid)
][
  ,KEY_ORDER := ifelse(Q_ORDER==REX_ORDER, Q_ORDER, "DIFF")
][
  ,KEY_SUPFAM := ifelse(Q_SUPFAM==REX_SUPFAM, Q_SUPFAM, "DIFF")
]

print("# select the best hit per query (each TE only keep one subject result )--------------")
df_report2 <- df_temp1[
  order(-LENGTH, -pident, evalue -bitscore), head(.SD, 1) ,by = .(qseqid)
][
  ,KEY_ORDER := ifelse(Q_ORDER==REX_ORDER, Q_ORDER, "DIFF")
][
  ,KEY_SUPFAM := ifelse(Q_SUPFAM==REX_SUPFAM, Q_SUPFAM, "DIFF")
]

glimpse(df_report)
glimpse(df_report2)

print("# per query is better. However, I save the result of per query and per subject in order to have more flexibility ================================================")
fwrite(
  df_report, 
  file = paste0("Summary_REXdb_besthit_", NAME, ".txt"), 
  col.names = T, sep = "\t"
)


print("# NEXT, using per-query df_report2 for downstream analysis ======================")
fam_list <- sort(unique(df_report2$REX_FAM))

name_list <- sort(unique(df_proc_rex$Q_FAM))


print("# prepare count table for allu plot, using family vs. family count ===========")
df_draw <- df_report2 %>%
  group_by(
    Q_FAM, REX_FAM
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  mutate(
    Q_FAM = factor(
      Q_FAM, 
      levels = name_list, 
      labels = name_list
    ), 
    REX_FAM = factor(
      REX_FAM, 
      levels = fam_list, 
      labels = fam_list
    )
  )

glimpse(df_draw)

print("# give color label to distinguish the single family and multiple families ==========")
df_allu <- df_draw %>%
  group_by(
    Q_FAM
  ) %>%
  summarise(
    COUNT = n()
  ) %>%
  ungroup() %>%
  mutate(
    GROUP = ifelse(
      COUNT == 1, 
      "Single", "Multiple"
    )
  ) %>%
  select(
    -COUNT
  ) %>%
  left_join(
    x =  df_draw, 
    y = ., 
    by = "Q_FAM"
  ) %>%
  arrange(
    Q_FAM, -COUNT
  ) %>%
  group_by(
    Q_FAM
  ) %>%
  mutate(
    PERC = COUNT/sum(COUNT)
  ) %>%
  ungroup() %>%
  group_by(
    Q_FAM
  ) %>%
  mutate(
    ANY = any(PERC > 0.95), 
    GROUP2 = ifelse(
      GROUP=="Multiple" && ANY==TRUE, 
      "Near_single", GROUP
    )
  ) %>%
  ungroup() %>%
  mutate(
    GROUP2 = factor(
      GROUP2,
      levels = c("Single", "Near_single", "Multiple"),
      labels = c("Single", "Near_single", "Multiple")
    )
  ) %>%
  select(
    -PERC, -ANY, -GROUP
  )

glimpse(df_allu)


ggplot(
  data = df_allu, 
  mapping = aes(
    axis1 = Q_FAM, axis2 = REX_FAM, y = COUNT
  )
)+
  scale_x_discrete(
    limits = c("Our\nfamily", "REXdb\nfamily"), 
    expand = expansion(mult = 0, add = c(0.1, 0.1))
  )+
  geom_alluvium(
    aes(
      fill = GROUP2
    )
  )+
  geom_stratum(width = 1/3) +
  geom_text(
    size = 2,
    stat = "stratum", 
    aes(label = after_stat(stratum))
  )+
  scale_y_continuous(
    expand = expansion(mult = 0, add = 0)
  )+
  scale_fill_manual(
    breaks = c("Single", "Near_single", "Multiple"), 
    values = c("#1b9e77", "#7570b3", "#d95f02")
  )+
  theme_cowplot()+
  theme(
    axis.line = element_blank(), 
    axis.ticks = element_blank(), 
    axis.title = element_blank(), 
    axis.text.y = element_blank(), 
    plot.margin = margin(0, 5, 0 ,0, unit = "pt"), 
    axis.text.x = element_text(
      size = 7
    ), 
    legend.text = element_text(
      size = 7
    ), 
    legend.position = "bottom"
  )+
  labs(
    fill = ""
  ) -> plot_p

plot_p
cowplot::ggsave2(
  plot = plot_p, 
  filename = paste0("alluvium_REXdb_comparison_", NAME, ".pdf"), 
  width = 3.5, 
  height = 6
)


#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================







# only 38 results....
# df_temp2 <- df_temp1[
#   order(-bitscore, -LENGTH, evalue, -pident), .N ,by = .(qseqid, S_SUPFAM)
# ][
#   order(qseqid, -N), head(.SD, 1), by = .(qseqid)
# ]


dim(df_temp2)
#[1] 1889    3

#[1] 1038    3


# unique(df_clean$Q_ORDER_REXDB)
# table(df_temp1$S_SUPFAM)
