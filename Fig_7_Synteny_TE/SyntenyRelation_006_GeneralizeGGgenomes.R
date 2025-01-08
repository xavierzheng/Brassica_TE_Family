#---------------------------------------------------------
# AIM
#     Use gggenome to draw decrease plot
#
#---------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(viridis, quietly = T))
suppressPackageStartupMessages(library(RColorBrewer, quietly = T))
suppressPackageStartupMessages(library(ggsci, quietly = T))
suppressPackageStartupMessages(library(cowplot))
library(gggenomes)

setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

REGION <- "/Fig_7_Synteny_TE/example_region_decreas_2.txt"
SYNTENY <- "/Fig_7_Synteny_TE/example_synteny_decrease_2.txt"

df_region <- read.table(
  file = REGION,
  header = F, sep = "\t"
)

print("# prepare gene list =======================================")
list_gene <- lapply(df_region$V1, FUN = function(NAM){
  
  FEAT <- df_region %>% filter(V1 %in% NAM) %>% select(V2) %>% pull()
  CHROM <- df_region %>% filter(V1 %in% NAM) %>% select(V3) %>% pull()
  START <- df_region %>% filter(V1 %in% NAM) %>% select(V4) %>% pull()
  END <- df_region %>% filter(V1 %in% NAM) %>% select(V5) %>% pull()
  
  df_temp <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/GeneAnnotate_",
                  NAM, ".CC.InChr.sort.gff.gz"),
    header = F, sep = "\t"
  ) %>%
    filter(
      V1 %in% CHROM,
      V3 %in% FEAT,
      V4 > START,
      V5 < END
    ) %>%
    mutate(
      PLANT_ID = NAM,
      feature_id = str_extract(V9, "ID=[A-Za-z0-9_]{1,}") %>% str_remove(., "ID=")
    )
  
  return(df_temp)
  
})

df_gene <- bind_rows(list_gene) %>%
  mutate(
    feature_id = ifelse(
      str_detect(feature_id, "NA"),
      NA, str_remove(feature_id, "^ ") %>% str_remove(., ",[ A-Za-z0-9]{1,}")
    )
  ) %>%
  dplyr::rename(
    seq_id = V1,
    start = V4,
    end = V5,
    strand = V7,
  )
glimpse(df_gene)

print("# prepare sequene file ----------------------")
list_seq <- lapply(df_region$V1, FUN = function(NAM){
  
  CHROM <- df_region %>% filter(V1 %in% NAM) %>% select(V3) %>% pull()
  START <- df_region %>% filter(V1 %in% NAM) %>% select(V4) %>% pull()
  END <- df_region %>% filter(V1 %in% NAM) %>% select(V5) %>% pull()
  
  df_temp <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/Genome_Fasta/",
                  NAM, ".CC.genome.fa.fai") # remember to generate it by samtools/1.13: samtools faidx XXX.fa
  ) %>%
    filter(
      V1 %in% CHROM
    ) %>%
    select(
      V1, V2
    ) %>%
    rename(
      seq_id = V1,
      length = V2
    ) %>%
    mutate(
      start = START,
      end = END,
      PLANT_ID = NAM
    )
  
  return(df_temp)
  
})

df_seq <- bind_rows(list_seq)
glimpse(df_seq)


print("# read structural TE ------------------")
list_TE <- lapply(df_region$V1, FUN = function(NAM){
  
  CHROM <- df_region %>% filter(V1 %in% NAM) %>% select(V3) %>% pull()
  START <- df_region %>% filter(V1 %in% NAM) %>% select(V4) %>% pull()
  END <- df_region %>% filter(V1 %in% NAM) %>% select(V5) %>% pull()
  
  df_temp <- fread(
    file = paste0("/Data/Fig_4_TE_Feature/GFF_Gene_TE/",
                  NAM, ".CC.TE.InChr.IntactSingleLine.sort.gff3"),
    header = F, sep = "\t"
  ) %>%
    filter(
      V1 %in% CHROM,
      V4 > START,
      V5 < END
    ) %>%
    mutate(
      feature_id = str_extract(V9, "Name=[A-Za-z0-9]{1,}") %>% str_remove(., "Name="),
      TE_CLASS = str_extract(V9, "Classification=[A-Za-z0-9\\/]{1,}") %>% str_remove(., "Classification="),
      PLANT_ID = NAM
    ) %>%
    dplyr::rename(
      seq_id = V1,
      start = V4,
      end = V5,
      strand = V7,
    )
  
  return(df_temp)
  
})
df_TE <- bind_rows(list_TE) %>%
  mutate(
    strand = "+"
  )
glimpse(df_TE)

print("# generate link/synteny -----------------------")
df_link <- fread(
  file = SYNTENY,
  header = T, sep = "\t"
)

print("# give TE as factor (because of color)------------------")
sort(unique(df_TE$TE_CLASS))

df_TE <- df_TE %>%
  mutate(
    TE_CLASS = factor(
      TE_CLASS, 
      levels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", "TIR/DTM"), 
      labels = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", "TIR/DTM")
    )
  )
glimpse(df_TE)
sort(unique(df_TE$TE_CLASS))

print("# drawing ---------------------")
gggenomes::gggenomes(
  genes = df_gene,
  seqs = df_seq,
  feats = list(TE = df_TE),
  links = df_link
) %>%
  sync() +
  gggenomes::geom_seq()+
  gggenomes::geom_bin_label()+
  gggenomes::geom_link(
    size = 1
  )+
  gggenomes::geom_gene()+
  gggenomes::geom_gene_tag(
    aes(label=feature_id),
    nudge_y=0.1,
    size = 5*0.35
    #check_overlap = TRUE
  )+
  gggenomes::geom_feat(
    data = feats(TE),
    mapping = aes(
      color = TE_CLASS
    ),
    position = "identity",
    linewidth = 4
  )+
  gggenomes::geom_feat_tag(
    data = feats(TE),
    mapping = aes(label = feature_id),
    nudge_y=0.1,
    size = 5*0.35
    #check_overlap = TRUE
  )+
  #ggsci::scale_color_npg()+
  scale_color_manual(
    breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron", "TIR/DTC", "TIR/DTM"), 
    values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8", "#6baed6", "#a1d99b")
  )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(
      size = 7
    ),
    legend.text = element_text(
      size = 6
    )
  )+
  labs(
    color = "TE superfamily"
  )
ggsave2(
  filename = "/Fig_7_Synteny_TE/example_region_decreas_2.pdf",
  width = 8.5,
  height = 5
)


# define color code -----------------------------
# scale_fill_manual(
#   breaks = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown", "DNA/Helitron",
#              "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT"),
#   values = c("#a50f15", "#fb6a4a", "#fee5d9", "#9e9ac8",
#              "#08519c", "#6baed6", "#31a354", "#a1d99b", "#edf8e9")
# )
  