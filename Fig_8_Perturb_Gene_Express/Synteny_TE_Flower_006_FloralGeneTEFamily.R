#==================================
# AIM
#	  Check how many different TE superfamily and family are inserted aroung floral gene
#	  Using Syntenic relationship
#
#==================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
library(ggupset)
library(patchwork)
library(treemapify)

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 讓arrow可以多核心讀寫
arrow::set_cpu_count(12)
arrow::set_io_thread_count(12)
arrow::cpu_count()

# read data============================================
print("# read floral gene list -----------")
df_flor_AA <- read.table(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/gene_florID_DA.txt.gz", 
  header = T, sep = "\t"
) %>%
  mutate(
    BoA = str_remove_all(BoA, "^DA_")
  )

df_flor_CC <- read.table(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/gene_florID_IJ.txt.gz", 
  header = T, sep = "\t"
) %>%
  mutate(
    BoA = str_split_fixed(BoA, "\\.",2)[,1],
    BoA = str_split_fixed(BoA, "_", 2)[,2]
  )

floral_AA_vec <- sort(unique(df_flor_AA$BoA))
floral_CC_vec <- sort(unique(df_flor_CC$BoA))

# statistic analysis-------------------------------------------------------------------------
print("# Specifically focus on syntenic core genes: Do floral genes associate with syntenic core gene?")

stat_core_floral <- function(PLANT, REF, FLORAL_VEC){
  
  df_ref <- fread(
    file = paste0("/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneAnnotate_", REF, "_Chr_gene.txt"),
    header = F, sep = "\t"
  ) %>%
    mutate(
      GENEID = str_extract(V9, "ID=[A-Za-z0-9]{1,}") %>% str_remove(., "ID="), 
    ) %>%
    dplyr::rename(
      CHROM = V1
    ) %>%
    select(
      CHROM, GENEID
    )
  
  df_core <- fread(
    file = paste0("/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/PANGENOME_cds_", PLANT, "_info.txt"), 
    header = T, sep = "\t"
  ) %>%
    filter(
      PANGENOME %in% "CORE", 
      PREFIX %in% REF
    ) %>%
    mutate(
      GENEID = str_split_fixed(GENEID, "_", 2)[,2]
    )
  
  df_stat <- left_join(
    x = df_ref, 
    y = df_core,
    by = "GENEID"
  ) %>%
    mutate(
      PANGENOME = ifelse(
        is.na(PANGENOME), 
        "NON_CORE", PANGENOME
      ), 
      FLORAL = ifelse(
        GENEID %in% FLORAL_VEC,
        "FLORAL", "NON_FLORAL"
      ),
      PANGENOME = factor(
        PANGENOME, 
        levels = c("CORE", "NON_CORE"), 
        labels = c("CORE", "NON_CORE")
      ), 
      FLORAL = factor(
        FLORAL,
        levels = c("FLORAL", "NON_FLORAL"), 
        labels = c("FLORAL", "NON_FLORAL")
      )
    )
  
  return(df_stat)
  
}

df_stat_AA <- stat_core_floral(PLANT = "AA", REF = "DA", FLORAL_VEC = floral_AA_vec)
df_stat_CC <- stat_core_floral(PLANT = "CC", REF = "IJ", FLORAL_VEC = floral_CC_vec)

table(df_stat_AA$PANGENOME, df_stat_AA$FLORAL)
fisher.test(table(df_stat_AA$PANGENOME, df_stat_AA$FLORAL))

table(df_stat_CC$PANGENOME, df_stat_CC$FLORAL)
fisher.test(table(df_stat_CC$PANGENOME, df_stat_CC$FLORAL))

# AA---------------------------
#           FLORAL NON_FLORAL
# CORE        578      25717
# NON_CORE    112      20065

# Fisher's Exact Test for Count Data
# 
# data:  table(df_stat$PANGENOME, df_stat$FLORAL)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.280589 4.979107
# sample estimates:
# odds ratio 
#   4.026344 

# CC---------------------------
#           FLORAL NON_FLORAL
# CORE        481      28454
# NON_CORE     64      29382

# Fisher's Exact Test for Count Data
# 
# data:  table(df_stat_CC$PANGENOME, df_stat_CC$FLORAL)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   5.966173 10.244602
# sample estimates:
# odds ratio 
#    7.76043 
print("# Conclusion: floral genes have significantly higher chance to be syntenic core genes in AA and CC ------")



# GO enrichment like analysis--------------------------------------------------
print("Do floral genes also enrich in TE-syntenic anchor genes ?-----------------------")

stat_ref_TE_floral <- function(PLANT, REF, TAB_LOC, FLORAL_VEC){
  
  df_ref <- fread(
    file = paste0("/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneAnnotate_", REF, "_Chr_gene.txt"),
    header = F, sep = "\t"
  ) %>%
    mutate(
      GENEID = str_extract(V9, "ID=[A-Za-z0-9]{1,}") %>% str_remove(., "ID="), 
    ) %>%
    dplyr::rename(
      CHROM = V1
    ) %>%
    select(
      CHROM, GENEID
    )
  
  df_TE_temp <- read.table(
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
  
  df_TE <- data.frame(
    GENEID = sort(unique(c(df_TE_temp$GENE1, df_TE_temp$GENE2))), 
    TE = "TE"
  )
  
  df_stat <- left_join(
    x = df_ref, 
    y = df_TE, 
    by = "GENEID"
  ) %>%
    mutate(
      TE = ifelse(
        is.na(TE), 
        "NON_TE", TE
      ), 
      FLORAL = ifelse(
        GENEID %in% FLORAL_VEC, 
        "FLORAL", "NON_FLORAL"
      ), 
      TE = factor(
        TE, 
        levels = c("TE", "NON_TE"),
        labels = c("TE", "NON_TE")
      ), 
      FLORAL = factor(
        FLORAL, 
        levels = c("FLORAL", "NON_FLORAL"), 
        labels = c("FLORAL", "NON_FLORAL")
      )
    )
  
  return(df_stat)
  
}

rm(df_stat_AA)
rm(df_stat_CC)

df_stat_AA <- stat_ref_TE_floral(
  PLANT = "AA", REF = "DA", 
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_AA.txt.gz", 
  FLORAL_VEC = floral_AA_vec)

df_stat_CC <- stat_ref_TE_floral(
  PLANT = "CC", REF = "IJ", 
  TAB_LOC = "/Data/Fig_7_Synteny_TE/Intact_Class_TECounts_betweenanchors_CC.txt.gz", 
  FLORAL_VEC = floral_CC_vec)

table(df_stat_AA$TE, df_stat_AA$FLORAL)
table(df_stat_CC$TE, df_stat_CC$FLORAL)
fisher.test(table(df_stat_AA$TE, df_stat_AA$FLORAL))
fisher.test(table(df_stat_CC$TE, df_stat_CC$FLORAL))

fwrite(
  df_stat_AA,
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/enrichment_TE_syntenic_AA.txt", 
  col.names = T, sep = "\t"
)

fwrite(
  df_stat_CC, 
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/enrichment_TE_syntenic_CC.txt",
  col.names = T, sep = "\t"
)

# AA------------------------------------
#         FLORAL NON_FLORAL
# TE        317      14825
# NON_TE    373      30957
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(df_stat_AA$TE, df_stat_AA$FLORAL)
# p-value = 2.245e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.521011 2.069536
# sample estimates:
# odds ratio 
#   1.774541 

# CC-------------------------------------
#         FLORAL NON_FLORAL
# TE        265      14433
# NON_TE    280      43403
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(df_stat_CC$TE, df_stat_CC$FLORAL)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.394161 3.382432
# sample estimates:
# odds ratio 
#   2.846164 
print("# Conclusion: Besides the GO enrichment results, we also found that floral genes are also enriched in TE-synteic gene pairs =======================")



# Draw: How many different TE family near floral gene ------------------------------------
print("# counting function for whole syntenic TE family ---------")
count_TE_fam <- function(TAB_FAM_LOC){
  
  df_TE_temp <- read.table(
    file = TAB_FAM_LOC,
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
    ) %>%
    mutate(
      TE_CLASS = ifelse(
        str_detect(TE_FAM, "^Copia"),
        "LTR/Copia", ifelse(
          str_detect(TE_FAM, "^Gypsy"),
          "LTR/Gypsy", ifelse(
            str_detect(TE_FAM, "^unknown"),
            "LTR/unknown", ifelse(
              str_detect(TE_FAM, "^Helitron"),
              "DNA/Helitron", ifelse(
                str_detect(TE_FAM, "^DTA"),
                "TIR/DTA", ifelse(
                  str_detect(TE_FAM, "^DTC"),
                  "TIR/DTC", ifelse(
                    str_detect(TE_FAM, "^DTH"),
                    "TIR/DTH", ifelse(
                      str_detect(TE_FAM, "^DTM"),
                      "TIR/DTM", ifelse(
                        str_detect(TE_FAM, "^DTT"),
                        "TIR/DTT", NA
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  
  df_TE_count <- df_TE_temp %>%
    summarise(
      .by = c(TE_CLASS, TE_FAM), 
      N_FAM = n()
    ) %>%
    mutate(
      .by = c(TE_CLASS), 
      N_CLAS = sum(N_FAM)
    ) %>%
    arrange(
      desc(N_CLAS)
    )
  
  return(df_TE_count)
}

# print("# counting function for TE family near floral genes ------")
# count_TE_fam_floral <- function(TAB_FAM_LOC, FLORAL_VEC){
  
#   df_TE_temp <- read.table(
#     file = TAB_FAM_LOC,
#     header = T, sep = "\t"
#   ) %>%
#     tibble::rownames_to_column(
#       var = "KEY"
#     ) %>%
#     mutate(
#       PAIR = str_split_fixed(KEY, " ", 2)[,1],
#       TE_FAM = str_split_fixed(KEY, " ", 2)[,2], 
#       GENE1 = str_split_fixed(PAIR, "-", 2)[,1], 
#       GENE2 = str_split_fixed(PAIR, "-", 2)[,2], 
#       GENE1 = str_remove(GENE1, "\\.3\\.1C$"),
#       GENE2 = str_remove(GENE2, "\\.3\\.1C$")
#     ) %>%
#     mutate(
#       TE_CLASS = ifelse(
#         str_detect(TE_FAM, "^Copia"),
#         "LTR/Copia", ifelse(
#           str_detect(TE_FAM, "^Gypsy"),
#           "LTR/Gypsy", ifelse(
#             str_detect(TE_FAM, "^unknown"),
#             "LTR/unknown", ifelse(
#               str_detect(TE_FAM, "^Helitron"),
#               "DNA/Helitron", ifelse(
#                 str_detect(TE_FAM, "^DTA"),
#                 "TIR/DTA", ifelse(
#                   str_detect(TE_FAM, "^DTC"),
#                   "TIR/DTC", ifelse(
#                     str_detect(TE_FAM, "^DTH"),
#                     "TIR/DTH", ifelse(
#                       str_detect(TE_FAM, "^DTM"),
#                       "TIR/DTM", ifelse(
#                         str_detect(TE_FAM, "^DTT"),
#                         "TIR/DTT", NA
#                       )
#                     )
#                   )
#                 )
#               )
#             )
#           )
#         )
#       )
#     )
  
#   df_TE_floral <- df_TE_temp %>%
#     filter(
#       GENE1 %in% FLORAL_VEC | GENE2 %in% FLORAL_VEC
#     ) %>%
#     summarise(
#       .by = c(TE_CLASS, TE_FAM), 
#       N_FAM = n()
#     ) %>%
#     mutate(
#       .by = c(TE_CLASS), 
#       N_CLAS = sum(N_FAM)
#     ) %>%
#     arrange(
#       desc(N_CLAS)
#     )
  
#   return(df_TE_floral)
# }

print("# drawing function for treemap -----------")
plot_TE_fam <- function(DF){
  
  DF %>%
    ggplot(
      aes(
        area = N_FAM,
        fill = N_FAM, 
        label = TE_FAM,
        subgroup = TE_CLASS
      )
    )+
    treemapify::geom_treemap(
      color = "white"
    )+
    treemapify::geom_treemap_text(
      place = "center", 
      grow = T, 
      min.size = 5*0.34, 
      color = "white", 
    )+
    treemapify::geom_treemap_subgroup_border(
      color = "white", size = 3
    )+
    treemapify::geom_treemap_subgroup_text(
      place = "center", 
      grow = T, 
      min.size = 5*0.34, 
      color = "#08306b", 
      alpha = 0.7, 
      fontface = "italic"
    )+
    scale_fill_distiller(
      palette = "RdPu",
      direction = 1
    )+
    theme(
      legend.title = element_text(
        size = 6, vjust = 1
      ),
      legend.text = element_text(
        size = 6
      ),
      legend.position = "bottom", 
      legend.key.width = unit(20, units = "pt"), 
      legend.key.height = unit(5, units = "pt")
    )+
    labs(
      fill = "TE family count"
    ) -> plot_ret
  
  return(plot_ret)
  
}

print("# using these function to obtain table and plot -----------------")
df_TE_all_AA <- count_TE_fam(
  TAB_FAM_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_AA.txt.gz"
)

df_TE_all_CC <- count_TE_fam(
  TAB_FAM_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_CC.txt.gz"
)

# df_TE_floral_AA <- count_TE_fam_floral(
#   TAB_FAM_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_AA.txt.gz", 
#   FLORAL_VEC = floral_AA_vec
# )

# df_TE_floral_CC <- count_TE_fam_floral(
#   TAB_FAM_LOC = "/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_CC.txt.gz", 
#   FLORAL_VEC = floral_CC_vec
# )

plot_all_AA <- plot_TE_fam(DF = df_TE_all_AA)
plot_all_CC <- plot_TE_fam(DF = df_TE_all_CC)
# plot_AA <- plot_TE_fam(DF = df_TE_floral_AA)
# plot_CC <- plot_TE_fam(DF = df_TE_floral_CC)

print("# which family is the most dominant ?---------------------------------")
glimpse(df_TE_all_AA)
glimpse(df_TE_all_CC)

df_TE_all_AA %>%
  arrange(
    desc(N_FAM)
  ) %>%
  glimpse()

df_TE_all_CC %>%
  arrange(
    desc(N_FAM)
  ) %>%
  glimpse()

print("# saving plot -------------------------------")
# plot_floral <- plot_AA + plot_CC + patchwork::plot_layout(
#   ncol = 2, widths = c(1,1))

# cowplot::ggsave2(
#   filename = "Treemap_syntenic_floral_TE_family_diversity.pdf", 
#   plot = plot_floral,
#   width = 7.2,
#   height = 3.2, 
#   units = "in"
# )

plot_all <- plot_all_AA + plot_all_CC + patchwork::plot_layout(
  ncol = 2, widths = c(1,1))

cowplot::ggsave2(
  filename = "/Supplement_Data/Treemap_syntenic_TE_family_diversity.pdf", 
  plot = plot_all,
  width = 7.2,
  height = 3.2, 
  units = "in"
)


