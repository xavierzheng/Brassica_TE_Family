#-----------------------------
# AIM
#	* Compare the expression variant between TE-anchor genes (they are core, anchor, within multiplicon) 
#	and all of core, anchor, within multiplicon genes
#
#	* Majorly focus on the part of all of core, anchor, within multiplicon genes
#
#-----------------------------


suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(ggsci, quietly = F))
library(arrow)

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 讓arrow可以多核心讀寫
arrow::set_cpu_count(12)
arrow::set_io_thread_count(12)
arrow::cpu_count()

#-----------------------------------
print("# read pangenome info ----")

define_TE_core_synteny <- function(SPECIES, PREFIX){
  
  SPECIES <- as.character(SPECIES)
  PREFIX <- as.character(PREFIX)
  
  # read pangenome table
  df_pan <- fread(
    file = paste0("/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/PANGENOME_cds_", SPECIES, "_info.txt")
  )
  
  # subset the CORE 
  df_core <- df_pan[PANGENOME=="CORE",]
  
  # read anchor-gene-anchor table
  df_TE_temp <- read.table(
    file = paste0("/Data/Fig_7_Synteny_TE/Intact_family_TECounts_betweenanchors_", SPECIES, ".txt.gz"), 
    header = T, sep= "\t"
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
  
  # extract the unique TE anchor gene
  df_TE_anchor <- data.frame(
    GENEID = sort(unique(c(df_TE_temp$GENE1, df_TE_temp$GENE2))),
    TE_SYNTENY = "TE"
  ) %>%
    mutate(
      GENEID = paste0(PREFIX, "_", GENEID)
    )
  
  print(paste0("# QC: are all of the TE-synteny anchor genes belong to core gene? ",
               dim(df_TE_anchor)[1] == sum(df_TE_anchor$GENEID %in% df_core$GENEID)))
  
  # joint result
  df_ret <- left_join(
    x = df_core,
    y = df_TE_anchor,
    by = "GENEID"
  ) %>%
    mutate(
      TE_SYNTENY = ifelse(
        is.na(TE_SYNTENY),
        "NON_TE", TE_SYNTENY
      ),
      TE_SYNTENY = factor(
        TE_SYNTENY, 
        levels = c("TE", "NON_TE"),
        labels = c("TE", "NON_TE")
      )
    ) %>%
    select(
      GENEID, PANGENOME, TE_SYNTENY
    )
  
  return(df_ret)
  
}

df_AA <- define_TE_core_synteny(SPECIES = "AA", PREFIX = "DA")
df_CC <- define_TE_core_synteny(SPECIES = "CC", PREFIX = "IJ")


print("# read expression ------")
# prepare function to add range and cv-------------------
add_cv_range <- function(SPECIES, DF_SYNTENY){
  
  SPECIES <- as.character(SPECIES)
  
  # read expression
  df_exp <- arrow::read_feather(
    file = paste0("/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/Result_kallistro_", SPECIES, "_AddInfo.arrow")
  )
  
  # focus on syntenic anchor gene, and estimate log2 TPM
  df_temp <- df_exp[
    GENEID %in% DF_SYNTENY$GENEID
  ][
    , "LOG2_TPM" := log2(TPM+0.001)
  ][
    ,.(SAMPLE, GENEID, LOG2_TPM, TPM)
  ]
  
  # based on log2 TPM to estimate CV and range, and filter
  df_key <- df_temp[
    , .(RANGE = (max(LOG2_TPM)-min(LOG2_TPM)), CV = (sd(LOG2_TPM)/mean(LOG2_TPM))), by = .(GENEID)
  ][
    RANGE > 5 & CV > 0.2,
  ][
    ,.(GENEID, CV, RANGE)
  ]
  
  # add gene info back, including CORE and TE_SYNTENY
  df_ret <- left_join(
    x = df_key,
    y = DF_SYNTENY, 
    by = "GENEID"
  )
  
  return(df_ret)
}

df_AA_cv <- add_cv_range(SPECIES = "AA", DF_SYNTENY = df_AA)
df_CC_cv <- add_cv_range(SPECIES = "CC", DF_SYNTENY = df_CC)

str(df_AA_cv)
str(df_CC_cv)

print("# No. of high variant gene ---------------------------------")
df_AA_cv[
  ,.(COUNT = .N), by = .(TE_SYNTENY)
]

#     TE_SYNTENY COUNT
# 1:     NON_TE  2366
# 2:         TE  2929

df_CC_cv[
  ,.(COUNT = .N), by = .(TE_SYNTENY)
]

#     TE_SYNTENY COUNT
# 1:         TE  2122
# 2:     NON_TE  1779

#----------------------------------------------------------------------
# function for prepare accumulate data frame ------------

draw_accum_df <- function(DF, FEATURE){
  
  df_ret <- DF %>%
    group_by(
      TE_SYNTENY
    ) %>%
    arrange(
      .by_group = T,
      desc(get(FEATURE))
    ) %>%
    mutate(
      RANK = 1:n(), 
      CUM = cumsum(get(FEATURE))
    ) %>%
    ungroup() 
  
  return(df_ret)
}

# function for plot ---------
draw_accum_plot <- function(DF_ACCUM, AXIS_Y){
  
  DF_ACCUM %>%
    ggplot(
      aes(
        x = RANK,
        y = CUM,
        color = TE_SYNTENY,
        group = TE_SYNTENY
      )
    )+
    geom_point(
      size = 1.2
    )+
    geom_line()+
    scale_color_manual(
      breaks = c("TE", "NON_TE"), 
      values = c("#FF7043FF", "#41A5F4FF")
    )+
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.1), add = c(0, 0))
    )+
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1), add = c(0, 0))
    )+
    theme_cowplot()+
    theme(
      axis.title = element_text(
        size = 7
      ),
      axis.text = element_text(
        size = 6
      ),
      legend.title = element_text(
        size = 7
      ),
      legend.text = element_text(
        size = 6
      ),
      legend.position.inside = c(0.7, 0.2)
    )+
    labs(
      x = "Ranked geneid",
      y = paste0("Accumulated ", AXIS_Y)
    ) -> plot_ret
  return(plot_ret)
  
}

# non-linear statistic analysis
draw_accum_stat <- function(DF_ACCUM){
  
  gam_model_ <- mgcv::gam(CUM ~ s(RANK) + TE_SYNTENY + s(RANK, by = TE_SYNTENY), data = DF_ACCUM)
  summary(gam_model_)
  
}


print("# drawing the accumulating plot ------------------")
df_draw_AA_CV <- draw_accum_df(DF = df_AA_cv, FEATURE = "CV")
plot_AA_CV <- draw_accum_plot(DF_ACCUM = df_draw_AA_CV, AXIS_Y = "CV")
draw_accum_stat(DF_ACCUM = df_draw_AA_CV)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE_SYNTENY + s(RANK, by = TE_SYNTENY)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      14649.775      5.447 2689.30   <2e-16 ***
#   TE_SYNTENYNON_TE -1647.205    162.467  -10.14   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F p-value    
# s(RANK)                  8.6615 8.6667 3134.77  <2e-16 ***
#   s(RANK):TE_SYNTENYTE     0.6667 0.6667   43.98  <2e-16 ***
#   s(RANK):TE_SYNTENYNON_TE 8.6501 8.6664  530.78  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =  0.983   Deviance explained = 98.3%
# GCV =  83559  Scale est. = 83243     n = 5295


df_draw_AA_RANGE <- draw_accum_df(DF = df_AA_cv, FEATURE = "RANGE")
plot_AA_RANGE <- draw_accum_plot(DF_ACCUM = df_draw_AA_RANGE, AXIS_Y = "range")
draw_accum_stat(DF_ACCUM = df_draw_AA_RANGE)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE_SYNTENY + s(RANK, by = TE_SYNTENY)
# 
# Parametric coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)      18909.670      0.466 40579.45   <2e-16 ***
#   TE_SYNTENYNON_TE  -661.973     13.890   -47.66   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(RANK)                  8.6653 8.6665 51775  <2e-16 ***
#   s(RANK):TE_SYNTENYTE     0.6668 0.6668  1474  <2e-16 ***
#   s(RANK):TE_SYNTENYNON_TE 8.6489 8.6664 32570  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 611.45  Scale est. = 609.14    n = 5295

df_draw_CC_CV <- draw_accum_df(DF = df_CC_cv, FEATURE = "CV")
plot_CC_CV <- draw_accum_plot(DF_ACCUM = df_draw_CC_CV, AXIS_Y = "CV")
draw_accum_stat(DF_ACCUM = df_draw_CC_CV)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE_SYNTENY + s(RANK, by = TE_SYNTENY)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       36630.763      3.507 10446.3   <2e-16 ***
#   TE_SYNTENYNON_TE -29315.949     40.670  -720.8   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(RANK)                  8.6529 8.6662 765.27  <2e-16 ***
#   s(RANK):TE_SYNTENYTE     8.6482 8.6661 375.80  <2e-16 ***
#   s(RANK):TE_SYNTENYNON_TE 0.6667 0.6667  59.48  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV =  25288  Scale est. = 25159     n = 3901

df_draw_CC_RANGE <- draw_accum_df(DF = df_CC_cv, FEATURE = "RANGE")
plot_CC_RANGE <- draw_accum_plot(DF_ACCUM = df_draw_CC_RANGE, AXIS_Y = "range")
draw_accum_stat(DF_ACCUM = df_draw_CC_RANGE)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE_SYNTENY + s(RANK, by = TE_SYNTENY)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      12542.9044     0.2526 49653.6   <2e-16 ***
#   TE_SYNTENYNON_TE  -814.2520     2.8779  -282.9   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(RANK)                  8.6651 8.6667 976924  <2e-16 ***
#   s(RANK):TE_SYNTENYTE     0.6667 0.6667   8281  <2e-16 ***
#   s(RANK):TE_SYNTENYNON_TE 8.5965 8.6622 109144  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 131.23  Scale est. = 130.56    n = 3901


plot_AA_CV + plot_AA_RANGE + plot_CC_CV + plot_CC_RANGE + patchwork::plot_layout(
  ncol = 2, nrow = 2, guides = "collect"
)
ggsave(
  filename = "/Supplement_Data/Accumulate_variaion_All_TE_anchor_gene.pdf",
  width = 7.2,
  height = 3
)
  
