#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Usage: draw_cforest_heatmap.R <TE class>", call.=FALSE)
}

#=================================================
# AIM
#	Summary the results of cforest:
#	* Importance and SHAP
#	* Draw feature vs. Expression
#
#=================================================


suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
suppressPackageStartupMessages(library(patchwork, quietly = F))
suppressPackageStartupMessages(library(ggsci, quietly = F))
suppressPackageStartupMessages(library(ComplexHeatmap, quietly = F))
suppressPackageStartupMessages(library(circlize, quietly = F))
library(ggrepel, quietly = T)


setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

INPUT <- as.character(x = args[1])
#INPUT <- "LTR_Copia"

print(paste0("#### RUN ", INPUT, "==========================================="))

print("# read important data -------------------")
df_imp <- fread(
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_importance_", INPUT, ".txt"), 
  header = T, sep = "\t"
)

df_imp_range <- df_imp[
  ,.("IMP05" = quantile(IMPORTANCE, prob = 0.05), 
     "IMP50" = quantile(IMPORTANCE, prob = 0.50),
     "IMP95" = quantile(IMPORTANCE, prob = 0.95)
     )
  , by = .(FEATURE)
][order(-IMP50)]

print("# read SHAP importance, SHapley Additive exPlanations ------------------")
df_SHAP <- fread(
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_SHAP_", INPUT, ".txt"),
  header = T, sep = "\t"
)

df_SHAP <- df_SHAP[,-5][order(-importance)]

colnames(df_SHAP) <- c("FEATURE", "IMP05", "IMP50", "IMP95")

print("# modify feature name, for easier reading ---------")

if(str_detect(INPUT, "LTR")){
 
  temp_feature_ori <- c("CLOSEST_GENE_DISTANCE_AVE", "CLOSEST_GENE_DISTANCE_SD", 
                        "CLOSEST_OTHERTE_DISTANCE_AVE", "CLOSEST_OTHERTE_DISTANCE_SD", 
                        "COUNT", "GC_AVE", "GC_SD", "GENE_CLASS1_TEcontainGENE", 
                        "GENE_CLASS2_TEwithinGENE", "GENE_CLASS3_partial", 
                        "GENE_CLASS4_ComplexEvent", "GENE_CLASS5_NotOverlapWithGene", 
                        "LEN_AVE", "LEN_SD", 
                        "LTR_AGE_AVE", "LTR_AGE_SD", 
                        "LTR_LEN_AVE", "LTR_LEN_SD", 
                        "METHYL_CG_AVE", "METHYL_CG_SD", 
                        "METHYL_CHG_AVE", "METHYL_CHG_SD", 
                        "METHYL_CHH_AVE", "METHYL_CHH_SD", 
                        "OUTER_GC_AVE", "OUTER_GC_SD", 
                        "PREVALENCE_AVE_ALL", "PREVALENCE_SD_ALL", 
                        "TE_CLASS1_TEcontainOTHERTE", "TE_CLASS2_TEwithinOTHERTE", 
                        "TE_CLASS3_partial", "TE_CLASS4_ComplexEvent", "TE_CLASS5_NotOverlapWithOTHERTE")
  
  temp_feature_mod <- c("Closest distance to gene (Mean)", "Closest distance to gene (sd)", 
                        "Closest distance to TE (Mean)",   "Closest distance to TE (sd)", 
                        "Count", "GC content (Mean)", "GC content (sd)", "TE contain gene", 
                        "TE within gene", "TE-gene partial overlap", 
                        "TE-gene complex event", "TE-gene no overlap",
                        "Length (Mean)", "Length (sd)", 
                        "LTR age (Mean)", "LTR age (sd)", 
                        "LTR length (Mean)", "LTR length (sd)", 
                        "CG motif (Mean)", "CG motif (sd)", 
                        "CHG motif (Mean)", "CHG motif (sd)", 
                        "CHH motif (Mean)", "CHH motif (sd)",
                        "Outer GC content (Mean)", "Outer GC content (sd)", 
                        "Prevelance (Mean)", "Prevelance (sd)",
                        "TE contain TE", "TE within TE", 
                        "TE-TE partial overlap", "TE-TE complex event", "TE-TE no overlap")
   
}else{
  
  temp_feature_ori <- c("CLOSEST_GENE_DISTANCE_AVE", "CLOSEST_GENE_DISTANCE_SD", 
                        "CLOSEST_OTHERTE_DISTANCE_AVE", "CLOSEST_OTHERTE_DISTANCE_SD", 
                        "COUNT", "GC_AVE", "GC_SD", "GENE_CLASS1_TEcontainGENE", 
                        "GENE_CLASS2_TEwithinGENE", "GENE_CLASS3_partial", 
                        "GENE_CLASS4_ComplexEvent", "GENE_CLASS5_NotOverlapWithGene", 
                        "LEN_AVE", "LEN_SD", 
                        "METHYL_CG_AVE", "METHYL_CG_SD", 
                        "METHYL_CHG_AVE", "METHYL_CHG_SD", 
                        "METHYL_CHH_AVE", "METHYL_CHH_SD", 
                        "OUTER_GC_AVE", "OUTER_GC_SD", 
                        "PREVALENCE_AVE_ALL", "PREVALENCE_SD_ALL", 
                        "TE_CLASS1_TEcontainOTHERTE", "TE_CLASS2_TEwithinOTHERTE", 
                        "TE_CLASS3_partial", "TE_CLASS4_ComplexEvent", "TE_CLASS5_NotOverlapWithOTHERTE")
  
  temp_feature_mod <- c("Closest distance to gene (Mean)", "Closest distance to gene (sd)", 
                        "Closest distance to TE (Mean)",   "Closest distance to TE (sd)", 
                        "Count", "GC content (Mean)", "GC content (sd)", "TE contain gene", 
                        "TE within gene", "TE-gene partial overlap", 
                        "TE-gene complex event", "TE-gene no overlap",
                        "Length (Mean)", "Length (sd)", 
                        "CG motif (Mean)", "CG motif (sd)", 
                        "CHG motif (Mean)", "CHG motif (sd)", 
                        "CHH motif (Mean)", "CHH motif (sd)",
                        "Outer GC content (Mean)", "Outer GC content (sd)", 
                        "Prevelance (Mean)", "Prevelance (sd)",
                        "TE contain TE", "TE within TE", 
                        "TE-TE partial overlap", "TE-TE complex event", "TE-TE no overlap")
  
}



df_imp_range <- df_imp_range %>%
  mutate(
    FEATURE = factor(
      FEATURE, 
      levels = temp_feature_ori,
      labels = temp_feature_mod
    ),
    FEATURE = as.character(FEATURE)
  )

df_SHAP <- df_SHAP %>%
  mutate(
    FEATURE = factor(
      FEATURE, 
      levels = temp_feature_ori,
      labels = temp_feature_mod
    ),
    FEATURE = as.character(FEATURE)
  )

print("# Check if top10 features is consistant --------")
df_key <- data.frame(
  FEATURE = c(head(df_imp_range$FEATURE, 10), head(df_SHAP$FEATURE, 10))
) %>%
  summarise(
    .by = FEATURE,
    COUNT = n()
  ) %>%
  filter(
    COUNT > 1
  ) %>%
  mutate(
    CONSIST = "Consist"
  ) %>%
  select(
    -COUNT
  )

# save this table
df_key %>%
  mutate(
    TE_CLASS = str_replace(INPUT, "_", "/"), 
  ) %>%
  fwrite(
    file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_ConsistTop10ImportantFeature_", INPUT, ".txt"),
    col.names = T, sep = "\t"
  )

# draw importance ------------------------
print("# prepare function for drawing --------")
draw_importance <- function(DF_IMP, DF_SHARE, AXIS_X){
  
  # joint top 10 features from importance table, and the feature-shared table
  df_draw_temp <- left_join(
    x = DF_IMP %>% slice_head(n = 10),
    y = DF_SHARE,
    by = "FEATURE"
  ) %>%
    mutate(
      CONSIST = ifelse(
        is.na(CONSIST),
        "No", CONSIST
      ),
      CONSIST = factor(
        CONSIST, 
        levels = c("Consist", "No"),
        labels = c("Consist", "No")
      )
    )
  
  # define y axis order
  name_list <- rev(df_draw_temp$FEATURE)
  df_draw <- df_draw_temp %>%
    mutate(
      FEATURE = factor(
        FEATURE,
        levels = name_list, 
        labels = name_list
      )
    )
  
  rm(df_draw_temp)
  
  # plot
  ggplot()+
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "black"
    )+
    geom_point(
      data = df_draw,
      mapping = aes(
        x = IMP50,
        y = FEATURE,
        color = CONSIST
      ),
      size = 1
    )+
    geom_errorbar(
      data = df_draw,
      mapping = aes(
        xmin = IMP05,
        xmax = IMP95,
        y = FEATURE,
        color = CONSIST
      ),
      width = 0.8,
      linewidth = 0.5
    )+
    scale_color_manual(
      breaks = c("Consist", "No"),
      values = c("#de2d26", "#bdbdbd")
    )+
    theme_cowplot()+
    theme(
      axis.title.x = element_text(
        size = 7
      ),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 6
      ),
      axis.text.y = element_text(
        size = 7
      ),
      legend.title = element_text(
        size = 7
      ),
      legend.text = element_text(
        size = 6
      ),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 0, r = 1, b = 0, l = 0, unit = "pt"),
      legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.box.spacing = ggplot2::unit(0, units = "pt")
    )+
    labs(
      x = paste0(AXIS_X, "\n", "importance"),
      color = ""
    ) -> plot_p
  
  return(plot_p)
}

print("# drawing importance ------------------------------------------------")
plot_s <- draw_importance(DF_IMP = df_SHAP, DF_SHARE = df_key, AXIS_X = "SHAP")
plot_p <- draw_importance(DF_IMP = df_imp_range, DF_SHARE = df_key, AXIS_X = "Permutated")

print("## this is importance plot ------------------------------------------")
# plot_p + plot_s + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")

print("# Visual the importance feature vs. expression ----------------------")

print("## read expression ------------------------")
df_exp_TEFAM <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_H37.txt",
  header = T, sep = "\t"
) 

df_exp_TEFAM_top <- df_exp_TEFAM %>%
  filter(
    TE_CLASS %in% str_replace(INPUT, "_", "/")
  ) %>%
  select(
    SAMPLE, TE_CLASS, TE_NAME, NORMALIZE_TPM
  ) %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  arrange(
    .by_group = T,
    desc(NORMALIZE_TPM)
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup()
#glimpse(df_exp_TEFAM_top)

print("# read TE family feature ------------------")
df_fam_feature <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix.txt",
  header = T, sep = "\t"
) %>%
  filter(
    TE_CLASS %in% str_replace(INPUT, "_", "/")
  )

df_feat_exp <- left_join(
  x = df_exp_TEFAM_top,
  y = df_fam_feature,
  by = c("TE_CLASS", "TE_NAME")
) 

rm(df_exp_TEFAM_top)
rm(df_fam_feature)


# define function for drawing features
draw_feature <- function(DF_FEA_EXP, FEATURE, AXIS_X){
  
  # choose to label top 5 TE family
  df_lab <- DF_FEA_EXP %>%
    arrange(
      desc(NORMALIZE_TPM)
    ) %>%
    slice_head(
      n = 5
    )
  
  ggplot()+
    geom_point(
      data = DF_FEA_EXP,
      mapping = aes(
        x = get(FEATURE), 
        y = NORMALIZE_TPM
      ),
      size = 0.5
    )+
    geom_point(
      data = df_lab,
      mapping = aes(
        x = get(FEATURE), 
        y = NORMALIZE_TPM
      ),
      color = "#de2d26",
      size = 0.5
    )+
    ggrepel::geom_text_repel(
      data = df_lab,
      mapping = aes(
        x = get(FEATURE), 
        y = NORMALIZE_TPM,
        label = TE_NAME
      ),
      color = "#de2d26",
      size = 6 * 0.352777778,
      max.time = 1,
      max.overlaps = 20
    )+
    theme_cowplot()+
    theme(
      axis.title = element_text(
        size = 7
      ),
      axis.text = element_text(
        size = 6
      )
    )+
    labs(
      x = AXIS_X,
      y = "Normalized TPM"
    ) -> plot_z
  
  return(plot_z)
}

plot_1 <- draw_feature(DF_FEA_EXP = df_feat_exp, 
                       FEATURE = temp_feature_ori[which(temp_feature_mod == df_imp_range$FEATURE[1])], 
                       AXIS_X = df_imp_range$FEATURE[1])
plot_2 <- draw_feature(DF_FEA_EXP = df_feat_exp, 
                       FEATURE = temp_feature_ori[which(temp_feature_mod == df_imp_range$FEATURE[2])], 
                       AXIS_X = df_imp_range$FEATURE[2])
plot_3 <- draw_feature(DF_FEA_EXP = df_feat_exp, 
                       FEATURE = temp_feature_ori[which(temp_feature_mod == df_imp_range$FEATURE[3])], 
                       AXIS_X = df_imp_range$FEATURE[3])


# print("# this is importance and correlation to expression ----------")
# plot_CForest <- (plot_p + plot_s + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom"))/
# (plot_CHG + plot_GC + plot_CG + plot_layout(ncol = 3, axes = "collect"))


print("# Re-draw TE family expression heatmap -------------")
df_thres <- fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Summary_kallistro_TEFam_NormalizedTPM_threshold_H37.txt",
  header = T, sep = "\t"
)

print("# prepare data frame for heatmap --------------------")
df_temp_top <- df_exp_TEFAM %>%
  filter(
    TE_CLASS %in% str_replace(INPUT, "_", "/")
  ) %>%
  arrange(
    desc(NORMALIZE_TPM)
  ) %>%
  left_join(
    x = .,
    y = df_thres,
    by = "TE_CLASS"
  ) %>%
  filter(
    NORMALIZE_TPM > T95
  ) 

# prepare the order of TE family
df_order <- df_temp_top %>%
  select(
    SAMPLE, TE_CLASS, TE_NAME, NORMALIZE_TPM
  ) %>%
  group_by(
    TE_CLASS, TE_NAME
  ) %>%
  arrange(
    .by_group = T,
    desc(NORMALIZE_TPM)
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  arrange(
    desc(NORMALIZE_TPM)
  )

# prepare data frame for heatmap, and make sure every order is correct
df_draw_heat <- df_temp_top %>%  
  mutate(
    SAMPLE = factor(
      SAMPLE,
      levels = c("H37_RTRoot", "H37_40CRoot", "H37_RTLeaf", "H37_40CLeaf", "H37_RTCurd"),
      labels = c("H37_RTRoot", "H37_40CRoot", "H37_RTLeaf", "H37_40CLeaf", "H37_RTCurd")
    )
  ) %>%
  as.data.table() %>%
  dcast.data.table(
    TE_CLASS + TE_NAME ~ SAMPLE,
    value.var = "NORMALIZE_TPM", fill = 0
  ) %>%
  mutate(
    TE_NAME = factor( # give the order to make sure TE family is order based on its maximum expression
      TE_NAME, 
      levels = df_order$TE_NAME, 
      labels = df_order$TE_NAME
    )
  ) %>%
  arrange( # order the TE family based on the provided order list
    TE_NAME
  )

#glimpse(df_draw_heat)

df_label <- data.frame(
  row.names = c("H37_RTRoot", "H37_40CRoot", "H37_RTLeaf", "H37_40CLeaf", "H37_RTCurd"),
  Tissue = c("Root", "Root", "Leaf", "Leaf", "Curd"),
  Treatment = c("Control", "Heat", "Control", "Heat", "Control")
)
color_list <- list(
  Tissue = c(Root = "#f4a259", Leaf = "#5b8e7d", Curd = "#f4e285"),
  Treatment = c(Control = "#118ab2", Heat = "#ef476f")
)

print("# maximum top 20 TE families were chosen for heatmap ---------")
if(dim(df_draw_heat)[1] >= 20){
  
  df_mat_temp <- df_draw_heat %>% 
    select(-TE_CLASS) %>%
    column_to_rownames(var = "TE_NAME") %>%
    slice_head(n = 20)
  
}else{
  
  df_mat_temp <- df_draw_heat %>% select(-TE_CLASS) %>% column_to_rownames(var = "TE_NAME")
  
}

df_mat1 <- as.matrix(df_mat_temp)

# define color
color_exp <- circlize::colorRamp2(
  breaks = c(min(df_mat1), (min(df_mat1)+ max(df_mat1))/2, max(df_mat1)), 
  colors = c("#FFF8E0" ,"#FFCA27" ,"#FF6E00")
)

df_mat2 <- matrix(data = df_label$Treatment, ncol = 5, nrow = 1)
rownames(df_mat2) <- "Treatment"
plot_heat_treat <- ComplexHeatmap::Heatmap(
  matrix = df_mat2,
  col = c("Control" = "#118ab2", "Heat" = "#ef476f"),
  name = "Treatment",
  cluster_rows = F, 
  cluster_columns = F, 
  show_column_names = F, 
  show_row_names = T,
  row_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))
)

df_mat3 <- matrix(data = df_label$Tissue, ncol = 5, nrow = 1)
rownames(df_mat3) <- "Tissue"
plot_heat_tissue <- ComplexHeatmap::Heatmap(
  matrix = df_mat3, #df_label$Tissue,
  col = c("Root" = "#f4a259", "Leaf" = "#5b8e7d", "Curd" = "#f4e285"),
  name = "Tissue",
  cluster_rows = F, 
  cluster_columns = F, 
  show_column_names = F, 
  show_row_names = T,
  row_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))
)

plot_heat_exp <- ComplexHeatmap::Heatmap(
  matrix = df_mat1, 
  name = "Normalized\nTPM", 
  col = color_exp, 
  cluster_rows = F, 
  cluster_columns = F, 
  show_column_names = F, 
  row_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6))
)

plot_heat <- plot_heat_treat %v% plot_heat_tissue %v% plot_heat_exp


# Convert pheatmap to a ggplot object
pheatmap_grob <- grid::grid.grabExpr(print(plot_heat))

patchwork::wrap_elements(full = pheatmap_grob) +
  (((plot_p + plot_s + plot_layout(ncol = 2, guides = "collect", widths = c(0.8, 1)) & theme(legend.position = "bottom"))/(plot_1 + plot_2 + plot_3 + plot_layout(ncol = 3, axes = "collect"))) +
     plot_layout(widths = c(1.1, 1))) +
  plot_layout(ncol = 2, widths = c(1, 2))

ggsave(
  filename = paste0("/Fig_6_TE_Express/TEExp_Cforest_", INPUT, ".pdf"),
  width = 7.2,
  height = 4
)
