#================================
# AIM
#	  
#	  - make heatmap more beautiful and readable. Mmmmmmm...
#
#================================


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(arrow)

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

# 讓arrow可以多核心讀寫
arrow::set_cpu_count(12)
arrow::set_io_thread_count(12)
arrow::cpu_count()

# read TE label------------------------------------
df_TE_AA <- fread(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/enrichment_TE_syntenic_AA.txt",
  header = T, sep = "\t"
)

df_TE_CC <- fread(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/enrichment_TE_syntenic_CC.txt",
  header = T, sep = "\t"
)

# read expression ---------------------------------
df_exp_AA <- arrow::read_feather(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/Result_kallistro_AA_AddInfo.arrow"
)

df_exp_CC <- arrow::read_feather(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/Result_kallistro_CC_AddInfo.arrow"
)

# read sample info ---------------------------------
df_sample_AA <- fread(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/SRA_describe_AA.txt",
  header = T, sep = "\t"
)

df_sample_CC <- fread(
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/SRA_describe_CC.txt",
  header = T, sep = "\t"
)

# purify data ----------------------------------------
purify_TE <- function(DF_TE, DF_EXP, PREFIX){

  df_temp_TE <- DF_TE %>%
    filter(
      FLORAL %in% "FLORAL"
    ) %>%
    mutate(
      GENEID = paste0(PREFIX, GENEID)
    ) %>%
    select(
      GENEID, TE
    )

  select_vec <- df_temp_TE %>%
    select(
      GENEID
    ) %>%
    pull()

  df_temp <- DF_EXP[
    PREFIX == PREFIX & PANGENOME == "CORE"
  ][
    GENEID %in% select_vec
  ][
    , "LOG2_TPM" := log2(TPM+0.001)
  ][
    ,.(SAMPLE, GENEID, LOG2_TPM, TPM)
  ]

  df_key <- df_temp[
    , .(RANGE = (max(LOG2_TPM)-min(LOG2_TPM)), CV = (sd(LOG2_TPM)/mean(LOG2_TPM))), by = .(GENEID)
  ][
    RANGE > 5 & CV > 0.2,
  ][
    ,.(GENEID, CV)
  ]

  df_ret <- df_temp[
    df_key, on = .(GENEID = GENEID)
  ][
    , "M_TPM" := ifelse(LOG2_TPM <1, 0, LOG2_TPM)
  ] %>%
    dcast.data.table(
      GENEID ~ SAMPLE, value.var = "M_TPM", fill = 0
    ) %>%
    left_join( # 整理格式，補回CV
      x = .,
      y = df_key,
      by = "GENEID"
    ) %>%
    left_join( # 整理GENEID順序，補回TE註解
      x = .,
      y = df_temp_TE,
      by = "GENEID"
    ) %>%
    arrange(
      desc(TE), desc(CV)
    ) %>%
    select(-TE, -CV) %>%
    tibble::column_to_rownames(
      var = "GENEID"
    )

  return(df_ret)
}

# get data----------------------------------------------------------------------
print("# remember the expression unit is log2 TPM ----------------------")
rm(df_select_AA)
rm(df_select_CC)
df_select_AA <- purify_TE(DF_TE = df_TE_AA, DF_EXP = df_exp_AA , PREFIX = "DA_")
df_select_CC <- purify_TE(DF_TE = df_TE_CC, DF_EXP = df_exp_CC , PREFIX = "IJ_")

str(df_select_AA)
str(df_select_CC)

fwrite(
  df_select_AA, 
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneExpLog2TPM_AA_filtered.txt",
  col.names = T, row.names = T, sep = "\t"
)

fwrite(
  df_select_CC, 
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneExpLog2TPM_CC_filtered.txt",
  col.names = T, row.names = T, sep = "\t"
)

rm(df_exp_AA)
rm(df_exp_CC)
gc()

print("# Estimate CV and Range -----------------")

estimate_diff <- function(DF){

  df_ret <- DF %>%
    tibble::rownames_to_column(
      var = "GENEID"
    ) %>%
    as.data.table() %>%
    data.table::melt.data.table(
      id.vars = "GENEID",
      variable.name = "SAMPLE",
      variable.factor = F, 
      value.name = "M_TPM"
    ) %>%
    summarise(
      .by = GENEID,
      CV = sd(M_TPM, na.rm = T)/mean(M_TPM, na.rm = T), 
      RANGE = (max(M_TPM, na.rm = T) - min(M_TPM, na.rm = T))
    ) %>%
    arrange(
      desc(CV), desc(RANGE)
    ) %>%
    tibble::column_to_rownames(
      var = "GENEID"
    )
  
  return(df_ret)
}

df_diff_AA <- estimate_diff(DF = df_select_AA)
df_diff_CC <- estimate_diff(DF = df_select_CC)

str(df_diff_AA)
str(df_diff_CC)

print("# perfom z scaling log2 TPM per gene -------------------------")

normaliz_Z <- function(DF){
  
  df_order <- data.frame(
    GENEID = rownames(DF)
  )
  
  df_melt <- DF %>%
    tibble::rownames_to_column(
      var = "GENEID"
    ) %>%
    as.data.table() %>%
    data.table::melt.data.table(
      id.vars = "GENEID", 
      variable.name = "SAMPLE", 
      variable.factor = F,
      value.name = "M_TPM"
    ) 
  
  df_Z<- df_melt[
    , "Z_SCAL" := (M_TPM - mean(M_TPM, na.rm = T))/sd(M_TPM, na.rm = T), by = .(GENEID)
  ]
  
  df_ret <- df_Z %>%
    mutate(
      Z_SCAL = ifelse(
        is.infinite(Z_SCAL), # make sure 1/0 --> become -3
        -4, Z_SCAL
      ), 
      Z_SCAL = ifelse(
        is.na(Z_SCAL), # make sure NA --> become -3
        -4, Z_SCAL
      )
    ) %>%
    select(
      -M_TPM # remove original log2 TPM
    ) %>%
    as.data.table() %>%
    data.table::dcast.data.table(
      GENEID ~ SAMPLE, value.var = "Z_SCAL", fill = -4
    ) %>%
    left_join(
      x = df_order, # 按照CV大到小排列
      y = .,
      by = "GENEID"
    ) %>%
    tibble::column_to_rownames(
      var = "GENEID"
    )
    
  return(df_ret)
  
}

df_select_Z_AA <- normaliz_Z(df_select_AA)
df_select_Z_CC <- normaliz_Z(df_select_CC)

glimpse(df_select_AA)
glimpse(df_select_Z_AA)
str(df_select_Z_CC)
str(df_select_CC)

fwrite(
  df_select_Z_AA, 
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneExpLog2TPM_ZScale_AA_filtered.txt",
  col.names = T, row.names = T, sep = "\t"
)

fwrite(
  df_select_Z_CC, 
  file = "/Data/Fig_8_Perturb_Gene_Express/Floral_Gene_Express/GeneExpLog2TPM_ZScale_CC_filtered.txt",
  col.names = T, row.names = T, sep = "\t"
)

# Start for drawing complexHeatMap ------------------
print("# Start for drawing complexHeatMap ---------------------------")

# Use kmean and complexHeatmap--------------------------------------------------------------------------
print("# define heatmap color code ---------------------------------")
packageVersion("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(ggsci)

print("# find color code")
ggsci::pal_material(palette = "blue")(10)[c(1,8,10)]
ggsci::pal_material(palette = "deep-orange")(10)[c(1,8,10)]


print("# find double color palette")
brewer.pal(5, "RdBu")
brewer.pal(5, "PuOr")

c("#bfd7ed", "#60a3d9","#0074b7","#003b73")
c("#fff3d9", "#dc4731", "#b8390e", "#3b0918")


N_col_fun <- circlize::colorRamp2(breaks = c(-3, 0, 1, 4), colors = c("#d4f1f4", "#60a3d9","#0074b7","#050a30"))
T_col_fun <- circlize::colorRamp2(breaks = c(-3, 0, 1, 4), colors = c("#fff3d9", "#dc4731", "#b8390e", "#3b0918"))

CV_col_fun <- circlize::colorRamp2(breaks = c(0, 1, 5), colors = c("#c3ceda", "#738fa7", "#0c4160"))
RANGE_col_fun <- circlize::colorRamp2(breaks = c(0, 5, 10), colors = c("#96D3D3", "#429E9D", "#193D3C"))

print("# prepare gene label ----------------------------------------")
df_row_C <- data.frame(
  GENEID = rownames(df_select_CC)
) %>%
  left_join(
    x =.,
    y = df_TE_CC %>% mutate(GENEID = paste0("IJ_", GENEID)),
    by = "GENEID"
  )
df_row_A <- data.frame(
  GENEID = rownames(df_select_AA)
) %>%
  left_join(
    x =.,
    y = df_TE_AA %>% mutate(GENEID = paste0("DA_", GENEID)),
    by = "GENEID"
  )

# glimpse(df_row_A)
# glimpse(df_row_C)

print("# define function to pick up gene closed to TE or without nearby TE ------")
select_TE <- function(DF, DF_ROW, KEY){
  DF %>%
    rownames_to_column(
      var = "GENEID"
    ) %>%
    left_join(
      x = .,
      y = DF_ROW %>% select(GENEID, TE),
      by = "GENEID"
    ) %>%
    filter(
      TE %in% KEY
    ) %>%
    column_to_rownames(
      var = "GENEID"
    ) %>%
    select(
      -TE
    ) -> df_ret
  return(df_ret)
}

print("# using function, remember to make sure using normalized expression or raw --------")
print("## first estimate the range and cv to determine the correct order")
df_diff_CC_TE <- select_TE(DF = df_diff_CC, DF_ROW = df_row_C, KEY = "TE") %>%
  arrange(desc(CV))

df_diff_CC_nonTE <- select_TE(DF = df_diff_CC, DF_ROW = df_row_C, KEY = "NON_TE") %>%
  arrange(desc(CV))

df_diff_AA_TE <- select_TE(DF = df_diff_AA, DF_ROW = df_row_A, KEY = "TE") %>%
  arrange(desc(CV))

df_diff_AA_nonTE <- select_TE(DF = df_diff_AA, DF_ROW = df_row_A, KEY = "NON_TE") %>%
  arrange(desc(CV))

print("## then, select TE and non-TE log2 TPM expression---------")
df_CC_TE <- select_TE(DF = df_select_Z_CC, DF_ROW = df_row_C, KEY =  "TE") %>%
  tibble::rownames_to_column(var = "GENEID") %>%
  left_join(
    x = data.frame(GENEID = rownames(df_diff_CC_TE)),
    y = .,
    by = "GENEID"
  ) %>%
  tibble::column_to_rownames(var = "GENEID")

df_CC_nonTE <- select_TE(DF = df_select_Z_CC, DF_ROW = df_row_C, KEY =  "NON_TE") %>%
  tibble::rownames_to_column(var = "GENEID") %>%
  left_join(
    x = data.frame(GENEID = rownames(df_diff_CC_nonTE)),
    y = .,
    by = "GENEID"
  ) %>%
  tibble::column_to_rownames(var = "GENEID")

df_AA_TE <- select_TE(DF = df_select_Z_AA, DF_ROW = df_row_A, KEY =  "TE") %>%
  tibble::rownames_to_column(var = "GENEID") %>%
  left_join(
    x = data.frame(GENEID = rownames(df_diff_AA_TE)),
    y = .,
    by = "GENEID"
  ) %>%
  tibble::column_to_rownames(var = "GENEID")

df_AA_nonTE <- select_TE(DF = df_select_Z_AA, DF_ROW = df_row_A, KEY =  "NON_TE") %>%
  tibble::rownames_to_column(var = "GENEID") %>%
  left_join(
    x = data.frame(GENEID = rownames(df_diff_AA_nonTE)),
    y = .,
    by = "GENEID"
  ) %>%
  tibble::column_to_rownames(var = "GENEID")
  


# Draw CC complexHeatmap and kmean--------------------------------------------------------
print("# prepare column label ----------------------------------------------")
df_col <- data.frame(
  SAMPLE = colnames(df_CC_TE)
) %>%
  left_join(
    x = .,
    y = df_sample_CC,
    by = "SAMPLE"
  )
glimpse(df_col)

matr_sample <- matrix(data = df_col$SimpleType, nrow = 1)
colnames(matr_sample) <- df_col$SAMPLE

print("# HINT: color should be vector, not list in complexHeatmap ----------")
color_sample <- c(
  "Flower" = "#ef476f",
  "Leafy" = "#118ab2",
  "Other" = "#ffd166"
)

print("# draw sample label ------------------------------------------------")
H_sample <- ComplexHeatmap::Heatmap(
  matrix = matr_sample, 
  col = color_sample, 
  name = "Sample",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H_sample

print("# draw gene closed to TE -------------------------------------------")
## left annotation
ha_CC_TE <- ComplexHeatmap::rowAnnotation(
  CV = as.matrix(df_diff_CC_TE %>% select(CV)),
  RANGE = as.matrix(df_diff_CC_TE %>% select(RANGE)),
  col = list(
    CV = CV_col_fun, 
    RANGE = RANGE_col_fun
  ),
  simple_anno_size = unit(5, "pt"), # Use simple_anno_size instead of width for controlling tile size
  show_annotation_name = FALSE,  # This disables the display of annotation names
  annotation_legend_param = list(
    # For the 'CV' annotation legend
    CV = list(
      title_gp = gpar(fontsize = 7),   # Legend title font size
      labels_gp = gpar(fontsize = 6)   # Legend labels font size
    ),
    # For the 'RANGE' annotation legend
    RANGE = list(
      title_gp = gpar(fontsize = 7),
      labels_gp = gpar(fontsize = 6)
    )
  )
)

ha_CC_nonTE <- ComplexHeatmap::rowAnnotation(
  CV = as.matrix(df_diff_CC_nonTE %>% select(CV)),
  RANGE = as.matrix(df_diff_CC_nonTE %>% select(RANGE)),
  col = list(
    CV = CV_col_fun, 
    RANGE = RANGE_col_fun
  ),
  simple_anno_size = unit(5, "pt"), # Use simple_anno_size instead of width for controlling tile size
  show_annotation_name = FALSE,  # This disables the display of annotation names
  annotation_legend_param = list(
    # For the 'CV' annotation legend
    CV = list(
      title_gp = gpar(fontsize = 7),   # Legend title font size
      labels_gp = gpar(fontsize = 6)   # Legend labels font size
    ),
    # For the 'RANGE' annotation legend
    RANGE = list(
      title_gp = gpar(fontsize = 7),   
      labels_gp = gpar(fontsize = 6)   
    )
  )
)


set.seed(1234)
H1 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(df_CC_TE),
  show_row_names = F,
  show_column_names = F,
  name = "Floral genes\nclose to TE\nZ-scaled TPM",
  col = T_col_fun,
  cluster_rows = F, 
  show_row_dend = F,
  left_annotation = ha_CC_TE,
  column_km = 5,
  column_km_repeats = 100,
  cluster_columns = T,
  show_column_dend = F,
  clustering_distance_columns = "spearman",
  column_gap = unit(2, "pt"),
  row_gap = unit(0, "pt"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H1

set.seed(1234)
H2 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(df_CC_nonTE),
  show_row_names = F,
  show_column_names = F,
  name = "Floral genes\nwithout TE\nZ-scaled TPM",
  col = N_col_fun,
  cluster_rows = F,
  show_row_dend = F,
  left_annotation = ha_CC_nonTE,
  cluster_columns = T,
  show_column_dend = F,
  clustering_distance_columns = "spearman",
  column_gap = unit(1, "pt"),
  row_gap = unit(0, "pt"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H2

set.seed(1234)
plot_list <- H_sample %v% H1 %v% H2

# draw and saving CC-------------------------------------
pdf(file = "/Fig_8_Perturb_Gene_Express/Heatmap_FloralTE_kmean_ZScale_CC.pdf", width = 4, height = 4)
set.seed(1234)
ComplexHeatmap::draw(plot_list, main_heatmap = "Floral genes\nclose to TE\nZ-scaled TPM")
dev.off()

#==================================================================
#==================================================================
#==================================================================
# THEN, AA ========================================================
#==================================================================
#==================================================================
#==================================================================

# Draw AA complexHeatmap and kmean--------------------------------------------------------
rm(df_col)
df_col <- data.frame(
  SAMPLE = colnames(df_AA_TE)
) %>%
  left_join(
    x = .,
    y = df_sample_AA,
    by = "SAMPLE"
  )
glimpse(df_col)

matr_sample <- matrix(data = df_col$SimpleType, nrow = 1)
colnames(matr_sample) <- df_col$SAMPLE

color_sample <- c(
  "Flower" = "#ef476f",
  "Leafy" = "#118ab2",
  "Others" = "#ffd166"
)

print("# draw sample label ------------------------------------------------")
H_sample <- ComplexHeatmap::Heatmap(
  matrix = matr_sample, 
  col = color_sample, 
  name = "Sample",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H_sample

print("# draw gene closed to TE -------------------------------------------")
## left annotation
ha_AA_TE <- ComplexHeatmap::rowAnnotation(
  CV = as.matrix(df_diff_AA_TE %>% select(CV)),
  RANGE = as.matrix(df_diff_AA_TE %>% select(RANGE)),
  col = list(
    CV = CV_col_fun, 
    RANGE = RANGE_col_fun
  ),
  simple_anno_size = unit(5, "pt"), # Use simple_anno_size instead of width for controlling tile size
  show_annotation_name = FALSE,  # This disables the display of annotation names
  annotation_legend_param = list(
    # For the 'CV' annotation legend
    CV = list(
      title_gp = gpar(fontsize = 7),   # Legend title font size
      labels_gp = gpar(fontsize = 6)   # Legend labels font size
    ),
    # For the 'RANGE' annotation legend
    RANGE = list(
      title_gp = gpar(fontsize = 7),   
      labels_gp = gpar(fontsize = 6)   
    )
  )
)

ha_AA_nonTE <- ComplexHeatmap::rowAnnotation(
  CV = as.matrix(df_diff_AA_nonTE %>% select(CV)),
  RANGE = as.matrix(df_diff_AA_nonTE %>% select(RANGE)),
  col = list(
    CV = CV_col_fun, 
    RANGE = RANGE_col_fun
  ),
  simple_anno_size = unit(5, "pt"), # Use simple_anno_size instead of width for controlling tile size
  show_annotation_name = FALSE,  # This disables the display of annotation names
  annotation_legend_param = list(
    # For the 'CV' annotation legend
    CV = list(
      title_gp = gpar(fontsize = 7),   # Legend title font size
      labels_gp = gpar(fontsize = 6)   # Legend labels font size
    ),
    # For the 'RANGE' annotation legend
    RANGE = list(
      title_gp = gpar(fontsize = 7),   
      labels_gp = gpar(fontsize = 6)   
    )
  )
)


set.seed(1234)
H1 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(df_AA_TE),
  show_row_names = F,
  show_column_names = F,
  name = "Floral genes\nclose to TE\nZ-scaled TPM",
  col = T_col_fun,
  cluster_rows = F, 
  show_row_dend = F,
  left_annotation = ha_AA_TE,
  column_km = 5,
  column_km_repeats = 100,
  cluster_columns = T,
  show_column_dend = F,
  clustering_distance_columns = "spearman",
  column_gap = unit(2, "pt"),
  row_gap = unit(0, "pt"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H1

set.seed(1234)
H2 <- ComplexHeatmap::Heatmap(
  matrix = as.matrix(df_AA_nonTE),
  show_row_names = F,
  show_column_names = F,
  name = "Floral genes\nwithout TE\nZ-scaled TPM",
  col = N_col_fun,
  cluster_rows = F,
  show_row_dend = F,
  left_annotation = ha_AA_nonTE,
  cluster_columns = T,
  show_column_dend = F,
  clustering_distance_columns = "spearman",
  column_gap = unit(1, "pt"),
  row_gap = unit(0, "pt"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7),  # Legend title font size
    labels_gp = gpar(fontsize = 6)  # Legend labels font size
  )
)
H2

set.seed(1234)
plot_list <- H_sample %v% H1 %v% H2

# draw and saving AA-------------------------------------
pdf(file = "/Fig_8_Perturb_Gene_Express/Heatmap_FloralTE_kmean_ZScale_AA.pdf", width = 4, height = 4)
set.seed(1234)
ComplexHeatmap::draw(plot_list, main_heatmap = "Floral genes\nclose to TE\nZ-scaled TPM")
dev.off()


#=================================================
#=================================================
#=================================================
# re-draw accumulation plot, decrease the width to 3 inch
#=================================================
#=================================================
#=================================================

glimpse(df_diff_AA_TE)

df_accu_AA <- rbind(df_diff_AA_TE %>% mutate(TE = "TE"), df_diff_AA_nonTE %>% mutate(TE = "Non_TE")) %>%
  mutate(
    TE = factor(
      TE, 
      levels = c("TE", "Non_TE"),
      labels = c("TE", "Non_TE")
    )
  )

df_accu_CC <- rbind(df_diff_CC_TE %>% mutate(TE = "TE"), df_diff_CC_nonTE %>% mutate(TE = "Non_TE")) %>%
  mutate(
    TE = factor(
      TE, 
      levels = c("TE", "Non_TE"),
      labels = c("TE", "Non_TE")
    )
  )

glimpse(df_accu_AA)
glimpse(df_accu_CC)

# prepare function for accumulate plot------------------------
plot_accumulate <- function(DF, VAR){
  
  DF %>%
    dplyr::group_by(
      TE
    ) %>%
    dplyr::arrange(
      desc(get(VAR))
    ) %>%
    dplyr::mutate(
      RANK = 1:n()
    ) %>%
    dplyr::mutate(
      CUM = cumsum(get(VAR))
    ) %>%
    ungroup() %>%
    ggplot(
      aes(
        x = RANK,
        y = CUM,
        color = TE,
        group = TE
      )
    )+
    geom_point(
      size = 1.7
    )+
    geom_line()+
    scale_color_manual(
      breaks = c("TE", "Non_TE"),
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
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.position.inside = c(0.7, 0.2)
    )+
    labs(
      x = "Ranked geneid",
      y = paste0("Accumulated ", VAR)
    ) -> plot_ret
  return(plot_ret)
}

plot_CV_CC <- plot_accumulate(DF = df_accu_CC, VAR = "CV")
plot_RANGE_CC <- plot_accumulate(DF = df_accu_CC, VAR = "RANGE")

plot_CV_AA <- plot_accumulate(DF = df_accu_AA, VAR = "CV")
plot_RANGE_AA <- plot_accumulate(DF = df_accu_AA, VAR = "RANGE")

library(patchwork)

plot_CV_CC + plot_RANGE_CC + patchwork::plot_layout(nrow = 2) + patchwork::plot_layout(guides = "collect")
ggsave2(
  filename = "/Fig_8_Perturb_Gene_Express/Accumulate_variaion_CC.pdf",
  width = 3.5,
  height = 3.7
)

plot_CV_AA + plot_RANGE_AA + patchwork::plot_layout(nrow = 2) + patchwork::plot_layout(guides = "collect")
ggsave2(
  filename = "/Fig_8_Perturb_Gene_Express/Accumulate_variaion_AA.pdf",
  width = 3.5,
  height = 3.7
)

# statistic analysis ========================
library(mgcv)

model_gam <- function(DF, VAR){
  
  df_temp_stat <- DF %>%
    dplyr::group_by(
      TE
    ) %>%
    dplyr::arrange(
      desc(get(VAR))
    ) %>%
    dplyr::mutate(
      RANK = 1:n()
    ) %>%
    dplyr::mutate(
      CUM = cumsum(get(VAR))
    ) %>%
    ungroup()
  
  gam_model_ <- mgcv::gam(CUM ~ s(RANK) + TE + s(RANK, by = TE), data = df_temp_stat)
  
  return(gam_model_)
}

gam_model_CV_CC <- model_gam(DF = df_accu_CC, VAR = "CV")
gam_model_RANGE_CC <- model_gam(DF = df_accu_CC, VAR = "RANGE")

summary(gam_model_CV_CC)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE + s(RANK, by = TE)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  18.73644    0.01598 1172.22   <2e-16 ***
#   TENon_TE    -10.63086    0.26297  -40.43   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F p-value    
# s(RANK)          7.8179 8.2065 272.885  <2e-16 ***
#   s(RANK):TETE     2.4948 3.0270 220.877  <2e-16 ***
#   s(RANK):TENon_TE 0.6667 0.6667   1.132    0.39    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 0.011433  Scale est. = 0.008735  n = 55


summary(gam_model_RANGE_CC)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE + s(RANK, by = TE)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  95.20611    0.01563 6092.05   <2e-16 ***
#   TENon_TE    -28.07002    1.39839  -20.07   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df       F p-value    
# s(RANK)          8.4818 8.6529 18517.7  <2e-16 ***
#   s(RANK):TETE     0.6667 0.6667   521.4  <2e-16 ***
#   s(RANK):TENon_TE 4.0826 4.4589   735.1  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 0.011539  Scale est. = 0.0083437  n = 55



gam_model_CV_AA <- model_gam(DF = df_accu_AA, VAR = "CV")
gam_model_RANGE_AA <- model_gam(DF = df_accu_AA, VAR = "RANGE")

summary(gam_model_CV_AA)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE + s(RANK, by = TE)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 23.090469   0.008474 2724.81   <2e-16 ***
#   TENon_TE    -1.373168   0.019711  -69.66   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(RANK)          7.5382 7.9733 5544.0  <2e-16 ***
#   s(RANK):TETE     8.2204 8.4553  485.4  <2e-16 ***
#   s(RANK):TENon_TE 0.6667 0.6667  108.7  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 0.0045428  Scale est. = 0.0036973  n = 99

summary(gam_model_RANGE_AA)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   CUM ~ s(RANK) + TE + s(RANK, by = TE)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 148.90535    0.03231  4608.5   <2e-16 ***
#   TENon_TE    -14.40565    0.09419  -152.9   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df      F p-value    
# s(RANK)          5.003  5.173 190.86  <2e-16 ***
#   s(RANK):TETE     7.769  7.875  61.98  <2e-16 ***
#   s(RANK):TENon_TE 4.156  4.359  19.14  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Rank: 28/29
# R-sq.(adj) =      1   Deviance explained =  100%
# GCV = 0.066455  Scale est. = 0.05375   n = 99