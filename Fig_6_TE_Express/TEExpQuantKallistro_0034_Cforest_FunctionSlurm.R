#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# AIM: cforest function 

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Usage: cforest.R <TE superfamily>", call.=FALSE)
}

CLAS <- as.character(x = args[1])


suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(viridis, quietly = F))
suppressPackageStartupMessages(library(RColorBrewer, quietly = F))
suppressPackageStartupMessages(library(cowplot, quietly = F))
library(party)
library(iml)
library(future) # SHAP 平行運算的基礎
library(future.callr)
library(Metrics) # mse, msa


setDTthreads(threads = 12)
options(dplyr.summarise.inform=FALSE)

df_for_tree <- data.table::fread(
  file = "/Data/Fig_6_TE_Express/TE_Expression/Whole_Family_Feature_Mix_AddExpRank.txt",
  header = T, sep = "\t"
)



# function: hyperparameters ---------------------

hyper_parameter <- function(DF){
  
  temp_list <- parallel::mclapply(seq(from = 10, to =  500, by = 10), 
                                  mc.set.seed = TRUE, 
                                  mc.cores = 12,
                                  FUN = function(NTREE){
    
    lapply(seq(from= 3, to = 10, by = 1), FUN = function(MTRY){
      
      RANDOM <- 4291
      set.seed(RANDOM)
      my_cforest_control <- cforest_control(teststat = "quad",
                                            testtype = "Univ", mincriterion = 0, ntree = NTREE, mtry = MTRY,
                                            replace = FALSE)
      set.seed(RANDOM)
      my_cforest <- cforest(NORMALIZE_EXP ~ ., data = DF, # Need to change
                            controls = my_cforest_control)
      
      df_draw <- DF %>%
        select(
          NORMALIZE_EXP  # Need to change
        ) %>%
        mutate(
          PREDICT = stats::predict(my_cforest, newdata = DF)[,1]
        )
      
      # true r2 value
      rss <- sum((df_draw$NORMALIZE_EXP - df_draw$PREDICT)^2)
      tss <- sum((df_draw$NORMALIZE_EXP - mean(df_draw$NORMALIZE_EXP))^2)
      R2 <- 1 - rss/tss
      
      # Calculate Mean Squared Error (MSE)
      mse_value <- Metrics::mse(df_draw$NORMALIZE_EXP, df_draw$PREDICT)
      
      # Calculate Root Mean Squared Error (RMSE)
      rmse_value <- Metrics::rmse(df_draw$NORMALIZE_EXP, df_draw$PREDICT)
      
      # Calculate Mean Absolute Error (MAE)
      mae_value <- Metrics::mae(df_draw$NORMALIZE_EXP, df_draw$PREDICT)
      
      df_ret <- data.frame(
        NTREE = NTREE,
        MTRY = MTRY,
        R2 = R2,
        MSE = mse_value,
        RMSE = rmse_value,
        MAE = mae_value,
        RANDOM = RANDOM
      )
      
      return(df_ret)
    })
  })
  
  df_hyp_para <- bind_rows(temp_list) %>%
    arrange(
      desc(R2), MSE, RMSE, MAE
    )
  
  return(df_hyp_para)
}

# function: default importance ------------------------------

analysis_importance <- function(INPUT_CFOREST_MODEL){
  
  temp_list <- parallel::mclapply(seq(1, 100), 
                                  mc.set.seed = TRUE, 
                                  mc.cores = 12, 
                                  FUN = function(RANDOM){
    
    set.seed(RANDOM)
    varimp_cforest <- varimp(INPUT_CFOREST_MODEL, nperm = 500, OOB = F, conditional = F)
    # scale
    varimp_cforest_scaled <- sort(varimp_cforest/sum(varimp_cforest)*100)
    
    df_ret <- data.frame(
      FEATURE = names(varimp_cforest_scaled),
      IMPORTANCE = varimp_cforest_scaled
    ) %>%
      mutate(
        RUN = RANDOM
      )
    
    return(df_ret)
    
  })
  
  df_imp <- bind_rows(temp_list)
  
  return(df_imp)
}

# function: SHAP importance ----------------------------------------
analysis_SHAP <- function(INPUT_CFOREST_MODEL, DF){
  
  # Define the predictor object using the fitted model
  predictor <- iml::Predictor$new(
    INPUT_CFOREST_MODEL, 
    data = DF[, -1], 
    y = DF$NORMALIZE_EXP,
    predict.fun = function(object, newdata) predict(object, newdata = newdata, OOB = FALSE) # MUST OOB = F
  )
  
  # Calculate feature importance using SHAP values for the whole dataset
  plan("callr", workers = 12)
  shapley_all <- iml::FeatureImp$new(predictor, loss = "mse", n.repetitions = 100) # n.repetitions = 10
  
  # data for ploting ----------
  return(shapley_all$results)
  
}


print("# usage: -------------------------------")
print(paste0("# ", Sys.time(), " START prepare input ", CLAS))
# rename TE class 
OUT <- str_replace_all(CLAS, "/", "_")

# make sure LTR and non LTR have right features
if(str_detect(CLAS, "LTR")){
  df_temp <- df_for_tree %>%
    filter(
      TE_CLASS %in% CLAS
    ) %>%
    select(-TE_CLASS, -TE_NAME, -NORMALIZE_TPM, -NORMALIZE_EXP_RANK) 
}else{
  df_temp <- df_for_tree %>%
    filter(
      TE_CLASS %in% CLAS
    ) %>%
    select(
      -TE_CLASS, -TE_NAME, -NORMALIZE_TPM, -NORMALIZE_EXP_RANK, -starts_with("LTR")
    ) 
}

print(paste0("# ", Sys.time(), " FINISH prepare input ", CLAS))

# hyper-parameter tuning
print(paste0("# ", Sys.time(), " START hyper-parameter ", CLAS))
df_hyp_para <- hyper_parameter(DF = df_temp)

df_hyp_para <- df_hyp_para %>%
  arrange(
    desc(R2), MSE, RMSE, MAE
  )

fwrite(
  df_hyp_para,
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_hyperparameter_", OUT, ".txt"),
  col.names = T, sep = "\t"
)

print(paste0("# ", Sys.time(), " FINISH hyper-parameter ", CLAS))

NTREE <- df_hyp_para[1, 1]
MTRY <-  df_hyp_para[1, 2]
RANDOM <- df_hyp_para[1, 7]

# cforest 
print(paste0("# ", Sys.time(), " START cforest ", CLAS))
set.seed(RANDOM)
my_cforest_control <- cforest_control(teststat = "quad",
                                      testtype = "Univ", mincriterion = 0, ntree = NTREE, mtry = MTRY,
                                      replace = FALSE)
set.seed(RANDOM)
my_cforest <- cforest(NORMALIZE_EXP ~ ., data = df_temp, # Need to change
                      controls = my_cforest_control)
print(paste0("# ", Sys.time(), " FINISH cforest ", CLAS))

# analysis importance
print(paste0("# ", Sys.time(), " START importance ", CLAS))
df_imp <- analysis_importance(INPUT_CFOREST_MODEL = my_cforest)

fwrite(
  df_imp,
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_importance_", OUT, ".txt"),
  col.names = T, sep = "\t"
)
print(paste0("# ", Sys.time(), " FINISH importance ", CLAS))

# analysis SHAP importance
print(paste0("# ", Sys.time(), " START SHAP ", CLAS))
df_SHAP <- analysis_SHAP(INPUT_CFOREST_MODEL = my_cforest, DF = df_temp)

fwrite(
  df_SHAP,
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_SHAP_", OUT, ".txt"),
  col.names = T, sep = "\t"
)

print(paste0("# ", Sys.time(), " FINISH SHAP ", CLAS))

# model performance
print(paste0("# ", Sys.time(), " START prediction ", CLAS))
df_draw <- df_for_tree %>%
  filter(
    TE_CLASS %in% CLAS
  ) %>%
  select(
    TE_CLASS, TE_NAME, NORMALIZE_EXP # Need to change
  ) %>%
  mutate(
    PREDICT = predict(my_cforest, newdata = df_temp)[,1]
  ) 

fwrite(
  df_draw, 
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_PredictResult_", OUT, ".txt"),
  col.names = T, sep = "\t"
)

print(paste0("# ", Sys.time(), " FINISH prediction ", CLAS))

# Calculate R-square
print(paste0("# ", Sys.time(), " START performance ", CLAS))
rss <- sum((df_draw$NORMALIZE_EXP - df_draw$PREDICT)^2)
tss <- sum((df_draw$NORMALIZE_EXP - mean(df_draw$NORMALIZE_EXP))^2)
r_squared <- 1 - rss/tss

df_performance <- data.frame(
  TE_CLASS = CLAS,
  R_SQUARE = r_squared,
  MSE = Metrics::mse(actual = df_draw$NORMALIZE_EXP, predicted = df_draw$PREDICT),
  RMSE = Metrics::rmse(actual = df_draw$NORMALIZE_EXP, predicted = df_draw$PREDICT),
  MAE = Metrics::mae(actual = df_draw$NORMALIZE_EXP, predicted = df_draw$PREDICT)
)

fwrite(
  df_performance,
  file = paste0("/Data/Fig_6_TE_Express/Cforest_Result/Cforest_PredictPerformance_", OUT, ".txt"),
  col.names = T, sep = "\t"
)

print(paste0("# ", Sys.time(), " FINISH performance ", CLAS))
print(paste0("# ", Sys.time(), " FINISH ", CLAS))

