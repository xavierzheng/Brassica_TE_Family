#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#===============================
# AIM
#   summarize CACTA and Helitron all-against-all blast result
#   Because the blast result is tooo big, I need to compress them into feather format by arrow
#   This script is conceptionally identical to TEAge_003_ParserBlastnResult.R
#
#===============================
# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
stop("Usage: TEAge_003_ParserBlastnResult_CACTA.R <PREFIX> <TE superfamily>", call.=FALSE)
}

#setwd("/nfs/project_ssd/project3/pxzhe/TEMP_BLAST_DTC")

print("# Load library ----------------------------------")
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
library(arrow)

setDTthreads(threads = 24)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

VERSION <- as.character(x = args[1])
INPUT_FILE <- as.character(x = args[2])

#VERSION <- "Relaxed"
#INPUT_FILE <- "DTC"

print("# read input ================================================")
temp_list <- list.files(
  path = paste0("/Data/Fig_3_Nucleotide_Identity/Blastn_Relaxed_", INPUT_FILE, "_Small_Table/"),  # Adjusted Location
  pattern = paste0("Blastn_", VERSION, "_IntackTE_", INPUT_FILE, ".fa_p.*") # files
)


# for each quary and subject, only keep one raw =============

print("# Take long time to read data and select the longest alignment---------------------------- ")
list2 <- lapply(temp_list, FUN = function(NAM){
  df_temp <- fread(
    file =  paste0("/Data/Fig_3_Nucleotide_Identity/Blastn_Relaxed_", INPUT_FILE, "_Small_Table/", NAM), # Adjusted Location, and file name .gz
    header = F, sep = "\t", 
    nThread = 20, 
    col.names = c("query_id", "subject_id", "query_length", "subject_length", "identity", "alignment_length")
  )
  
  df_return <- df_temp[
    order(-alignment_length), head(.SD, n = 1), by = .(query_id, subject_id)
  ]
  
  return(df_return)
})

df_temp <- bind_rows(list2)


# create NAME, FAM, CLASS ==================================
df_temp <- as.data.table(df_temp)
df_temp[,
        c("Q_ID", "Q_FAM", "Q_CLASS"):=tstrsplit(query_id, ";",fixed = TRUE)][,
                                                                              c("S_ID", "S_FAM", "S_CLASS"):=tstrsplit(subject_id, ";",fixed = TRUE)]

# only select 3 columns, including c("query_id", "subject_id", "identity") ====================
df_temp2 <- df_temp[order(Q_FAM, Q_ID), c("query_id", "subject_id", "identity")]

rm(df_temp)

# dcast, in order to fill all of comparison as 0 ==============
df_dcast <- dcast.data.table(
  df_temp2, 
  query_id ~ subject_id, value.var = "identity", fill = 0
)

# re-melt the result ==========================================
df_melt <- melt.data.table(
  df_dcast, 
  id.vars = c("query_id"), 
  variable.name = "subject_id", 
  variable.factor = F, 
  value.name = "identity"
)

rm(df_temp2)

# add column name ============================================
df_melt[,
        c("Q_ID", "Q_FAM", "Q_CLASS"):=tstrsplit(query_id, ";",fixed = TRUE)][,
                                                                              c("S_ID", "S_FAM", "S_CLASS"):=tstrsplit(subject_id, ";",fixed = TRUE)]


# mean identity between family ============================
df_temp2 <- df_melt[
  , .(mean(identity)), by = c("Q_FAM", "S_FAM")
][
  ,.(Q_FAM = str_remove_all(Q_FAM, "^Name="), S_FAM = str_remove_all(S_FAM, "^Name="), AVE = V1)
][
  order(Q_FAM, S_FAM),
]

fwrite(
  df_temp2, 
  file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEFamily_", VERSION, "_",INPUT_FILE, ".txt"), 
  col.names = T, sep = "\t", quote = F
)

# name_list <- sort(unique(df_temp2$Q_FAM))

rm(df_temp2)
gc()

print("# estimate mean ID ===========================")

df_temp4 <- df_melt[
  ,.(AVE = mean(identity)), by = c("query_id", "Q_FAM", "S_FAM")
][
  ,.(ifelse(Q_FAM==S_FAM, "Y", "N"), AVE, query_id, Q_FAM=str_remove_all(Q_FAM, "^Name="), S_FAM=str_remove_all(S_FAM, "^Name="))
][
  V1=="Y", c("query_id", "Q_FAM", "S_FAM", "AVE")
]

#glimpse(df_temp4)

fwrite(
  df_temp4, 
  file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEWithinFamily_", VERSION, "_", INPUT_FILE, ".txt"), 
  col.names = T, sep = "\t", quote = F
)

gc()
