#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#=============================
# AIM
#	  Parser the blastn result to get the HSP with hightest alignment length
#
#=============================

# test if there is at least one argument: if not, return an error
if (length(args)!=2) {
  stop("Usage: TEAge_003_ParserBlastnResult.R <VERSION> <INPUT>", call.=FALSE)
}

gc()

#setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
#library(csvread)


setDTthreads(threads = 24)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

#VERSION <- "MoreOut"
VERSION <- as.character(x = args[1])
INPUT_FILE <- as.character(x = args[2])

#VERSION <- "Relaxed"
#INPUT_FILE <- "DTA"

df <- fread(
  cmd = paste0("grep -v '^#' /Data/Fig_3_Nucleotide_Identity/Blastn_", VERSION, "_", INPUT_FILE, ".tab"), # HINT: remeber to unzip the file firstly !!!!!!!!!!
  header = F, sep = "\t"
)


print("# add column name ========================================")
colnames(df) <- c("query_id", "subject_id", "query_length", "subject_length", "subject_strand", 
                  "identity", "alignment_length", "mismatches", "gap_opens", 
                  "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")


print("# for each quary and subject, only keep one raw =============")
df_temp <- df[order(-alignment_length),
              head(.SD,1),
              by = c("query_id", "subject_id")]

rm(df)

print("# create NAME, FAM, CLASS ==================================")
df_temp[,
   c("Q_ID", "Q_FAM", "Q_CLASS"):=tstrsplit(query_id, ";",fixed = TRUE)][,
                                                                         c("S_ID", "S_FAM", "S_CLASS"):=tstrsplit(subject_id, ";",fixed = TRUE)]


#str(df_temp)

df_temp2 <- df_temp[order(Q_FAM, Q_ID), c("query_id", "subject_id", "identity")]

rm(df_temp)
gc()
#str(df_temp2)

df_dcast <- dcast.data.table(
  df_temp2, 
  query_id ~ subject_id, value.var = "identity", fill = 0
)

rm(df_temp2)
gc()
#str(df_dcast)

df_melt <- melt.data.table(
  df_dcast, 
  id.vars = c("query_id"), 
  variable.name = "subject_id", 
  variable.factor = F, 
  value.name = "identity"
)

#str(df_melt)
df_melt[,
        c("Q_ID", "Q_FAM", "Q_CLASS"):=tstrsplit(query_id, ";",fixed = TRUE)][,
                                                                              c("S_ID", "S_FAM", "S_CLASS"):=tstrsplit(subject_id, ";",fixed = TRUE)]

#str(df_melt)
gc()

print("# mean identity between family ============================")
df_temp2 <- df_melt[
  , .(mean(identity)), by = c("Q_FAM", "S_FAM")
][
  ,.(Q_FAM = str_remove_all(Q_FAM, "^Name="), S_FAM = str_remove_all(S_FAM, "^Name="), AVE = V1)
][
  order(Q_FAM, S_FAM),
]

#str(df_temp2)
gc()

fwrite(
  df_temp2, 
  file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEFamily_", VERSION, "_", INPUT_FILE, ".txt"), 
  col.names = T, sep = "\t", quote = F
)


# estimate mean ID ===========================
#rm(df_temp4)
df_temp4 <- df_melt[
  ,.(AVE = mean(identity)), by = c("query_id", "Q_FAM", "S_FAM")
][
  ,.(ifelse(Q_FAM==S_FAM, "Y", "N"), AVE, query_id, Q_FAM=str_remove_all(Q_FAM, "^Name="), S_FAM=str_remove_all(S_FAM, "^Name="))
][
  V1=="Y", c("query_id", "Q_FAM", "S_FAM", "AVE")
]

gc()
#glimpse(df_temp4)

fwrite(
  df_temp4, 
  file = paste0("/Data/Fig_3_Nucleotide_Identity/PerTEWithinFamily_", VERSION, "_", INPUT_FILE, ".txt"), 
  col.names = T, sep = "\t", quote = F
)
