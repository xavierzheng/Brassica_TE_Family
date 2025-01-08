#==========================
# AIM
#	  1. Summarize FIMO table output, 
#	  2. add FDR-adjusted p value to table, 
#	  3. select proper TFBS by adjusted p value
#
#==========================

gc()
setwd("/nfs/projects/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/02_IntactSingleLine/TE_Age_BLASTN/TF_FIMO")

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))

setDTthreads(threads = 12)
options(tidyverse.quiet=TRUE)
options(dplyr.summarise.inform=FALSE)

INPUT_FILE <- "unknown"

temp_list <- list.files(
  path = ".", 
  pattern = paste0("FIMO_OUT_IntackTE_", INPUT_FILE,".*")
)

list2 <- lapply(temp_list, FUN = function(NAM){
  df <- fread(
    file = NAM, 
    header = T, sep = "\t"
  )
  return(df)
})

df <- bind_rows(list2)
rm(list2)
df <- as.data.table(df)
str(df)

colnames(df) <- c("motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score", "pvalue", "qvalue", "matched_sequence")

# sequence_name, start, stop, strand, score, matched_sequence
#keep_name <- c("sequence_name", "start", "stop", "strand", "score", "matched_sequence")


Sys.time()
df_temp1 <- df[
  ,-c("qvalue")
][
  , qvalue := p.adjust(pvalue, method = "fdr"), 
  by = .(motif_id, motif_alt_id)
][
  qvalue < 1e-5,
][
  , c("INFO", "POS") := tstrsplit(sequence_name, "::", fixed = T)
][
  , c("CHROM", "INFO2") := tstrsplit(POS, ":", fixed = T)
][
  , c("TE_START", "TE_END") := tstrsplit(INFO2, "-", fixed = T)
][
  , "TE_START2" := as.integer(TE_START)
][
  , "TFBS_GFF_START" := (TE_START2 + start )
][
  , "TFBS_GFF_END" := (TE_START2 + stop)
][
  , -c("INFO", "POS", "INFO2","TE_START", "TE_END", "TE_START2")
]
Sys.time()


fwrite(
  df_temp1, 
  file = paste0("FIMO_Qfiltered_", INPUT_FILE, ".txt"), 
  col.names = T, sep = "\t"
)




# df_temp2 <- df_temp1[
#   qvalue < 1e-5,
# ]
# 
# str(df_temp1)
# str(df_temp2)
# 
# length(sort(unique(df_temp1$motif_alt_id)))
# length(sort(unique(df_temp2$motif_alt_id)))


# make sure data table result is identical to dplyr ===============
# df_temp3 <- df %>%
#   group_by(
#     motif_id, motif_alt_id
#   ) %>%
#   mutate(
#     qvalue = p.adjust(pvalue, method = "fdr")
#   ) %>%
#   ungroup() %>%
#   glimpse()
# 
# df_temp4 <- df_temp3 %>%
#   filter(
#     qvalue < 1e-5
#   )
# length(sort(unique(df_temp3$motif_alt_id)))
# length(sort(unique(df_temp4$motif_alt_id)))
# 
# 
# identical(sort(unique(df_temp4$motif_alt_id)) , sort(unique(df_temp2$motif_alt_id)))
#                                            


