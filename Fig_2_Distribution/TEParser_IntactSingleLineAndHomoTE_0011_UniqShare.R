#=====================
# AIM
#	Generate unique or share family info
#
#=====================

# setwd("/nfs/project1/Brassica_TE/1000_plot_TotalTE/03_ParserTE/04_IntactSingleLineAndHomoTE")

suppressPackageStartupMessages(library(tidyverse, quietly = F))
suppressPackageStartupMessages(library(data.table, quietly = F))
suppressPackageStartupMessages(library(UpSetR))
setDTthreads(threads = 4)
options(dplyr.summarise.inform=FALSE)

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements  #只有多這一行
  return(data)
}


df <- fread(
  file = "/Data/Fig_2_Distribution/TE_INFO_AA_CC.txt.gz", 
  header = T, sep = "\t"
) %>%
  filter(
    # remove NA
    !is.na(TE_NAME)
  ) %>%
  filter(
    # remove strange TE family, they should be NA
    nchar(TE_NAME) > 1
  )

fwrite(
  df,
  file = "/Data/Fig_2_Distribution/TE_INFO_NoTEFamilyNA_AA_CC.txt",
  col.names = T, row.names = F, sep = "\t", quote = F
)


print("# generate unique list, and save out ==================")

name_list <- c("DNA/Helitron", "LTR/Copia", "LTR/Gypsy", "LTR/unknown", "TIR/DTA", "TIR/DTC", "TIR/DTH", "TIR/DTM", "TIR/DTT")

lapply(name_list, FUN = function(NAM){

  OUT_PREFIX <- str_replace(NAM, "/", "_")

  df_temp <-df %>%
    filter(
      TE_CLASS %in% NAM
    )

  ## family found in AA
  df_temp %>%
    filter(
      str_detect(PLANT_ID, "AA$")
    ) %>%
    select(
      TE_NAME
    ) %>%
    unique() %>%
    fwrite(
      file = paste0("TEMP_", OUT_PREFIX, "_AA.txt"),
      col.names = F
    )

  ## family found in CC
  df_temp %>%
    filter(
      str_detect(PLANT_ID, "CC$")
    ) %>%
    select(
      TE_NAME
    ) %>%
    unique() %>%
    fwrite(
      file = paste0("TEMP_", OUT_PREFIX, "_CC.txt"),
      col.names = F
    )

})

print("# using fromList to modify table=================")

temp_list <- lapply(name_list, FUN = function(NAM){
  
  OUT_PREFIX <- str_replace(NAM, "/", "_")
  
  AA<- fread(
    file = paste0("TEMP_", OUT_PREFIX, "_AA.txt"), header = F, sep = "\t"
  ) %>% select(V1) %>% pull
  
  CC <- fread(
    file = paste0("TEMP_", OUT_PREFIX, "_CC.txt"), header = F, sep = "\t"
  ) %>% select(V1) %>% pull
  
  temp_list <- list(AA, CC)
  
  names(temp_list) <- c("AA", "CC")
  
  df_temp <- fromList(temp_list)
  
  df_temp <- rownames_to_column(df_temp, var = "TE_Family")
  
  df_out <- df_temp %>%
    mutate(
      Share_or_Unique = ifelse(
        AA==1&CC==1, 
        "Share_in_AA_CC", 
        ifelse(
          AA==1&CC==0,
          "Unique_in_AA", 
          ifelse(
            AA==0&CC==1, 
            "Unique_in_CC", 
            "ERROR"
          )
        )
      )
    ) %>%
    mutate(
      TE_superfamily = NAM
    ) %>%
    select(
      TE_superfamily, TE_Family, Share_or_Unique
    )
  
  return(df_out)
  
})

df_ShareUniq <- bind_rows(temp_list)

glimpse(df_ShareUniq)

fwrite(
  df_ShareUniq, 
  file = "/Data/Fig_2_Distribution/TE_family_uniqueness.txt", 
  col.names = T, row.names = F, sep = "\t", quote = F
)


# print("# Overall descriptive table")
# df_ShareUniq %>%
#   group_by(
#     TE_superfamily, Share_or_Unique
#   ) %>%
#   summarise(
#     COUNT = n()
#   ) %>%
#   ungroup() %>%
#   as.data.table() %>%
#   data.table::dcast.data.table(
#     TE_superfamily ~ Share_or_Unique, value.var = "COUNT"
#   )
