library(dplyr)
library(stringr)
library(readr)
args <- commandArgs(trailingOnly = TRUE)
input_file <- commandArgs(TRUE)[1]  # Input file containing list of filepaths

process_file <- function(filepath) {
  cat("Processing file:", filepath, "\n")
  tbl<- read.table(filepath, comment.char="", fill = TRUE) %>% filter(str_detect(V1, "#", negate = TRUE)) %>% 
    select(-c(23:30))
  colnames(tbl) <- c("target_name", "blank", "tlen", "query_name", "pfam", "qlen", "evalue", "score1", "bias1", 
                     "#", "of", "c-evalue", "i-evalue", "score2", "bias2", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from",
                     "env_to", "number", "accession")
  tbl <- tbl %>% select(-c("blank", "tlen", "qlen", "number"))
  filtered_info <- tbl %>%
    filter(grepl("PF11999", pfam))
  return(filtered_info)
  }

# Usage example:
results <- process_file(input_file)
output_name <- commandArgs(TRUE)[2]
write.csv(results, file=output_name)
# Print results or perform additional operations
cat("Processing completed.\n")