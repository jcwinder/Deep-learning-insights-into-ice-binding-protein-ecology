# trim alignment, realign, filter for conserved regions
library(tidyverse)
library(Biostrings)
setwd("")

count_start_gaps <- function(alignment){
  split_align <- strsplit(as.character(alignment), "[A-Z]+[-]+")[[1]]
  split_align_1 <- split_align[nzchar(split_align)]
  length_first_gaps <- ifelse(length(split_align_1)==0, 0, nchar(split_align_1))
  length_df <- cbind(length_first_gaps, names(alignment))
  return(length_df)
}
count_end_gaps <- function(alignment){
  split_align <- strsplit(as.character(alignment), "[-]+[A-Z]+")[[1]]
  split_align_1 <- split_align[nzchar(split_align)]
  length_end_gaps <- ifelse(length(split_align_1)==0, 0, nchar(split_align_1))
  length_df <- cbind(length_end_gaps, names(alignment))
  return(length_df)
}

dynamic_alignment_trim <- function(unaligned_seqs_fp, global_align_seqs_fp){
  unaligned_seqs <- readAAStringSet(unaligned_seqs_fp)
  global_aligned_seqs <- global_align_seqs_fp 
  print("Counting start and end gaps")
  #Count start and end gaps: 
  pb = txtProgressBar(min = 0, max = length(sed_align), initial = 0) 
  gaps_df <- data.frame()
  for(i in 1:length(global_aligned_seqs)){
    alignment <- global_aligned_seqs[i]
    start_gaps <- count_start_gaps(alignment)
    end_gaps <- count_end_gaps(alignment)
    out <- as.data.frame(cbind(start_gaps, end_gaps))
    gaps_df <- rbind(gaps_df, out)
    setTxtProgressBar(pb,i)
    close(pb)
  }
  gaps_df$length_first_gaps <- as.numeric(gaps_df$length_first_gaps)
  gaps_df$length_end_gaps <- as.numeric(gaps_df$length_end_gaps)
  freq_start_gaps <- as.data.frame(table(gaps_df$length_first_gaps)) %>% arrange(desc(Freq))
  freq_end_gaps <- as.data.frame(table(gaps_df$length_end_gaps)) %>% arrange(desc(Freq))
  
  #Determine trim value
  print(paste("Most frequent start gaps"))
  print(head(freq_start_gaps))
  trim_start_value <- readline(prompt="Enter start trim value: ")
  trim_start_value <- as.numeric(trim_start_value)
  print(paste("Most frequent end gaps"))
  print(head(freq_end_gaps))
  trim_end_value <- readline(prompt="Enter end trim value: ")
  trim_end_value <- as.numeric(trim_end_value)
  #Trimming
  print("Dynamically trimming sequences based on alignment gappyness")
  df_totrim <- gaps_df %>% filter(length_first_gaps <= trim_start_value | length_end_gaps <= trim_end_value)
  df_nottrim <- gaps_df %>% 
    filter((length_first_gaps > trim_start_value & length_first_gaps < 300) & 
             (length_end_gaps > trim_end_value & length_end_gaps < 100))
  #print(df_totrim$V4)
  seqs_totrim <- unaligned_seqs[which(names(unaligned_seqs) %in% df_totrim$V2)]
  seqs_nottrim <- unaligned_seqs[which(names(unaligned_seqs) %in% df_nottrim$V2)]
  trimmed_seqs <- AAStringSet()
  for(i in 1:length(seqs_totrim)){
    sequence <- seqs_totrim[i]
    sequence_cr <- as.character(sequence)
    length_first_gaps <- df_totrim$length_first_gaps[i]
    length_end_gaps <- df_totrim$length_end_gaps[i]
    sub_sequence_st <- str_sub(sequence_cr, start = (trim_start_value-length_first_gaps),end = nchar(sequence_cr))
    sub_sequence_st <- AAStringSet(sub_sequence_st)
    sequence_pe <- as.character(sub_sequence_st)
    sub_sequence_end <- str_sub(sequence_pe, start = 0, end = (nchar(sequence_pe)-(trim_end_value-length_end_gaps)))
    names(sub_sequence_end) <- names(sequence)
    trimmed_seqs <- c(trimmed_seqs, sub_sequence_end)
  }
  trimmed_full_seqs <- AAStringSet(c(trimmed_seqs, seqs_nottrim))
  #Remove those which are longer than 300 AA and shorter than 90
  trimmed_filt_seq <- trimmed_full_seqs[which(width(trimmed_full_seqs)<300 & width(trimmed_full_seqs)>80)]
  return(trimmed_filt_seq)
}

# pre trim to remove positions which are > 80% gaps
tst_4c <- readAAStringSet("/path/to/folder/pf11999_magusehmm_alignment_masked.fasta")
tst_4c <- tst_4c[which(names(tst_4c) %in% pf11999_lbl7$unique_name)]


process_column2 <- function(i, align) {
  aa <- as.character(subseq(align, start = i, end = i))
  tbl <- as.data.frame(table(aa))
  sumfreq <- sum(tbl$Freq)
  prop_gap <- ifelse(any(tbl$aa == "-"), tbl$Freq[which(tbl$aa == "-")]/sumfreq, 0)
  tbl_ng <- tbl %>% filter(aa != "-")
  if(nrow(tbl_ng)==0)
  {
    ab_val = "-"
    div_val = 1
  }
  else {
  ab_val <- tbl_ng %>% filter(Freq == max(Freq)) %>% mutate(val = ifelse(n() == 1, as.character(aa), paste(aa, sep = "|"))) %>% 
    select(val) %>% unique()
  div_val <- diversity(t(tbl_ng$Freq), index = "shannon")
  }
  
  div_row <- data.frame(value = as.character(ab_val), shannon_div = div_val, prop_gap = prop_gap, pos = i)
  return(div_row)
}

# Apply the function to each column and combine results
div_table_list <- lapply(1:max(width(tst_4c)), process_column2, align = tst_4c)
div_table4 <- do.call(rbind, div_table_list)
div_table5 <- div_table4 %>% filter(prop_gap <0.5)

range_aa <- IRanges(div_table5$pos)

align_filt <- extractAt(tst_4c, at= range_aa)

align_paste <- unlist(lapply(align_filt, function(x) {
  aastr <- paste(unlist(x), sep = "")
  name <- names(x)
  names(aastr) <- name
  return(aastr)
}))
align_paste2 <- AAStringSet(align_paste)
gaps_df <- as.data.frame(names(align_paste2)) %>% mutate(ngaps = str_count(as.character(align_paste2), "-")) %>%
  mutate(prop_gaps = ngaps/as.numeric(unique(width(align_paste2)))) %>% filter(prop_gaps < 0.3) # columns (positions) which are less than 30% gaps
align_tfilt <- align_paste2[which(names(align_paste2) %in% gaps_df$`names(align_paste2)`)]

writeXStringSet(align_tfilt, "filtalign.fasta")