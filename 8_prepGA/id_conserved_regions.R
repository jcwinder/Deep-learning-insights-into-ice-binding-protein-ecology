# ID conserved regions of trimmed alignment
library(Biostrings)
library(vegan)
library(ggplot2)

align_tfilt <- readAAStringSet("/path/to/folder/filtalign.fasta") # from the clean_align.R script

process_column <- function(i, align) {
  aa <- as.character(subseq(align, start = i, end = i))
  tbl <- as.data.frame(table(aa))
  sumfreq <- sum(tbl$Freq)
  prop_gap <- ifelse(any(tbl$aa == "-"), tbl$Freq[which(tbl$aa == "-")]/sumfreq, 0)
  tbl_ng <- tbl %>% filter(aa != "-")
  ab_val <- tbl_ng %>% filter(Freq == max(Freq)) %>% mutate(val = ifelse(n() == 1, as.character(aa), paste(aa, sep = "|"))) %>% 
    select(val) %>% unique()
  ab_val <- as.character(ab_val[[1]])
  div_val <- diversity(t(tbl_ng$Freq), index = "shannon")
  div_row <- data.frame(value = ab_val, shannon_div = div_val, prop_gap = prop_gap, pos = i)
  return(div_row)
}

div_table_list <- lapply(1:max(width(align_tfilt)), process_column, align = align_tfilt)
div_table3 <- do.call(rbind, div_table_list)
div_tbl <- div_table3 %>% mutate(cut_n = cut_number(pos, n=43, labels = FALSE)) %>% group_by(cut_n) %>% 
  filter(prop_gap < 0.5) %>% mutate(sd_mean = mean(shannon_div)) # can change threshold: sequences composed of < 50% gaps

# for every window of n AAs, by steps of 1, measure the average shannon diversity
rolling_av <- data.frame()
for(i in 1:nrow(div_table3)){
  positions <- div_table3 %>% filter(pos >= i & pos <= i+8) # vary +8 to whatever window size you want
  mean_sd <- mean(positions$shannon_div)
  full_row <- as.data.frame(t(c(mean_sd, i)))
  rolling_av<- rbind(full_row, rolling_av)
}
colnames(rolling_av) <- c("shannon_div", "pos")
#ggplot(data = rolling_av, aes(x=pos, y = shannon_div)) + geom_point() + theme_classic() + geom_hline(yintercept=1.5) # plot diversity 

#Export 2 spreadsheets: alignment of only the conserved regions, alignment of all the regions, both need labels
seq_labels <- read.csv("/path/to/environment/labels.csv") #environmental labels, based on metadata information

# full alignment
full_align <- as.data.frame(as.character(align_tfilt)) %>% rename("seq" = 1) %>% mutate(names = rownames(.)) %>% 
  mutate(v=str_split(seq,"(?=.)")) %>% unnest_wider(v, names_sep="_") %>% select(-seq) %>% select(where(~!all(is.na(.)))) %>% 
  left_join(seq_labels, by = join_by(x$names==y$names)) %>% 
  mutate(env_lbl2 = case_when(env_id == "0" ~ "frozen_sediment",
                              env_id == "1" ~ "rock",
                              env_id == "2" ~ "subsurface",
                              env_id == "3" ~ "polar_marine",
                              env_id == "4" ~ "glacier_ice")) 


# conserved regions, ID'd by plotting & setting a shannon diversity threshold

seq_tosub <- c(5:10, 17:20, 26:30, 71:83, 94:98, 105:109, 124:127, 132:135, 142:145, 151:154, 169:173, 176:179, 185) 
  
sub_align <- IRanges(seq_tosub)
align_filt <- extractAt(align_tfilt, at= sub_align)
align_paste <- unlist(lapply(align_filt, function(x) {
  aastr <- paste(unlist(x), sep = "")
  name <- names(x)
  names(aastr) <- name
  return(aastr)
}))
align_paste2 <- AAStringSet(align_paste)

cons_align <- as.data.frame(as.character(align_paste2)) %>% rename("seq" = 1) %>% mutate(names = rownames(.)) %>% 
  mutate(v=str_split(seq,"(?=.)")) %>% unnest_wider(v, names_sep="_") %>% select(-seq) %>% select(where(~!all(is.na(.)))) %>% 
  left_join(seq_labels, by = join_by(x$names==y$names)) %>%
  mutate(env_lbl2 = case_when(env_id == "0" ~ "frozen_sediment",
                              env_id == "1" ~ "rock",
                              env_id == "2" ~ "subsurface",
                              env_id == "3" ~ "polar_marine",
                              env_id == "4" ~ "glacier_ice",))


write.csv(full_align, "/path/to/folder/filtered_full_alignment.csv")
write.csv(cons_align, "/path/to/folder/conserved_regions_only.csv")
  