setwd("~/Desktop/PhD/mg_duf3494/8a_pfam_ab_2/data/")
# e-value filtering & renaming
prok_data <- read.csv("/path/to/folder/prok_duf3494.csv") %>% mutate(unique_name = paste("seq_", (row_number()), sep ="")) %>% 
  mutate(domain = "prokaryotic") # rename sequences, repeat for uc 
prok_seqs <- readAAStringSet("prok_duf3494_genes.fasta") # do for uc

#Get sequence subsets using the list of names in prok data and repeatedly calling prok seqs
prok_subseq <- AAStringSet()
for(i in 1:nrow(prok_data)){
  seq_name <- prok_data$target_name[i]
  start_ind <- prok_data$env_from[i]
  fin_ind <- prok_data$env_to[i]
  find_seq <- prok_seqs[which(grepl(seq_name, names(prok_seqs))==TRUE)]
  new_seq <- AAStringSet(substring(find_seq, start_ind, fin_ind))
  names(new_seq) <- prok_data$unique_name[i]
  prok_subseq <- c(prok_subseq, new_seq)
}
#When euk data is read in: 
euk_data <- euk_data %>% 
  mutate(target_name2 = str_replace_all(target_name,
                                        pattern = c('\\|'= '_', '\\+'= '_', '\\:'= '_', '\\['= '(', '\\]'= ')', '\\-'= '_'))) %>%
  mutate(target_name3 = str_extract(target_name2, "^[^ ]*" ), 
         target_name3 = str_remove(target_name3, "___.*"))
nes <- as.data.frame(names(euk_seqs)) %>% 
  mutate(target_name=str_replace_all(names(euk_seqs), 
                                     pattern = c('\\|'= '_', '\\+'= '_', '\\:'= '_', '\\['= '(', '\\]'= ')', '\\-'= '_'))) %>%
  mutate(target_name2 = str_extract(target_name, "^[^ ]*" ), 
         target_name2 = str_remove(target_name2, "___.*"))
names(euk_seqs) <- nes$target_name2
euk_subseq <- AAStringSet()
for(i in 1:nrow(euk_data)){
  seq_name <- euk_data$target_name3[i]
  start_ind <- euk_data$env_from[i]
  fin_ind <- euk_data$env_to[i]
  find_seq <- euk_seqs[which(grepl(seq_name, names(euk_seqs))==TRUE)]
  new_seq <- AAStringSet(substring(find_seq, start_ind, fin_ind))
  names(new_seq) <- euk_data$unique_name[i]
  euk_subseq <- c(euk_subseq, new_seq)
}
euk_data <- euk_data %>% select(-c(target_name2, target_name3))

# Once euk, unclassified and prokaryotic data have been read in 
full_data <- as.data.frame(rbind(prok_data, euk_data, uc_data)) %>% mutate(pf_len = env_to-env_from) %>% 
  filter(evalue<1e-10, pf_len>150) # Filter for hit quality and length


full_seqs <- c(prok_subseq, euk_subseq, uc_subseq)

filter_seqs <- full_seqs[which(names(full_seqs) %in% full_data$unique_name)]
writeXStringSet(filter_seqs, "filtered_duf3494.fasta")
write.csv(full_data, "filtered_metadata.csv")

tp_md <- read.csv("metadata.csv") %>% select(c("label_env", "prj_accession"))

full_data2 <- full_data %>% left_join(tp_md, by=join_by(x$accession==y$prj_accession))
