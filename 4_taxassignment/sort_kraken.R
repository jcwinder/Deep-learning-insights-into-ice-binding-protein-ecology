#Sorts kraken output into 3 streams: unclassified, prokaryotic (bact + arch), and eukaryotic
library(Biostrings)
# Set the working directory
args <- commandArgs(trailingOnly = TRUE)

# Define the file path
file_id <- commandArgs(TRUE)[1]  # The first argument will be the file path
sample_directory <- commandArgs(TRUE)[2]
output_directory <- commandArgs(TRUE)[3]

# Helper function to retrieve taxonomic information for a given taxid
get_tax_info <- function(row) {
  classify <- classification(row$taxid, "ncbi")
  return(classify)
}

#Helper function for obtaining classifiable information from the unclassified reads
proc_unclass <- function(data_frame) {
  for (row in 1:nrow(data_frame)) {
    if (grepl(" ", data_frame[row, "kmer_mapping"])) {
      tmp_class <- as.data.frame(strsplit(data_frame[row, "kmer_mapping"], " "))
      tmp_class2 <- apply(tmp_class, 2, function(x) {
        sub(":.*", "", x)
      })
      paste_tmp <- as.data.frame(tmp_class2[which(tmp_class2[, 1] != "0"), ])
      apply_paste <- paste(paste_tmp[,1], collapse = ",")
      data_frame[row, "kmer_mapping"] <- apply_paste
      data_frame[row, "classified"] <- "M"
    }
  }
  return(data_frame)
}

# Main function to process kraken output
process_kraken <- function(filepath, verbose = TRUE) {
  # Read the kraken output data from the file
  data <- read.table(filepath, sep = "\t", header = FALSE)
  colnames(data) <- c("classified", "contig", "taxid", "length", "kmer_mapping")
  
  # Filter classified and unclassified data
  classified <- data[data$classified == "C", ]
  unclassified <- data[data$classified == "U", ]
  
  # Split unclassified into classifiable
  print("classifying ambiguous sequences")
  tmp_unclass <- proc_unclass(unclassified)
  print("sorting ambiguous sequences")
  for (i in 1:nrow(tmp_unclass)) {
    # Extract the kmer_mapping value for the current row
    kmer_mapping <- tmp_unclass[i, "kmer_mapping"]
    
    # Check if kmer_mapping contains ","
    if (grepl(",", kmer_mapping)) {
      temp_split <- as.data.frame(strsplit(kmer_mapping, ","))
      temp_split <- as.data.frame(temp_split[-which(grepl("A", temp_split[, 1])), ])
      
      if (is.na(temp_split[1, 1]) == TRUE) {
        tmp_unclass[i, "classified"] <- "U"
      } else if (length(unique(temp_split[, 1])) > 1) {
        print("checking for domain congruence")
        
        # Use your Python script here to obtain the taxonomic information
        # Assuming temp_split is a matrix or data frame where each row corresponds to a sequence
        # You can use apply to process each row
        
        # Define a function to process each row of temp_split
        process_taxonomic_info <- function(row) {
          # Use your Python script here to obtain the taxonomic information
          cmd <- "python"
          args <- c("/path/to/file/search_dict.py", row[1]) 
          result <- system2(cmd, args, stdout=TRUE)
          # Process the result as needed (assuming it returns the desired taxonomic information)
          result <- as.data.frame(strsplit(result[[2]], ";"))[2, ]
          true_dom <- result
          # Check if there is only one unique domain in the taxonomic information
          if (length(unique(true_dom[, 1])) == 1) {
            return(row[1])  # Return the consistent taxonomic identifier
          } else {
            return("U")  # Return "U" for unclassified if taxonomic information is incongruent
          }
        }
        # Apply the function to each row of temp_split
        processed_taxids <- apply(temp_split, 1, process_taxonomic_info)
        # Assign the processed taxids back to tmp_unclass
        tmp_unclass[, "taxid"] <- processed_taxids
        
        if (length(unique(true_dom[, 1])) == 1) {
          tmp_unclass[i, "taxid"] <- temp_split[1, 1]
        } else {
          tmp_unclass[i, "classified"] <- "U"
        }
      } else {
        tmp_unclass[i, "taxid"] <- temp_split[1,1]
      }
    } else if ((tmp_unclass[i, "kmer_mapping"]!= "0" && !grepl("0:", tmp_unclass[i, "kmer_mapping"])) && !grepl("A", tmp_unclass[i, "kmer_mapping"])) {
      tmp_unclass[i, "taxid"] <- tmp_unclass[i, "kmer_mapping"]
    } else if (grepl("A", tmp_unclass[i, "kmer_mapping"])) {
      tmp_unclass[i, "classified"] <- "U"
    }
  }
  
  
  classified <- as.data.frame(rbind(classified, tmp_unclass[which(tmp_unclass$classified == "M"), ]))
  unclassified <- tmp_unclass[which(tmp_unclass$classified == "U"), ]
  # Get unique taxids from the classified data
  unique_taxids <- unique(classified$taxid)
  # Create a dictionary-like structure to store taxid-domain mapping
  taxid_to_domain <- list()
  
  # Loop through unique taxids to retrieve taxonomic information
  # Number of taxids to process at a time
  for (taxid in unique_taxids) {
    cmd = "python"
    args = c("/path/to/file/search_dict.py", taxid)
    result <- system2(cmd, args, stdout=TRUE)
    result <- result[[2]]
    result <- as.data.frame(strsplit(result, ": "))[2, ]
    taxid_to_domain[[as.character(taxid)]] <- result
  }
  # Add taxonomic domain information to the classified data
  classified$domain <- sapply(classified$taxid, function(taxid) {
    if (as.character(taxid) %in% names(taxid_to_domain)) {
      return(taxid_to_domain[[as.character(taxid)]])
    } else {
      return(NA)
    }
  })
  
  # Filter data based on taxonomic domain and create output list
  prokaryotic <- classified[classified$domain %in% c("2", "2157"), ]
  eukaryotic <- classified[classified$domain %in% c("2759"), ]
  virus <- classified[classified$domain %in% c("10239"), ]
  unclass.class <- classified[which(classified$domain == "NA"), ]
  unclassified <- as.data.frame(rbind(unclassified, unclass.class[,-which(names(unclass.class)=="domain")]))
  output_list <- list(prokaryotic, eukaryotic, virus, unclassified)
  
  # Return the output list
  return(output_list)
}

#Output
result <- process_kraken(file_id, verbose = TRUE)
prokaryotic <- result[[1]]
eukaryotic <- result[[2]]
viral <- result[[3]]
unclass <- result[[4]]


# Split fasta files based on this output

# Read the fasta file
fasta_path <- paste0("/path/to/file/", sample_directory, "/scaffolds.fasta")
fasta_gz_path <- paste0(fasta_path, ".gz")

# Check if the .fasta or .fasta.gz file exists and read the correct one
if (file.exists(fasta_path)) {
  print(".fasta found")
  fasta_file <- readDNAStringSet(fasta_path)
} else if (file.exists(fasta_gz_path)) {
  print(".fasta.gz found")
  fasta_file <- readDNAStringSet(fasta_gz_path)
} else {
  stop("Neither scaffolds.fasta nor scaffolds.fasta.gz was found in the directory.")
}
# Create a function to write a fasta file given a data frame and a filename
write_fasta <- function(data_frame, filename, output_dir) {
  selected_contigs <- data_frame$contig
  selected_sequences <- fasta_file[selected_contigs]
  writeXStringSet(selected_sequences, file.path(output_dir, filename), format = "fasta")
}

# Specify the output filenames
prokaryotic_output <- "prokaryotic.fasta"
eukaryotic_output <- "eukaryotic.fasta"
virus_output <- "viral.fasta"
unclass_output <- "unclass.fasta"

# Use the function to write the fasta files
write_fasta(prokaryotic, prokaryotic_output, output_directory)
write_fasta(eukaryotic, eukaryotic_output, output_directory)
write_fasta(viral, virus_output, output_directory)
write_fasta(unclass, unclass_output, output_directory)
