#!/bin/bash -e
#SBATCH slurm/HPC parameters
#SBATCH --mem=128G

module add hmmer/3.3
cd /path/to/folder
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/folder/pf_annotate_md.txt)
sample_id=$(echo "$line")
sample_dir=$(echo "$sample_id" | awk -F'/' '{print $9}')
domain1=$(echo "$sample_id" | awk -F'm_' '{print $2}')
domain=$(echo "$domain1" | awk -F'_' '{print $1}')
echo "directory is $sample_dir and domain is $domain"
if [ -f "$sample_id" ]; then
    echo "file $sample_id exists"
    echo "processing $sample_id"
    cp "/path/to/db/Pfam-A.hmm" "/tmp"
	hmmsearch --domtblout "/path/to/folder/$sample_dir/${domain}_hmmsearch3.pfam" "/tmp/Pfam-A.hmm" "$sample_id"
	rm -f "/tmp/Pfam-A.hmm"
else 
    echo "$sample_id is not a file"
fi