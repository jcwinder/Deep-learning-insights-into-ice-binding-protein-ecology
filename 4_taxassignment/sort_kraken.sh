#!/bin/bash -e
#SBATCH #HPC/slurm parameters here
#SBATCH --mem=450G

#sort kraken files into prokaryotic, eukaryotic, and unclassified
module add R/4.3.1
module add python/anaconda/2020.11/3.8

cd /path/to/file
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/sort_kraken_md.txt) # filepaths, if running as array job
sample_id=$(echo "$line")
base_name=$(echo "$sample_id" | awk -F'/' '{print $9}') # $9 will vary depending on filepath structure
if [ -d "$sample_id" ]; then
	if [ -f "$sample_id/kraken_output.txt" ]; then
		echo "$sample_id/kraken_output.txt output found, proceeding to sort"
		echo "$base_name"
		Rscript /path/to/file/sort_kraken.R "$sample_id/kraken_output.txt" "$base_name" "$sample_id"
	else
		echo "$sample_id/kraken_output.txt not found"
	fi
else
	echo "$sample_id directory not found"
fi