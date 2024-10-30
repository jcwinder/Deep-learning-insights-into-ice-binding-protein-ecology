#!/bin/bash -e
#SBATCH # HPC/slurm parameters here
#SBATCH --mem=498G

module add kraken2
cd /path/to/folder
line=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/kraken_md.txt) # containing filepaths to spades output folders
sample_id=$(echo $line)

if [ -f "${sample_id}/scaffolds.fasta" ] && [ ! -f "${sample_id}/kraken_output.txt" ]; then
	echo "scaffolds found for ${sample_id}, proceeding"
	rm -f "${sample_id}/scaffolds.fasta.gz"
	cd ${sample_id}
	echo "gzipping"
	gzip -5 scaffolds.fasta
	echo "running kraken2"
	kraken2 --db /path/to/dbs/kraken_db \
	${sample_id}/scaffolds.fasta.gz \
	 --output ${sample_id}/kraken_output.txt \
	 --threads 64
elif [ -f "${sample_id}/scaffolds.fasta.gz" ] && [ ! -f "${sample_id}/kraken_output.txt" ]; then
	echo "gzipped scaffolds found for ${sample_id}, proceeding"
	cd ${sample_id}
	echo "running kraken2"
	kraken2 --db /path/to/dbs/kraken_db/kraken_db \
	${sample_id}/scaffolds.fasta.gz \
	 --output ${sample_id}/kraken_output.txt \
	 --threads 64
else
	echo "scaffolds not found for ${sample_id}, exiting"
fi