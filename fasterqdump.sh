#!/bin/bash -e
#SBATCH # SLURM/HPC commands here
#SBATCH --mem=495G

module add sra # this is the SRAtools
module add python/anaconda/2020.11/3.8
module add biopython
cd /path/to/files/sra_files
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/fasterq_md.txt) # if running as an array job (recommended for large # of samples)
sample_path=$(echo "$line") #filepath
fasterq-dump $sample_path -O $sample_path
