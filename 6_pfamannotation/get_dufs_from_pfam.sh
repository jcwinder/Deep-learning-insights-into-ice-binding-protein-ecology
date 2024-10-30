#!/bin/bash -e
#SBATCH slurm/HPC parameters
#SBATCH --mem=240G

module add R/4.3.1
cd /path/to/folder
Rscript /path/to/folder/get_dufs_from_pfam.R "/path/to/folder/prok_pfam_comb.pfam" "/path/to/folder/prok_duf3494.csv" #provide paths to concatenated pfam file & output name
# need to run for uc and euk as well