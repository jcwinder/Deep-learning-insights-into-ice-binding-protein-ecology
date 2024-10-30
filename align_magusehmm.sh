#!/bin/bash -e
#SBATCH -p hmem-754   # Partition
#SBATCH -J clustalo       # Job name
#SBATCH -o clustalo_%j.out          # Standard output messages
#SBATCH -e clustalo_%j.err          # Standard error messages
#SBATCH --time=7-00:00:00 #7 day limit
#SBATCH --mem=200G

module add FastTree/2.1.11

cd /path/to/file
/path/to/installs/.local/bin/magus -i /path/to/file/subsample_duf3494.fasta -o /path/to/file/subsample_duf3494_backbone.fasta  --recurse false

FastTreeMP /path/to/file/subsample_duf3494_align.fasta > /path/to/file/subsample_duf3494.tree
python3 /path/to/installed/run_upp.py -s /path/to/file/duf3494_seqs.fasta -m amino \
-a /path/to/file/subsample_duf3494_backbone.fasta \
-t /path/to/file/subsample_duf3494.tree -d /path/to/file/magus_ehmm \
-o pf11999_magusehmm

