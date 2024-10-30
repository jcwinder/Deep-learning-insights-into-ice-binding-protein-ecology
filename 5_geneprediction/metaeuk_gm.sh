#!/bin/bash -e
#SBATCH slurm/HPC parameters
#SBATCH --mem=88G

module add metaeuk/3-8dc7e0b

cd /path/to/folder
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/run_kraken_md.txt) # contains accessions
sample_id=$(echo "$line")
if [ -d "$sample_id" ]; then
	if [ -f "$sample_id/eukaryotic.fasta" ]; then
		cd "$sample_id"
        echo "modelling eukaryotic genes for $sample_id"
        metaeuk easy-predict "$sample_id/eukaryotic.fasta" \
        /path/to/db/uniprot_trembl \
        m_euk_genes "$sample_id" --threads 24
    fi
fi
