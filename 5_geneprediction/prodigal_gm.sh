#!/bin/bash -e
#SBATCH slurm/HPC parameters here
#SBATCH --mem=48G

module add Prodigal/2.6.3
cd /path/to/folder
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/run_kraken_md.txt) # contains accessions
sample_id=$(echo "$line")
if [ -d "$sample_id" ]; then
    if [ -f "$sample_id/prokaryotic.fasta" ] && [ ! -f "$sample_id/gm_prok_nucl.faa" ]; then
        echo "modelling prokaryotic genes for $sample_id"
        prodigal -i "$sample_id/prokaryotic.fasta" -d "$sample_id/gm_prok_nucl.faa" -a "$sample_id/gm_prok_aa.faa" \
        -p meta -o "$sample_id/gm_prok_coords.gbk"
    fi
    
    if [ -f "$sample_id/unclass.fasta" ] && [ ! -f "$sample_id/gm_unclass_nucl.faa" ]; then
        echo "modelling unclassified genes as prokaryotic for $sample_id"
        prodigal -i "$sample_id/unclass.fasta" -d "$sample_id/gm_unclass_nucl.faa" -a "$sample_id/gm_unclass_aa.faa" \
        -p meta -o "$sample_id/gm_unclass_coords.gbk"
    fi
    
    if [ ! -f "$sample_id/unclass.fasta" ] || [ ! -f "$sample_id/prokaryotic.fasta" ]; then
    	echo "$sample_id has no kraken output"
    fi
    
    if [ -f "$sample_id/gm_prok_nucl.faa" ] || [ -f "$sample_id/gm_unclass_nucl.faa" ]; then
    	echo "gene models already exist for $sample_id"
    fi
    
elif [ ! -d "$sample_id" ]; then
	echo "$sample_id is not a directory"
fi
