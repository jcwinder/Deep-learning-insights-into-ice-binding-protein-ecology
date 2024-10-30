#!/bin/bash -e
#SBATCH # slurm/HPC commands here
#SBATCH --mem=498G

module add SPAdes/3.15.5
module add python/anaconda/2020.11/3.8

cd /path/to/folder
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/contig_assemble_md.txt) # if running as an array job, contains full filenames
sample_id=$(echo $line)
sample_name=$(echo "$sample_id" | awk -F'/' '{print $9}') # may be a number other than $9, depends on your filepath
   
    # Check if the sample_id directory exists
   if [ -d "$sample_id" ]; then
        if [ -f "$sample_id/scaffolds.fasta" ]; then
            echo "contigs found for $sample_id, exiting"
            continue  # Move to the next sample
		elif [ -f "$sample_id/trimmed/metaspades/spades.log" ]; then
            cd "$sample_id/trimmed/metaspades"
            if [[ -f "$sample_id/trimmed/trimmed_${sample_name}_1.fastq.gz" ]] && [[ -f "$sample_id/trimmed/trimmed_${sample_name}_2.fastq.gz" ]]; then
                metaspades.py \
                --restart-from last \
                --threads=64 \
                --memory=498 \
                -o "$sample_id/trimmed/metaspades"
                echo "Processed $sample_id, gzipping"
                mv "$sample_id/trimmed/metaspades/scaffolds.fasta" "$sample_id"
                mv "$sample_id/trimmed/trimmed_${sample_name}_1.fastq.gz" "fastq_dump"
                mv "$sample_id/trimmed/trimmed_${sample_name}_2.fastq.gz" "fastq_dump"
                tar -czvf "$sample_id/archive.tar.gz" "$sample_id/trimmed"
                rm -r "$sample_id/trimmed"          
            else
                echo "$sample_id only has single-end data"
                continue
            fi
            cd "/path/to/folder"  # Go back to the main directory
       else
            echo "Directory $sample_id/trimmed/metaspades is empty or doesn't exist"
            mkdir -p "$sample_id/trimmed/metaspades"
            cd "$sample_id/trimmed/metaspades"
           if [[ -f "$sample_id/trimmed/trimmed_${sample_name}_1.fastq.gz" ]] && [[ -f "$sample_id/trimmed/trimmed_${sample_name}_2.fastq.gz" ]]; then
                 metaspades.py \
                -1 "$sample_id/trimmed/trimmed_${sample_name}_1.fastq.gz" \
                -2 "$sample_id/trimmed/trimmed_${sample_name}_2.fastq.gz" \
                --threads=64 \
                --memory=498 \
                --only-assembler \
                -o "$sample_id/trimmed/metaspades"
                echo "Processed $sample_id, gzipping"
                mv "$sample_id/trimmed/metaspades/scaffolds.fasta" "$sample_id"
                mv "$sample_id/trimmed/trimmed_${sample_name}_1.fastq.gz" "fastq_dump"
                mv "$sample_id/trimmed/trimmed_${sample_name}_2.fastq.gz" "fastq_dump"
                tar -czvf "$sample_id/archive.tar.gz" "$sample_id/trimmed"
                rm -r "$sample_id/trimmed"
            else
                echo "$sample_id only has single-end data"
                continue
            fi
            cd "/path/to/folder"  # Go back to the main directory
        fi
    else
        echo "Directory $sample_id doesn't exist"
    fi
