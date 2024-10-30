#!/bin/bash -e
#SBATCH slurm/HPC commands here 
#SBATCH --mem=498G

module add BBMAP/38.86
cd /gpfs/home/hwe21ndu/scratch/sra_files
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/metadata/files/trim_zip_md.txt) # if running as an array
# md should be a file containing ids and file paths

sample_name=$(echo $line)
sample_id=$(echo "$sample_name" | awk -F'/' '{print $7}') # adapt according to your specific filepath structure - may not be $7
echo "$sample_id"
shopt -s globstar  # Enable globstar option
source_directory="cd /path/to/folder/sra_files"
# Check if the directory already exists
if [[ -d "${sample_id}" ]]; then
    echo "Sample directory ${sample_id} already exists."
    # Check if the directory contains a file with "trimmed" in its name
    if ls "${sample_id}"/trimmed/*trimmed* &>/dev/null; then
        echo "Trimmed sequences found: zipping and moving originals."
        gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}_1.fastq"
        gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}_2.fastq"
        mv "/gpfs/path/to/folder/sra_files/${sample_id}/${sample_id}_1.fastq" /path/to/folder/sra_files/dump_fastq
        mv "/path/to/folder/sra_files/${sample_id}/${sample_id}_2.fastq" /path/to/folder/sra_files/dump_fastq
        gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}.fastq"
        mv "/path/to/folder/sra_files/${sample_id}/${sample_id}.fastq" /path/to/folder/sra_files/dump_fastq
        exit 0
    else
        echo "No trimmed sequences found, proceeding with adapter trimming."
        mkdir -p "${sample_id}/trimmed"
        echo "creating trimmed directory"
    fi
else
	echo "${sample_id} is not a directory"
fi

cd "/path/to/folder/sra_files/${sample_id}"


# Adapter trimming: need to account for paired or single reads

if [[ -f "${sample_id}_1.fastq" || -f "${sample_id}_1.fastq.gz" ]] &&
   [[ -f "${sample_id}_2.fastq" || -f "${sample_id}_2.fastq.gz" ]]; then
    # Paired-end processing
    echo "Paired-end processing for sample $sample_id"
    # Add your paired-end processing commands here
    in1="${sample_id}_1.fastq"
    if [[ -f "${sample_id}_1.fastq.gz" ]]; then
        in1="${sample_id}_1.fastq.gz"
    fi
    
    in2="${sample_id}_2.fastq"
    if [[ -f "${sample_id}_2.fastq.gz" ]]; then
        in2="${sample_id}_2.fastq.gz"
    fi

    bbduk.sh \
        in1="$in1" \
        in2="$in2" \
        out1="trimmed/trimmed_${sample_id}_1.fastq" \
        out2="trimmed/trimmed_${sample_id}_2.fastq" \
        k=25 \
        mink=8 \
        ktrim=r \
        ref=/path/to/bbmap/resources/adapters.fa \
        hdist=1 > bbduk_log.txt
    echo "zipping and moving files $out1 and $out2"
    gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}_1.fastq"
    gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}_2.fastq"
    mv "/path/to/folder/sra_files/${sample_id}/${sample_id}_1.fastq" /path/to/folder/sra_files/dump_fastq
    mv "/path/to/folder/sra_files/${sample_id}/${sample_id}_2.fastq" /path/to/folder/sra_files/dump_fastq


elif [[ -f "${sample_id}.fastq" || -f "${sample_id}.fastq.gz" ]]; then
    # Single-end processing
    echo "Single-end processing for sample $sample_id"
    # Add your single-end processing commands here
    in="${sample_id}.fastq"
    if [[ -f "${sample_id}.fastq.gz" ]]; then
        in="${sample_id}.fastq.gz"
    fi

    bbduk.sh \
        in="$in" \
        out="trimmed/trimmed_${sample_id}.fastq" \
        k=25 \
        mink=8 \
        ktrim=r \
        ref=/path/to/bbmap/resources/adapters.fa \
        hdist=1 > bbduk_log.txt
    echo "zipping and moving files"
    gzip -5 "/path/to/folder/sra_files/${sample_id}/trimmed/trimmed_${sample_id}.fastq"
    mv "/path/to/folder/sra_files/${sample_id}/${sample_id}.fastq" /path/to/folder/sra_files/dump_fastq

else
    # No matching files found
    echo "No matching files found for sample $sample_id"
fi

#####
# gzip fastq

