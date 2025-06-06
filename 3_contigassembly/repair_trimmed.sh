#!/bin/bash -e
#SBATCH HPC/SLURM commands here
#SBATCH --mem=500G

module add bbmap/37.28
module add SPAdes/3.15.5
module add python/anaconda/2020.11/3.8
cd /path/to/folder
line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /path/to/folder/repair_trimmed_md.txt) # list of spades.log files to identify broken sequences
sample_id=$(echo $line)
sample_name=$(echo "$sample_id" | awk -F'/' '{print $9}') # $9 may vary depending on your filepath structure
if  grep -q "Unequal number of read-pairs detected" "$sample_id"; then
	echo "unequal read pairs detected in $sample_id, repairing"
	repair.sh in1="/path/to/folder/$sample_name/trimmed/trimmed_${sample_name}_1.fastq.gz" \
	 in2="/path/to/folder/$sample_name/trimmed/trimmed_${sample_name}_2.fastq.gz" \
	 out=stdout.fq \
	 outsingle="/path/to/folder/$sample_name/trimmed/broken.fastq.gz" | reformat.sh \
	 in=stdin.fq out1="/path/to/folder/$sample_name/trimmed/${sample_name}_1_fixed.fastq" \
	 out2="/path/to/folder/$sample_name/trimmed/${sample_name}_2_fixed.fastq" interleaved addcolon
	 rm -r "/path/to/folder/$sample_name/trimmed/metaspades"
	mkdir "/path/to/folder/$sample_name/trimmed/metaspades"
	cd "$sample_name/trimmed/metaspades"
	metaspades.py \
                -1 "/path/to/folder/$sample_name/trimmed/${sample_name}_1_fixed.fastq" \
                -2 "/path/to/folder/$sample_name/trimmed/${sample_name}_2_fixed.fastq" \
                --memory=500 \
                --only-assembler \
                -o "/path/to/folder/$sample_name/trimmed/metaspades"
                echo "Processed $sample_id, gzipping"
                mv "/path/to/folder/$sample_name/trimmed/metaspades/scaffolds.fasta" "$sample_id"
                mv "/path/to/folder/$sample_name/trimmed/${sample_name}_1_fixed" "/path/to/folder"
                mv "/path/to/folder/$sample_name/trimmed/${sample_name}_2_fixed" "/path/to/folder"
                tar -czvf "/path/to/folder/$sample_name/archive.tar.gz" "$sample_id/trimmed"
                #rm -r "/path/to/folder/$sample_id/trimmed"  # ONLY IF DESIRED       
else
	echo "$sample_id is broken for another reason"
fi
