import csv
from Bio import SeqIO
import sys
import os
import pandas as pd

sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)
data = pd.read_csv('/gpfs/home/hwe21ndu/mg_duf3494/metadata/1/test_process_2/uc_index.csv')
accession_ext = os.environ["sample_id"].strip('"')
sub_acc=data.loc[data['filepath']==accession_ext]
gene_numbers=sub_acc["gene_number"]
fasta_fp=accession_ext
sequences = []
with open(fasta_fp, 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        found = False
        for gene_number in gene_numbers:
            if record.id == gene_number:
                found = True
                break  # exit the loop as soon as a match is found
        if found:
            print(f"Sequence found: {record.id}")
            sequences.append(record)   # add it to the list of sequences
                

output_filepath = f"/gpfs/home/hwe21ndu/mg_duf3494/metadata/1/test_process_2/output_uc/{os.getenv('SLURM_ARRAY_TASK_ID')}.fasta"
with open(output_filepath, 'w') as output_file:
    SeqIO.write(sequences, output_file, 'fasta')