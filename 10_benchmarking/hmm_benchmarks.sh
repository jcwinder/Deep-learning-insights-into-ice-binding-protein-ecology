#!/bin/bash -e
#SBATCH slurm/HPC parameters
#SBATCH --mem=488G

cd /path/to/hmm_ecotype
module add hmmer/3.3
module load clustal-omega/1.2.4
clustalo -i env0_train.fasta --auto -o env0_train_aln.fasta -v 
hmmbuild ecotype_0.hmm env0_train_aln.fasta
clustalo -i env1_train.fasta --auto -o env1_train_aln.fasta -v 
hmmbuild ecotype_1.hmm env1_train_aln.fasta
clustalo -i env2_train.fasta --auto -o env2_train_aln.fasta -v 
hmmbuild ecotype_2.hmm env2_train_aln.fasta
clustalo -i env3_train.fasta --auto -o env3_train_aln.fasta -v 
hmmbuild ecotype_3.hmm env3_train_aln.fasta
clustalo -i env4_train.fasta --auto -o env4_train_aln.fasta -v 
hmmbuild ecotype_4.hmm env4_train_aln.fasta

# Compare all v all: 
hmmsearch --tblout ecotype_0_results.txt all_eco.hmm env0_test.fasta
hmmsearch --tblout ecotype_1_results.txt all_eco.hmm env1_test.fasta
hmmsearch --tblout ecotype_2_results.txt all_eco.hmm env2_test.fasta
hmmsearch --tblout ecotype_3_results.txt all_eco.hmm env3_test.fasta
hmmsearch --tblout ecotype_4_results.txt all_eco.hmm env4_test.fas