#!/bin/bash -e
#SBATCH #HPC/slurm parameters here
#SBATCH --mem=450G

module load MMSeqs2/14
mmseqs databases UniRef90 /path/to/uniref90 tmp
mmseqs createtaxdb UniRef90 tmp
mmseqs createdb  /path/to/all_dufgenes.fasta queryDB
mmseqs taxonomy queryDB uniref90 taxonomyResult tmp

# Additionally annotate marine sequences using MARFerret: 
module load diamond/2.0.10
diamond blastp -b 100 -c 1 -d ' /path/to/MarFERReT.v1.1.1.dmnd' -e 1e-5 --top 10 -f 102 -q ' /path/to/ms_dufgenes.fasta' -o ' /path/to/mar_out.lca.tab'