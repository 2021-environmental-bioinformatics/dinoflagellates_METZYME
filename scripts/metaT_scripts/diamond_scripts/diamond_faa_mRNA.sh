#!/bin/bash
#SBATCH --partition=compute            # Queue selection
#SBATCH --job-name=diamond_faa_k1    # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=12               # Number of CPU cores per task
#SBATCH --mem=50gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=diamond_faa_k1%j.log  # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate diamond

diamond makedb --in /vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/PhyloDB/phylodb_1.076.pep.fa.gz -d PhyloDB

#Flags
#The --in flag indicates the input database file. MUST BE A FASTA FILE
#The -d flag give that indexed database a name to refer to in the blastp scripts below

#40m
diamond blastp -q 30B8Z_S11_001_40m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B8Z_S11_001_40m_mRNA_diamond_faa_out.tsv 
#70m
diamond blastp -q 30B90_S12_001_70m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B90_S12_001_70m_mRNA_diamond_faa_out.tsv 
#380m
diamond blastp -q 30B91_S28_001_380m_mRNA_annotated.faa -d PhyloDB --evalue 0.001 -k 1 -o 30B91_S28_001_380m_mRNA_diamond_faa_out.tsv

#Flags
#The -q flag indicates the query sequence
#The -d flag points to the diamond-formatted indexed database, made using the diamond makedb function
#The --evalue specifies the false positive rate
#The -k flag specifies the number of hits to show. We chose one, since we only wanted the top-scoring hit. 
#The -o flag specifies the name of the output file
