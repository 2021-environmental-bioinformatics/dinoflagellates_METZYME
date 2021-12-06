#!/bin/bash
#SBATCH --partition=scavenger            # Queue selection
#SBATCH --job-name=diamond_coassembly_faa_k1    # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=diamond_coassembly_faa_k1%j.log  # Standard output/error
#export OMP_NUM_THREADS=12

module load anaconda/5.1
source activate diamond

dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA/dino_metzyme_annotated_coassembly

diamond blastp -q ${dir}/dino_metzyme_annotated_coassembly.faa -d PhyloDB --evalue 0.001 -k 1 -o dino_metzyme_annotated_coassembly_diamond_out.tsv

#Flags
#The -q	flag indicates the query sequence
#The -d	flag points to the diamond-formatted indexed database, made using the diamond makedb function
#The --evalue specifies	the false positive rate
#The -k	flag specifies the number of hits to show. We chose one, since we only wanted the top-scoring hit. 
#The -o	flag specifies the name of the output file
