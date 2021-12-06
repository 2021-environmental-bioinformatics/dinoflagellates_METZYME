#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch_pf3             # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=selkassas@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch_pf3%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate hmmer

hmmsearch --domtblout=dino_metzyme_annotated_coassembly_proteorhodopsin_out -E 0.001 proteorhodopsin.hmm dino_metzyme_annotated_coassembly.faa

hmmsearch --domtblout=dino_metzyme_annotated_coassembly_ISIP2A_out -E 0.001 ISIP2A.hmm dino_metzyme_annotated_coassembly.faa

hmmsearch --domtblout=dino_metzyme_annotated_coassembly_ISIP3_out -E 0.001 ISIP3.hmm dino_metzyme_annotated_coassembly.faa
