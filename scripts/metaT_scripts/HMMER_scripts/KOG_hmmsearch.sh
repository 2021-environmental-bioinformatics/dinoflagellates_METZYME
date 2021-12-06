#!/bin/bash
#SBATCH --partition=compute              # Queue selection
#SBATCH --job-name=hmmsearch_kog             # Job name
#SBATCH --mail-type=ALL                  # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kaabbott@whoi.edu   # Where to send mail
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=36               # Number of CPU cores per task
#SBATCH --mem=180gb                      # Job memory request
#SBATCH --time=24:00:00                  # Time limit hrs:min:sec
#SBATCH --output=hmmsearch_kog3%j.log         # Standard output/error
#export OMP_NUM_THREADS=16

module load anaconda/5.1
source activate hmmer

kog_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/databases/KOG
data_dir=/vortexfs1/omics/env-bio/collaboration/dinoflagellates_METZYME/data/metaT_trimmed_reads/fasta_files/paired/mRNA

cd $data_dir

echo $data_dir


hmmsearch -o hmmer_output/dino_metzyme_annotated_coassembly_KOG2348_out -E 0.001 ${kog_dir}/KOG2348.hmm dino_metzyme_annotated_coassembly.faa
